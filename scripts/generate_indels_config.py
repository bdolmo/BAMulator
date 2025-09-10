#!/usr/bin/env python3

import argparse
import glob
import os
import random
import subprocess
from natsort import natsorted
import numpy as np
from intervaltree import Interval, IntervalTree
import matplotlib.pyplot as plt
import seaborn as sns
import csv

def main():
    parser = argparse.ArgumentParser(description="Generate config file for Indels")
    parser.add_argument('-i', '--indir', type=str, help='Input BAM directory')
    parser.add_argument('-l', '--list', type=str, help='Input BAM list file')
    parser.add_argument('-b', '--bed', type=str, required=True, help='ROI bed')
    parser.add_argument('-g', '--genome', type=str, required=True, help='Reference genome in FASTA format')
    parser.add_argument('--num_variants', type=int, required=True, help='Total number of variants to simulate')
    parser.add_argument('--autosomes', action='store_true', help='Restrict to autosomes')
    parser.add_argument('--insertion_rate', type=float, default=0.45, help='Proportion of insertions')
    parser.add_argument('--deletion_rate', type=float, default=0.45, help='Proportion of deletions')
    parser.add_argument('--delins_rate', type=float, default=0.1, help='Proportion of delins')
    parser.add_argument('--chromosome', type=str, help='Specific chromosome to introduce variants')
    parser.add_argument('--proportions', type=str, required=True, help='Comma-separated list of proportions for clones')
    parser.add_argument('--hom_rate', type=float, default=0.1, help='Rate of homozygous variants (default: 0.1)')
    parser.add_argument('--output', type=str, required=True, help='Output TSV file for variant list')
    args = parser.parse_args()

    if not (args.indir or args.list):
        parser.error("Please provide either --indir or --list.")

    total_rate = args.insertion_rate + args.deletion_rate + args.delins_rate
    if not np.isclose(total_rate, 1.0):
        parser.error("The sum of insertion_rate, deletion_rate, and delins_rate must be 1.0")

    if not os.path.exists(args.bed):
        parser.error(f"ROI BED file not found: {args.bed}")

    if not os.path.exists(args.genome):
        parser.error(f"Reference genome not found: {args.genome}")

    samtools = subprocess.getoutput('which samtools')
    if not samtools:
        parser.error("samtools not found in PATH")

    bams = []
    if args.indir:
        bams = glob.glob(f"{args.indir}/*.bam")
        if not bams:
            parser.error(f"No BAM files found in directory: {args.indir}")
    elif args.list:
        if not os.path.exists(args.list):
            parser.error(f"BAM list file not found: {args.list}")
        with open(args.list) as f:
            bams = f.read().splitlines()
        if not bams:
            parser.error(f"No BAMs listed in: {args.list}")

    with open(args.bed) as f:
        lines = [l for l in f if not l.startswith(('browser', 'track'))]

    if args.autosomes:
        lines = [l for l in lines if not l.startswith(('X', 'Y'))]

    if args.chromosome:
        lines = [l for l in lines if l.startswith(args.chromosome + "\t")]

    if not lines:
        parser.error("No valid regions in BED after filtering")

    regions = parse_bed_data(lines)

    proportions = [float(x) for x in args.proportions.split(',')]
    if not np.isclose(sum(proportions), 1.0):
        parser.error("Proportions must sum to 1.0")

    all_variants = []
    for bam in natsorted(bams):
        sample_name = os.path.basename(bam)
        variants = generate_variants(regions, args.num_variants, args.genome, samtools,
                                     args.insertion_rate, args.deletion_rate, args.delins_rate,
                                     proportions, args.hom_rate)
        for v in variants:
            v['sample'] = bam
        all_variants.extend(variants)

    output_fields = ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'CLONE', 'CLONAL_PROPORTION', 'GENOTYPE', 'ADDITIONAL_INFO']
    with open(args.output, 'w') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(output_fields)
        for v in all_variants:
            writer.writerow([
                v['sample'], v['chrom'], v['pos'], v['ref'], v['alt'], v['type'],
                f"clone{v['clone']}", f"{v['proportion']:.4f}", v['genotype'], "."
            ])

def parse_bed_data(lines):
    regions = []
    for line in lines:
        parts = line.strip().split('\t')
        chr, start, end = parts[0], int(parts[1]), int(parts[2])
        regions.append((chr, start, end))
    return regions

def select_random_region(regions):
    return random.choice(regions)

def get_reference_sequence(chr, position, length, genome, samtools):
    region = f"{chr}:{position}-{position + length - 1}"
    result = subprocess.getoutput(f"{samtools} faidx {genome} {region}")
    lines = result.split('\n')
    return ''.join(lines[1:]).upper() if len(lines) > 1 else ''

def generate_inserted_sequence(chr, position, length, genome, samtools, ref_seq=None):
    if random.random() < 1/3:
        start = max(1, position - random.randint(1, 1000))
        dup_seq = get_reference_sequence(chr, start, length, genome, samtools)
        inserted = ''.join(random.choice('ACGT') if random.random() < 0.1 else b for b in dup_seq)
    else:
        inserted = ''.join(random.choice('ACGT') for _ in range(length))
    if ref_seq and (ref_seq + inserted).upper() == ref_seq.upper():
        inserted = ''.join(random.choice('ACGT') for _ in range(length))
    return inserted.upper()

def generate_indel_size():
    size = np.random.geometric(0.1)
    return min(max(size, 1), 50)

def generate_delins(chr, position, genome, samtools):
    while True:
        deletion_size = generate_indel_size()
        if deletion_size > 1:
            break
    insertion_size = random.randint(1, deletion_size - 1)
    if insertion_size == 1:
        insertion_size = 2
    ref_seq = get_reference_sequence(chr, position, deletion_size, genome, samtools)
    inserted_seq = generate_inserted_sequence(chr, position, insertion_size, genome, samtools, ref_seq)
    alt_seq = (ref_seq[:1] + inserted_seq).upper()
    return ref_seq, alt_seq

def generate_variants(regions, num_variants, genome, samtools,
                      insertion_rate, deletion_rate, delins_rate,
                      proportions, hom_rate):
    variants = []
    used_intervals = {r[0]: IntervalTree() for r in regions}
    num_insertions = int(num_variants * insertion_rate)
    num_deletions = int(num_variants * deletion_rate)
    num_delins = num_variants - num_insertions - num_deletions
    clone_counts = [int(num_variants * p) for p in proportions]

    for idx, proportion in enumerate(proportions):
        for _ in range(clone_counts[idx]):
            if len(variants) >= num_variants:
                break
            region = select_random_region(regions)
            position = random.randint(region[1], region[2])
            indel_size = generate_indel_size()
            interval = Interval(position, position + indel_size)
            if used_intervals[region[0]].overlaps(interval):
                continue

            variant_type = None
            ref_seq, alt_seq = "", ""

            if len(variants) < num_deletions:
                if position <= 1:
                    continue
                ref_seq = get_reference_sequence(region[0], position, indel_size, genome, samtools)
                if indel_size == 1:
                    pre_base = get_reference_sequence(region[0], position - 1, 1, genome, samtools)
                    if not pre_base:
                        continue
                    ref_seq = pre_base + ref_seq
                    alt_seq = pre_base
                    position -= 1
                else:
                    alt_seq = ref_seq[0]
                variant_type = 'deletion'

            elif len(variants) < num_deletions + num_insertions:
                ref_seq = get_reference_sequence(region[0], position, 1, genome, samtools)
                inserted = generate_inserted_sequence(region[0], position, indel_size, genome, samtools, ref_seq)
                alt_seq = (ref_seq + inserted).upper()
                variant_type = 'insertion'

            else:
                ref_seq, alt_seq = generate_delins(region[0], position, genome, samtools)
                variant_type = 'delins'

            if not ref_seq or not alt_seq or ref_seq == alt_seq or 'N' in ref_seq or 'N' in alt_seq:
                continue

            genotype = '1/1' if random.random() < hom_rate else '0/1'
            variants.append({
                'chrom': region[0],
                'pos': position,
                'ref': ref_seq,
                'alt': alt_seq,
                'type': variant_type,
                'clone': idx + 1,
                'proportion': proportion,
                'genotype': genotype
            })
            used_intervals[region[0]].add(interval)

    return variants

if __name__ == "__main__":
    main()
