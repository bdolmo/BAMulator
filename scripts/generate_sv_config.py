#!/usr/bin/env python3

import argparse
import glob
import os
import random
import csv
import re
import numpy as np
import pandas as pd
from natsort import natsorted
from pycirclize import Circos
from pycirclize.utils import ColorCycler, load_eukaryote_example_dataset

def parse_bed_file(bed_path, autosomes=False, chromosome=None):
    bed_records = []
    with open(bed_path) as f:
        for line in f:
            if not line.startswith(('browser', 'track')):
                if chromosome and not line.startswith(chromosome + "\t"):
                    continue
                if autosomes and line.startswith(('X\t', 'Y\t')):
                    continue
                parts = line.strip().split('\t')
                record = (parts[0], int(parts[1]), int(parts[2]), parts[3])
                bed_records.append(record)
    return bed_records

def load_gap_regions(genome):
    gap_file = os.path.join(os.path.dirname(__file__), "annotations", f"{genome}.gap.txt.gz")
    df = pd.read_csv(gap_file, sep='\t', comment='#', header=None)
    df.columns = ['bin', 'chrom', 'start', 'end', 'ix', 'type', 'size', 'bridge', 'centromere']
    df = df[df['type'] == 'N']
    gap_records = []
    for _, row in df.iterrows():
        gap_records.append((row['chrom'], row['start'], row['end']))
    return gap_records

def region_overlaps(existing, chr, start, end):
    for c, s, e in existing:
        if c == chr and not (end < s or start > e):
            return True
    return False

def simulate_mechanism(mech, ntmin, ntmax):
    if mech == 'nhej':
        return random_nt_insertion(1, 10), '.'
    elif mech == 'mmbir':
        hom = random_nt_insertion(2, 15)
        return hom + random_nt_insertion(ntmin, ntmax), hom
    elif mech == 'alt_ebr':
        return random_nt_insertion(10, 50), '.'
    else:
        return random_nt_insertion(ntmin, ntmax), '.'

def random_nt_insertion(min_nt, max_nt):
    nucleotides = ['A', 'C', 'G', 'T']
    length = random.randint(min_nt, max_nt)
    sequence = ''
    for _ in range(length):
        sequence += random.choice(nucleotides)
    return sequence

def plot_sv_circos(df, sample_name, chrom_sizes, output_path='circos_sv_plot.png', genome_version="hg19"):
    if genome_version == "hg19":
        chr_bed_file, cytoband_file, _ = load_eukaryote_example_dataset("hg38")
        label = f"{sample_name}\n(hg19 assumed as hg38 layout)"
    else:
        chr_bed_file, cytoband_file, _ = load_eukaryote_example_dataset("hg38")
        label = f"{sample_name}\n(hg38)"

    circos = Circos.initialize_from_bed(chr_bed_file, space=3)
    circos.text(label, deg=315, r=150, size=12)
    circos.add_cytoband_tracks((95, 100), cytoband_file)

    ColorCycler.set_cmap("tab10")
    chr_names = []
    for sector in circos.sectors:
        chr_names.append(sector.name)
    colors = ColorCycler.get_color_list(len(chr_names))
    chr_color = {}
    for name, color in zip(chr_names, colors):
        chr_color[name] = color

    for sector in circos.sectors:
        sector.text(sector.name, r=120, size=10, color="black")
        sector.get_track("cytoband").xticks_by_interval(
            40000000,
            label_size=8,
            label_orientation="vertical",
            label_formatter=lambda v: f"{v / 1000000:.0f} Mb",
        )

    sv_colors = {
        'INS': 'red',
        'DEL': 'blue',
        'DUP': 'green',
        'INV': 'orange',
        'BALANCED_TRANSLOCATION': 'purple',
        'WHOLE_ARM_TRANSLOCATION': 'brown',
        'TRA': 'purple',
    }

    for _, row in df.iterrows():
        sv_type = row['SVTYPE']
        color = sv_colors.get(sv_type, "gray")
        chr1 = row['CHR1']
        start1 = int(row['START1'])
        end1 = int(row['END1'])
        chr2 = row['CHR2']
        start2 = int(row['START2'])
        end2 = int(row['END2'])

        if chr2 == '.' or chr2 == "":
            chr2 = chr1
            start2 = end1
            end2 = end1 + 1

        circos.link((chr1, start1, end1), (chr2, start2, end2), color=color, alpha=0.7, linewidth=0.7)

    fig = circos.plotfig()
    _ = circos.ax.legend(loc="upper right", fontsize=10)
    fig.savefig(output_path, dpi=130)

def main():
    parser = argparse.ArgumentParser(description="Create SV simulations")
    parser.add_argument('-i', '--indir', type=str, required=True)
    parser.add_argument('-b', '--bed', type=str, required=True)
    parser.add_argument('--num_variants', type=int, required=True)
    parser.add_argument('--insertion_rate', type=float, default=0.2)
    parser.add_argument('--deletion_rate', type=float, default=0.2)
    parser.add_argument('--duplication_rate', type=float, default=0.2)
    parser.add_argument('--inversion_rate', type=float, default=0.2)
    parser.add_argument('--translocation_rate', type=float, default=0.2,
                        help="Total rate for translocations (split evenly if the two specialized rates below are not given)")
    parser.add_argument('--balanced_translocation_rate', type=float, default=None,
                        help="Optional. If set, overrides how much of the total goes to balanced translocations")
    parser.add_argument('--whole_arm_rate', type=float, default=None,
                        help="Optional. If set, overrides how much of the total goes to whole-arm (single-junction) translocations")
    parser.add_argument('--chromosome', type=str)
    parser.add_argument('--autosomes', action='store_true')
    parser.add_argument('--minsize', type=int, default=100)
    parser.add_argument('--maxsize', type=int, default=5000)
    parser.add_argument('--minins', type=int, default=100)
    parser.add_argument('--maxins', type=int, default=5000)
    parser.add_argument('--mindist', type=int, default=1_000_000)
    parser.add_argument('--ntmin', type=int, default=1)
    parser.add_argument('--ntmax', type=int, default=20)
    parser.add_argument('--chromfile', type=str)
    parser.add_argument('--mechanism', type=str, default='random')
    parser.add_argument('--proportions', type=str, required=True)
    parser.add_argument('--output', type=str, default='sv_config.tsv')
    parser.add_argument('--seed', type=int, default=None)
    parser.add_argument('--genome_version', choices=['hg19', 'hg38'], default='hg19')
    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)
        np.random.seed(args.seed)

    proportions = []
    for p in args.proportions.split(','):
        proportions.append(float(p))
    if not np.isclose(sum(proportions), 1.0):
        raise ValueError("Clone proportions must sum to 1.0")

    base_types = ['INS', 'DEL', 'DUP', 'INV']
    trans_types = ['BALANCED_TRANSLOCATION', 'WHOLE_ARM_TRANSLOCATION']
    sv_types = base_types + trans_types

    if args.balanced_translocation_rate is None and args.whole_arm_rate is None:
        tra_bal = args.translocation_rate * 0.5
        tra_wat = args.translocation_rate * 0.5
    else:
        tra_bal = float(args.balanced_translocation_rate or 0.0)
        tra_wat = float(args.whole_arm_rate or 0.0)

    sv_weights = [
        float(args.insertion_rate),
        float(args.deletion_rate),
        float(args.duplication_rate),
        float(args.inversion_rate),
        tra_bal,
        tra_wat,
    ]

    total = sum(sv_weights)
    if total <= 0:
        raise ValueError("All SV rates are zero; at least one must be > 0.")
    if not np.isclose(total, 1.0):
        print(f" [warn] SV rates sum to {total:.4f}, auto-normalizing to 1.0")
        normalized_weights = []
        for w in sv_weights:
            normalized_weights.append(w / total)
        sv_weights = normalized_weights

    print(" [info] SV weights =>", dict(zip(sv_types, [f"{w:.4f}" for w in sv_weights])))

    bams = glob.glob(f"{args.indir}/*.bam")
    if not bams:
        raise ValueError("No BAM files found in input directory")

    region_info = parse_bed_file(args.bed, args.autosomes, args.chromosome)
    gap_regions = load_gap_regions(args.genome_version)

    def in_gap(chr, start, end):
        for c, s, e in gap_regions:
            if c == chr and not (end <= s or start >= e):
                return True
        return False

    chrom_sizes = []
    chrom_pattern = re.compile(r"^chr([1-9]|1[0-9]|2[0-2]|X|Y)$")
    if args.chromfile:
        with open(args.chromfile) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split('\t')
                chrom = parts[0]
                if chrom_pattern.match(chrom):
                    chrom_sizes.append((chrom, int(parts[1])))

    def random_position_outside_bed():
        for _ in range(1000):
            chrom, size = random.choice(chrom_sizes)
            pos = random.randint(1000, size - 1000)
            overlaps = False
            for r in region_info:
                if chrom == r[0] and r[1] <= pos < r[2]:
                    overlaps = True
                    break
            if not overlaps and not in_gap(chrom, pos, pos + 1):
                return chrom, pos
        return None, None

    def random_position_outside_bed_fallback(region_info, in_gap_fn, max_tries=2000):
        by_chr = {}
        for c, s, e, _ in region_info:
            if c not in by_chr:
                by_chr[c] = []
            by_chr[c].append((s, e))

        for c in by_chr:
            ivs = sorted(by_chr[c])
            merged = []
            for s, e in ivs:
                if not merged or s > merged[-1][1]:
                    merged.append([s, e])
                else:
                    merged[-1][1] = max(merged[-1][1], e)
            by_chr[c] = merged

        chroms = list(by_chr.keys())
        if not chroms:
            return None, None

        for _ in range(max_tries):
            c = random.choice(chroms)
            ivs = by_chr[c]
            candidates = []
            prev_end = 0
            for s, e in ivs:
                if s - prev_end > 2000:
                    candidates.append((prev_end + 1000, s - 1000))
                prev_end = e
            tail_start = prev_end + 1000
            tail_end = prev_end + 5_000_000
            candidates.append((tail_start, tail_end))

            if candidates:
                a, b = random.choice(candidates)
                if b > a:
                    pos = random.randint(a, b)
                    if not in_gap_fn(c, pos, pos + 1):
                        return c, pos
        return None, None

    def pick_random_region_span(minlen, maxlen):
        for _ in range(1000):
            region = random.choice(region_info)
            chrA, rstart, rend, geneA = region
            svlen = random.randint(minlen, maxlen)
            startA = rstart + max(0, (rend - rstart - svlen) // 2)
            endA = startA + svlen
            if not in_gap(chrA, startA, endA):
                return chrA, startA, endA, geneA, svlen
        return None

    all_records = []
    for bam in natsorted(bams):
        used_regions = []
        clone_counts = []
        for p in proportions:
            clone_counts.append(int(round(args.num_variants * p)))
        total_requested = sum(clone_counts)
        if total_requested > args.num_variants:
            clone_counts[-1] -= total_requested - args.num_variants

        variant_total = 0
        for clone_id, proportion in enumerate(proportions, 1):
            count = clone_counts[clone_id - 1]
            written = 0
            while written < count and variant_total < args.num_variants:
                sv_type = random.choices(sv_types, weights=sv_weights, k=1)[0]

                if sv_type in ['DEL', 'DUP', 'INV']:
                    region = random.choice(region_info)
                    chr, start, _, gene = region
                    svlen = random.randint(args.minsize, args.maxsize)
                    end = start + svlen
                    if region_overlaps(used_regions, chr, start, end) or in_gap(chr, start, end):
                        continue
                    used_regions.append((chr, start, end))
                    chr2, start2, end2 = ".", ".", "."

                elif sv_type == 'INS':
                    region = random.choice(region_info)
                    chr, start, _, gene = region
                    end = start + 1
                    if in_gap(chr, start, end):
                        continue
                    chr2, start2, end2 = ".", ".", "."

                elif sv_type == 'BALANCED_TRANSLOCATION':
                    pickA = pick_random_region_span(args.minsize, args.maxsize)
                    pickB = pick_random_region_span(args.minsize, args.maxsize)
                    if not pickA or not pickB:
                        continue
                    chr, start, end, geneA, svlenA = pickA
                    chr2, start2, end2, geneB, svlenB = pickB

                    if chr == chr2 and abs(((start + end)//2) - ((start2 + end2)//2)) < args.mindist:
                        continue
                    if in_gap(chr, start, end) or in_gap(chr2, start2, end2):
                        continue

                    gene = f"{geneA}__{geneB}"

                elif sv_type == 'WHOLE_ARM_TRANSLOCATION':
                    region = random.choice(region_info)
                    chr, region_start, region_end, gene = region
                    start = (region_start + region_end) // 2
                    end = start + 1
                    chr2, start2 = random_position_outside_bed()
                    if not chr2:
                        chr2, start2 = random_position_outside_bed_fallback(region_info, in_gap)
                        if not chr2:
                            continue
                    end2 = start2 + 1

                ntins, hom = simulate_mechanism(args.mechanism, args.ntmin, args.ntmax)

                ploidy_map = {
                    'DEL': '1',
                    'DUP': '3',
                    'INV': '2',
                    'INS': '.',
                    'BALANCED_TRANSLOCATION': '.',
                    'WHOLE_ARM_TRANSLOCATION': '.'
                }
                ploidy = ploidy_map.get(sv_type, '.')

                svid = f"{gene}_{sv_type.lower()}"

                if sv_type == 'WHOLE_ARM_TRANSLOCATION':
                    add_info = "WHOLE_ARM_SINGLE_JUNCTION"
                elif sv_type == 'BALANCED_TRANSLOCATION':
                    add_info = "BALANCED_TWO_ENDS"
                else:
                    add_info = "."

                if ploidy in ('1', '3', '.'):
                    hap = random.choice(['1', '2'])
                elif ploidy == '2':
                    hap = '1,2'
                else:
                    hap = random.choice(['1', '2'])

                record = [
                    bam, chr, start, end,
                    chr2, start2, end2, sv_type, svid,
                    ploidy, args.mechanism, ntins or '.', hom or '.',
                    f"clone{clone_id}", f"{proportion:.4f}", add_info, hap
                ]
                all_records.append(record)
                written += 1
                variant_total += 1

    with open(args.output, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow([
            'SAMPLE', 'CHR1', 'START1', 'END1', 'CHR2', 'START2', 'END2',
            'SVTYPE', 'SVID', 'PLOIDY', 'SV_MECHANISM', 'NTINS', 'HOMOLOGY',
            'CLONE', 'CLONAL_PROPORTION', 'ADDITIONAL_INFO', 'HAPLOTYPE'
        ])
        for record in all_records:
            writer.writerow(record)

    df = pd.DataFrame(all_records, columns=[
        'SAMPLE', 'CHR1', 'START1', 'END1', 'CHR2', 'START2', 'END2',
        'SVTYPE', 'SVID', 'PLOIDY', 'SV_MECHANISM', 'NTINS', 'HOMOLOGY',
        'CLONE', 'CLONAL_PROPORTION', 'ADDITIONAL_INFO', 'HAPLOTYPE'
    ])

    for sample_name, group_df in df.groupby("SAMPLE"):
        plot_path = f"{sample_name}_circos.png"
        plot_sv_circos(group_df, sample_name, chrom_sizes, output_path=plot_path, genome_version=args.genome_version)

if __name__ == "__main__":
    main()