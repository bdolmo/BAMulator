#!/usr/bin/env python3

import argparse
import glob
import os
import random
import subprocess
import csv
from natsort import natsorted

def main():
    parser = argparse.ArgumentParser(description="Generate config file for panelCNVsim")
    parser.add_argument('-i', '--indir', type=str, required=True, help='Input BAM directory')
    parser.add_argument('-b', '--bed', type=str, required=True, help='ROI bed')
    parser.add_argument('-m', '--mode', type=str, required=True, choices=['partial', 'single', 'multiple', 'whole', 'simul'], help='Mode (partial, single, multiple, whole, simul)')
    parser.add_argument('-c', '--var_type', type=str, required=True, choices=['del', 'dup'], help='CNV class (del, dup) - not used in simul mode')
    parser.add_argument('--minsize', type=int, default=20, help='Minimum size for partial CNV')
    parser.add_argument('--maxsize', type=int, default=200, help='Maximum size for partial CNV')
    parser.add_argument('--minexons', type=int, default=2, help='Minimum number of exons for multiple CNVs')
    parser.add_argument('--maxexons', type=int, default=8, help='Maximum number of exons for multiple CNVs')
    parser.add_argument('--autosomes', action='store_true', help='Simulate only over autosomes')
    parser.add_argument('--ncnv', type=int, default=2, help='Number of CNVs to simulate per sample (simul mode only)')
    parser.add_argument('--proportions', type=str, required=False, help='Comma-separated list of clone proportions (must sum to 1.0)')
    parser.add_argument('--output', type=str, default='cnv_config.tsv', help='Output config TSV')
    args = parser.parse_args()

    if not os.path.exists(args.bed):
        print(f" ERROR: No input ROI bed file found at {args.bed}")
        exit()

    samtools = subprocess.getoutput('which samtools')
    if not samtools:
        print(" ERROR: samtools was not found on path")
        exit()

    bams = glob.glob(f"{args.indir}/*.bam")
    if not bams:
        print(f" ERROR: no BAM files were found in {args.indir} directory!")
        exit()

    with open(args.bed, 'r') as file:
        lines = file.read().splitlines()

    if args.autosomes:
        lines = [line for line in lines if not (line.startswith('X') or line.startswith('Y'))]

    gene_exon_info = parse_bed_data(lines)

    if not args.proportions:
        proportions = [1]
    else:
        proportions = list(map(float, args.proportions.split(',')))
    if not abs(sum(proportions) - 1.0) < 1e-5:
        raise ValueError("Clone proportions must sum to 1.0")

    all_records = []

    for bam in natsorted(bams):
        clone_counts = [int(round(args.ncnv * p)) for p in proportions]
        total_requested = sum(clone_counts)
        if total_requested > args.ncnv:
            clone_counts[-1] -= total_requested - args.ncnv

        for clone_id, (proportion, count) in enumerate(zip(proportions, clone_counts), 1):
            written = 0
            while written < count:
                sample_name = bam

                if args.mode == 'partial':
                    result = get_partial_exon(gene_exon_info, args.minsize, args.maxsize)
                elif args.mode == 'single':
                    result = get_single_exon(gene_exon_info)
                elif args.mode == 'multiple':
                    result = get_multiple_exon(gene_exon_info, args.minexons, args.maxexons)
                elif args.mode == 'simul':
                    events = simulate_multiple_cnvs(gene_exon_info, count, args.minexons, args.maxexons)
                    for e in events:
                        var_type = random.choice(['del', 'dup'])
                        ploidy = '1' if var_type == 'del' else '3'
                        all_records.append([
                            sample_name, e[0], e[1], e[2], '.', '.', '.', var_type.upper(), e[3], ploidy, '.', '.', '.', f"clone{clone_id}", f"{proportion:.4f}", args.mode
                        ])
                    break
                else:
                    continue

                if result is None:
                    continue

                var_type = args.var_type
                ploidy = '1' if var_type == 'del' else '3'
                all_records.append([
                    sample_name, result[0], result[1], result[2], '.', '.', '.', var_type.upper(), result[3], ploidy, '.', '.', '.', f"clone{clone_id}", f"{proportion:.4f}", args.mode
                ])
                written += 1

    with open(args.output, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(['SAMPLE', 'CHR1', 'START1', 'END1', 'CHR2', 'START2', 'END2', 'SVTYPE', 'ID', 'PLOIDY', 'SV_MECHANISM', 'NTINS', 'HOMOLOGY', 'CLONE', 'CLONAL_PROPORTION', 'ADDITIONAL_INFO'])
        writer.writerows(all_records)

def parse_bed_data(lines):
    gene_exon_info = []
    for line in lines:
        parts = line.split('\t')
        if line.startswith("browser") or line.startswith("track name"):
            continue
        chr, start, end, gene = parts[0], int(parts[1]), int(parts[2]), parts[3]
        gene_exon_info.append((chr, start, end, gene))
    return gene_exon_info


def get_partial_exon(gene_exon_info, min_size, max_size):
    index = random.randint(0, len(gene_exon_info) - 1)
    chr, start, end, info = gene_exon_info[index]
    new_start = start + 10
    new_end = new_start + random.randint(min_size, max_size)
    return (chr, new_start, new_end, info)

def get_single_exon(gene_exon_info):
    seen = set()
    while True:
        index = random.randint(0, len(gene_exon_info) - 1)
        chr, start, end, info = gene_exon_info[index]
        
        # Check the neighbors
        valid = True
        
        # Check the previous exon if it's not the first exon
        if index > 0:
            chr_prev, start_prev, end_prev, _ = gene_exon_info[index - 1]
            if chr == chr_prev and start <= end_prev + 500 and (end-start) > 20:
                valid = False
        
        # Check the next exon if it's not the last exon
        if index < len(gene_exon_info) - 1:
            chr_next, start_next, end_next, _ = gene_exon_info[index + 1]
            if chr == chr_next and start_next <= end + 500  and (end-start) > 20:
                valid = False
        
        if valid:
            # Add buffer if the exon is valid and hasn't been seen before
            if (chr, start, end, info) not in seen:
                seen.add((chr, start, end, info))
                start = max(0, start - 500)  # Ensure start does not go below 0
                end += 500
                return (chr, start, end, info)

def get_multiple_exon(gene_exon_info, min_exons, max_exons):
    gene_to_exons = {}
    # Organize exons by gene
    for chr, start, end, info in gene_exon_info:
        if info not in gene_to_exons:
            gene_to_exons[info] = []
        gene_to_exons[info].append((chr, start, end))

    # Filter genes with enough exons
    valid_genes = {gene: exons for gene, exons in gene_to_exons.items() if len(exons) >= min_exons}
    
    attempts = 0
    max_attempts = 1000  # to prevent infinite loops
    seen = set()

    while attempts < max_attempts:
        attempts += 1
        # Randomly select a gene
        gene, exons = random.choice(list(valid_genes.items()))
        exons = natsorted(exons, key=lambda x: x[1])  # Sort by start position
        
        # Randomly choose number of exons to select
        num_exons = random.randint(min_exons, min(max_exons, len(exons)))

        # Randomly select starting exon
        if len(exons) - num_exons > 0:
            start_index = random.randint(0, len(exons) - num_exons)
        else:
            continue

        selected_exons = exons[start_index:start_index + num_exons]
        
        # Check space conditions
        if start_index > 0:
            prev_exon = exons[start_index - 1]
            if selected_exons[0][1] <= prev_exon[2] + 2000:  # Checking space before the first selected exon
                continue
        if start_index + num_exons < len(exons):
            next_exon = exons[start_index + num_exons]
            if next_exon[1] <= selected_exons[-1][2] + 2000:  # Checking space after the last selected exon
                continue
        
        # Extract chromosome, start, end, and gene info
        chr = selected_exons[0][0]
        start = selected_exons[0][1] - 2000  # Add buffer space
        end = selected_exons[-1][2] + 2000   # Add buffer space
        gene_info = f"{gene}_{start_index+1}_to_{start_index+num_exons}"

        # Check if this selection has been seen
        if (chr, start, end, gene_info) not in seen:
            seen.add((chr, start, end, gene_info))
            return (chr, start, end, gene_info)

    # Fallback if no suitable selection is found after max attempts
    print("Warning: Could not find a suitable set of exons after multiple attempts.")
    return None

def simulate_multiple_cnvs(gene_exon_info, num_cnvs, min_exons, max_exons):
    """
    Simulate num_cnvs non-overlapping CNVs per sample.
    Each event is randomly chosen as either a single exon or a multiple exon CNV.
    """
    events = []
    attempts = 0
    max_attempts = 1000
    while len(events) < num_cnvs and attempts < max_attempts:
        attempts += 1
        # 50/50 chance for single vs. multiple exon event.
        if random.random() < 0.5:
            event = get_single_exon(gene_exon_info)
        else:
            event = get_multiple_exon(gene_exon_info, min_exons, max_exons)
        if event is None:
            continue
        new_chr, new_start, new_end, new_info = event
        # Check for overlaps with previously selected events (on the same chromosome)
        overlap = False
        for ev in events:
            ev_chr, ev_start, ev_end, _ = ev
            if new_chr == ev_chr and not (new_end < ev_start or new_start > ev_end):
                overlap = True
                break
        if overlap:
            continue
        events.append(event)
    if len(events) < num_cnvs:
        print("Warning: Could not simulate the requested number of non-overlapping CNVs.")
    return events

if __name__ == "__main__":
    main()

