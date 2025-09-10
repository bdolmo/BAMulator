# BAMulator — Fast variant simulation on BAMs

Simulate variants directly on existing BAM files

## Install

```bash
pip3 install -r requirements.txt
```

## Usage

```text
Usage: bamulator [OPTIONS]

Options:
  -v, --variants TEXT     Tab-separated file listing BAM files and variants
                          [required]
  -t, --threads INTEGER   Total number of threads
  -r, --reference TEXT    Genome reference in FASTA format  [required]
  -o, --output TEXT       Output directory for simulated BAMs and phased VCFs
                          [required]
  -s, --suffix TEXT       Suffix for simulated BAM files
  --input_vcf_list TEXT   File listing input VCFs or directory containing VCFs
  --genome_version TEXT   Genome version for phasing (e.g., hg38, hg19)
  --map_dir TEXT          Directory containing PLINK .map files for phasing
  --ref_dir TEXT          Directory containing .bref3 reference files for
                          phasing
  --help                  Show this message and exit.
```

### Quick start

```bash
bamulator \
  --variants variants.txt \
  --reference /path/to/ref.fasta \
  --threads 4 \
  --suffix .simulated \
  --output /path/to/outdir
```

Output BAMs are named `<sample><suffix>.bam` (default: `.simulated`).

## Variant config (`variants.txt`)

One variant per line, tab-separated.

### SNVs / small indels
```
/home/user/input/sample1.bam    chr7     55242467    GAATTAAGAGAAGCAACA   GTTGCT
/home/user/input/sample2.bam    chr18    29104352    A                    T
```
Columns: BAM, chrom, pos (1-based), REF, ALT

### CNVs
```
/home/user/input/sample3.bam    chr1    156082722   156089883   DEL   SOS1_7_to_8   1   multiple
/home/user/input/sample4.bam    chr7    150652435   150661275   DUP   TAZ_1_to_4    1   multiple
```
Columns: BAM, chrom, start, end, type (`DEL|DUP`), label, copy-number/state, mode

## Phasing (optional)

Provide `--input_vcf_list` plus:
- `--genome_version` (e.g., `hg38`, `hg19`)
- `--map_dir` with PLINK `.map`
- `--ref_dir` with `.bref3`

Phased VCFs are written to `--output`.

## Features

- SNVs, indels (including complex), CNVs
- VAFs adjusted for SNVs overlapping simulated CNVs
- Multithreaded execution

## Limitations

- No SVs requiring sequenced breakpoints
- No somatic variant simulation
- Copy-number changes >3 not supported; no homo/hemizygous deletions

## Helper script

Generate CNV configs:

```bash
python3 scripts/generate_cnv_config.py \
  --indir /path/to/bam_dir \
  -c <del|dup> \
  -b /path/to/gene_panel.bed \
  --mode <single|multiple> > variants.txt
```

