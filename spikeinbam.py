import click
import subprocess
import os
import time
import glob
from modules.utils import process_all_bam_files, get_file_prefix, get_input_bams, get_input_vcfs, match_bam_vcf, link_bam_vcf
from modules.phasing.phasing import Beagle
import sys

def phase_vcfs(vcf_files: list[str], output_dir: str, genome_version: str, map_dir: str, ref_dir: str | None, threads: int) -> list[str]:
    phased_vcfs = []
    for vcf in vcf_files:
        prefix = get_file_prefix(vcf)
        print(f" INFO: Phasing VCF for sample {prefix}")
        try:
            beagle = Beagle(
                input_vcf=vcf,
                out_dir=output_dir,
                genome_version=genome_version,
                threads=threads,
                map_dir=map_dir,
                ref_dir=ref_dir
            )
            outputs = beagle.process_chromosomes()
            if not outputs:
                print(f" ERROR: No phased VCFs produced for {vcf}")
                continue
            phased_vcf = beagle.concatenate_phased_vcfs(outputs)
            print(f" INFO: Phased VCF generated: {phased_vcf}")
            phased_vcfs.append(phased_vcf)
        except Exception as e:
            print(f" ERROR: Failed to phase VCF {vcf}: {str(e)}")
            continue
    return phased_vcfs


@click.command()
@click.option('--variants', '-v', required=True, help='Tab-separated file listing BAM files and variants')
@click.option('--threads', '-t', type=int, default=1, help='Total number of threads')
@click.option('--reference', '-r', required=True, help='Genome reference in FASTA format')
@click.option('--output', '-o', required=True, help='Output directory for simulated BAMs and phased VCFs')
@click.option('--suffix', '-s', default=".simulated", help='Suffix for simulated BAM files')
@click.option('--input_vcf_list', default=None, help='File listing input VCFs or directory containing VCFs')
@click.option('--genome_version', default='hg19', help='Genome version for phasing (e.g., hg38, hg19)')
@click.option('--map_dir', default=None, help='Directory containing PLINK .map files for phasing')
@click.option('--ref_dir', default=None, help='Directory containing .bref3 reference files for phasing')
def main(variants, threads, reference, output, suffix, input_vcf_list, genome_version, map_dir, ref_dir):
    # parser = argparse.ArgumentParser(description='Spike: simulate variants in existing BAMs')
    # parser.add_argument('--variants', '-v', required=True, help='Directory containing input BAM files')
    # parser.add_argument('--threads', '-t', type=int, default=1, help='Total number of threads')
    # parser.add_argument('--reference', '-r', required=True, help='Genome reference in FASTA format')
    # parser.add_argument('--output', '-o', required=True, help='Output directory for simulated BAMs and phased VCFs')
    # parser.add_argument('--suffix', '-s', default=".simulated", help='Suffix for simulated BAM files')
    # parser.add_argument('--input_vcf_list', default=None, help='File listing input VCFs or directory containing VCFs')
    # parser.add_argument('--genome_version', default='hg19', help='Genome version for phasing (e.g., hg38, hg19)')
    # parser.add_argument('--map_dir', default=None, help='Directory containing PLINK .map files for phasing')
    # parser.add_argument('--ref_dir', default=None, help='Directory containing .bref3 reference files for phasing')
    # args = parser.parse_args()

    start_time = time.time()
    script_dir = os.path.dirname(os.path.abspath(__file__))
    spikeinbam_exe = os.path.join(script_dir, "src", "spikeinbam")

    if not os.path.isdir(output):
        os.makedirs(output, exist_ok=True)

    # Get input BAMs and VCFs
    bam_files = get_input_bams(variants)
    if not bam_files:
        raise ValueError(" ERROR: No BAM files found in --variants file")
    print(f" INFO: Found {len(bam_files)} BAM files")

    vcf_files = get_input_vcfs(input_vcf_list) if input_vcf_list else []
    if input_vcf_list and not vcf_files:
        print(f" ERROR: No VCF files found in --input_vcf_list: {input_vcf_list}")
    else:
        print(f" INFO: Found {len(vcf_files)} VCF files")

    # Match BAMs and VCFs by prefix
    matched_files = match_bam_vcf(bam_files, vcf_files)
    print(f" INFO: Matched {len(matched_files)} samples")
 

    # Phase VCFs if provided
    phased_vcfs = []
    if vcf_files:
        if not map_dir:
            raise ValueError(" ERROR: --map_dir required for phasing VCFs")
        phased_vcfs = phase_vcfs(
            vcf_files=vcf_files,
            output_dir=output,
            genome_version=genome_version,
            map_dir=map_dir,
            ref_dir=ref_dir,
            threads=threads
        )
        print(f" INFO: Generated {len(phased_vcfs)} phased VCFs")

    bam_vcf_dict = link_bam_vcf(bam_files, phased_vcfs)
    
    bam_vcf_txt = os.path.join(output, "bam_vcf.txt")
    o = open(bam_vcf_txt, "w") 
    for sample in bam_vcf_dict:
        o.write(bam_vcf_dict[sample]['bam']+"\t"+bam_vcf_dict[sample]['vcf'])
    o.close()

    # Create BAM-VCF-haplotype mapping file
    # haplotype_map_file = create_bam_vcf_haplotype_map(variants, matched_files, phased_vcfs, output)
    # if not haplotype_map_file:
    #     raise RuntimeError(" ERROR: Failed to create BAM-VCF-haplotype mapping file")
    # sys.exit()

    # Run SpikeInBAM
    command = [
        spikeinbam_exe,
        variants,
        bam_vcf_txt,
        reference,
        output,
        suffix,
    ]
    try:
        print(f" INFO: Simulating variants:")
        print(f' INFO: {" ".join(command)}')
        subprocess.run(command, check=True)
        print(f" INFO: SpikeInBAM completed successfully")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f" ERROR: Could not execute SpikeInBAM: {str(e)}")

    # Process BAM files
    print(f" INFO: Processing simulated BAM files")
    process_all_bam_files(output, reference, threads)
    print(f" INFO: BAM file processing completed")

    # end_time = time.time()
    # total_time = end_time - start_time
    # print(f" INFO: Total execution time: {total_time:.2f} seconds")


if __name__ == "__main__":
    main()