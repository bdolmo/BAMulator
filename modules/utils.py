import pysam
import glob
import os
import subprocess
import sys
import csv
import pandas as pd
from typing import List, Dict

def get_file_prefix(filename: str) -> str:
    return os.path.basename(filename).split('.')[0]


def get_input_bams(variants_file: str) -> list[str]:
    if not os.path.isfile(variants_file):
        raise FileNotFoundError(f" ERROR: Variants file not found: {variants_file}")
    bam_files = set()
    with open(variants_file, 'r') as f:
        header = next(f).strip().split('\t')
        if not header or header[0].upper() != 'SAMPLE':
            raise ValueError(f" ERROR: Variants file must have 'SAMPLE' as first column")
        for line in f:
            fields = line.strip().split('\t')
            if fields and os.path.isfile(fields[0]):
                bam_files.add(fields[0])
    if not bam_files:
        raise ValueError(f" ERROR: No valid BAM files found in {variants_file}")
    return sorted(bam_files)


def get_input_vcfs(input_path: str | None) -> list[str]:
    if not input_path:
        return []
    if os.path.isdir(input_path):
        patterns = [os.path.join(input_path, "*.vcf"), os.path.join(input_path, "*.vcf.gz")]
        vcf_files = []
        for pattern in patterns:
            vcf_files.extend(glob.glob(pattern))
        if not vcf_files:
            print(f" ERROR: No VCF files found in directory: {input_path}")
        return sorted(vcf_files)
    elif os.path.isfile(input_path):
        if input_path.endswith('.vcf') or input_path.endswith('.vcf.gz'):
            return [input_path]
        try:
            with open(input_path, 'r', encoding='utf-8') as f:
                files = [line.strip() for line in f if line.strip() and (line.strip().endswith('.vcf') or line.strip().endswith('.vcf.gz'))]
            valid_files = [f for f in files if os.path.isfile(f)]
            if not valid_files:
                print(f" ERROR: No valid VCF files found in {input_path}")
            return valid_files
        except UnicodeDecodeError:
            print(f" ERROR: Cannot read {input_path}: Invalid text file (not UTF-8 encoded)")
            return []
    print(f" ERROR: Invalid --input_vcf_list: {input_path} is neither a valid VCF file nor a directory nor a text file")
    return []


def match_bam_vcf(bam_files: list[str], vcf_files: list[str]) -> dict[str, tuple[str, str | None]]:
    bam_prefixes = {get_file_prefix(bam): bam for bam in bam_files}
    vcf_prefixes = {get_file_prefix(vcf): vcf for vcf in vcf_files}
    matched = {}
    for prefix in bam_prefixes:
        matched[prefix] = (bam_prefixes[prefix], vcf_prefixes.get(prefix))
    return matched


def link_bam_vcf(bams: List[str], vcfs: List[str]) -> Dict[str, Dict[str, str]]:
    sample_map = {}
    
    def extract_sample_name(filepath: str) -> str:
        return os.path.basename(filepath).split('.')[0]
    
    for bam in bams:
        sample = extract_sample_name(bam)
        sample_map[sample] = sample_map.get(sample, {'bam': '', 'vcf': ''})
        sample_map[sample]['bam'] = bam
    
    for vcf in vcfs:
        sample = extract_sample_name(vcf)
        sample_map[sample] = sample_map.get(sample, {'bam': '', 'vcf': ''})
        sample_map[sample]['vcf'] = vcf
    
    return {s: f for s, f in sample_map.items() if f['bam'] and f['vcf']}



def sort_bam(input_bam, output_bam=None, n_cpus=4):
    """
    Sort a BAM file.

    :param input_bam: Path to the input BAM file.
    :param output_bam: Path to the output sorted BAM file. If None, replaces the input BAM file name's suffix with '.sorted.bam'.
    :return: Path to the sorted BAM file.
    """
    filedir = os.path.dirname(os.path.abspath(__file__))

    samtools_binary = filedir + "/../samtools/samtools"
    if not os.path.isfile(samtools_binary):
        msg = " ERROR: Missing samtools binary"
        print(msg)
        sys.exit()

    if output_bam is None:
        output_bam = input_bam.rsplit(".", 1)[0] + ".sorted.bam"
    
    # Constructing the samtools sort command
    cmd = [samtools_binary, "sort", "-@", str(n_cpus),  "-T test", "-o", output_bam, input_bam]
    
    try:
        msg = f" INFO: Executing Samtools sorting for {input_bam}"
        print(msg)
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        msg = " ERROR: Could not sort bam:", str(e)
        print(msg)

    # os.remove(input_bam)
    # os.rename(output_bam, input_bam)

    # index_bam(input_bam)
    index_bam(output_bam)


def index_bam(bam_file):
    """
    Indexes a BAM file.

    :param bam_file: Path to the BAM file to be indexed.

    **Example**::

        bam_file = "path/to/your/sorted.bam"
        index_bam(bam_file)
        print("BAM file indexed.")
    """
    pysam.index(bam_file)



def process_all_bam_files(directory, reference, threads):
    """
    For each sample with simulated FASTQs and corresponding unsimulated BAM, 
    map reads, merge simulated + unsimulated BAMs, and index the result.
    """
    filedir = os.path.dirname(os.path.abspath(__file__))
    samtools = os.path.join(filedir, "../samtools/samtools")
    bwa = os.path.join(filedir, "../bwa/bwa")

    if not os.path.isfile(samtools) or not os.path.isfile(bwa):
        print(" ERROR: Missing bwa or samtools binaries.")
        sys.exit(1)

    # reference = os.path.join(directory, "ref.fa")
    if not os.path.isfile(reference):
        print(f" ERROR: Missing reference: {reference}")
        sys.exit(1)

    # Check BWA index files
    for ext in [".bwt", ".pac", ".ann", ".amb", ".sa"]:
        if not os.path.isfile(reference + ext):
            print(f" ERROR: Missing BWA index file: {reference + ext}")
            sys.exit(1)

    fq1_files = glob.glob(os.path.join(directory, "*_simulated_R1.fq"))

    empty_fq_files = False

    for fq in fq1_files:
        file_size = os.path.getsize(fq)
        if file_size == 0:
            empty_fq_files = True
    
    if not empty_fq_files:

        for fq1 in fq1_files:
            base_name = os.path.basename(fq1).replace("_simulated_R1.fq", "")
            fq2 = os.path.join(directory, base_name + "_simulated_R2.fq")
            unsim_bam = os.path.join(directory, base_name + "_unsimulated.bam")
            sorted_unsim_bam = unsim_bam.replace(".bam", ".sorted.bam")


            simulated_bam = os.path.join(directory, base_name + "_simulated.bam")
            merged_bam = os.path.join(directory, base_name + "_merged.bam")

            if not os.path.isfile(fq2):
                print(f" WARNING: Missing mate FASTQ: {fq2}")
                continue
            if not os.path.isfile(unsim_bam):
                print(f" WARNING: Missing unsimulated BAM: {unsim_bam}")
                continue

            # print(f" INFO: Sorting {unsim_bam}")
            # sort_bam(input_bam=unsim_bam, n_cpus=threads)
            # index_bam(sorted_unsim_bam)


            print(f" INFO: Mapping {fq1} and {fq2} with BWA MEM")
            cmd = (
                f"{bwa} mem -t {threads} {reference} {fq1} {fq2} | "
                f"{samtools} view -@ {threads} -b - | "
                f"{samtools} sort -@ {threads} -o {simulated_bam} -"
            )
            try:
                subprocess.run(cmd, shell=True, check=True)
            except subprocess.CalledProcessError as e:
                print(f" ERROR: Mapping failed for {base_name}: {e}")
                continue

            print(f" INFO: Indexing simulated BAM: {simulated_bam}")
            index_bam(simulated_bam)

            print(f" INFO: Merging {simulated_bam} + {unsim_bam} → {merged_bam}")
            merge_cmd = [
                samtools, "merge", "-@", str(threads),
                "-f", merged_bam, simulated_bam, unsim_bam
            ]
            merge_cmd_str = " ".join(merge_cmd)
            print(merge_cmd_str)

            try:
                subprocess.run(merge_cmd, check=True)
            except subprocess.CalledProcessError as e:
                print(f" ERROR: Merge failed: {e}")
                continue

            print(f" INFO: Indexing merged BAM: {merged_bam}")
            index_bam(merged_bam)

            print(f" INFO: Done processing sample: {base_name}")

    else:
        # Fetch all bam files 
        bam_files = glob.glob(os.path.join(directory, "*.bam"))

        for bam_file in bam_files:
            msg = f" INFO: Sorting BAM {bam_file}"
            print(msg)
            sort_bam(input_bam=bam_file, n_cpus=threads)


# def process_all_bam_files(directory, threads):
#     """
#     Sorts and indexes all BAM files in the specified directory.

#     :param directory: The directory to search for BAM files.
    
#     **Example**::

#         directory = "/path/to/your/directory/with/bamfiles"
#         process_all_bam_files(directory)
#         print("All BAM files processed.")

#     """

#     # first, map fastq files...




#     # Fetch all bam files 
#     bam_files = glob.glob(os.path.join(directory, "*.bam"))

#     for bam_file in bam_files:
#         msg = f" INFO: Sorting BAM {bam_file}"
#         print(msg)
#         sort_bam(input_bam=bam_file, n_cpus=threads)


def map_fastq_with_bwa(fq1, fq2, reference_fasta, output_bam, n_cpus=4):
    """
    Maps paired-end FASTQ files to a reference using BWA MEM and pipes the output to generate a sorted and indexed BAM.

    :param fq1: Path to the FASTQ R1 file.
    :param fq2: Path to the FASTQ R2 file.
    :param reference_fasta: Path to the reference FASTA file.
    :param output_bam: Path to the final sorted BAM output.
    :param n_cpus: Number of threads to use.
    """
    filedir = os.path.dirname(os.path.abspath(__file__))
    samtools = os.path.join(filedir, "../samtools/samtools")
    bwa = os.path.join(filedir, "../bwa/bwa")

    # Check binaries
    if not os.path.isfile(bwa):
        print(" ERROR: Missing bwa binary at", bwa)
        sys.exit(1)

    if not os.path.isfile(samtools):
        print(" ERROR: Missing samtools binary at", samtools)
        sys.exit(1)

    # Check FASTA index files
    for ext in [".bwt", ".pac", ".ann", ".amb", ".sa"]:
        if not os.path.isfile(reference_fasta + ext):
            print(f" ERROR: Missing BWA index file: {reference_fasta + ext}")
            sys.exit(1)

    # Construct pipeline command
    cmd = (
        f"{bwa} mem -t {n_cpus} {reference_fasta} {fq1} {fq2} | "
        f"{samtools} view -@ {n_cpus} -b - | "
        f"{samtools} sort -@ {n_cpus} -o {output_bam} -"
    )

    try:
        print(f" INFO: Mapping {fq1} and {fq2} to {reference_fasta}")
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f" ERROR: Mapping failed: {e}")
        sys.exit(1)

    # Index the final BAM
    index_bam(output_bam)
    print(f" INFO: Finished mapping and indexing: {output_bam}")
