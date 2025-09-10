import os
import re
import glob
import subprocess
import shutil
from typing import Dict, List, Optional
import sys


script_path = os.path.dirname(os.path.realpath(__file__))


class Beagle:
    def __init__(self,
        input_vcf: str,
        out_dir: str,
        genome_version: str,
        threads: int = 4,
        java_mem_gb: int = 8,
        map_dir: Optional[str] = None,
        ref_dir: Optional[str] = None,
        run_chroms: Optional[List[str]] = None,
        use_imputation: bool = True):

        self.input_vcf = input_vcf
        self.out_dir = out_dir
        self.threads = threads
        self.java_mem_gb = java_mem_gb
        self.genome_version = genome_version
        self.map_dir = map_dir
        self.ref_dir = ref_dir 
        self.use_imputation = use_imputation

        self.beagle_jar = os.path.join(script_path, "beagle", "beagle.27Feb25.75f.jar")

        if not os.path.isdir(self.out_dir):
            os.makedirs(self.out_dir, exist_ok=True)

        tools = ["bgzip", "tabix", "bcftools", "java"]

        for tool in tools:
            if not shutil.which(tool):
                raise FileNotFoundError(f" ERROR: '{tool}' not found in PATH")

        if not os.path.isfile(self.beagle_jar):
            raise FileNotFoundError(f" ERROR: Beagle jar not found: {self.beagle_jar}")

        mapping = {
            "hg19": "GRCh37", 
            "b37": "GRCh37", 
            "hg38": "GRCh38", 
            "b38": "GRCh38"
        }
        if genome_version in mapping:
            self.build_tag = mapping[genome_version]
        else:
            raise KeyError(f" ERROR: Unallowed genome_version: {genome_version}")

        self.input_vcf = self.compress_and_index_vcf(self.input_vcf)
        self.input_vcf = self.normalize_vcf_chromosomes(self.input_vcf)
        self.use_chr_prefix = self.check_chromosome_prefix()
        self.set_chromosomes(run_chroms)
        self.map_by_chr = {}
        self.ref_by_chr = {}
        self.find_plink_maps()
        if self.ref_dir is not None:
            self.find_bref3_references()
        self.check_plink_maps()
        self.check_bref3_references()
        num_samples = self.count_vcf_samples(self.input_vcf)
        if num_samples == 1 and self.ref_dir is None:
            raise RuntimeError(" ERROR: For single-sample analysis it is required a reference panel (.bref3)")

    def set_chromosomes(self, run_chroms=None):
        if run_chroms is not None:
            self.chroms = run_chroms
        else:
            base = [str(i) for i in range(1, 23)]
            self.chroms = [f"chr{c}" if self.use_chr_prefix else c for c in base]

    def execute_shell_command(self, command_args: list[str]) -> subprocess.CompletedProcess:
        if not command_args:
            raise ValueError(" ERROR: Command arguments list cannot be empty")
        if not all(isinstance(arg, str) for arg in command_args):
            raise ValueError(" ERROR: All command arguments must be strings")

        command_str = " ".join(command_args)
        # print(f" INFO: Executing command: {command_str}")

        try:
            result = subprocess.run(
                command_args,
                capture_output=True,
                text=True,
                check=False
            )
        except FileNotFoundError:
            raise FileNotFoundError(f" ERROR: Command not found: {command_args[0]}")

        if result.returncode != 0:
            error_list = result.stderr.split("\n")
            error_message = ""
            for line in error_list:
                error_message+= "\t" +  line +"\n"
            raise RuntimeError()

        return result

    def compress_and_index_vcf(self, path: str) -> str:
        if not os.path.isfile(path):
            raise FileNotFoundError(f" ERROR: Input file not found: {path}")

        file_ext = path.lower()

        if file_ext.endswith(".vcf"):
            gz_path = f"{path}.gz"
            if not os.path.isfile(gz_path):
                try:
                    with open(gz_path, "wb") as out_f:
                        self.execute_shell_command(["bgzip", "-c", path])
                except subprocess.CalledProcessError:
                    raise RuntimeError(" ERROR: Failed to compress VCF file with bgzip")
            path = gz_path
            file_ext = path.lower()

        if file_ext.endswith(".vcf.gz"):
            try:
                self.execute_shell_command(["tabix", "-f", "-p", "vcf", path])
            except subprocess.CalledProcessError:
                raise RuntimeError(" ERROR: Failed to index VCF.GZ file with tabix")
            return path

        if file_ext.endswith(".bcf"):
            try:
                self.execute_shell_command(["bcftools", "index", "-f", path])
            except subprocess.CalledProcessError:
                raise RuntimeError(" ERROR: Failed to index BCF file with bcftools")
            return path

        raise ValueError(" ERROR: Input must be .vcf, .vcf.gz, or .bcf")

    def check_chromosome_prefix(self) -> bool:
        try:
            proc = subprocess.run(
                ["bcftools", "view", "-H", "-n", "1", self.input_vcf],
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL,
                text=True,
                check=True
            )
            txt = proc.stdout.strip()
            if txt:
                chrom = txt.split("\t", 1)[0]
                return chrom.startswith("chr")
        except (subprocess.CalledProcessError, IndexError):
            pass

        try:
            proc = subprocess.run(
                ["bcftools", "view", "-h", self.input_vcf],
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL,
                text=True,
                check=True
            )
            return "ID=chr" in proc.stdout
        except subprocess.CalledProcessError:
            return False

    def normalize_vcf_chromosomes(self, input_vcf: str) -> str:

        input_vcf_name = os.path.basename(self.input_vcf).replace(".vcf.gz", "")
        
        output_vcf = os.path.join(self.out_dir, f"{input_vcf_name}.normalized.vcf.gz")
        map_file = os.path.join(self.out_dir, "chrom_map.txt")
        with open(map_file, 'w') as f:
            for i in range(1, 23):
                f.write(f"chr{i}\t{i}\n")
            f.write("chrX\tX\n")
            f.write("chrY\tY\n")

        try:
            self.execute_shell_command(
                ["bcftools", "annotate", "--rename-chrs", map_file, input_vcf, "-Oz", "-o", output_vcf]
            )
            print(f" INFO: Normalizing chromosome names for {input_vcf}")
        except subprocess.CalledProcessError as e:
            error_msg = e.stderr.strip() or "Failed to normalize VCF chromosome names"
            raise RuntimeError(f" ERROR: {error_msg}")

        try:
            self.execute_shell_command(["tabix", "-p", "vcf", output_vcf])
        except subprocess.CalledProcessError as e:
            error_msg = e.stderr.strip() or "Failed to index normalized VCF"
            raise RuntimeError(f" ERROR: {error_msg}")

        return output_vcf

    def canonical_chromosome_label(self, label: str) -> str:
        pattern = r'^(?:chr)?([0-9]{1,2}|X|Y)$'
        match = re.match(pattern, label, re.IGNORECASE)
        if match:
            return match.group(1).upper()
        return label.upper()

    def count_vcf_samples(self, vcf_gz: str) -> int:
        if not os.path.exists(vcf_gz):
            raise FileNotFoundError(f" ERROR: VCF file not found: {vcf_gz}")

        try:
            result = self.execute_shell_command(
                ["bcftools", "query", "-l", vcf_gz]
            )
            sample_lines = result.stdout.strip().splitlines()
            return len(sample_lines)
        except subprocess.CalledProcessError as e:
            error_message = e.stderr.strip()
            raise RuntimeError(f" ERROR: bcftools query failed for {vcf_gz}: {error_message}")
        except FileNotFoundError:
            raise RuntimeError(f" ERROR: bcftools not found in system PATH")

    def find_plink_maps(self):
        if not self.map_dir:
            raise FileNotFoundError(f" ERROR: Map directory not provided")
        if not os.path.isdir(self.map_dir):
            raise FileNotFoundError(f" ERROR: Map directory not found: {self.map_dir}")

        pattern = os.path.join(self.map_dir, "**", "*.map")
        paths = glob.glob(pattern, recursive=True)
        if not paths:
            raise FileNotFoundError(f" ERROR: No .map files found under {self.map_dir}")

        self.map_by_chr = {}
        chromosome_pattern = re.compile(r'(?:^|[^A-Za-z])(?:chr)?([0-9]{1,2}|X|Y)(?:[^0-9A-Za-z]|$)', re.IGNORECASE)
        preferred_pattern = re.compile(rf'plink\.chr([0-9]{{1,2}}|X|Y)\.{self.build_tag}\.map', re.IGNORECASE)

        for path in paths:
            filename = os.path.basename(path)
            match = chromosome_pattern.search(filename)
            if match:
                chromosome = match.group(1).upper()
                if preferred_pattern.fullmatch(filename):
                    self.map_by_chr[chromosome] = path
                elif chromosome not in self.map_by_chr:
                    self.map_by_chr[chromosome] = path

    def find_bref3_references(self):
        if not os.path.isdir(self.ref_dir):
            raise FileNotFoundError(f" ERROR: Reference directory not found: {self.ref_dir}")

        patterns = [
            os.path.join(self.ref_dir, "**", "*.bref3"),
            os.path.join(self.ref_dir, "**", "*.vcf"),
            os.path.join(self.ref_dir, "**", "*.vcf.gz"),
        ]

        all_paths = []
        for pattern in patterns:
            all_paths.extend(glob.glob(pattern, recursive=True))

        if not all_paths:
            raise FileNotFoundError(f" ERROR: No reference files found under {self.ref_dir}")

        non_bref3 = []
        self.ref_by_chr = {}
        chromosome_pattern = re.compile(r'(?:^|[^A-Za-z])(?:chr)?([0-9]{1,2}|X|Y)(?:[^0-9A-Za-z]|$)', re.IGNORECASE)

        for path in all_paths:
            if path.endswith(".bref3"):
                filename = os.path.basename(path)
                match = chromosome_pattern.search(filename)
                if match:
                    chromosome = match.group(1).upper()
                    self.ref_by_chr[chromosome] = path
            else:
                non_bref3.append(path)

        if non_bref3:
            raise ValueError(
                f" ERROR: Only .bref3 references are supported. Found non-bref3 files:\n  " +
                "\n  ".join(non_bref3)
            )

        if not self.ref_by_chr:
            raise FileNotFoundError(f" ERROR: No .bref3 reference files mapped to chromosomes at {self.ref_dir}")

    def check_plink_maps(self):
        target_chromosomes = []
        for chrom in self.chroms:
            target_chromosomes.append(self.canonical_chromosome_label(chrom))

        missing_chromosomes = []
        for chrom in target_chromosomes:
            if chrom not in self.map_by_chr:
                missing_chromosomes.append(chrom)

        if missing_chromosomes:
            readable_list = ", ".join(missing_chromosomes)
            raise FileNotFoundError(f" ERROR: Missing PLINK maps for: {readable_list} at {self.map_dir}")

    def check_bref3_references(self):
        is_single_sample = self.count_vcf_samples(self.input_vcf) < 2
        references_needed = self.ref_dir is not None or is_single_sample
        if not references_needed:
            return

        if self.ref_dir is None:
            raise FileNotFoundError(" ERROR: Reference directory not provided but required for .bref3 files")

        target_chromosomes = []
        for chrom in self.chroms:
            target_chromosomes.append(self.canonical_chromosome_label(chrom))

        missing_chromosomes = []
        for chrom in target_chromosomes:
            if chrom not in self.ref_by_chr:
                missing_chromosomes.append(chrom)

        if missing_chromosomes:
            readable_list = ", ".join(missing_chromosomes)
            raise FileNotFoundError(f" ERROR: Missing .bref3 references for chromosomes: {readable_list} (searched in {self.ref_dir})")

    def extract_chromosome_vcf(self, chrom_label: str) -> Optional[str]:
        out_path = os.path.join(self.out_dir, f"sample.{chrom_label}.vcf.gz")

        try:
            probe = subprocess.run(
                f"bcftools view -H -r {chrom_label} {self.input_vcf} | head -n 1",
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
            if not probe.stdout.strip():
                print(f" INFO: No variants found for {chrom_label}")
                return None
        except subprocess.CalledProcessError as e:
            print(f" INFO: Probe failed for {chrom_label}: {e.stderr.strip()}")
            return None

        if not os.path.isfile(out_path):
            try:
                self.execute_shell_command(
                    ["bcftools", "view", "-r", chrom_label, "-Oz", "-o", out_path, self.input_vcf]
                )
                self.execute_shell_command(["tabix", "-p", "vcf", out_path])
            except subprocess.CalledProcessError as e:
                error_msg = e.stderr.strip() or "Failed to extract chromosome region"
                raise RuntimeError(f" ERROR: {error_msg}")

        return out_path

    def concatenate_phased_vcfs(self, vcf_files: List[str]) -> str:
        if not vcf_files:
            raise ValueError(" ERROR: No phased VCF files provided for concatenation")

        for vcf in vcf_files:
            if not os.path.isfile(vcf):
                raise FileNotFoundError(f" ERROR: Phased VCF file not found: {vcf}")
            if not os.path.isfile(f"{vcf}.tbi"):
                raise FileNotFoundError(f" ERROR: Index file not found for {vcf}")

        input_vcf_name = os.path.basename(self.input_vcf).replace(".vcf.gz", "")
       
        output_vcf = os.path.join(self.out_dir, f"{input_vcf_name}.phased.vcf.gz")
        try:
            command_str = " ".join(["bcftools", "concat", "-Oz", "-o", output_vcf] + vcf_files)
            max_line_length = 120
            if len(command_str) > max_line_length:
                words = ["bcftools", "concat", "-Oz", "-o", output_vcf] + vcf_files
                lines = []
                current_line = []
                current_length = 0
                for word in words:
                    word_length = len(word) + 1
                    if current_length + word_length > max_line_length and current_line:
                        lines.append(" ".join(current_line) + " \\")
                        current_line = ["\t" + word]
                        current_length = word_length
                    else:
                        current_line.append("\t" + word)
                        current_length += word_length
                if current_line:
                    lines.append(" ".join(current_line))
                command_display = "\n".join(lines)
            else:
                command_display = command_str

            print(f" INFO: Concatenating phased VCF files => {output_vcf}")
            # print(f"{command_display}")
            self.execute_shell_command(["bcftools", "concat", "-Oz", "-o", output_vcf] + vcf_files)
            print(f" INFO: Indexing concatenated VCF {output_vcf}")
            self.execute_shell_command(["tabix", "-p", "vcf", output_vcf])
        except subprocess.CalledProcessError as e:
            error_msg = e.stderr.strip() or "Failed to concatenate or index VCF files"
            raise RuntimeError(f" ERROR: {error_msg}")

        return output_vcf

    def process_chromosomes(self) -> List[str]:
        if not os.path.isfile(self.input_vcf):
            raise FileNotFoundError(f" ERROR: Input VCF not found: {self.input_vcf}")

        outputs = []
        window_sizes = [100, 150, 200]

        for chrom_label in self.chroms:
            canon = self.canonical_chromosome_label(chrom_label)
            print(f" INFO: Processing chromosome {chrom_label}")

            # Extract input slice for this chromosome
            try:
                cohort_chr = self.extract_chromosome_vcf(canon)
                if cohort_chr is None:
                    print(f" INFO: {canon} not found or empty in input; skipping")
                    continue
            except Exception as e:
                print(f" ERROR: Failed to extract VCF for {canon}: {str(e)}")
                continue

            # Genetic map
            try:
                map_path = self.map_by_chr[canon]
                print(f" INFO: genetic_map: {map_path}")
            except KeyError:
                print(f" ERROR: No map found for {canon}")
                continue

            out_prefix = os.path.join(self.out_dir, f"beagle.{canon}")
            out_vcf = f"{out_prefix}.vcf.gz"
            success = False

            # Try multiple window sizes
            for window in window_sizes:
                args = [
                    "java", f"-Xmx{self.java_mem_gb}g", "-jar", self.beagle_jar,
                    f"gt={cohort_chr}",
                    f"map={map_path}",
                    f"out={out_prefix}",
                    f"nthreads={self.threads}",
                    f"window={window}"
                ]

                if self.ref_dir is not None:
                    try:
                        ref_path = self.ref_by_chr[canon]
                        args.append(f"ref={ref_path}")
                        args.append("impute=true" if self.use_imputation else "impute=false")
                        print(f" INFO: reference: {ref_path}")
                    except KeyError:
                        print(f" ERROR: No reference file found for {canon}")
                        break

                try:
                    # (optional) pretty print
                    command_str = " ".join(args)
                    max_line_length = 120
                    if len(command_str) > max_line_length:
                        words = args
                        lines, current_line, current_length = [], [], 0
                        for w in words:
                            wl = len(w) + 1
                            if current_length + wl > max_line_length and current_line:
                                lines.append(" ".join(current_line) + " \\")
                                current_line = ["\t" + w]
                                current_length = wl
                            else:
                                current_line.append("\t" + w)
                                current_length += wl
                        if current_line:
                            lines.append(" ".join(current_line))
                        print("\n".join(lines))
                    else:
                        print(command_str)

                    print(f" INFO: Running Beagle for chromosome {canon} (window={window})")
                    self.execute_shell_command(args)
                    success = True
                    print(" INFO: Success")
                    break
                except Exception as e:
                    print(f" WARNING: Failed for chromosome {canon} (window={window}): {str(e)}")
                    continue

            if not success:
                print(f" ERROR: All window sizes failed for {canon}; skipping")
                continue

            # Ensure Beagle output is indexed
            if os.path.isfile(out_vcf) and not os.path.isfile(f"{out_vcf}.tbi"):
                try:
                    self.execute_shell_command(["tabix", "-p", "vcf", out_vcf])
                except Exception as e:
                    print(f" ERROR: Failed to index VCF for {canon}: {str(e)}")
                    continue

            # ---- Keep ONLY positions present in cohort_chr (position-only, no rename map) ----
            try:
                # Build a site list (CHROM\tPOS). Include both with- and without-'chr'
                sites_txt = os.path.join(self.out_dir, f"typed.sites.{canon}.txt")
                res = self.execute_shell_command(["bcftools", "query", "-f", "%CHROM\t%POS\n", cohort_chr])

                seen = set()
                with open(sites_txt, "w") as f:
                    for line in res.stdout.splitlines():
                        if not line:
                            continue
                        chrom, pos = line.split("\t")
                        # no-chr
                        chrom_nochr = chrom[3:] if chrom.startswith("chr") else chrom
                        key1 = (chrom_nochr, pos)
                        if key1 not in seen:
                            f.write(f"{chrom_nochr}\t{pos}\n")
                            seen.add(key1)
                        # with-chr
                        chrom_chr = chrom if chrom.startswith("chr") else f"chr{chrom}"
                        key2 = (chrom_chr, pos)
                        if key2 not in seen:
                            f.write(f"{chrom_chr}\t{pos}\n")
                            seen.add(key2)

                # Subset Beagle VCF to those positions
                pos_only_vcf = os.path.join(self.out_dir, f"beagle.{canon}.only_input_sites.vcf.gz")
                self.execute_shell_command([
                    "bcftools", "view", "-T", sites_txt, out_vcf, "-Oz", "-o", pos_only_vcf
                ])
                self.execute_shell_command(["tabix", "-f", "-p", "vcf", pos_only_vcf])

                print(f" INFO: Position-only filter applied for {canon} using {sites_txt} => {pos_only_vcf}")
                out_vcf = pos_only_vcf
            except Exception as e:
                print(f" ERROR: Position-only subsetting failed for {canon}: {e}")
                continue
            # ----------------------------------------------------------------------------------

            outputs.append(out_vcf)

        return outputs


if __name__ == "__main__":

    input_vcf="/raw-data/projectes/SPIKE_IN_BAM/BAM_PANELS/RB16187.gatk.vcf" 
    output_vcf="/raw-data/projectes/SPIKE_IN_BAM/BAM_PANELS/test.phased.vcf"
    output_dir = "/raw-data/projectes/SPIKE_IN_BAM/BAM_PANELS/test"
    genome_version="b37"
    threads= 4
    map_dir="/home/udmmp/Desktop/SpikeInBAM2/modules/phasing/beagle/genetic_maps"
    ref_dir = "/raw-data/ANNOTATIONS/ANN_DIR/1000Genomes/hg19/beagle"

    beagle = Beagle(
        input_vcf=input_vcf,
        out_dir=output_dir,
        genome_version=genome_version,
        threads=threads,
        map_dir=map_dir,
        ref_dir=ref_dir,         # can be None if no reference panel
        run_chroms=None,         # or e.g. ["1", "2", "X"] if you want only certain chroms
        use_imputation=False      # set False to phase without imputing missing variants
    )

    # Run phasing
    phased_vcfs = beagle.process_chromosomes()
    beagle.concatenate_phased_vcfs(phased_vcfs)
