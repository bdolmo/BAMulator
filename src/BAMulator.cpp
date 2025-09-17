#include "BamReader.h"
#include "BamRecord.h"
#include "BamWriter.h"
#include "RefFasta.h"
#include <iostream>
#include <fstream>
#include <filesystem>
#include <htslib/sam.h>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include "sw.cpp"
#include <random>
#include "SnvCaller.cpp"
#include "Simulations.cpp"
#include <chrono>
#include <iomanip>
#include <libgen.h>
#include <unordered_set>
#include <set>
#include "VcfRecord.h"

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

std::string getBasename(const std::string& path) {
    size_t lastSlash = path.find_last_of("/\\");
    if (lastSlash != std::string::npos)
        return path.substr(lastSlash + 1);
    return path;
}

std::string replaceSuffix(const std::string& filename, const std::string& newSuffix) {
    size_t dotPos = filename.rfind(".bam");
    if (dotPos != std::string::npos) {
        return filename.substr(0, dotPos) + newSuffix;
    }
    return filename;
}


void replace_all(std::string& str, const std::string& from, const std::string& to) {
    if (from.empty()) return;
    size_t start_pos = 0;
    while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // Move past the replacement
    }
}


std::vector<std::string> split(const std::string &line, char delimiter) {
    std::vector<std::string> fields;
    size_t start = 0, end;

    while ((end = line.find(delimiter, start)) != std::string::npos) {
        fields.emplace_back(line.substr(start, end - start));
        start = end + 1;
    }
    fields.emplace_back(line.substr(start)); // last field
    return fields;
}
// Remove leading "chr" if present
inline std::string stripChrPrefix(const std::string& input) {
    if (input.rfind("chr", 0) == 0) { // starts with "chr"
        return input.substr(3);
    }
    return input;
}

std::map<std::string, std::vector<Variant>> parseBed(const std::string& bedFile) {

    std::map<std::string, std::vector<Variant>> variantsMap;
    std::ifstream file(bedFile);
    std::string line;

    // Read and parse header
    std::getline(file, line);
    std::vector<std::string> headers;
    std::stringstream headerStream(line);
    std::string col;
    while (std::getline(headerStream, col, '\t')) {
        headers.push_back(col);
    }


    bool isSVFormat   = std::find(headers.begin(), headers.end(), "SVTYPE") != headers.end();
    bool isSNVFormat  = std::find(headers.begin(), headers.end(), "REF")    != headers.end() &&
                        std::find(headers.begin(), headers.end(), "SAMPLE") != headers.end();

    if (!isSVFormat && !isSNVFormat) {
        std::cerr << "ERROR: Unknown or unsupported format in file: " << bedFile << std::endl;
        return variantsMap;
    }

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::stringstream ss(line);
        std::vector<std::string> fields;
        std::string value;
        while (std::getline(ss, value, '\t')) {
            fields.push_back(value);
        }

        Variant variant;

        try {
            if (isSVFormat) {
                if (fields.size() < 17) {
                    std::cerr << "WARNING: Skipping malformed SV line: " << line << std::endl;
                    continue;
                }

                variant.bamFile        = fields[0];
                variant.chr            = fields[1];
                variant.start          = std::stoll(fields[2]) - 1;
                variant.end            = std::stoll(fields[3]);
                variant.chr2           = fields[4];
                variant.start2         = (fields[5] == "." ? -1 : std::stoll(fields[5]));
                variant.end2           = (fields[6] == "." ? -1 : std::stoll(fields[6]));
                variant.varType        = fields[7];
                variant.svId           = fields[8];
                variant.svMechanism    = fields[10];
                variant.homology       = fields[11];
                variant.clone          = fields[13];
                variant.vaf            = std::stof(fields[14]);
                variant.additionalInfo = fields[15];
                variant.position       = variant.start;
                variant.haplotype      = fields[16];


                int haplotype;
                //validate haplotype must be 0,1 or 2
                try {
                    haplotype = stoi(variant.haplotype);
                }
                catch( const std::exception& e) {
                    std::cerr << " ERROR: Variant haplotype: " << variant.haplotype << " is not an integer" << std::endl;
                } 
                if (haplotype > 2) {
                    std::cerr << " ERROR: Variant haplotype: " << variant.haplotype << " must be equal or smaller than 2" << std::endl;
                }

                std::cout << " INFO: Parsed SV " << variant.varType << " (" << variant.svId << ") from " << variant.bamFile << std::endl;

            } else if (isSNVFormat) {
                if (fields.size() < 10) {
                    std::cerr << "WARNING: Skipping malformed SNV/indel line: " << line << std::endl;
                    continue;
                }

                variant.bamFile        = fields[0];
                variant.chr            = fields[1];
                variant.start          = std::stoll(fields[2]) - 1;
                variant.position       = variant.start;
                variant.refSeq         = fields[3];
                variant.altSeq         = fields[4];
                variant.varType        = fields[5];
                variant.clone          = fields[6];
                variant.vaf            = std::stof(fields[7]);
                variant.genotype       = fields[8];
                variant.additionalInfo = fields[9];
                variant.end            = variant.start + variant.refSeq.length() - 1;

                std::cout << "INFO: Parsed SNV/INDEL " << variant.varType
                          << " at " << variant.chr << ":" << variant.start
                          << " (VAF: " << variant.vaf << ") from " << variant.bamFile << std::endl;
            }

            variantsMap[variant.bamFile].push_back(variant);
        } catch (const std::exception& e) {
            std::cerr << "ERROR: Failed to parse line:\n" << line << "\n" << e.what() << std::endl;
        }
    }

    return variantsMap;
}


struct Sample {
    std::string bamFile;
    std::string vcfFile;
};


std::vector<Sample> linkBamWithPhasedVcf(const std::string& bamVcfFile) {
    std::ifstream ifile(bamVcfFile);
    if (!ifile.is_open()) {
        std::cerr << " ERROR: cannot open input file: " << bamVcfFile << '\n';
        std::runtime_error(" ERROR FOUND!");
    }

    std::vector<Sample> samples;
    std::string line;
    size_t nLine = 0;

    while (std::getline(ifile, line)) {
        ++nLine;

        // Remove trailing CR if present (Windows files)
        if (!line.empty() && line.back() == '\r') line.pop_back();

        // Skip empty lines and comment lines
        if (line.empty() || line[0] == '#') continue;

        auto fields = split(line, '\t');
        if (fields.size() < 2) {
            std::cerr << " ERROR: line " << nLine
                      << " must contain at least two tab-separated fields (BAM\tVCF)."
                      << " Got: \"" << line << "\"\n";
            std::runtime_error(" ERROR FOUND!");
        }
        Sample s;
        s.bamFile = fields[0];
        s.vcfFile = fields[1];

        // Optional: warn if files don’t exist
        // (Comment out if you don’t want filesystem checks.)
        namespace fs = std::filesystem;
        if (!s.bamFile.empty() && !fs::exists(s.bamFile)) {
            std::cerr << " WARNING: line " << nLine << " BAM not found: " << s.bamFile << '\n';
        }
        if (!s.vcfFile.empty() && !fs::exists(s.vcfFile)) {
            std::cerr << " WARNING: line " << nLine << " VCF not found: " << s.vcfFile << '\n';
        }

        samples.emplace_back(std::move(s));
    }

    return samples;
}



std::mt19937 rng{std::random_device{}()};

int main(int argc, char *argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <BAM file> <BED file> <REF FASTA> <output BAM file>" << std::endl;
        return 1;
    }

    std::string variantsFile = argv[1];
    std::string bamVcfFile = argv[2];
    std::string reference = argv[3];
    std::string outputDir = argv[4];
    std::string suffix = argv[5];

    std::vector<Sample> sampleFiles = linkBamWithPhasedVcf(bamVcfFile);

    // for (auto& sample : sampleFiles) {
    //     std::cout << sample.bamFile << " " << sample.vcfFile << std::endl;
    //     std::string region = "1:1167790-23328358";
    //     std::vector<VcfRecord> listOfSnv = fetchSnvs(sample.vcfFile, region);
    //     std::cout << listOfSnv.size() << std::endl;
    //     for (auto& snv : listOfSnv) {
    //         std::cout << snv.getGenotype() << std::endl;
    //     }
    // }
    // return 0;

    int maxReads = 1000000;
    auto variantsMap = parseBed(variantsFile);
    RefFasta ref(reference);
        std::cout << "DONE" << std::endl;

    for (const auto& [bamFile, variantList] : variantsMap) {
        std::string vcfFile;
        for (auto& sample : sampleFiles) {
            if (sample.bamFile == bamFile) {
                vcfFile = sample.vcfFile;
            }
        }

        std::unordered_set<std::string> skipReads;
        std::unordered_map<std::string, FastqRead> fqBuffer;

        // here we sill store those reads to be skipped (for large deletions)

        // here we will store all the modified reads // simulated
        std::unordered_map<std::string, std::vector<FastqRead>> modifiedReads;

        std::cout << " INFO: Processing BAM file: " << bamFile << "\n";

        std::ifstream ifile(bamFile);
        if (!ifile) {
            std::cerr << " ERROR: BAM file does not exist: " << bamFile << std::endl;
            continue;
        }

        std::cout << " INFO: Estimating insert size stats from: " << bamFile << std::endl;
        BamProfile bamProfile = calculateBamProfile(bamFile, maxReads);

        std::cout << " INFO: Mean insert size = " << bamProfile.meanInsertSize
                << ", StdDev = " << bamProfile.stdInsertSize << "\n";

        // Make a copy and sort variants by start
        std::vector<Variant> sortedVariants = variantList;
        std::sort(sortedVariants.begin(), sortedVariants.end(),
            [](const Variant& a, const Variant& b) {
                return a.start < b.start;
            });
        std::string bamName = getBasename(bamFile);
        std::string sampleName = bamName;
        replace_all(sampleName, ".bam", "");

        std::string bamNonSimulated = outputDir + "/" +  sampleName + "_unsimulated.bam";

        std::string fq1_path = outputDir + "/" + sampleName + "_simulated_R1.fq";
        std::string fq2_path = outputDir + "/" + sampleName + "_simulated_R2.fq";
        std::ofstream fastq1(fq1_path);
        std::ofstream fastq2(fq2_path);

        // Read raw bam file
        BamReader reader(bamFile);
        BamRecord record;

        // write nons imualted reads
        // BamWriter writer(bamOut, reader.getHeader());
        BamWriter writer(bamNonSimulated, reader.getHeader());

        for (const auto& variant : variantList) { 

            std::cout << " INFO: Simulating variant: " << variant.varType << " at " << variant.chr << ":" << variant.start << "-" << variant.end << "\n";
            // std::vector<SNV> snvs;

            std::vector<std::string> regions;
            std::string region1;
            std::vector<VcfRecord> snvs;

            if (variant.varType == "SNV" || variant.varType == "INDEL") {
                region1 = variant.chr + ":" + std::to_string(std::max(static_cast<int64_t>(0), variant.start - 50)) + "-" + std::to_string(variant.end + 50);
            } 
            else {
                region1 = variant.chr + ":" + std::to_string(std::max(static_cast<int64_t>(0), variant.start-2000)) + "-" + std::to_string(variant.end+2000);

                std::string regiontest = "1:80000-90000";
                // snvs = identifySNVs(bamFile, ref, region1, variant.chr);

                std::string regionVcf = stripChrPrefix(variant.chr) + ":" + std::to_string(std::max(static_cast<int64_t>(0), variant.start-2000)) + "-" + std::to_string(variant.end+2000);
                snvs = fetchSnvs(vcfFile, regionVcf);
                std::cout << snvs.size() << std::endl;
                
            }
            for (auto& variant : snvs) {
                std::cout << variant.getChrom() << ":" << variant.getPosition() << " " << variant.getRef() << " " << variant.getAlt() << " " << std::endl;
            }


            regions.push_back(region1);
            if (variant.varType == "BALANCED_TRANSLOCATION" || variant.varType == "WHOLE_ARM_TRANSLOCATION" ) {
                std::string region2 = variant.chr2 + ":" + std::to_string(std::max(static_cast<int64_t>(0), variant.start2-2000)) + "-" + std::to_string(variant.end2+2000);
                regions.push_back(region2);
            }
            int adjacency = 0;
            for (auto& region : regions) {
                adjacency++;

                std::cout << "Region: " << region << std::endl;
                reader.SetRegion(region);

                std::unordered_map<std::string, ReadPair> ReadPairs;
                while (reader.GetNextRecord(record)) {
                    if (!record.IsPaired()) {
                        continue;
                    }
                    if (record.IsUnmapped()) {
                        continue;
                    }
                    if (record.IsSecondary()) {
                        continue;
                    } 
                    if (record.IsSupplementary()) {
                        continue;
                    }

                    // if (variant.varType == "SNV") {
                    //     if (record.Position() <= variant.start &&
                    //         (record.Position() + record.Seq().length()) >= variant.start) {
                    //         simulateSNV(record, variant, ref, fastq1, fastq2, fqBuffer);
                    //     }
                    // }
                    // else if (variant.varType == "INDEL") {
                    //     if (record.Position() <= variant.start && (record.Position() + record.Seq().length()) >= variant.start) {
                    //         bool is_edit = simulateIndel(record, variant, ref, fqBuffer);
                    //     }
                    // }
                    if (variant.varType == "DEL" || variant.varType == "DUP"|| variant.varType == "INV" 
                        || variant.varType == "BALANCED_TRANSLOCATION" || variant.varType == "WHOLE_ARM_TRANSLOCATION") {
                        createReadPairs(record, ReadPairs);
                    }
                }

                for (auto& item : ReadPairs) {
                    std::string qname = item.first;
                    ReadPair& pair = item.second;
                    std::string olap = classifyOverlap(pair, variant, adjacency);
                    if (olap == "None") {
                        continue;
                    }

                    if (variant.varType == "DEL") {          
                        ReadPair newpair = simulateDeletion(pair, variant, olap, ref, snvs);
                        skipReads.insert(qname);
                        if (!newpair.isSkipped) {
                            fastq1 << "@"+newpair.readGroup << "\n"
                                    << newpair.read1Seq << "\n+\n"
                                    << newpair.read1Qual << "\n";

                            fastq2 << "@"+newpair.readGroup << "\n"
                                    << newpair.read2Seq << "\n+\n"
                                    << newpair.read2Qual << "\n";
                        }
                    }
                    if (variant.varType == "DUP") {
                        std::vector<ReadPair> newPairs = simulateTandemDuplication(pair, variant, olap, ref);
                        for (auto& newpair : newPairs) {

                            // std::cout << newpair.read1Seq << std::endl;
                            // std::cout << newpair.read2Seq << std::endl;

                            if (newpair.read1Seq.length()!= newpair.read1Qual.length()) {
                                std::cout << " error: unequal sizes" << std::endl;
                            }
                           if (newpair.read2Seq.length()!= newpair.read2Qual.length()) {
                                std::cout << " error: unequal sizes" << std::endl;
                            }

                            
                            skipReads.insert(newpair.readGroup);
                            fastq1 << "@" + newpair.readGroup << "\n"
                                    << newpair.read1Seq << "\n+\n"
                                    << newpair.read1Qual << "\n";

                            fastq2 << "@" + newpair.readGroup << "\n"
                                    << newpair.read2Seq << "\n+\n"
                                    << newpair.read2Qual << "\n";
                            skipReads.insert(newpair.readGroup);
                        }
                    }
                    if(variant.varType == "INV") {
                        if (olap != "7" && olap != "9" && olap != "10") {  
                            skipReads.insert(pair.readGroup);
                        }
                        else {
                            continue;
                        }
                            ReadPair newpair = simulateInversion(pair, variant, olap, ref);
                            fastq1 << "@" + newpair.readGroup << "\n"
                                    << newpair.read1Seq << "\n+\n"
                                    << newpair.read1Qual << "\n";

                            fastq2 << "@" + newpair.readGroup << "\n"
                                    << newpair.read2Seq << "\n+\n"
                                    << newpair.read2Qual << "\n";
                                    
                    }
                    if (variant.varType == "BALANCED_TRANSLOCATION") {
                        if (olap == "1" || olap == "2" || olap == "3" || olap == "4" || olap == "5" || olap == "6"|| olap == "11" || olap == "12" || olap == "13") {

                            ReadPair newpair = simulateBalancedTranslocation(pair, variant, olap, ref, adjacency);                         

                            if (newpair.read1Seq.length()!= newpair.read1Qual.length()) {
                                std::cout << pair.read1Pos << "-" << pair.read1End  << "Olap:" << olap << " " << newpair.readGroup << " error: unequal sizes" << std::endl;
                            }
                            if (newpair.read2Seq.length()!= newpair.read2Qual.length()) {
                                std::cout << pair.read2Pos << "-" << pair.read2End << " Olap:" << olap << " " << newpair.readGroup <<  " error: unequal sizes" << std::endl;
                            }

                            fastq1 << "@" + newpair.readGroup << "\n"
                                    << newpair.read1Seq << "\n+\n"
                                    << newpair.read1Qual << "\n";

                            fastq2 << "@" + newpair.readGroup << "\n"
                                    << newpair.read2Seq << "\n+\n"
                                    << newpair.read2Qual << "\n";
                            skipReads.insert(pair.readGroup);
                        }
                    }
                    if (variant.varType == "WHOLE_ARM_TRANSLOCATION") {
                        if (olap == "1" || olap == "2" || olap == "11" || olap == "12" || olap == "13") {

                            ReadPair newpair = simulateWholeArmTranslocation(pair, variant, olap, ref, adjacency);
                            fastq1 << "@" + newpair.readGroup << "\n"
                                    << newpair.read1Seq << "\n+\n"
                                    << newpair.read1Qual << "\n";

                            fastq2 << "@" + newpair.readGroup << "\n"
                                    << newpair.read2Seq << "\n+\n"
                                    << newpair.read2Qual << "\n";
                            skipReads.insert(pair.readGroup);

                        }

                    }
                }
            
            }
        }
    
        BamReader fullReader(bamFile);
        BamRecord newRecord;

        while (fullReader.GetNextRecord(newRecord)) {
            std::string qname = newRecord.Qname();
            // Avoid re-writing if already done
            if (skipReads.find(qname) != skipReads.end()) {
                continue;
            }

            // If not in modifiedReads: write original record
            writer.WriteRawRecord(newRecord);
        }
    
    }
    return 0;
}
