#ifndef VCFRECORD_H
#define VCFRECORD_H
#include <iostream>

#include <htslib/vcf.h>
#include <string>
#include <vector>

class VcfRecord {
public:
    // Constructor from bcf1_t record, header, and optional ALT allele index
    VcfRecord(bcf_hdr_t* hdr, bcf1_t* rec, int alt_index = 1);
    // Default constructor for empty record
    VcfRecord();
    // Destructor
    ~VcfRecord() = default;

    // Check if the record is a valid SNV (single-nucleotide REF and ALT)
    bool isValidSNV() const;

    // Getters
    const std::string& getChrom() const { return chr_; }
    int64_t getPosition() const { return position_; }
    const std::string& getRef() const { return ref_; }
    const std::string& getAlt() const { return alt_; }
    int getDepth() const { return depth_; }
    double getVaf() const { return vaf_; }
    const std::vector<std::string>& getQnameVector() const { return qnameVector_; }
    const std::string& getGenotype() const { return genotype_; }

private:
    std::string chr_;                    // Chromosome name
    int64_t position_;                   // 1-based position
    std::string ref_;                    // Reference allele
    std::string alt_;                    // Alternate allele
    int depth_;                          // Depth (default -1)
    double vaf_;                         // Variant allele frequency (default NaN)
    std::vector<std::string> qnameVector_; // QNAME vector (default empty)
    std::string genotype_;               // Genotype (e.g., "1|0", "1|1")
};

#endif
