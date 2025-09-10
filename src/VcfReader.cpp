#include "VcfReader.h"
#include <sstream>
#include <string>
#include <iostream>
#include <cstring>

// Constructor
VcfReader::VcfReader(const std::string& path) : path_(path), itr_(nullptr), tbx_(nullptr), 
                                               in_(nullptr), header_(nullptr), buffer_(nullptr), eof_(false) {
    // Initialize kstring_t
    kstr_ = {0, 0, nullptr};

    // Open VCF.gz file
    in_ = hts_open(path_.c_str(), "r");
    if (!in_) {
        throw std::runtime_error("ERROR: Failed to open VCF.gz file: " + path_);
    }

    // Read header
    header_ = bcf_hdr_read(in_);
    if (!header_) {
        hts_close(in_);
        throw std::runtime_error("ERROR: Failed to read VCF header: " + path_);
    }

    // Initialize buffer
    buffer_ = bcf_init();
    if (!buffer_) {
        bcf_hdr_destroy(header_);
        hts_close(in_);
        throw std::runtime_error("ERROR: Failed to initialize BCF buffer for: " + path_);
    }

    // Load Tabix index
    tbx_ = tbx_index_load2(path_.c_str(), nullptr);
    if (!tbx_) {
        std::string idx = path_ + ".tbi";
        tbx_ = tbx_index_load2(path_.c_str(), idx.c_str());
        if (!tbx_) {
            bcf_destroy(buffer_);
            bcf_hdr_destroy(header_);
            hts_close(in_);
            throw std::runtime_error("ERROR: Failed to load Tabix index (.tbi) for: " + path_);
        }
    }
}

// Destructor
VcfReader::~VcfReader() {
    if (itr_) hts_itr_destroy(itr_);
    if (tbx_) tbx_destroy(tbx_);
    if (buffer_) bcf_destroy(buffer_);
    if (header_) bcf_hdr_destroy(header_);
    if (in_) hts_close(in_);
    if (kstr_.s) free(kstr_.s);
}

// Get header
bcf_hdr_t* VcfReader::getHeader() { 
    return header_; 
}

// Set region for querying
void VcfReader::SetRegion(const std::string& region) {
    if (itr_) {
        hts_itr_destroy(itr_);
        itr_ = nullptr;
    }
    eof_ = false;

    if (!tbx_) {
        throw std::runtime_error(" ERROR: No Tabix index loaded for " + path_);
    }

    itr_ = tbx_itr_querys(tbx_, region.c_str());
    if (!itr_) {
        throw std::runtime_error(" ERROR: Region not found in VCF.gz/.tbi index: " + region +
                                 " (check contig names, e.g., '1' vs 'chr1')");
    }
}

// Check if more records are available
bool VcfReader::HasNext() { return !eof_; }

// Get next record(s) as VcfRecord objects for valid SNVs
bool VcfReader::GetNextRecord(std::vector<VcfRecord>& records) {
    records.clear();
    if (eof_) return false;

    // Sequential reading (no region specified)
    if (!itr_) {
        int ret = bcf_read(in_, header_, buffer_);
        if (ret < 0) {
            eof_ = true;
            return false;
        }

        // Unpack the record to access REF/ALT strings
        if (bcf_unpack(buffer_, BCF_UN_STR) < 0) {
            std::cerr << " WARNING: Failed to unpack record at position " << buffer_->pos + 1 << ", skipping" << std::endl;
            return GetNextRecord(records); // Skip invalid record
        }

        // Check for valid alleles
        if (!buffer_->d.allele || buffer_->n_allele < 2) {
            std::cerr << " WARNING: Record at position " << buffer_->pos + 1 << " has insufficient alleles, skipping" << std::endl;
            return GetNextRecord(records); // Skip invalid record
        }

        // Check if REF is a single nucleotide
        const char* ref = buffer_->d.allele[0];
        if (!ref || std::strlen(ref) != 1) {
            return GetNextRecord(records); // Skip non-SNV records
        }

        // Process each ALT allele
        for (int ai = 1; ai < buffer_->n_allele; ++ai) {
            const char* alt = buffer_->d.allele[ai];
            if (!alt || std::strlen(alt) != 1) {
                continue; // Skip non-SNV ALT alleles
            }
            VcfRecord vcf_rec(header_, buffer_, ai);
            if (vcf_rec.isValidSNV()) {
                records.emplace_back(std::move(vcf_rec));
            }
        }
        return !records.empty();
    }

    // Tabix index-based reading
    kstr_.l = 0;
    int ret = hts_itr_next(in_->fp.bgzf, itr_, &kstr_, tbx_);
    if (ret < 0) {
        eof_ = true;
        return false;
    }
    if (kstr_.l == 0 || kstr_.s[0] == '#') {
        return GetNextRecord(records); // Skip headers/empty lines
    }
    if (vcf_parse1(&kstr_, header_, buffer_) < 0) {
        std::cerr << " WARNING: Skipping malformed VCF line: " << kstr_.s << std::endl;
        return GetNextRecord(records); // Skip malformed lines
    }

    // Unpack the record to access REF/ALT strings
    if (bcf_unpack(buffer_, BCF_UN_STR) < 0) {
        std::cerr << " WARNING: Failed to unpack record at position " << buffer_->pos + 1 << ", skipping" << std::endl;
        return GetNextRecord(records); // Skip invalid record
    }

    // Check for valid alleles
    if (!buffer_->d.allele || buffer_->n_allele < 2) {
        std::cerr << " WARNING: Record at position " << buffer_->pos + 1 << " has insufficient alleles, skipping" << std::endl;
        return GetNextRecord(records); // Skip invalid record
    }

    // Check if REF is a single nucleotide
    const char* ref = buffer_->d.allele[0];
    if (!ref || std::strlen(ref) != 1) {
        return GetNextRecord(records); // Skip non-SNV records
    }

    // Process each ALT allele
    for (int ai = 1; ai < buffer_->n_allele; ++ai) {
        const char* alt = buffer_->d.allele[ai];
        if (!alt || std::strlen(alt) != 1) {
            continue; // Skip non-SNV ALT alleles
        }
        VcfRecord vcf_rec(header_, buffer_, ai);
        if (vcf_rec.isValidSNV()) {
            records.emplace_back(std::move(vcf_rec));
        }
    }
    return !records.empty();
}