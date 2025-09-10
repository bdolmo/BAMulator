#include "VcfRecord.h"
#include <iostream>
#include <limits>

VcfRecord::VcfRecord(bcf_hdr_t* hdr, bcf1_t* rec) 
    : depth_(-1), vaf_(std::numeric_limits<double>::quiet_NaN()) {
    if (!hdr || !rec) {
        std::cerr << "WARNING: Null header or record provided to VcfRecord constructor" << std::endl;
        chr_ = ".";
        position_ = 0;
        ref_ = "";
        alt_ = "";
        return;
    }

    // Unpack the record to access REF/ALT
    if (bcf_unpack(rec, BCF_UN_STR) < 0) {
        std::cerr << "WARNING: Failed to unpack record at position " << rec->pos + 1 << std::endl;
        chr_ = ".";
        position_ = 0;
        ref_ = "";
        alt_ = "";
        return;
    }

    // Get chromosome
    const char* chrom = bcf_hdr_id2name(hdr, rec->rid);
    if (chrom) {
        chr_ = chrom;
    } else {
        std::cerr << "WARNING: Invalid chromosome ID " << rec->rid << " at position " << rec->pos + 1 << std::endl;
        chr_ = ".";
    }

    // Set position (1-based)
    position_ = static_cast<int64_t>(rec->pos) + 1;

    // Get REF and first ALT allele
    if (rec->d.allele && rec->n_allele >= 2) {
        ref_ = rec->d.allele[0] ? rec->d.allele[0] : "";
        alt_ = rec->d.allele[1] ? rec->d.allele[1] : "";
    } else {
        std::cerr << "WARNING: Record at position " << position_ << " has insufficient alleles" << std::endl;
        ref_ = "";
        alt_ = "";
    }
}

VcfRecord::VcfRecord() 
    : chr_("."), position_(0), ref_(""), alt_(""), depth_(-1), 
      vaf_(std::numeric_limits<double>::quiet_NaN()) {
}

bool VcfRecord::isValidSNV() const {
    return !chr_.empty() && chr_ != "." && 
           position_ > 0 && 
           ref_.length() == 1 && 
           alt_.length() == 1;
}