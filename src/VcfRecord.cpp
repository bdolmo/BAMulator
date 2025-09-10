#include "VcfRecord.h"
#include <iostream>
#include <limits>
#include <cstring>
#include <sstream>

VcfRecord::VcfRecord(bcf_hdr_t* hdr, bcf1_t* rec, int alt_index) 
    : depth_(-1), vaf_(std::numeric_limits<double>::quiet_NaN()), genotype_(".") {
    if (!hdr || !rec) {
        std::cerr << "WARNING: Null header or record provided to VcfRecord constructor" << std::endl;
        chr_ = ".";
        position_ = 0;
        ref_ = "";
        alt_ = "";
        return;
    }

    // Unpack the record to access REF/ALT and FORMAT
    if (bcf_unpack(rec, BCF_UN_STR | BCF_UN_FMT) < 0) {
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

    // Get REF and specified ALT allele
    if (rec->d.allele && rec->n_allele > alt_index && alt_index >= 1) {
        ref_ = rec->d.allele[0] ? rec->d.allele[0] : "";
        alt_ = rec->d.allele[alt_index] ? rec->d.allele[alt_index] : "";
    } else {
        std::cerr << "WARNING: Record at position " << position_ << " has insufficient alleles or invalid alt_index " << alt_index << std::endl;
        ref_ = "";
        alt_ = "";
    }

    // Parse GT field for the first sample
    int ngt = 0;
    int32_t* gt_arr = nullptr;
    int gt_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "GT");
    if (gt_id >= 0 && bcf_get_format_int32(hdr, rec, "GT", &gt_arr, &ngt) > 0 && ngt >= 1) {
        // Assuming single sample (ngt >= 1 for at least one genotype)
        int32_t g1 = gt_arr[0];
        int32_t g2 = ngt > 1 ? gt_arr[1] : -1; // Second allele if diploid
        std::stringstream ss;
        if (bcf_gt_is_missing(g1)) {
            ss << ".";
        } else {
            ss << bcf_gt_allele(g1);
        }
        if (g2 >= 0) {
            ss << (bcf_gt_is_phased(g2) ? "|" : "/") << (bcf_gt_is_missing(g2) ? "." : std::to_string(bcf_gt_allele(g2)));
        }
        genotype_ = ss.str();
        free(gt_arr);
    } else {
        std::cerr << "WARNING: No GT field or invalid GT data at position " << position_ << std::endl;
        genotype_ = ".";
    }
}

VcfRecord::VcfRecord() 
    : chr_("."), position_(0), ref_(""), alt_(""), depth_(-1), 
      vaf_(std::numeric_limits<double>::quiet_NaN()), genotype_(".") {
}

bool VcfRecord::isValidSNV() const {
    return !chr_.empty() && chr_ != "." && 
           position_ > 0 && 
           ref_.length() == 1 && 
           alt_.length() == 1;
}
