#ifndef VCFREADER_H
#define VCFREADER_H

#include <htslib/vcf.h>
#include <htslib/tbx.h>
#include <htslib/kstring.h>
#include <string>
#include <vector>
#include "VcfRecord.h"

class VcfReader {
public:
    // Constructor
    VcfReader(const std::string& path);
    // Destructor
    ~VcfReader();
    // Get header
    bcf_hdr_t* getHeader();
    // Set region for querying
    void SetRegion(const std::string& region);
    // Check if more records are available
    bool HasNext();
    // Get next record(s) as VcfRecord objects for valid SNVs
    bool GetNextRecord(std::vector<VcfRecord>& records);

private:
    std::string path_;          // Path to VCF.gz file
    htsFile* in_;               // File handle
    bcf_hdr_t* header_;         // VCF header
    bcf1_t* buffer_;            // Record buffer
    tbx_t* tbx_;                // Tabix index
    hts_itr_t* itr_;            // Region iterator
    kstring_t kstr_;            // Buffer for Tabix line reading
    bool eof_;                  // End-of-file flag
};
#endif