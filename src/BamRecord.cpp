#include "BamRecord.h"
#include "RefFasta.h"
#include <stdexcept>
#include <string>
#include <cstring>  // For memcpy and memset
#include <cstdlib>  // For malloc, free
#include <stdexcept> // For std::runtime_error
#include <htslib/sam.h>
#include <regex>
#include <iostream>

BamRecord::BamRecord() {
}

BamRecord::BamRecord(bam1_t *b, bam_hdr_t *h) : b_(b), h_(h) {
    if (b_ == nullptr) {
        throw std::runtime_error("Null pointer received for BAM record");
    }
    if (h_ == nullptr) {
        throw std::runtime_error("Null pointer received from BAM header");
    }

    _position = b->core.pos;
    _chrID = b->core.tid;
    _chrMateID = b->core.mtid;

    if (b->core.tid >= 0) {
        _chrName = std::string(h->target_name[b->core.tid]);
    } 
    else {
        _chrName = "*";
    }
    if (b->core.mtid >= 0) {
        _chrMateName = std::string(h->target_name[b->core.mtid]);
    } 
    else {
        _chrMateName = "*";
    }
    _qname = bam_get_qname(b);  // Extract QNAME from the bam1_t structure
    _flag = b->core.flag; // Derive the flag from the bam1_t structure
    const uint32_t *cigar = bam_get_cigar(b);

    std::string merdatest;
    for (unsigned i = 0; i < b->core.n_cigar; ++i) {
        int cigarLen = bam_cigar_oplen(cigar[i]);
        char cigarOp = bam_cigar_opchr(cigar[i]);
        _cigarString += std::to_string(cigarLen) + cigarOp;
        merdatest += std::to_string(cigarLen) + cigarOp;
        
    }
    _insertSize = b->core.isize;
    // Borrowed from SeqLib
    uint8_t* p = bam_aux_get(b, "MD");
    if (!p) {
        _mdString = "";
    }
    else {
        char* pp = bam_aux2Z(p);
        if (!pp) {
            _mdString = "";
        }
        else {
            _mdString = std::string(pp);
        }       
    }

    _mapQual = b->core.qual;
    _matePos = b->core.mpos;


    uint8_t *seq = bam_get_seq(b); // Get compressed sequence
    int32_t len = b->core.l_qseq; // Get length of the sequence
    _seq.reserve(len); // Reserve space to avoid multiple reallocations

    for (int i = 0; i < len; ++i) {
        uint8_t base = bam_seqi(seq, i); // Get base in 4-bit integer format
        switch (base) {
            case 1: _seq += 'A'; break;
            case 2: _seq += 'C'; break;
            case 4: _seq += 'G'; break;
            case 8: _seq += 'T'; break;
            case 15: _seq += 'N'; break; // 15 represents 'N' in BAM format
            default: _seq += '?'; break; // Just in case there is an unknown base
        }
    }
    _end = _position + len;

    // Set the quality
    uint8_t* _qual = bam_get_qual(b);

    // std::string _qualString;
    _qualString.reserve(len);
    for (int i = 0; i < len; ++i) {
        // Convert the Phred quality score to the ASCII representation (Phred+33)
        _qualString.push_back(_qual[i] + 33);
    }


    int calculated_length = 0;
    std::string number;


    for (char ch : _cigarString) {
        if (std::isdigit(ch)) {
            number += ch; // Build the number as a string
        } else {
            int length = std::stoi(number); // Convert the number string to an integer
            number.clear(); // Clear the number string for the next operation

            if (ch == 'M' || ch == 'I' || ch == 'S' || ch == '=' || ch == 'X') {
                calculated_length += length;
            }
            // Note: 'D', 'N', 'H', 'P' are skipped as they do not consume query sequence bases
        }
    }

    // std::cout << calculated_length << std::endl;
    if (calculated_length != len) {
        if (!IsUnmapped()) {
            throw std::runtime_error("CIGAR string and query sequence are of different lengths1");
        }
    }

}


BamRecord BamRecord::DeepCopy() const {
    if (!b_ || !h_)
        throw std::runtime_error("Cannot deep copy: null BAM record or header");

    // Clone the bam1_t structure
    bam1_t* b_new = bam_dup1(b_);
    if (!b_new)
        throw std::runtime_error("Failed to duplicate bam1_t structure");

    // Header can be shared safely, as it’s immutable
    BamRecord copy(b_new, h_);

    // Manually copy all cached string fields and values
    copy._qname       = _qname;
    copy._seq         = _seq;
    copy._qualString  = _qualString;
    copy._cigarString = _cigarString;
    copy._chrName     = _chrName;
    copy._chrMateName = _chrMateName;
    copy._mdString    = _mdString;

    // Copy numeric attributes
    copy._position    = _position;
    copy._matePos     = _matePos;
    copy._chrID       = _chrID;
    copy._chrMateID   = _chrMateID;
    copy._mapQual     = _mapQual;
    copy._insertSize  = _insertSize;
    copy._flag        = _flag;

    return copy;
}


bool BamRecord::IsPaired() const {
    return (b_->core.flag & 0x1) != 0;
}

bool BamRecord::IsProperlyPaired() const {
    return (b_->core.flag & 0x2) != 0;
}

bool BamRecord::IsFirstMate() const {
    return b_ && (b_->core.flag & BAM_FREAD1);
}

std::string BamRecord::Qual() const {
    if (b_ == nullptr) {
        throw std::runtime_error("Null pointer encountered in BamRecord");
    }
    return _qualString;
}

bool BamRecord::IsSecondary() const {
    return (b_->core.flag & 0x100) != 0;
}

bool BamRecord::IsSupplementary() const {
    return (b_->core.flag & 0x800) != 0;
}

std::string BamRecord::CalculateMDTag(const std::string& seq, const std::string& ref, const std::string& cigar) {
    std::string md_tag;
    int matches = 0;
    int ref_pos = 0;  // Position in reference sequence
    int seq_pos = 0;  // Position in read sequence
    std::string number;

    bool last_was_mismatch_or_deletion = false;  // Track if the last base was a mismatch or deletion

    // Loop through the CIGAR string to process each operation
    for (size_t i = 0; i < cigar.length(); ++i) {
        char ch = cigar[i];
        if (isdigit(ch)) {
            number += ch;  // Build the length of the CIGAR operation
        } else {
            int len = std::stoi(number);  // Convert accumulated number to an integer
            number.clear();  // Clear the number for the next CIGAR operation

            switch (ch) {
                case 'M':  // Match or mismatch
                    for (int j = 0; j < len; ++j) {
                        if (seq[seq_pos] == ref[ref_pos]) {
                            if (last_was_mismatch_or_deletion) {
                                last_was_mismatch_or_deletion = false;  // Reset flag
                            }
                            matches++;  // Record a match
                        } else {
                            if (matches > 0) {
                                md_tag += std::to_string(matches);  // AppFposition match count
                                matches = 0;
                            }
                            if (last_was_mismatch_or_deletion) {
                                md_tag += "0";  // Add '0' between mismatches or after deletions
                            }
                            md_tag += ref[ref_pos];  // Append the mismatched reference base
                            last_was_mismatch_or_deletion = true;  // Mark as mismatch
                        }
                        seq_pos++;  // Move to the next base in the sequence
                        ref_pos++;  // Move to the next base in the reference
                    }
                    break;

                case 'I':  // Insertion (skip read bases, no reference comparison)
                    seq_pos += len;  // Insertion consumes read bases, but not reference bases
                    break;

                case 'D':  // Deletion (append deleted reference bases)
                    if (matches > 0) {
                        md_tag += std::to_string(matches);  // Append match count before deletion
                        matches = 0;
                    }
                    if (last_was_mismatch_or_deletion) {
                        md_tag += "0";  // Insert '0' between a deletion and previous mismatch or deletion
                    }
                    md_tag += "^";  // Indicate deletion in MD tag
                    for (int j = 0; j < len; ++j) {
                        md_tag += ref[ref_pos];  // Append the deleted reference base
                        ref_pos++;  // Move to the next base in the reference
                    }
                    last_was_mismatch_or_deletion = true;  // Mark as deletion
                    break;

                case 'S':  // Soft clipping (skip read bases, no reference comparison)
                    seq_pos += len;  // Soft clipping consumes read bases, but not reference bases
                    break;
            }
        }
    }

    // Append any remaining match count at the end of the sequence
    if (matches > 0) {
        md_tag += std::to_string(matches);  // Append remaining match count
    }

    return md_tag;  // Return the final MD tag
}

std::vector<uint32_t> BamRecord::getCigarVector() const {
    std::vector<uint32_t> cigarVector;
    std::regex cigarRegex("(\\d+)([MIDNSHP=X])"); // Regular expression to parse CIGAR string
    std::smatch match;
    std::string tempCigar = _cigarString;

    // Use regular expression to extract each CIGAR operation
    while (std::regex_search(tempCigar, match, cigarRegex)) {
        uint32_t len = std::stoi(match[1].str()); // Length of the CIGAR operation
        char opChar = match[2].str()[0]; // Character representing the CIGAR operation
        uint32_t op; // Numeric code for the CIGAR operation

        // Convert the CIGAR operation character to BAM CIGAR operation code
        switch (opChar) {
            case 'M': op = BAM_CMATCH; break;
            case 'I': op = BAM_CINS; break;
            case 'D': op = BAM_CDEL; break;
            case 'N': op = BAM_CREF_SKIP; break;
            case 'S': op = BAM_CSOFT_CLIP; break;
            case 'H': op = BAM_CHARD_CLIP; break;
            case 'P': op = BAM_CPAD; break;
            case '=': op = BAM_CEQUAL; break;
            case 'X': op = BAM_CDIFF; break;
            default: throw std::runtime_error("Unknown CIGAR operation: " + std::string(1, opChar));
        }

        // Encode the CIGAR operation into the 32-bit format used by BAM files
        uint32_t encodedCigar = (len << BAM_CIGAR_SHIFT) | op;
        cigarVector.push_back(encodedCigar);

        // Continue searching from the end of the last match
        tempCigar = match.suffix().str();
    }

    return cigarVector;
}

// Helper function for converting base characters to their 4-bit encoded values
uint8_t base2val(char base) {
    switch (base) {
        case 'A': case 'a': return 1;
        case 'C': case 'c': return 2;
        case 'G': case 'g': return 4;
        case 'T': case 't': return 8;
        default: return 15; // 'N' or unknown base
    }
}

bam1_t* BamRecord::GetBam1_t() const {
    return b_;
}

bam1_t* BamRecord::ToBam1_t() const {
    // Allocate a new bam1_t structure
    bam1_t* b_new = bam_init1();
    if (!b_new) throw std::runtime_error("Failed to allocate bam1_t structure");

    // Convert the CIGAR string to a vector of uint32_t
    std::vector<uint32_t> cigarVector = getCigarVector();
    size_t qname_len = _qname.length() + 1; // +1 for null terminator
    size_t seq_len = (_seq.length() + 1) / 2; // Sequence is 4-bit encoded
    size_t qual_len = _seq.length(); // Quality score length equals sequence length
    size_t cigar_bytes = cigarVector.size() * 4; // 4 bytes per CIGAR op

    // Calculate total data size needed
    size_t total_data_size = qname_len + cigar_bytes + seq_len + qual_len;
    // std::cout << "Total data size: " << total_data_size << std::endl;
    // Ensure space for all data components is allocated
    if (b_new->m_data < total_data_size) {
        b_new->data = (uint8_t*)realloc(b_new->data, total_data_size);
        if (!b_new->data) {
            bam_destroy1(b_new); // Free memory if reallocation failed
            throw std::runtime_error("Failed to reallocate memory for bam1_t structure");
        }
        b_new->m_data = total_data_size;
    }

    // Set basic fields from BamRecord to bam1_t
    b_new->core.tid = _chrID;
    b_new->core.pos = _position;
    b_new->core.qual = _mapQual;
    b_new->core.mtid = _chrMateID;
    b_new->core.l_qseq = _seq.length();
    b_new->core.n_cigar = cigarVector.size();
    b_new->core.flag = _flag;
    b_new->l_data = total_data_size;
    b_new->core.isize = _insertSize;
    b_new->core.mpos = _matePos;
    b_new->core.l_qname = qname_len;

    if (!_mdString.empty()) {
        // Append the MD tag as a null-terminated string
        if (bam_aux_append(b_new, "MD", 'Z', _mdString.length() + 1, 
                        reinterpret_cast<const uint8_t*>(_mdString.c_str())) < 0) {
            throw std::runtime_error("Failed to append MD tag to bam1_t structure");
        }
    }

    // Set the query name (QNAME)
    if (qname_len > 0) {
        std::memcpy(b_new->data, _qname.c_str(), qname_len);
    }

    // Set the CIGAR operations
    uint32_t* cigar_ptr = bam_get_cigar(b_new);
    for (size_t i = 0; i < cigarVector.size(); ++i) {
        cigar_ptr[i] = cigarVector[i];
    }

    // Set the sequence
    uint8_t* seq = bam_get_seq(b_new);
    for (size_t i = 0; i < _seq.length(); ++i) {
        seq[i >> 1] = (i & 1) ? (seq[i >> 1] | base2val(_seq[i])) : (base2val(_seq[i]) << 4);
    }

    // Set the quality
    uint8_t* quality = bam_get_qual(b_new);
    for (size_t i = 0; i < qual_len; ++i) {
        quality[i] = _qualString[i]- 33;
    }

    memcpy(bam_get_qual(b_new), quality, qual_len); // Copy converted quality scores back into the BAM structure

    return b_new;
}

void BamRecord::SetPosition(int64_t newPos) {
    if (b_ == nullptr) {
        throw std::runtime_error("Null pointer encountered in BamRecord");
    }
    _position = newPos; // Update the internal representation
    b_->core.pos = newPos; // Update the position in the bam1_t struct
}

void BamRecord::SetMatePosition(int64_t pos) {
    if (b_ == nullptr) {
        throw std::runtime_error("Null pointer encountered in BamRecord");
    }
    _matePos = pos;
    b_->core.mpos = pos;
}

void BamRecord::SetFlag(uint16_t newFlag) {
    if (b_ == nullptr) {
        throw std::runtime_error("Null pointer encountered in BamRecord");
    }
    b_->core.flag = newFlag;
    _flag = newFlag;
}


void BamRecord::UpdateSeq(const std::string& seq, const std::string& new_cigar, RefFasta &refFasta) {
    // Step 1: Skip update for unmapped reads
    if (b_->core.flag & BAM_FUNMAP) {
        std::cout << "Skipping update for unmapped read." << std::endl;
        return;  // No need to update sequence or MD tag for unmapped reads
    }

    // Step 2: Verify that new CIGAR string length matches sequence length
    int calculated_length = 0;
    std::string number;
    for (char ch : new_cigar) {
        if (std::isdigit(ch)) {
            number += ch;  // Build the number as a string
        } else {
            int length = std::stoi(number);  // Convert the number string to an integer
            number.clear();  // Clear the number string for the next operation

            if (ch == 'M' || ch == 'I' || ch == 'S' || ch == '=' || ch == 'X') {
                calculated_length += length;
            }
        }
    }

    if (calculated_length != seq.length()) {
        std::cout << "ERROR: CIGAR: " << new_cigar << " Sequence: " << seq << " CIGAR length: " << calculated_length << " Sequence length: " << seq.length() << std::endl;
        std::cout << _qname << " " << _position << std::endl << std::endl;
        throw std::runtime_error("CIGAR string and query sequence lengths differ.");
    }

    // Step 3: Calculate the MD tag before memory operations
    int ref_span = CalculateReferenceSpan(new_cigar);
    int read_end_pos = _position + ref_span;

    std::string ref_seq = refFasta.fetchSequence(_chrName, _position + 1, read_end_pos);
    std::transform(ref_seq.begin(), ref_seq.end(), ref_seq.begin(), ::toupper);

    if (ref_seq.empty()) {
        throw std::runtime_error("Failed to fetch reference sequence.");
    }

    std::string md_tag = CalculateMDTag(seq, ref_seq, new_cigar);
    size_t md_tag_size = md_tag.length() + 1;  // Add 1 for null terminator

    _mdString = md_tag;
    // Step 4: Calculate the new required size for BAM data, including the new CIGAR
    int qname_len = b_->core.l_qname;
    int old_cigar_len = b_->core.n_cigar * 4;  // Original CIGAR size (4 bytes per CIGAR operation)
    
    // New CIGAR length calculation
    int new_cigar_len = 0;
    number.clear();
    for (char ch : new_cigar) {
        if (std::isdigit(ch)) {
            number += ch;
        } else {
            int length = std::stoi(number);  // Length of the new CIGAR operation
            number.clear();
            new_cigar_len++;  // Increment for each operation
        }
    }
    new_cigar_len *= 4;  // Each CIGAR operation is 4 bytes
    _cigarString = new_cigar;

    int seq_len = (seq.length() + 1) / 2;  // 4-bit encoded sequence
    int qual_len = seq.length();  // 1 byte per base for quality scores
    int old_aux_len = bam_get_l_aux(b_);
    int new_aux_len = old_aux_len + md_tag_size;  // Additional space for new MD tag

    // The required size is now based on the new CIGAR length instead of the old one
    size_t required_size = qname_len + new_cigar_len + seq_len + qual_len + new_aux_len;

    // Backup old data
    uint8_t* oldd = (uint8_t*)malloc(b_->l_data);
    if (!oldd) {
        throw std::runtime_error("Failed to allocate memory for backing up old data.");
    }
    memcpy(oldd, b_->data, b_->l_data);

    // Check if we need to reallocate BAM data to accommodate new MD tag, sequence, and CIGAR
    // if (required_size > b_->m_data) {
        uint8_t* new_data = (uint8_t*)calloc(required_size, sizeof(uint8_t));
        if (!new_data) {
            free(oldd);
            throw std::runtime_error("Failed to allocate memory for BAM record.");
        }

        // Copy old data (before sequence and auxiliary modifications) into the newly allocated space
        memcpy(new_data, b_->data, qname_len);
        free(b_->data);  // Free old memory
        b_->data = new_data;  // Point to the newly allocated data
        b_->m_data = required_size;  // Update max data size
    // }

    // Step 5: Add back qname and new cigar
    memcpy(b_->data, oldd, qname_len);

    // Update BAM record sizes
    b_->l_data = required_size;
    b_->core.n_cigar = new_cigar_len / 4;  // Update CIGAR length in BAM core
    b_->core.l_qseq = seq.length();

    // Step 6: Update the sequence (SEQ) and quality (QUAL) fields
    int slen = seq.length();

    _seq = "";
    _seq.reserve(slen);
    uint8_t* m_bases = b_->data + qname_len + new_cigar_len;
    for (int i = 0; i < seq.length(); ++i) {
        uint8_t base = 15;  // Default to 'N'

        if (seq.at(i) == 'A') {
            base = 1;
            _seq += 'A';
        }
        else if (seq.at(i) == 'C') {
            base = 2;
            _seq += 'C';
        }
        else if (seq.at(i) == 'G') {
            base = 4;
            _seq += 'G';
        }
        else if (seq.at(i) == 'T') {
            base = 8;
            _seq += 'T';
        }
        else {
            _seq += 'N';
        }

        if (i % 2 == 0) {
            m_bases[i / 2] = base << 4;  // First base into the higher 4 bits
        } else {
            m_bases[i / 2] |= base;  // Second base into the lower 4 bits
        }
    }

    // // Update quality scores
    // uint8_t* quality = bam_get_qual(b_);
    // for (int i = 0; i < qual_len; ++i) {
    //     quality[i] = _qualString[i] - 33;  // Update quality scores
    // }

    // _qualString = "";
    // _qualString.reserve(seq_len);
    // for (int i = 0; i < seq_len; ++i) {
    //     // Convert the Phred quality score to the ASCII representation (Phred+33)
    //     _qualString.push_back(quality[i] + 33);
    // }

    // Update quality scores
    // uint8_t* quality = bam_get_qual(b_);
    // int original_qual_len = qual_len;  // Store original length of quality scores

    // // If the sequence length is different, we may need to modify the quality scores
    // if (seq_len > qual_len) {
    //     // Sequence is longer than the original: add artificial quality scores
    //     _qualString.reserve(seq_len);  // Adjust the size of the string buffer
    //     for (int i = 0; i < seq_len; ++i) {
    //         if (i < original_qual_len) {
    //             // Convert existing Phred quality scores to ASCII representation (Phred+33)
    //             _qualString.push_back(quality[i] + 33);
    //         } else {
    //             // Generate a default Phred quality score for new bases
    //             uint8_t artificial_quality = 30;  // Use default Phred score 30
    //             _qualString.push_back(artificial_quality + 33);  // Convert to ASCII
    //         }
    //     }
    // } else if (seq_len < qual_len) {
    //     // Sequence is shorter than the original: truncate the quality scores
    //     _qualString.reserve(seq_len);  // Adjust the size of the string buffer
    //     for (int i = 0; i < seq_len; ++i) {
    //         // Convert the existing Phred quality scores to ASCII representation (Phred+33)
    //         _qualString.push_back(quality[i] + 33);
    //     }
    // } else {
    //     // Sequence length matches quality score length: no need for adjustment
    //     _qualString = "";
    //     _qualString.reserve(seq_len);
    //     for (int i = 0; i < seq_len; ++i) {
    //         _qualString.push_back(quality[i] + 33);  // Convert Phred+33 back to ASCII
    //     }
    // }


        // uint8_t* quality = bam_get_qual(b_);
        // _qualString.clear();
        // _qualString.reserve(seq.length());

        // for (int i = 0; i < seq.length(); ++i) {
        //     uint8_t q = 0;  // Set Phred score = 30 for all bases, adjust as needed
        //     quality[i] = q;
        //     _qualString.push_back(q + 33);  // Store printable quality string (Phred+33)
        // }

        // Location of original quality scores in the old buffer
        uint8_t* old_quality = oldd + qname_len + old_cigar_len + ((b_->core.l_qseq + 1) / 2);

        // Location of new quality scores
        uint8_t* quality = bam_get_qual(b_);

        // Copy over qualities
        _qualString.clear();
        _qualString.reserve(seq.length());

        for (int i = 0; i < seq.length(); ++i) {
            quality[i] = old_quality[i];  // Copy original Phred score
            _qualString.push_back(static_cast<char>(quality[i] + 33));  // Phred+33 for printable string
        }




    // Step 7: Append the MD tag as an auxiliary field
    bam_aux_append(b_, "MD", 'Z', md_tag_size, (const uint8_t*)md_tag.c_str());

    // Step 8: Restore old auxiliary data after appending MD tag
    uint8_t* t = bam_get_aux(b_);
    int old_aux_spot = qname_len + new_cigar_len + (b_->core.l_qseq + 1) / 2 + b_->core.l_qseq;
    memcpy(t, oldd + old_aux_spot, old_aux_len);  // Copy old auxiliary data

    // Update auxiliary field size
    b_->l_data = qname_len + new_cigar_len + seq_len + qual_len + new_aux_len;

    // Step 9: Clean up
    free(oldd);
}


// Function to calculate the reference span (reference start to end position) based on CIGAR string
int BamRecord::CalculateReferenceSpan(const std::string& cigar) {
    int ref_span = 0;
    std::string number;
    
    for (char ch : cigar) {
        if (isdigit(ch)) {
            number += ch;  // Build the number as a string
        } else {
            int length = std::stoi(number);
            number.clear();

            // Only count operations that consume reference positions
            if (ch == 'M' || ch == 'D' || ch == 'N' || ch == '=' || ch == 'X') {
                ref_span += length;  // These operations consume reference bases
            }
            // Insertions ('I') and soft clips ('S') do not consume reference bases
        }
    }
    return ref_span;
}

void BamRecord::SetQname(const std::string& n) {
    if (!b_) throw std::runtime_error("BAM record is null");

    // Compute the lengths
    size_t new_qname_len = n.length() + 1;
    size_t old_qname_len = b_->core.l_qname;
    size_t old_data_len = b_->l_data;
    size_t nonq_len = old_data_len - old_qname_len;

    // Allocate new data block
    uint8_t* new_data = (uint8_t*)calloc(new_qname_len + nonq_len, 1);
    if (!new_data) throw std::runtime_error("Failed to allocate memory for new BAM data");

    // Copy new qname
    memcpy(new_data, n.c_str(), n.length());
    new_data[n.length()] = '\0';  // null terminator

    // Copy rest of data after old qname
    memcpy(new_data + new_qname_len, b_->data + old_qname_len, nonq_len);

    // Replace old data
    free(b_->data);
    b_->data = new_data;

    // Update fields
    b_->core.l_qname = new_qname_len;
    b_->l_data = new_qname_len + nonq_len;
    b_->m_data = b_->l_data;

    // Update cached string if applicable
    _qname = n;
}


// void BamRecord::SetQname(const std::string& n) {
//     // Copy out the non-qname data
//     size_t nonq_len = b_->l_data - b_->core.l_qname;
//     uint8_t* nonq = (uint8_t*)malloc(nonq_len);
//     if (!nonq) throw std::runtime_error("Failed to allocate memory for non-qname data");
//     memcpy(nonq, b_->data + b_->core.l_qname, nonq_len);

//     // Clear the old data and allocate the new amount
//     free(b_->data);
//     b_->data = (uint8_t*)calloc(nonq_len + n.length() + 1, 1);
//     if (!b_->data) {
//         free(nonq); // Avoid memory leak
//         throw std::runtime_error("Failed to allocate memory for new data");
//     }
    
//     // Add in the new qname
//     memcpy(b_->data, (uint8_t*)n.c_str(), n.length() + 1); // +1 for '\0' terminator

//     // Update the sizes
//     b_->l_data = nonq_len + n.length() + 1;
//     b_->core.l_qname = n.length() + 1;
    
//     // Copy over the old data after the new qname
//     memcpy(b_->data + b_->core.l_qname, nonq, nonq_len);
//     free(nonq); // Free the temporary buffer

//     // Reset the max size
//     b_->m_data = b_->l_data;
// }

void BamRecord::SetQualities(const std::string& n, int offset) {
    if (!n.empty() && n.length() != static_cast<size_t>(b_->core.l_qseq))
        throw std::invalid_argument("New quality score string should be the same length as sequence length");
    
    // Length of qual is always same as seq. If empty qual, just set first bit of qual to 0 (equivalent to '!')
    if (n.empty()) {
        uint8_t* r = bam_get_qual(b_);
        r[0] = 0; // In BAM format, a quality score of 0 is used to represent '!'
        return;
    }

    // Convert quality string to numeric values
    char * q = strdup(n.c_str());
    if (!q) throw std::runtime_error("Failed to duplicate quality string");
    for (size_t i = 0; i < n.length(); ++i) {
        q[i]; // Adjust quality scores based on the offset (typically 33 or 64)
    }
    memcpy(bam_get_qual(b_), q, n.length()); // Copy converted quality scores back into the BAM structure
    free(q); // Free the duplicated quality string
}

void BamRecord::SetInsertSize(int32_t isize) {
    _insertSize = isize;
    if (b_) {
        b_->core.isize = isize;
    }
}

int BamRecord::chrID() const {
    return _chrID;
}

int BamRecord::chrMateID() const {
    return _chrMateID;
}

int64_t BamRecord::Position() const{
    return _position;
}


int64_t BamRecord::End() const{
    return _end;
}

std::string BamRecord::chrName() const {
    return _chrName;
}

std::string BamRecord::chrMateName() const {
    return _chrMateName;
}

bool BamRecord::IsReverseStrand() const {
        return (b_->core.flag & 0x10) != 0;
}


std::string BamRecord::cigarString() const {
    return _cigarString;
}

int BamRecord::mapQual() const {
    return _mapQual;
}

std::string BamRecord::MDtag() const {
    return _mdString;
}

std::string BamRecord::Seq() const {
    return _seq;
}

std::string BamRecord::Qname() const {
    // Check if the bam1_t structure is valid
    if (!b_ || !(b_->data)) {
        throw std::runtime_error("Invalid BAM record");
    }

    // The QNAME is stored at the beginning of the data section of bam1_t structure
    // and is NULL-terminated.
    const char* qname = bam_get_qname(b_);

    // Convert C string to C++ string and return
    return std::string(qname);
}

int32_t BamRecord::InsertSize() const {
    return _insertSize;
}

int64_t BamRecord::matePos() const {
    return _matePos;
}

// Member function to check if the read is unmapped
bool BamRecord::IsUnmapped() const {
    if (b_ == nullptr) {
        throw std::runtime_error("Null pointer encountered in BamRecord");
    }
    return (b_->core.flag & BAM_FUNMAP) != 0;
}

// Member function to check if the mate is unmapped
bool BamRecord::IsMateUnmapped() const {
    if (b_ == nullptr) {
        throw std::runtime_error("Null pointer encountered in BamRecord");
    }
    return (b_->core.flag & BAM_FMUNMAP) != 0;
}

bool BamRecord::IsSecondMate() const {
    if (b_ == nullptr) {
        throw std::runtime_error("Null pointer encountered in BamRecord");
    }
    return (b_->core.flag & BAM_FREAD2) != 0;
}

bool BamRecord::IsMateReverseStrand() const {
    if (b_ == nullptr) {
        throw std::runtime_error("Null pointer encountered in BamRecord");
    }
    return (b_->core.flag & BAM_FMREVERSE) != 0;
}

