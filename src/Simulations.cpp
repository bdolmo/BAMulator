#include "BamReader.h"
#include "BamRecord.h"
#include "BamWriter.h"
#include "RefFasta.h"
#include <iostream>
#include <fstream>
#include <htslib/sam.h>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <random>
#include "ssw_cpp.h"
#include <unordered_set>
#include <cctype> // for std::toupper
#include <tuple>

struct Variant
{
    std::string bamFile;        // BAM file this variant comes from
    std::string chr;            // Primary chromosome
    int64_t position = -1;      // Useful for SNVs/indels, usually == start

    int64_t start = -1;         // For CNVs/SVs/SNVs
    int64_t end = -1;

    std::string refSeq;         // For SNVs/indels
    std::string altSeq;

    std::string varType;        // SNV, INDEL, DEL, DUP, INS, INV, TRA, etc.
    float vaf = 0.0f;           // Variant Allele Frequency (Clonal Proportion)

    // Optional SV-specific fields
    std::string chr2 = ".";     // For translocations or complex SVs
    int64_t start2 = -1;        // Partner SV coordinates
    int64_t end2 = -1;

    std::string svId;           // SVID (e.g., TTN_ins)
    std::string svMechanism;    // SV mechanism (e.g., random, templated)
    std::string homology;       // Sequence homology if applicable
    std::string clone;          // Clone name or label (e.g., clone1)
    std::string genotype;       // e.g., 0/1
    std::string additionalInfo; // Anything else
    std::string haplotype;
    
    // Optional helper methods (if needed later)
    bool isStructuralVariant() const {
        return varType == "DEL" || varType == "DUP" || varType == "INS" || varType == "INV" || varType == "TRA";
    }

    bool isSnvOrIndel() const {
        return varType == "SNV" || varType == "INDEL" || varType == "delins";
    }
};


struct AlignmentResult2
{
    int query_start;
    int query_end;
    int ref_start;
    int ref_end;
    std::string compactedCigar;
    std::string extendedCigar;
};

AlignmentResult2 alignReadToContig(const std::string &read, const std::string &contig)
{
    // Initialize a default scoring matrix for DNA sequences
    const int8_t match = 2;     // Match score
    const int8_t mismatch = 8;  // Mismatch penalty
    const int8_t gapOpen = 10;  // Gap open penalty
    const int8_t gapExtend = 1; // Gap extend penalty

    // Initialize the SSW aligner
    StripedSmithWaterman::Aligner aligner(match, mismatch, gapOpen, gapExtend);

    // Initialize the SSW filter (optional settings, can be adjusted as needed)
    StripedSmithWaterman::Filter filter;

    // Create an alignment object to store the result
    StripedSmithWaterman::Alignment alignment;

    // Perform the alignment
    aligner.Align(read.c_str(), contig.c_str(), contig.size(), filter, &alignment);

    // Extract the compacted CIGAR string
    std::string compactedCigar = alignment.cigar_string;

    // Generate the extended CIGAR string and adjust query_start and query_end to skip soft-clips
    std::string extendedCigar;
    int query_start = alignment.query_begin;
    int query_end = alignment.query_end;
    int ref_start = alignment.ref_begin;
    int ref_end = alignment.ref_end;

    int readPos = alignment.query_begin;
    int contigPos = alignment.ref_begin;

    // Adjust query_start by skipping leading soft-clips
    if (compactedCigar[0] == 'S')
    {
        size_t i = 0;
        int len = 0;
        while (isdigit(compactedCigar[i]))
        {
            len = len * 10 + (compactedCigar[i] - '0');
            ++i;
        }
        if (compactedCigar[i] == 'S')
        {
            query_start += len; // Skip the soft-clipped positions
        }
    }

    // Adjust query_end by skipping trailing soft-clips
    if (compactedCigar.back() == 'S')
    {
        size_t i = compactedCigar.size() - 1;
        int len = 0;
        while (isdigit(compactedCigar[i]))
        {
            len = len + (compactedCigar[i] - '0') * pow(10, compactedCigar.size() - i - 1);
            --i;
        }
        if (compactedCigar[i] == 'S')
        {
            query_end -= len; // Skip the soft-clipped positions
        }
    }

    for (size_t i = 0; i < compactedCigar.length(); ++i)
    {
        char op = compactedCigar[i];
        if (isdigit(op))
        {
            int len = 0;
            while (isdigit(compactedCigar[i]))
            {
                len = len * 10 + (compactedCigar[i] - '0');
                ++i;
            }
            op = compactedCigar[i];

            for (int j = 0; j < len; ++j)
            {
                extendedCigar += op;
                if (op == 'M')
                {
                    readPos++;
                    contigPos++;
                }
                else if (op == 'I')
                {
                    readPos++;
                }
                else if (op == 'D')
                {
                    contigPos++;
                }
            }
        }
    }

    // Create the AlignmentResult struct
    AlignmentResult2 result;
    result.query_start = query_start;
    result.query_end = query_end;
    result.ref_start = ref_start;
    result.ref_end = ref_end;
    result.compactedCigar = compactedCigar;
    result.extendedCigar = extendedCigar;

    return result;
}


std::string mutateWithQuality(const std::string& seq, const std::vector<float>& qualityProfile, std::mt19937& gen) {
    std::uniform_real_distribution<> dis(0.0, 1.0);
    std::string mutated = seq;

    for (size_t i = 0; i < seq.size(); ++i) {
        float q = (i < qualityProfile.size()) ? qualityProfile[i] : 30.0f;  // Default Phred if not known
        float errorProb = std::pow(10.0f, -q / 10.0f);  // Convert Phred to error probability

        if (dis(gen) < errorProb) {
            // Mutate this base
            char orig = seq[i];
            std::string bases = "ACGT";
            bases.erase(std::remove(bases.begin(), bases.end(), orig), bases.end());

            std::uniform_int_distribution<> base_dis(0, 2);
            mutated[i] = bases[base_dis(gen)];
        }
    }
    return mutated;
}


std::string str_toupper(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(),
                   // static_cast<int(*)(int)>(std::toupper)         // wrong
                   // [](int c){ return std::toupper(c); }           // wrong
                   // [](char c){ return std::toupper(c); }          // wrong
                   [](unsigned char c)
                   { return std::toupper(c); } // correct
    );
    return s;
}

struct FastqRead
{
    std::string seq;
    std::string qual;
    bool isFirstMate;
    bool isSecondMate;
    bool isReverseStrand;
    bool isMateReverseStrand;
};

struct BamProfile
{
    float meanInsertSize;
    float stdInsertSize;
    std::vector<float> qualityProfile;
};

BamProfile calculateBamProfile(const std::string &bamFile, int &maxReads) {

    // Read raw bam file
    BamReader reader(bamFile);
    BamRecord record;

    std::vector<int> insertSizes;
    std::vector<int> qualSum;
    std::vector<int> qualCount;

    insertSizes.reserve(100000);

    std::string chrM = "M";

    int longestRead = 0;
    int nReads = 0;

    while (reader.GetNextRecord(record))
    {
        // skipping mitochondrial DNA
        bool isFound = record.chrName().find(chrM) != std::string::npos;
        if (isFound)
        {
            continue;
        }
        if (!record.IsProperlyPaired())
        {
            continue;
        }
        if (record.IsUnmapped())
        {
            continue;
        }
        if (record.IsSecondary())
        {
            continue;
        }
        if (record.IsSupplementary())
        {
            continue;
        }

        nReads++;
        if (nReads > maxReads)
        {
            break;
        }

        int isize = std::abs(record.InsertSize());
        insertSizes.push_back(isize);

        const std::string &qual = record.Qual();
        int len = qual.length();
        if (len == 0)
        {
            continue;
        }

        if (len > static_cast<int>(qualSum.size()))
        {
            qualSum.resize(len, 0);
            qualCount.resize(len, 0);
        }

        for (int i = 0; i < len; ++i)
        {
            qualSum[i] += static_cast<int>(qual[i]) - 33;
            qualCount[i] += 1;
        }

        if (len > longestRead)
        {
            longestRead = len;
        }
    }
    BamProfile stats = {
        0.0f,
        0.0f,
        {}};

    // Mean
    double sum = 0.0;
    for (int val : insertSizes)
        sum += val;
    stats.meanInsertSize = static_cast<float>(sum / insertSizes.size());

    // Stddev
    double variance = 0.0;
    for (int val : insertSizes)
        variance += (val - stats.meanInsertSize) * (val - stats.meanInsertSize);
    stats.stdInsertSize = static_cast<float>(std::sqrt(variance / insertSizes.size()));

    // Quality profile
    stats.qualityProfile.resize(longestRead, 0.0f);
    for (size_t i = 0; i < stats.qualityProfile.size(); ++i)
    {
        if (qualCount[i] > 0)
        {
            stats.qualityProfile[i] = static_cast<float>(qualSum[i]) / qualCount[i];
        }
    }

    return stats;
}

std::string reverseComplement(const std::string &seq)
{
    std::string rc;
    rc.reserve(seq.size());
    for (auto it = seq.rbegin(); it != seq.rend(); ++it)
    {
        switch (std::toupper(*it))
        {
        case 'A':
            rc += 'T';
            break;
        case 'T':
            rc += 'A';
            break;
        case 'C':
            rc += 'G';
            break;
        case 'G':
            rc += 'C';
            break;
        case 'N':
            rc += 'N';
            break;
        default:
            rc += 'N';
            break;
        }
    }
    return rc;
}

void writeFastq(const std::string &qname,
                bool isFirstMate,
                const std::string &seq,
                const std::string &qual,
                bool isReverseStrand,
                std::ofstream &fq1,
                std::ofstream &fq2)
{
    std::string finalSeq = seq;
    std::string finalQual = qual;
    if (isReverseStrand)
    {
        finalSeq = reverseComplement(seq);
        finalQual = std::string(qual.rbegin(), qual.rend());
    }

    std::string name = "@" + qname;
    if (isFirstMate)
    {
        fq1 << name << "\n"
            << finalSeq << "\n+\n"
            << finalQual << "\n";
    }
    else
    {
        fq2 << name << "\n"
            << finalSeq << "\n+\n"
            << finalQual << "\n";
    }
}

void bufferOrWritePaired(const std::string &qname,
    bool isFirstMate, const std::string &seq, const std::string &qual, bool isReverseStrand,
                         std::ofstream &fq1,
                         std::ofstream &fq2,
                         std::unordered_map<std::string, FastqRead> &buffer)
{

    if (isFirstMate) {
        writeFastq(qname, true, seq, qual, isReverseStrand, fq1, fq2);
    }
    else {
        writeFastq(qname, false, seq, qual, isReverseStrand, fq1, fq2);
    }
}


std::pair<FastqRead, FastqRead> generateBreakpointPairedEnd(const Variant &variant, RefFasta &ref,
    float meanIsize, float stdIsize, int readLength = 101) {
    // Compute flanks around the breakpoint
    int flankSize = static_cast<int>(meanIsize);
    int leftStart = std::max(1L, variant.start - flankSize);
    int leftEnd = variant.start - 1;
    int rightStart = variant.end;
    int rightEnd = variant.end + flankSize;

    std::string leftFlank = ref.fetchSequence(variant.chr, leftStart, leftEnd);
    std::string rightFlank = ref.fetchSequence(variant.chr, rightStart, rightEnd);

    std::string contig;
    if (variant.varType == "DEL") {
        // For deletion, skip the deleted segment
        contig = leftFlank + rightFlank;
    }
    else if (variant.varType == "DUP") {
        // For duplication, include the duplicated segment twice
        std::string duplicatedSegment = ref.fetchSequence(variant.chr, variant.start, variant.end);
        contig = leftFlank + duplicatedSegment + duplicatedSegment + rightFlank;
    }
    else {
        throw std::runtime_error("Unsupported variant type for paired-end simulation: " + variant.varType);
    }

    // std::string contig = leftFlank + rightFlank;

    if (contig.length() < 2 * readLength) {
        throw std::runtime_error("Contig too short to simulate paired-end reads.");
    }

    // Sample insert size
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dist(meanIsize, stdIsize);
    int insertSize = std::clamp(static_cast<int>(dist(gen)), 2 * readLength, static_cast<int>(contig.size()));

    std::string fragment;
    if (variant.varType == "DEL") {
        // Sample fragment start
        std::uniform_int_distribution<> startDist(0, contig.length() - insertSize);
        int fragmentStart = startDist(gen);
        fragment = contig.substr(fragmentStart, insertSize);
    }
    if (variant.varType == "DUP") {
        // Sample fragment start
        std::uniform_int_distribution<> startDist((contig.length() / 2) - insertSize, (contig.length() / 2) + insertSize);
        int fragmentStart = startDist(gen);
        fragment = contig.substr(fragmentStart, insertSize);
    }
    // std::cout << " FRAGMENT:" << fragment << std::endl;

    // Generate read sequences
    std::string r1_seq = fragment.substr(0, readLength);
    std::string r2_seq = fragment.substr(fragment.length() - readLength);
    std::string r2_seq_rc = reverseComplement(r2_seq);

    // Simulate quality
    std::string r1_qual(readLength, 'I');
    std::string r2_qual(readLength, 'I');

    // Build FastqRead pair
    FastqRead r1{
        .seq = r1_seq,
        .qual = r1_qual,
        .isFirstMate = true,
        .isSecondMate = false,
        .isReverseStrand = false,
        .isMateReverseStrand = true};

    FastqRead r2{
        .seq = r2_seq,
        .qual = r2_qual,
        .isFirstMate = false,
        .isSecondMate = true,
        .isReverseStrand = true,
        .isMateReverseStrand = false};

    // FastqRead r2 {
    //     .seq = r2_seq_rc,
    //     .qual = std::string(r2_qual.rbegin(), r2_qual.rend()),
    //     .isFirstMate = false,
    //     .isSecondMate = true,
    //     .isReverseStrand = true,
    //     .isMateReverseStrand = false
    // };
    return std::make_pair(r1, r2);
}



std::pair<BamRecord, BamRecord> generatePairedEnd(BamRecord& record, int offset, RefFasta& ref, std::vector<SNV>& snvs) {

    std::string newQname = "SIMULATED_" + std::to_string(offset) + record.Qname();
    std::string modifiedSeq = record.Seq();
    bool modified = false;

    // === Apply SNVs to original read sequence ===
    for (const auto &snv : snvs)
    {
        if (record.Position() <= snv.position &&
            (record.Position() + modifiedSeq.length()) > snv.position)
        {
            int snvPosInRead = snv.position - record.Position();
            if (snvPosInRead >= 0 && snvPosInRead < static_cast<int>(modifiedSeq.length()))
            {
                modifiedSeq[snvPosInRead] = snv.alt[0];
                modified = true;
            }
        }
    }

    // === DUP1: Same as original read, modified if SNVs apply ===
    int insertSize = 0;

    BamRecord dup1 = record.DeepCopy();
    BamRecord dup2 = record.DeepCopy();

    dup1.SetQname(newQname);

    if (modified){
        dup1.UpdateSeq(modifiedSeq, record.cigarString(), ref);
    }

    int readLength = dup1.Seq().length();
    insertSize = abs(dup1.Position() - dup1.matePos())+offset;


    int dup2_pos = dup1.Position() + insertSize;
    dup2.SetPosition(dup2_pos);

    std::string dup2_seq = ref.fetchSequence(dup1.chrName(), dup2_pos + 1, dup2_pos + readLength);

    // convert the sequence to uppercase
    std::transform(dup2_seq.begin(), dup2_seq.end(), dup2_seq.begin(),
                    [](unsigned char c)
                    { return std::toupper(c); });

    std::string dup2_cigar = std::to_string(readLength) + "M";

    if (dup1.IsReverseStrand()) {
        newQname = "SIMULATED_R1_REV_" + std::to_string(offset) + record.Qname();
        dup2.SetQname(newQname);

        dup1.SetQname(newQname);
        int64_t mate_pos = dup1.Position() - insertSize;
        
        dup2.SetPosition(mate_pos);
        dup1.SetInsertSize(-insertSize);
        dup2.SetInsertSize(insertSize);

        std::string dup2_seq = ref.fetchSequence(dup1.chrName(), mate_pos + 1, mate_pos + dup1.Seq().length());

        dup2.UpdateSeq(dup2_seq, dup1.cigarString(), ref);
        dup2.SetFlag(BAM_FPAIRED | BAM_FPROPER_PAIR | BAM_FREAD2 | BAM_FMREVERSE);
    }
    else {
        newQname = "SIMULATED_R1_FWD_"+ std::to_string(offset) + record.Qname();

        dup2.SetQname(newQname);
        dup1.SetQname(newQname);
        dup2.SetPosition(dup2_pos);
        dup2.UpdateSeq(dup2_seq, dup2_cigar, ref);

        dup1.SetInsertSize(insertSize);
        dup2.SetInsertSize(-insertSize);
        dup2.SetFlag(BAM_FPAIRED | BAM_FPROPER_PAIR | BAM_FREAD2 | BAM_FREVERSE);
    }
    dup2.SetMatePosition(dup1.Position());
    dup1.SetMatePosition(dup2.Position());

    return {dup1, dup2};

}


bool simulateCNV(BamRecord &record, const std::string &region, const Variant &variant,
                 std::vector<SNV> &snvs, RefFasta &ref, BamWriter &writer,
                 const std::string &action, double probability, BamProfile bamStats,
                 std::unordered_map<std::string, FastqRead> &fqBuffer, std::ofstream &fq1, std::ofstream &fq2, 
                 std::unordered_map<std::string, std::vector<FastqRead>>& modifiedReads,
                 std::unordered_map<std::string, int>& skipReads) {

    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    double randomValue = dis(gen);

    if (variant.varType == "DEL") probability = 0.5;
    else if (variant.varType == "DUP") probability = 0.66;

    bool overlapsBreakpoint5prime = (record.Position() <= variant.start &&
                                     (record.Position() + record.Seq().length()) >= variant.start);

    bool overlapsBreakpoint3prime = (record.Position() <= variant.end &&
                                     (record.Position() + record.Seq().length()) >= variant.end);


    if (randomValue < probability) {
        if (overlapsBreakpoint5prime || overlapsBreakpoint3prime) {
            std::string newQname = record.Qname() + "_sim";

            // Generate synthetic paired-end read crossing the breakpoint
            auto [r1, r2] = generateBreakpointPairedEnd(variant, ref, bamStats.meanInsertSize, bamStats.stdInsertSize);

            writeFastq(newQname, true, r1.seq, r1.qual, r1.isReverseStrand, fq1, fq2);
            writeFastq(newQname, false, r2.seq, r2.qual, r2.isReverseStrand, fq1, fq2);
            return 1;
        }
    }


     // ---------- Handle DEL once per qname ----------
    if (variant.varType == "DEL") {
        std::string qname = record.Qname();
        // First time seeing this qname? Decide fate
        if (skipReads.count(qname) == 0) {
            double r = dis(gen);
            if (r < probability) {
                skipReads[qname] = 1;  // Mark to skip (deleted)
                return false;
            } 
            else {
                skipReads[qname] = -1; // Mark to keep
            }
        }

        // If marked to keep, apply SNVs and store modified read
        if (skipReads[qname] == -1) {
            std::string modifiedSeq = record.Seq();
            bool modified = false;

            for (const auto& snv : snvs) {
                if (record.Position() <= snv.position &&
                    (record.Position() + modifiedSeq.length()) > snv.position) {
                    int snvPosInRead = snv.position - record.Position();
                    if (snvPosInRead >= 0 && snvPosInRead < static_cast<int>(modifiedSeq.length())) {
                        modifiedSeq[snvPosInRead] = snv.alt[0];
                        modified = true;
                    }
                }
            }

            FastqRead fq;
            fq.seq = modified ? modifiedSeq : record.Seq();
            fq.qual = record.Qual();
            fq.isFirstMate = record.IsFirstMate();
            fq.isSecondMate = record.IsSecondMate();
            fq.isReverseStrand = record.IsReverseStrand();
            fq.isMateReverseStrand = record.IsMateReverseStrand();

            if (modified) {
                record.UpdateSeq(modifiedSeq, record.cigarString(), ref);
            }

            modifiedReads[qname].push_back(fq);
            return true;
        }

        // Otherwise: it's marked to be skipped
        return false;
    }

    if (randomValue < probability) {
        if (action == "DUP") {
            if (!record.IsFirstMate()) {
                return 0;
            }
            std::string newQname = "SIMULATED_" + record.Qname();
            writeFastq(newQname, true, record.Seq(), record.Qual(), record.IsReverseStrand(), fq1, fq2);

            int read2_pos = record.Position()+record.InsertSize();
            int readLength = record.Seq().length();

            std::string read2_seq = ref.fetchSequence(record.chrName(), read2_pos + 1, read2_pos + readLength);

            // convert the sequence to uppercase
            std::transform(read2_seq.begin(), read2_seq.end(), read2_seq.begin(),
                            [](unsigned char c)
                            { return std::toupper(c); });

            // 🔁 Introduce simulated sequencing errors in read2
            std::string mutated = mutateWithQuality(read2_seq, bamStats.qualityProfile, gen);

            bool isRev = record.IsReverseStrand() ? false : true;  // Flip strand for R2
            writeFastq(newQname, false, mutated, record.Qual(), isRev, fq1, fq2);
            // if (record.IsReverseStrand()) {
            //     writeFastq(newQname, false, read2_seq, record.Qual(), false, fq1, fq2);
            // }
            // else {
            //     writeFastq(newQname, false, read2_seq, record.Qual(), true, fq1, fq2);
            // }

            return true;
        }
    }

}

bool simulateIndel(const BamRecord &record, const Variant &indel, RefFasta &ref,
                   std::unordered_map<std::string, FastqRead> &fqBuffer) {

    const int64_t recordStart = record.Position();
    const int64_t recordEnd = recordStart + record.Seq().length();
    const int64_t indelStart = indel.start;
    const int64_t indelEnd = indel.start + indel.refSeq.length();

    std::string modifiedSeq = record.Seq();
    std::string originalQual = record.Qual();
    std::string modifiedQual = originalQual; // will fix length after edit

    int64_t overlapStart = std::max(recordStart, indelStart);
    int64_t overlapEnd = std::min(recordEnd, indelEnd);
    int64_t indelPosInRead = overlapStart - recordStart;
    int64_t indelLengthInRead = overlapEnd - overlapStart;

    if (indelStart >= recordStart && indelEnd <= recordEnd) {
        modifiedSeq.replace(indelPosInRead, indel.refSeq.length(), indel.altSeq);
    }
    else {
        std::string partialIndelSeq = indel.altSeq.substr(0, indelLengthInRead);
        modifiedSeq.replace(indelPosInRead, indelLengthInRead, partialIndelSeq);
    }

    // Adjust quality string to match new length
    if (modifiedSeq.length() < originalQual.length()) {
        modifiedQual = originalQual.substr(0, modifiedSeq.length());
    }
    else if (modifiedSeq.length() > originalQual.length()){
        // Extend using last quality score char
        char pad = originalQual.empty() ? 'I' : originalQual.back();
        modifiedQual += std::string(modifiedSeq.length() - originalQual.length(), pad);
    }

    std::random_device rd;                        // Obtain a random number from hardware
    std::mt19937 gen(rd());                       // Seed the generator
    std::uniform_int_distribution<> distr(0, 99); // Define the range

    int prob = indel.vaf * 100;

    if (distr(gen) < prob) {
        auto it = fqBuffer.find(record.Qname());
        if (it == fqBuffer.end())
        {
            fqBuffer[record.Qname()] = FastqRead{
                modifiedSeq,
                modifiedQual,
                record.IsFirstMate(),
                record.IsSecondMate(),
                record.IsReverseStrand(),
                record.IsMateReverseStrand()};
        }
        return true;
    }
    else
    {
        return false;
    }
}

void simulateSNV(const BamRecord &record, const Variant &variant, RefFasta &ref,
                 std::ofstream &fq1, std::ofstream &fq2,
                 std::unordered_map<std::string, FastqRead> &fqBuffer)
{
    if (variant.refSeq.length() == 1 && variant.altSeq.length() == 1)
    {
        int64_t variantPosInRead = variant.start - record.Position();
        if (variantPosInRead >= 0 && variantPosInRead < static_cast<int64_t>(record.Seq().length()))
        {
            std::string modifiedSeq = record.Seq();
            modifiedSeq[variantPosInRead] = variant.altSeq[0];
            bool isReverseStrand = record.IsReverseStrand();
            bufferOrWritePaired(record.Qname(), record.IsFirstMate(), modifiedSeq, record.Qual(), isReverseStrand, fq1, fq2, fqBuffer);
        }
    }
}


struct ReadPair {

    std::string readGroup;
    std::string read1Chr;
    int64_t read1Pos = -1;
    int64_t read1End = -1;
    std::string read1Seq;
    std::string read1Qual;

    std::string read2Chr;
    int64_t read2Pos = -1;
    int64_t read2End = -1;
    std::string read2Seq;
    std::string read2Qual;

    std::string tag;
    bool isRead1Reverse;
    bool isRead2Reverse;
    bool isSkipped = false;
};


void createReadPairs(const BamRecord &record, std::unordered_map<std::string, ReadPair>& ReadPairs) {
    if (record.IsSecondary() || record.IsSupplementary() || record.IsUnmapped())
        return;  // Skip non-primary or unmapped reads

    const std::string& qname = record.Qname();
    ReadPair& pair = ReadPairs[qname];  // Get or create the entry

    // Set read group only once (optional: implement record.ReadGroup() if needed)
    if (pair.readGroup.empty()) {
        pair.readGroup = record.Qname(); // or fetch from BAM aux if you implement it
    }

    if (record.IsFirstMate()) {
        pair.read1Chr = record.chrName();
        pair.read1Pos = record.Position();
        pair.read1End = record.End();
        pair.read1Seq = record.Seq();
        pair.read1Qual = record.Qual();
        pair.isRead1Reverse = record.IsReverseStrand();

        if (record.IsReverseStrand()) {
            pair.tag = "F2R1";
        }
        else{
            pair.tag = "F1R2";
        }
    }
    else if (record.IsSecondMate()) {
        pair.read2Chr = record.chrName();
        pair.read2Pos = record.Position();
        pair.read2End = record.End();
        pair.read2Seq = record.Seq();
        pair.read2Qual = record.Qual();
        pair.isRead2Reverse = record.IsReverseStrand();
        if (record.IsReverseStrand()) {
            pair.tag = "F1R2";
        }
        else{
            pair.tag = "F2R1";
        }
    }   
}


std::string classifyOverlap(const ReadPair &pair, const Variant& variant, int adjacency) {
    int variantStart, variantEnd;

    if (adjacency == 1) {
        variantStart = variant.start;
        variantEnd = variant.end;
    } else {
        variantStart = variant.start2;
        variantEnd = variant.end2;
    }

    std::string overlap_class = "None";

    // Standard left/right read assignment
    int64_t leftReadPos, leftReadEnd;
    int64_t rightReadPos, rightReadEnd;

    if (pair.read2Pos > pair.read1Pos) {
        leftReadPos = pair.read1Pos;
        leftReadEnd = pair.read1End;
        rightReadPos = pair.read2Pos;
        rightReadEnd = pair.read2End;
    } else {
        leftReadPos = pair.read2Pos;
        leftReadEnd = pair.read2End;
        rightReadPos = pair.read1Pos;
        rightReadEnd = pair.read1End;
    }

    // Quick sanity check
    if (leftReadPos < 0 || leftReadEnd < 0 || rightReadPos < 0 || rightReadEnd < 0) {
        return overlap_class;
    }

    // --- Overlap Types ---

    // Type 1: Right read overlaps 5' breakpoint, left is upstream
    if (leftReadEnd < variantStart &&
        rightReadPos < variantStart && rightReadEnd > variantStart) {
        overlap_class = "1";
    }

    // Type 2: Left read overlaps 3' breakpoint, right is downstream
    else if (leftReadPos <= variantEnd && leftReadEnd > variantEnd &&
             rightReadPos >= variantEnd) {
        overlap_class = "2";
    }

    // Type 3: Left read overlaps 5', right read inside SV
    else if (leftReadPos < variantStart && leftReadEnd > variantStart &&
             rightReadPos > variantStart && rightReadEnd < variantEnd) {
        overlap_class = "3";
    }

    // Type 4: Right read overlaps 3', left read inside SV
    else if (leftReadPos > variantStart && leftReadEnd < variantEnd &&
             rightReadPos < variantEnd && rightReadEnd > variantEnd) {
        overlap_class = "4";
    }

    // Type 5: Discordant reads — outside SV but no overlap
    else if (leftReadEnd <= variantStart &&
             rightReadPos >= variantStart && rightReadEnd <= variantEnd) {
        overlap_class = "5";
    }

    // Type 6: Discordant reads — other side
    else if (leftReadPos >= variantStart && leftReadEnd <= variantEnd &&
             rightReadPos >= variantEnd) {
        overlap_class = "6";
    }

    // Type 7: Both reads inside SV
    else if (leftReadPos >= variantStart && leftReadEnd <= variantEnd &&
             rightReadPos >= variantStart && rightReadEnd <= variantEnd) {
        overlap_class = "7";
    }

    // Type 8: Both reads span both breakpoints
    else if (leftReadPos < variantStart && leftReadEnd > variantStart &&
             rightReadPos < variantEnd && rightReadEnd > variantEnd) {
        overlap_class = "8";
    }

    // Type 9: Both reads upstream of SV
    else if (leftReadEnd < variantStart && rightReadEnd < variantStart) {
        overlap_class = "9";
    }

    // Type 10: Both reads downstream of SV
    else if (leftReadPos > variantEnd && rightReadPos > variantEnd) {
        overlap_class = "10";
    }

    // Type 11: Both reads overlap 5' breakpoint
    else if (leftReadPos <= variantStart && leftReadEnd >= variantStart &&
             rightReadPos <= variantStart && rightReadEnd >= variantStart) {
        overlap_class = "11";
    }

    // Type 12: Both reads overlap 3' breakpoint
    else if (leftReadPos <= variantEnd && leftReadEnd >= variantEnd &&
             rightReadPos <= variantEnd && rightReadEnd >= variantEnd) {
        overlap_class = "12";
    }

    // Type 13: Left read upstream of SV, right read downstream (spanning)
    else if (leftReadPos <= variantStart && leftReadEnd <= variantStart &&
             rightReadPos >= variantEnd && rightReadEnd >= variantEnd) {
        overlap_class = "13";
    }

    return overlap_class;
}

// Complement a DNA base (uppercase). Non-ACGT are returned unchanged.
inline char comp(char b) {
    switch (std::toupper(static_cast<unsigned char>(b))) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default:  return b;
    }
}
void applySnvsToRead(int64_t readPos,
                     int64_t readEnd,
                     std::string& readSeq,
                     bool isReverse,
                     const std::vector<VcfRecord>& snvs)
{
    for (const auto& snv : snvs) {
        const int64_t pos = snv.getPosition()-1;
        // overlap test (inclusive)
        if (pos < readPos || pos > readEnd) continue;

        // 0-based offset within the read sequence
        const int64_t offset = pos - readPos; // 1-based -> 0-based
        if (offset < 0 || static_cast<size_t>(offset) >= readSeq.size()) continue;

        // only handle strict SNVs
        if (snv.getRef().size() != 1 || snv.getAlt().size() != 1) continue;

        char alt = snv.getAlt()[0];
        // If read is reverse-strand, write the complement of ALT into the read string
        if (isReverse) alt = comp(alt);

        readSeq[static_cast<size_t>(offset)] = alt;
        // (Optionally adjust quality here, e.g., set to a high Q)
        // if (readQual.size() == readSeq.size()) readQual[offset] = 'I';
    }
}




ReadPair simulateDeletion( ReadPair& pair, const Variant& variant, const std::string& overlap_type, RefFasta& ref, std::vector<VcfRecord>& snvs) {

    double probability;
    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    double randomValue = dis(gen);
    
    ReadPair newPair;

    int64_t leftReadPos;
    int64_t leftReadEnd;
    std::string leftReadSeq;
    std::string leftReadQual;
    bool isLeftReverse;
    bool isLeftFirst;

    int64_t rightReadPos;
    int64_t rightReadEnd;
    std::string rightReadSeq;
    std::string rightReadQual;
    bool isRightReverse;
    bool isRightFirst;

    if (pair.read2Pos > pair.read1Pos) {
        leftReadPos = pair.read1Pos;
        leftReadEnd = pair.read1End;
        leftReadSeq = pair.read1Seq;
        leftReadQual = pair.read1Qual;

        isLeftReverse = pair.isRead1Reverse;
        isLeftFirst = true;

        rightReadPos = pair.read2Pos;
        rightReadEnd = pair.read2End;
        rightReadSeq = pair.read2Seq;
        rightReadQual = pair.read2Qual;
        isRightReverse = pair.isRead2Reverse;
        isRightFirst = false;
    }
    else {
        leftReadPos = pair.read2Pos;
        leftReadEnd = pair.read2End;
        leftReadSeq = pair.read2Seq;
        leftReadQual = pair.read2Qual;
        isLeftReverse = pair.isRead2Reverse;
        isLeftFirst = false;

        rightReadPos = pair.read1Pos;
        rightReadEnd = pair.read1End;
        rightReadSeq = pair.read1Seq;
        rightReadQual = pair.read1Qual;

        isRightReverse = pair.isRead1Reverse;
        isRightFirst = true;
    }

    
    // if leftRead or rightRead overlap with the SNV,
    // for (auto &snv : snvs) {
    //     if (snv.getPosition() <= leftReadPos && snv.getPosition() <= leftReadEnd) {
    //         // modify the leftReadSeq according with the offset of the snv, and ALT allele of the snv
    //     }
    //     if (snv.getPosition() <= rightReadPos && snv.getPosition() <= rightReadEnd) {
    //         // modify the leftReadSeq according with the offset of the snv, and ALT allele of the snv

    //     }

    // }
    applySnvsToRead(leftReadPos,  leftReadEnd,  leftReadSeq,  isLeftReverse,  snvs);
    applySnvsToRead(rightReadPos, rightReadEnd, rightReadSeq, isRightReverse, snvs);
    
    newPair.readGroup = "OLAP" + overlap_type + pair.tag + pair.readGroup;

    if (overlap_type == "9" || overlap_type == "10") {
        if (isLeftReverse) leftReadSeq = reverseComplement(leftReadSeq);
        if (isRightReverse) rightReadSeq = reverseComplement(rightReadSeq);

        //  Assign to correct struct fields
        if (isLeftFirst) {
            newPair.read1Seq = leftReadSeq;
            newPair.read1Qual = leftReadQual;
            newPair.read2Seq = rightReadSeq;
            newPair.read2Qual = rightReadQual;
        } else {
            newPair.read2Seq = leftReadSeq;
            newPair.read2Qual = leftReadQual;
            newPair.read1Seq = rightReadSeq;
            newPair.read1Qual = rightReadQual;
        }
        return newPair;
    }
    if (variant.varType == "DEL") {
        probability = 0.5;
    }
    if (randomValue > probability) {

        if (isLeftReverse) leftReadSeq = reverseComplement(leftReadSeq);
        if (isRightReverse) rightReadSeq = reverseComplement(rightReadSeq);

        //  Assign to correct struct fields
        if (isLeftFirst) {
            newPair.read1Seq = leftReadSeq;
            newPair.read1Qual = leftReadQual;
            newPair.read2Seq = rightReadSeq;
            newPair.read2Qual = rightReadQual;
        } else {
            newPair.read2Seq = leftReadSeq;
            newPair.read2Qual = leftReadQual;
            newPair.read1Seq = rightReadSeq;
            newPair.read1Qual = rightReadQual;
        }
        return newPair;
    }
        

    // Disallowed types for ALT allele (deletion-supporting)
    if (overlap_type == "3" || overlap_type == "4" || overlap_type == "7") {

        if (isLeftReverse) leftReadSeq = reverseComplement(leftReadSeq);
        if (isRightReverse) rightReadSeq = reverseComplement(rightReadSeq);

        //  Assign to correct struct fields
        if (isLeftFirst) {
            newPair.read1Seq = leftReadSeq;
            newPair.read1Qual = leftReadQual;
            newPair.read2Seq = rightReadSeq;
            newPair.read2Qual = rightReadQual;
        } 
        else {
            newPair.read2Seq = leftReadSeq;
            newPair.read2Qual = leftReadQual;
            newPair.read1Seq = rightReadSeq;
            newPair.read1Qual = rightReadQual;
        }
        newPair.isSkipped = true;
        return newPair;
    }


    if (overlap_type == "1") {

        std::string leftSeq = leftReadSeq;
        int nRightBases = rightReadEnd - variant.start;
        int matchBases = rightReadSeq.length() - nRightBases; 
        std::string rightSoftClip = ref.fetchSequence(variant.chr, variant.end + 1, variant.end+ nRightBases);
        std::string leftMatch = rightReadSeq.substr(0, matchBases);

        std::string rightSeq = leftMatch+rightSoftClip;

        //  Apply reverseComplement if needed
        if (isLeftReverse) leftSeq = reverseComplement(leftSeq);
        if (isRightReverse) rightSeq = reverseComplement(rightSeq);

        //  Assign to correct struct fields
        if (isLeftFirst) {
            newPair.read1Seq = leftSeq;
            newPair.read1Qual = leftReadQual;
            newPair.read2Seq = rightSeq;
            newPair.read2Qual = rightReadQual;
        } else {
            newPair.read2Seq = leftSeq;
            newPair.read2Qual = leftReadQual;
            newPair.read1Seq = rightSeq;
            newPair.read1Qual = rightReadQual;
        }
        return newPair;

    }
    if (overlap_type == "2") {
        std::string rightSeq = rightReadSeq;

        int nLeftBases = variant.end-leftReadPos;
        int matchBases = leftReadSeq.length()-nLeftBases;
        std::string leftSoftClip = ref.fetchSequence(variant.chr, variant.start-nLeftBases+1, variant.start);
        std::string rightMatch = leftReadSeq.substr(nLeftBases);
        std::string leftSeq = leftSoftClip+rightMatch;

        //  Apply reverseComplement if needed
        if (isLeftReverse) leftSeq = reverseComplement(leftSeq);
        if (isRightReverse) rightSeq = reverseComplement(rightSeq);

        //  Assign to correct struct fields
        if (isLeftFirst) {
            newPair.read1Seq = leftSeq;
            newPair.read1Qual = leftReadQual;
            newPair.read2Seq = rightSeq;
            newPair.read2Qual = rightReadQual;
        } else {
            newPair.read2Seq = leftSeq;
            newPair.read2Qual = leftReadQual;
            newPair.read1Seq = rightSeq;
            newPair.read1Qual = rightReadQual;
        }
        return newPair;

    }

    if (overlap_type == "5") {
        std::string leftSeq = leftReadSeq;
        int diff = rightReadPos-variant.start;
        std::string rightSeq = ref.fetchSequence(variant.chr, variant.end+diff+1, variant.end+diff+rightReadSeq.length());
        //  Apply reverseComplement if needed
        if (isLeftReverse) leftSeq = reverseComplement(leftSeq);
        if (isRightReverse) rightSeq = reverseComplement(rightSeq);

        //  Assign to correct struct fields
        if (isLeftFirst) {
            newPair.read1Seq = leftSeq;
            newPair.read1Qual = leftReadQual;
            newPair.read2Seq = rightSeq;
            newPair.read2Qual = rightReadQual;
        } else {
            newPair.read2Seq = leftSeq;
            newPair.read2Qual = leftReadQual;
            newPair.read1Seq = rightSeq;
            newPair.read1Qual = rightReadQual;
        }
        return newPair;

    }
    if (overlap_type == "6") {
        std::string rightSeq = rightReadSeq;
        int diff = variant.end-leftReadEnd;
        std::string leftSeq = ref.fetchSequence(variant.chr, variant.start-diff-leftReadSeq.length()+1, variant.start-diff);

        if (isLeftReverse) leftSeq = reverseComplement(leftSeq);
        if (isRightReverse) rightSeq = reverseComplement(rightSeq);

        //  Assign to correct struct fields
        if (isLeftFirst) {
            newPair.read1Seq = leftSeq;
            newPair.read1Qual = leftReadQual;
            newPair.read2Seq = rightSeq;
            newPair.read2Qual = rightReadQual;
        } 
        else {
            newPair.read2Seq = leftSeq;
            newPair.read2Qual = leftReadQual;
            newPair.read1Seq = rightSeq;
            newPair.read1Qual = rightReadQual;
        }
        return newPair;
    }

    if (overlap_type == "8") {

        int nRightBases = leftReadEnd-variant.start;
        int matchBases = variant.start-nRightBases;

        std::string leftMatch = leftReadSeq.substr(0, matchBases);
        std::string rightSoftClip = ref.fetchSequence(variant.chr, variant.end+1, variant.end+nRightBases);
        std::string leftSeq = leftMatch+rightSoftClip;

        matchBases = rightReadEnd-variant.end;
        int nLeftBases = rightReadSeq.length()-matchBases;

        std::string leftSoftClip = ref.fetchSequence(variant.chr, variant.start-nLeftBases+1, variant.start);
        std::string rightMatch = rightReadSeq.substr(nLeftBases);
        std::string rightSeq = leftSoftClip+rightMatch;

        if (isLeftReverse) leftSeq = reverseComplement(leftSeq);
        if (isRightReverse) rightSeq = reverseComplement(rightSeq);

        //  Assign to correct struct fields
        if (isLeftFirst) {
            newPair.read1Seq = leftSeq;
            newPair.read1Qual = leftReadQual;
            newPair.read2Seq = rightSeq;
            newPair.read2Qual = rightReadQual;
        } 
        else {
            newPair.read2Seq = leftSeq;
            newPair.read2Qual = leftReadQual;
            newPair.read1Seq = rightSeq;
            newPair.read1Qual = rightReadQual;
        }
        return newPair;
    }
    if (overlap_type == "11") {
        int nRightBases = leftReadEnd-variant.start;
        int matchBases = leftReadSeq.length()-nRightBases;

        std::string leftMatch = leftReadSeq.substr(0, matchBases);
        std::string rightSoftClip = ref.fetchSequence(variant.chr, variant.end+1, variant.end+nRightBases);
        std::string leftSeq = leftMatch+rightSoftClip;


        nRightBases = rightReadEnd-variant.start;
        matchBases = rightReadSeq.length()-nRightBases;

        leftMatch = rightReadSeq.substr(0, matchBases); 
        rightSoftClip = ref.fetchSequence(variant.chr, variant.end+1, variant.end+nRightBases);
        std::string rightSeq = leftMatch+rightSoftClip;


        if (isLeftReverse) leftSeq = reverseComplement(leftSeq);
        if (isRightReverse) rightSeq = reverseComplement(rightSeq);

        //  Assign to correct struct fields
        if (isLeftFirst) {
            newPair.read1Seq = leftSeq;
            newPair.read1Qual = leftReadQual;
            newPair.read2Seq = rightSeq;
            newPair.read2Qual = rightReadQual;
        } 
        else {
            newPair.read2Seq = leftSeq;
            newPair.read2Qual = leftReadQual;
            newPair.read1Seq = rightSeq;
            newPair.read1Qual = rightReadQual;
        }
        return newPair;

    }

    if (overlap_type == "12") {

        int matchBases = leftReadEnd-variant.end;
        int nLeftBases = leftReadSeq.length()-matchBases;

        std::string leftSoftClip = ref.fetchSequence(variant.chr, variant.start-nLeftBases+1, variant.start);
        std::string rightMatch = leftReadSeq.substr(nLeftBases);
        std::string leftSeq = leftSoftClip+rightMatch;

        matchBases = rightReadEnd-variant.end;
        nLeftBases =rightReadSeq.length()-matchBases;

        leftSoftClip = ref.fetchSequence(variant.chr, variant.start-nLeftBases+1, variant.start);
        rightMatch = rightReadSeq.substr(nLeftBases);
        std::string rightSeq = leftSoftClip+rightMatch;
        
        if (isLeftReverse) leftSeq = reverseComplement(leftSeq);
        if (isRightReverse) rightSeq = reverseComplement(rightSeq);

        //  Assign to correct struct fields
        if (isLeftFirst) {
            newPair.read1Seq = leftSeq;
            newPair.read1Qual = leftReadQual;
            newPair.read2Seq = rightSeq;
            newPair.read2Qual = rightReadQual;
        } 
        else {
            newPair.read2Seq = leftSeq;
            newPair.read2Qual = leftReadQual;
            newPair.read1Seq = rightSeq;
            newPair.read1Qual = rightReadQual;
        }
        return newPair;
    }
  
    return newPair;
}


std::vector<ReadPair> simulateTandemDuplication( ReadPair& pair, const Variant& variant, const std::string& overlap_type, RefFasta& ref) {
    
    std::vector<ReadPair> simulatedReads;

    double probability;
    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    double randomValue = dis(gen);

    int64_t leftReadPos;
    int64_t leftReadEnd;
    std::string leftReadSeq;
    std::string leftReadQual;
    bool isLeftReverse;
    bool isLeftFirst;

    int64_t rightReadPos;
    int64_t rightReadEnd;
    std::string rightReadSeq;
    std::string rightReadQual;
    bool isRightReverse;
    bool isRightFirst;

    if (pair.read2Pos > pair.read1Pos) {
        leftReadPos = pair.read1Pos;
        leftReadEnd = pair.read1End;
        leftReadSeq = pair.read1Seq;
        leftReadQual = pair.read1Qual;

        isLeftReverse = pair.isRead1Reverse;
        isLeftFirst = true;

        rightReadPos = pair.read2Pos;
        rightReadEnd = pair.read2End;
        rightReadSeq = pair.read2Seq;
        rightReadQual = pair.read2Qual;
        isRightReverse = pair.isRead2Reverse;
        isRightFirst = false;
    }
    else {
        leftReadPos = pair.read2Pos;
        leftReadEnd = pair.read2End;
        leftReadSeq = pair.read2Seq;
        leftReadQual = pair.read2Qual;
        isLeftReverse = pair.isRead2Reverse;
        isLeftFirst = false;

        rightReadPos = pair.read1Pos;
        rightReadEnd = pair.read1End;
        rightReadSeq = pair.read1Seq;
        rightReadQual = pair.read1Qual;

        isRightReverse = pair.isRead1Reverse;
        isRightFirst = true;
    }
    ReadPair newPair;
    newPair.readGroup = pair.readGroup;


    if (variant.varType == "DUP") {
        probability = 0.5;
    }
    probability = 1;
    if (randomValue > probability) {

        if (isLeftReverse) leftReadSeq = reverseComplement(leftReadSeq);
        if (isRightReverse) rightReadSeq = reverseComplement(rightReadSeq);

        //  Assign to correct struct fields
        if (isLeftFirst) {
            pair.read1Seq = leftReadSeq;
            pair.read1Qual = leftReadQual;
            pair.read2Seq = rightReadSeq;
            pair.read2Qual = rightReadQual;
        } else {
            pair.read2Seq = leftReadSeq;
            pair.read2Qual = leftReadQual;
            pair.read1Seq = rightReadSeq;
            pair.read1Qual = rightReadQual;
        }
        simulatedReads.push_back(pair);
        return simulatedReads;
    }

    if (overlap_type == "3") {
        int nLeftBases = variant.start-leftReadPos;
        int matchBases = leftReadEnd-variant.start;

        std::string leftSoftClip = ref.fetchSequence(variant.chr, variant.end-nLeftBases+1, variant.end);
        std::string rightMatch = leftReadSeq.substr(nLeftBases);
        std::string leftSeq = leftSoftClip+rightMatch;
        std::string rightSeq = rightReadSeq;

        if (isLeftReverse) leftSeq = reverseComplement(leftSeq);
        if (isRightReverse) rightSeq = reverseComplement(rightSeq);

        //  Assign to correct struct fields

        if (isLeftFirst) {
            newPair.read1Seq = rightSeq;
            newPair.read1Qual = rightReadQual;
            newPair.read2Seq = leftSeq;
            newPair.read2Qual = leftReadQual;
        } 
        else {
            newPair.read2Seq = rightSeq;
            newPair.read2Qual = rightReadQual;
            newPair.read1Seq = leftSeq;
            newPair.read1Qual = leftReadQual;
        }
        simulatedReads.push_back(newPair);
    }

    else if (overlap_type == "4") {
        int nRightBases = rightReadEnd-variant.end;
        int matchBases = variant.end-rightReadPos;

        std::string rightSoftClip = ref.fetchSequence(variant.chr, variant.start+1, variant.start+nRightBases);
        std::string leftMatch = rightReadSeq.substr(0,matchBases);
        std::string rightSeq = leftMatch+rightSoftClip;
        std::string leftSeq = leftReadSeq;
        if (isLeftReverse) leftSeq = reverseComplement(leftSeq);
        if (isRightReverse) rightSeq = reverseComplement(rightSeq);

        //  Assign to correct struct fields
        if (isLeftFirst) {
            newPair.read1Seq = rightSeq;
            newPair.read1Qual = rightReadQual;
            newPair.read2Seq = leftSeq;
            newPair.read2Qual = leftReadQual;
        } 
        else {
            newPair.read2Seq = rightSeq;
            newPair.read2Qual = rightReadQual;
            newPair.read1Seq = leftSeq;
            newPair.read1Qual = leftReadQual;
        }
        simulatedReads.push_back(newPair);
    }
    else if (overlap_type == "5") {

        int diff = variant.start-leftReadEnd;
        std::string leftSeq = ref.fetchSequence(variant.chr, variant.end-diff-leftReadSeq.length()+1, variant.end-diff);
        std::string rightSeq = rightReadSeq;
        
        if (isLeftReverse) leftSeq = reverseComplement(leftSeq);
        if (isRightReverse) rightSeq = reverseComplement(rightSeq);

        if (isLeftFirst) {
            newPair.read1Seq = rightSeq;
            newPair.read1Qual = rightReadQual;
            newPair.read2Seq = leftSeq;
            newPair.read2Qual = leftReadQual;
        } 
        else {
            newPair.read2Seq = rightSeq;
            newPair.read2Qual = rightReadQual;
            newPair.read1Seq = leftSeq;
            newPair.read1Qual = leftReadQual;
        }

        simulatedReads.push_back(newPair);
    }
    else if (overlap_type == "6") {
        int diff = rightReadPos-variant.end;

        std::string rightSeq = ref.fetchSequence(variant.chr, variant.start+diff+1, variant.start+diff+rightReadSeq.length());
        std::string leftSeq = leftReadSeq;

        if (isLeftReverse) leftSeq = reverseComplement(leftSeq);
        if (isRightReverse) rightSeq = reverseComplement(rightSeq);
        if (isLeftFirst) {
            newPair.read1Seq = rightSeq;
            newPair.read1Qual = rightReadQual;
            newPair.read2Seq = leftSeq;
            newPair.read2Qual = leftReadQual;
        } 
        else {
            newPair.read2Seq = rightSeq;
            newPair.read2Qual = rightReadQual;
            newPair.read1Seq = leftSeq;
            newPair.read1Qual = leftReadQual;
        }
        simulatedReads.push_back(newPair);
    }

    else if (overlap_type == "7") {
        ReadPair newPair2 = pair;
        newPair2.readGroup = "SIMULATED"+ pair.readGroup;

        int64_t newRead1Pos = leftReadPos+1;
        int64_t newRead1End = leftReadEnd+1;
        std::string newRead1Seq = ref.fetchSequence(variant.chr, newRead1Pos+1, newRead1End);
        newPair2.read1Pos = newRead1Pos;
        newPair2.read1End = newRead1End;

        int64_t newRead2Pos = rightReadPos+1;
        int64_t newRead2End = rightReadEnd+1;
        std::string newRead2Seq = ref.fetchSequence(variant.chr, newRead2Pos+1, newRead2End);

        newPair2.read2Pos = newRead2Pos;
        newPair2.read2End = newRead2End;

        if (isLeftReverse) leftReadSeq = reverseComplement(leftReadSeq);
        if (isRightReverse) rightReadSeq = reverseComplement(rightReadSeq);


        if (isLeftReverse) newRead1Seq = reverseComplement(newRead1Seq);
        if (isRightReverse) newRead2Seq = reverseComplement(newRead2Seq);

        if (isLeftFirst) {
            newPair.read1Seq = leftReadSeq;
            newPair.read1Qual = leftReadQual;
            newPair.read2Seq = rightReadSeq;
            newPair.read2Qual = rightReadQual;

            newPair2.read1Seq = newRead1Seq;
            newPair2.read1Qual = leftReadQual;
            newPair2.read2Seq = newRead2Seq;
            newPair2.read2Qual = rightReadQual;
        } 
        else {
            newPair.read2Seq = leftReadSeq;
            newPair.read2Qual = leftReadQual;
            newPair.read1Seq = rightReadSeq;
            newPair.read1Qual = rightReadQual;

            newPair2.read2Seq = newRead1Seq;
            newPair2.read2Qual = leftReadQual;
            newPair2.read1Seq = newRead2Seq;
            newPair2.read1Qual = rightReadQual;
        }
        simulatedReads.push_back(newPair);
        simulatedReads.push_back(newPair2);
    }
    else if (overlap_type == "8") {
        // if (pair.tag == "F1R2") {
        int nLeftBases = variant.start-leftReadPos;
        int matchBases = leftReadEnd-variant.start;

        std::cout << variant.chr << ":" << variant.start << "-" << variant.end << std::endl;
        std::cout << leftReadPos << "-" << leftReadEnd << "\t" << rightReadPos << "-" << rightReadEnd << std::endl;

        std::string leftSoftClip = ref.fetchSequence(variant.chr, variant.end-nLeftBases+1, variant.end);
        std::string rightMatch = pair.read1Seq.substr(nLeftBases);
        std::string leftSeq = leftSoftClip+rightMatch;
        // std::cout << "LeftRead:" << " nLeftBases:" << nLeftBases << " matchBases:" << matchBases << " " << leftSeq << std::endl;

        int nRightBases = rightReadEnd-variant.end;
        matchBases = variant.end-rightReadPos;

        std::string rightSoftClip = ref.fetchSequence(variant.chr, variant.start+1,  variant.start+nRightBases);
        std::string leftMatch = rightReadSeq.substr(0, matchBases);
        std::string rightSeq = leftMatch+rightSoftClip;
        // std::cout << "RightRead:" << " matchBases:" << matchBases << " nRightBases:" << nRightBases  << " " << rightSeq << std::endl << std::endl;
        if (isLeftReverse) leftSeq = reverseComplement(leftSeq);
        if (isRightReverse) rightSeq = reverseComplement(rightSeq);

        if (isLeftFirst) {
            newPair.read1Seq = rightSeq;
            newPair.read1Qual = rightReadQual;
            newPair.read2Seq = leftSeq;
            newPair.read2Qual = leftReadQual;
        } 
        else {
            newPair.read2Seq = rightSeq;
            newPair.read2Qual = rightReadQual;
            newPair.read1Seq = leftSeq;
            newPair.read1Qual = leftReadQual;
        }
        simulatedReads.push_back(newPair);
    }

    else if (overlap_type == "9" || overlap_type == "11") {

        std::string leftSeq = leftReadSeq;
        std::string rightSeq = rightReadSeq;

        if (isLeftReverse) leftSeq = reverseComplement(leftSeq);
        if (isRightReverse) rightSeq = reverseComplement(rightSeq);

        //  Assign to correct struct fields
        if (isLeftFirst) {
            newPair.read1Seq = leftSeq;
            newPair.read1Qual = leftReadQual;
            newPair.read2Seq = rightSeq;
            newPair.read2Qual = rightReadQual;
        } 
        else {
            newPair.read2Seq = leftSeq;
            newPair.read2Qual = leftReadQual;
            newPair.read1Seq = rightSeq;
            newPair.read1Qual = rightReadQual;
        }
        simulatedReads.push_back(newPair);

    }
    return simulatedReads;
}


ReadPair simulateInversion(ReadPair& pair, const Variant& variant, const std::string& overlap_type, RefFasta& ref) {

    // std::cout << overlap_type << std::endl;
    // ReadPair newPair = pair;
    ReadPair newPair;


    int64_t leftReadPos;
    int64_t leftReadEnd;
    std::string leftReadSeq;
    std::string leftReadQual;
    bool isLeftReverse;
    bool isLeftFirst;

    int64_t rightReadPos;
    int64_t rightReadEnd;
    std::string rightReadSeq;
    std::string rightReadQual;
    bool isRightReverse;
    bool isRightFirst;

    if (pair.read2Pos > pair.read1Pos) {
        leftReadPos = pair.read1Pos;
        leftReadEnd = pair.read1End;
        leftReadSeq = pair.read1Seq;
        leftReadQual = pair.read1Qual;

        isLeftReverse = pair.isRead1Reverse;
        isLeftFirst = true;

        rightReadPos = pair.read2Pos;
        rightReadEnd = pair.read2End;
        rightReadSeq = pair.read2Seq;
        rightReadQual = pair.read2Qual;
        isRightReverse = pair.isRead2Reverse;
        isRightFirst = false;
    }
    else {
        leftReadPos = pair.read2Pos;
        leftReadEnd = pair.read2End;
        leftReadSeq = pair.read2Seq;
        leftReadQual = pair.read2Qual;
        isLeftReverse = pair.isRead2Reverse;
        isLeftFirst = false;

        rightReadPos = pair.read1Pos;
        rightReadEnd = pair.read1End;
        rightReadSeq = pair.read1Seq;
        rightReadQual = pair.read1Qual;

        isRightReverse = pair.isRead1Reverse;
        isRightFirst = true;
    }
    
    newPair.readGroup = "OLAP" + overlap_type + pair.tag + pair.readGroup;

    if (overlap_type == "1") {
        int nRightBases = rightReadEnd - variant.start;
        int matchBases = rightReadSeq.length() - nRightBases;
        std::string leftMatch = rightReadSeq.substr(0, matchBases);

        std::string rightSoftClip = ref.fetchSequence(variant.chr, variant.end - nRightBases + 1, variant.end);
        rightSoftClip = reverseComplement(rightSoftClip);

        std::string leftSeq = leftReadSeq;
        std::string rightSeq = leftMatch + rightSoftClip;

        //  Apply reverseComplement if needed
        if (isLeftReverse) leftSeq = reverseComplement(leftSeq);
        if (isRightReverse) rightSeq = reverseComplement(rightSeq);

        //  Assign to correct struct fields
        if (isLeftFirst) {
            newPair.read1Seq = leftSeq;
            newPair.read1Qual = leftReadQual;
            newPair.read2Seq = rightSeq;
            newPair.read2Qual = rightReadQual;
        } else {
            newPair.read2Seq = leftSeq;
            newPair.read2Qual = leftReadQual;
            newPair.read1Seq = rightSeq;
            newPair.read1Qual = rightReadQual;
        }
    }


    else if (overlap_type == "2") {
        int matchBases = leftReadEnd - variant.end;
        int nLeftBases = leftReadSeq.length() - matchBases;

        std::string leftSoftClip = ref.fetchSequence(variant.chr, variant.start + 1, variant.start + nLeftBases);
        leftSoftClip = reverseComplement(leftSoftClip);

        std::string rightMatch = leftReadSeq.substr(nLeftBases);
        std::string leftSeq = leftSoftClip + rightMatch;

        std::string rightSeq = rightReadSeq;

        // Apply strand-based corrections
        if (isLeftReverse) leftSeq = reverseComplement(leftSeq);
        if (isRightReverse) rightSeq = reverseComplement(rightSeq);

        if (isLeftFirst) {
            newPair.read1Seq = leftSeq;
            newPair.read1Qual = leftReadQual;
            newPair.read2Seq = rightSeq;
            newPair.read2Qual = rightReadQual;
        } else {
            newPair.read2Seq = leftSeq;
            newPair.read2Qual = leftReadQual;
            newPair.read1Seq = rightSeq;
            newPair.read1Qual = rightReadQual;
        }
    }

   
    else if (overlap_type == "3") {
        int nRightBases = leftReadEnd-variant.start;
        int matchBases = leftReadSeq.length()-nRightBases;

        std::string rightSoftClip = ref.fetchSequence(variant.chr, variant.end-nRightBases+1, variant.end);
        rightSoftClip = reverseComplement(rightSoftClip);
        std::string leftMatch = leftReadSeq.substr(0, matchBases);

        int diff = abs(rightReadPos-leftReadEnd);
        std::string rightSeq = ref.fetchSequence(variant.chr, variant.end-nRightBases-diff-rightReadSeq.length()+1, variant.end-nRightBases-diff);
        rightSeq = reverseComplement(rightSeq);

        std::string leftSeq = leftMatch+rightSoftClip;

        // Apply strand-based corrections
        if (isLeftReverse) leftSeq = reverseComplement(leftSeq);
        if (isRightReverse) rightSeq = reverseComplement(rightSeq);

        if (isLeftFirst) {
            newPair.read1Seq = leftSeq;
            newPair.read1Qual = leftReadQual;
            newPair.read2Seq = rightSeq;
            newPair.read2Qual = rightReadQual;
        }
        else {
            newPair.read2Seq = leftSeq;
            newPair.read2Qual = leftReadQual;
            newPair.read1Seq = rightSeq;
            newPair.read1Qual = rightReadQual;
        }
    }
    else if (overlap_type == "4") {

        int matchBases = rightReadEnd-variant.end;
        int nLeftBases = rightReadSeq.length()-matchBases;

        std::string leftSoftClip = ref.fetchSequence(variant.chr, variant.start+1, variant.start+nLeftBases);
        leftSoftClip = reverseComplement(leftSoftClip);

        std::string rightMatch = rightReadSeq.substr(nLeftBases);

        int diff = abs(rightReadPos-leftReadEnd);
        std::string leftSeq = ref.fetchSequence(variant.chr, variant.start+nLeftBases+diff+1, variant.start+nLeftBases+diff+leftReadSeq.length());
        std::string rightSeq = leftSoftClip+rightMatch;

        // std::string rightSeq = reverseComplement(leftSoftClip+rightMatch);
        leftSeq = reverseComplement(leftSeq);


        // Apply strand-based corrections
        if (isLeftReverse) leftSeq = reverseComplement(leftSeq);
        if (isRightReverse) rightSeq = reverseComplement(rightSeq);

        if (isLeftFirst) {
            newPair.read1Seq = leftSeq;
            newPair.read1Qual = leftReadQual;

            newPair.read2Seq = rightSeq;
            newPair.read2Qual = rightReadQual;
        }
        else {
            newPair.read2Seq = leftSeq;
            newPair.read2Qual = leftReadQual;
            newPair.read1Seq = rightSeq;
            newPair.read1Qual = rightReadQual;
        }
        
    }
    else if (overlap_type == "5") {
            // newPair.readGroup = "SIMULATED5F1R2" + newPair.readGroup;          
            int diff = abs(rightReadEnd-variant.start);
            std::string rightSeq = ref.fetchSequence(variant.chr, variant.end-diff-rightReadSeq.length()+1, variant.end-diff);
            std::string leftSeq = leftReadSeq;
            rightSeq = reverseComplement(rightSeq);

            // Apply strand-based corrections
            if (isLeftReverse) leftSeq = reverseComplement(leftSeq);
            if (isRightReverse) rightSeq = reverseComplement(rightSeq);

            if (isLeftFirst) {
                newPair.read1Seq = leftSeq;
                newPair.read1Qual = leftReadQual;

                newPair.read2Seq = rightSeq;
                newPair.read2Qual = rightReadQual;
            }
            else {
                newPair.read2Seq = leftSeq;
                newPair.read2Qual = leftReadQual;
                newPair.read1Seq = rightSeq;
                newPair.read1Qual = rightReadQual;
            }

    }
    else if (overlap_type == "6") {
            int diff = abs(variant.end-leftReadEnd);
            std::string leftSeq = ref.fetchSequence(variant.chr, variant.start+diff+1, variant.start+diff+leftReadSeq.length());
            leftSeq = reverseComplement(leftSeq);
            std::string rightSeq = rightReadSeq;

            // Apply strand-based corrections
            if (isLeftReverse) leftSeq = reverseComplement(leftSeq);
            if (isRightReverse) rightSeq = reverseComplement(rightSeq);

            if (isLeftFirst) {
                newPair.read1Seq = leftSeq;
                newPair.read1Qual = leftReadQual;

                newPair.read2Seq = rightSeq;
                newPair.read2Qual = rightReadQual;
            }
            else {
                newPair.read2Seq = leftSeq;
                newPair.read2Qual = leftReadQual;
                newPair.read1Seq = rightSeq;
                newPair.read1Qual = rightReadQual;
            } 

    }
    else if (overlap_type == "11") {

        int matchBases = variant.start-leftReadPos;
        int nRightBases = leftReadEnd-variant.start;

        std::string firstReadClipped = ref.fetchSequence(variant.chr, variant.end-nRightBases+1, variant.end);
        std::string firstReadMatched = leftReadSeq.substr(0, matchBases);

        std::string leftSeq = firstReadMatched+reverseComplement(firstReadClipped);

        // Now for read2
        matchBases = variant.start-rightReadPos;
        nRightBases = rightReadEnd-variant.start;

        std::string secondReadClipped = ref.fetchSequence(variant.chr, variant.end-nRightBases+1, variant.end);
        std::string secondReadMatched = rightReadSeq.substr(0, matchBases);
        std::string rightSeq =  secondReadMatched+reverseComplement(secondReadClipped);


        //  Apply reverseComplement if needed
        if (isLeftReverse) leftSeq = reverseComplement(leftSeq);
        if (isRightReverse) rightSeq = reverseComplement(rightSeq);

        //  Assign to correct struct fields
        if (isLeftFirst) {
            newPair.read1Seq = leftSeq;
            newPair.read1Qual = leftReadQual;
            newPair.read2Seq = rightSeq;
            newPair.read2Qual = rightReadQual;
        } else {
            newPair.read2Seq = leftSeq;
            newPair.read2Qual = leftReadQual;
            newPair.read1Seq = rightSeq;
            newPair.read1Qual = rightReadQual;
        }
    }
    else if (overlap_type == "12") {
    

        int nLeftBases = variant.end-leftReadPos;
        int matchBases = leftReadEnd-variant.end;

        std::string firstReadClipped = ref.fetchSequence(variant.chr, variant.start+1, variant.start+nLeftBases);
        std::string firstReadMatched = leftReadSeq.substr(nLeftBases);
        std::string leftSeq = reverseComplement(firstReadClipped) + firstReadMatched;

        nLeftBases = variant.end-rightReadPos;
        matchBases = rightReadEnd-variant.end;

        std::string secondReadClipped = ref.fetchSequence(variant.chr, variant.start+1, variant.start+nLeftBases);
        std::string secondReadMatched = rightReadSeq.substr(nLeftBases);
        
        std::string rightSeq = reverseComplement(secondReadClipped) + secondReadMatched;
        // Apply strand-based corrections
        if (isLeftReverse) leftSeq = reverseComplement(leftSeq);
        if (isRightReverse) rightSeq = reverseComplement(rightSeq);

        if (isLeftFirst) {
            newPair.read1Seq = leftSeq;
            newPair.read1Qual = leftReadQual;

            newPair.read2Seq = rightSeq;
            newPair.read2Qual = rightReadQual;
        }
        else {
            newPair.read2Seq = leftSeq;
            newPair.read2Qual = leftReadQual;
            newPair.read1Seq = rightSeq;
            newPair.read1Qual = rightReadQual;
        } 
    }
    return newPair;
}


ReadPair simulateBalancedTranslocation(ReadPair& pair, const Variant& variant, const std::string& overlap_type, RefFasta& ref, int adjacency) {

    int64_t leftReadPos;
    int64_t leftReadEnd;
    std::string leftReadSeq;
    std::string leftReadQual;
    bool isLeftReverse;
    bool isLeftFirst;

    int64_t rightReadPos;
    int64_t rightReadEnd;
    std::string rightReadSeq;
    std::string rightReadQual;
    bool isRightReverse;
    bool isRightFirst;

    std::string chr2;
    int64_t variant2Start;
    int64_t variant2End;

    std::string chr1;
    int64_t variant1Start;
    int64_t variant1End;

    if (adjacency == 1) {
        chr1 = variant.chr;
        variant1Start = variant.start;
        variant1End = variant.end;

        chr2 = variant.chr2;
        variant2Start = variant.start2;
        variant2End = variant.end2;
    }
    if (adjacency == 2) {
        chr1 = variant.chr2;
        variant1Start = variant.start2;
        variant1End = variant.end2;

        chr2 = variant.chr;
        variant2Start = variant.start;
        variant2End = variant.end;
    }
    if (pair.read2Pos > pair.read1Pos) {
        leftReadPos = pair.read1Pos;
        leftReadEnd = pair.read1End;
        leftReadSeq = pair.read1Seq;
        leftReadQual = pair.read1Qual;

        isLeftReverse = pair.isRead1Reverse;
        isLeftFirst = true;

        rightReadPos = pair.read2Pos;
        rightReadEnd = pair.read2End;
        rightReadSeq = pair.read2Seq;
        rightReadQual = pair.read2Qual;
        isRightReverse = pair.isRead2Reverse;
        isRightFirst = false;
    }
    else {
        leftReadPos = pair.read2Pos;
        leftReadEnd = pair.read2End;
        leftReadSeq = pair.read2Seq;
        leftReadQual = pair.read2Qual;
        isLeftReverse = pair.isRead2Reverse;
        isLeftFirst = false;

        rightReadPos = pair.read1Pos;
        rightReadEnd = pair.read1End;
        rightReadSeq = pair.read1Seq;
        rightReadQual = pair.read1Qual;

        isRightReverse = pair.isRead1Reverse;
        isRightFirst = true;
    }
    ReadPair newPair;
    newPair.readGroup = "OLAP" + overlap_type + pair.tag  + "ADJ" + std::to_string(adjacency) + pair.readGroup;

    if (overlap_type == "1" ) {
        std::string leftSeq = leftReadSeq;
        int nRightBases = rightReadEnd-variant1Start;
        int matchBases = variant1Start-rightReadPos;


        std::string leftMatch = rightReadSeq.substr(0,matchBases);
        std::string rightSoftClip = ref.fetchSequence(chr2, variant2Start+1, variant2Start+nRightBases);
        std::string rightSeq = leftMatch+rightSoftClip;


        if (isLeftReverse) leftSeq = reverseComplement(leftSeq);
        if (isRightReverse) rightSeq = reverseComplement(rightSeq);
        if (isLeftFirst) {
            newPair.read1Seq = leftSeq;
            newPair.read1Qual = leftReadQual;
            newPair.read2Seq = rightSeq;
            newPair.read2Qual = rightReadQual;
        }
        else {
            newPair.read2Seq = leftSeq;
            newPair.read2Qual = leftReadQual;
            newPair.read1Seq = rightSeq;
            newPair.read1Qual = rightReadQual;
        } 
    }
    else if (overlap_type == "2") {
        // Variant => 7583259:7584260 
        
        // LeftRead: 7584248 7584349	
        // RightRead: 7584264 7584365
        int diff = rightReadPos - leftReadEnd;        
        if (diff < 0) diff = 0;

        std::string rightSeq = ref.fetchSequence(chr1, variant1End+diff+1, variant1End+diff+rightReadSeq.length());
        
        int nRightBases = variant1End-leftReadPos;
        int matchBases = leftReadEnd-variant1End;

        std::string rightMatch = leftReadSeq.substr(nRightBases);        
        std::string leftSoftClip = ref.fetchSequence(chr2, variant2End-nRightBases+1, variant2End);
        std::string leftSeq =  leftSoftClip+rightMatch;
        if (newPair.readGroup == "OLAP2F2R1ADJ1VH01226:3:AAFG23VM5:1:2603:56973:9254") {
            std::cout << variant1Start << ":" << variant1End << " " << leftReadPos << " " << leftReadEnd << "\t" << rightReadPos << " " << rightReadEnd << std::endl;
            std::cout << matchBases << " " << nRightBases << " " <<diff << std::endl;
        }
        if (isLeftReverse) leftSeq = reverseComplement(leftSeq);
        if (isRightReverse) rightSeq = reverseComplement(rightSeq);
        if (isLeftFirst) {
            newPair.read1Seq = leftSeq;
            newPair.read1Qual = leftReadQual;
            newPair.read2Seq = rightSeq;
            newPair.read2Qual = rightReadQual;
        }
        else {
            newPair.read2Seq = leftSeq;
            newPair.read2Qual = leftReadQual;
            newPair.read1Seq = rightSeq;
            newPair.read1Qual = rightReadQual;
        }
    }
    else if (overlap_type == "3") {
        int matchBases = variant1Start-leftReadPos;
        int nRightBases = leftReadSeq.length()-matchBases;

        std::string leftMatch = leftReadSeq.substr(0, matchBases);
        std::string rightSoftClip = ref.fetchSequence(chr2, variant2Start+1, variant2Start+nRightBases);
        std::string leftSeq = leftMatch+rightSoftClip;


        int diff = abs(rightReadPos-variant1Start);
        std::string rightSeq = ref.fetchSequence(chr2, variant2Start+diff+1, variant2Start+diff+rightReadSeq.length());
        if (isLeftReverse) leftSeq = reverseComplement(leftSeq);
        if (isRightReverse) rightSeq = reverseComplement(rightSeq);
        if (isLeftFirst) {
            newPair.read1Seq = leftSeq;
            newPair.read1Qual = leftReadQual;
            newPair.read2Seq = rightSeq;
            newPair.read2Qual = rightReadQual;
        }
        else {
            newPair.read2Seq = leftSeq;
            newPair.read2Qual = leftReadQual;
            newPair.read1Seq = rightSeq;
            newPair.read1Qual = rightReadQual;
        }
    }
    else if (overlap_type == "4") {

        int matchBases = rightReadEnd-variant1End;
        int nLeftBases = rightReadSeq.length()-matchBases;

        std::string rightMatch = rightReadSeq.substr(nLeftBases);
        std::string leftSoftClip = ref.fetchSequence(chr2, variant2End-nLeftBases+1, variant2End);
        std::string rightSeq = leftSoftClip+rightMatch;
        int diff = abs(variant1End-leftReadEnd);

        std::string leftSeq = ref.fetchSequence(chr2, variant2End-diff-leftReadSeq.length()+1, variant2End-diff);
        if (isLeftReverse) leftSeq = reverseComplement(leftSeq);
        if (isRightReverse) rightSeq = reverseComplement(rightSeq);
        if (isLeftFirst) {
            newPair.read1Seq = leftSeq;
            newPair.read1Qual = leftReadQual;
            newPair.read2Seq = rightSeq;
            newPair.read2Qual = rightReadQual;
        }
        else {
            newPair.read2Seq = leftSeq;
            newPair.read2Qual = leftReadQual;
            newPair.read1Seq = rightSeq;
            newPair.read1Qual = rightReadQual;
        }

    }
    else if (overlap_type == "5" ) {

        std::string leftSeq = leftReadSeq;
        int diff = rightReadPos-variant1Start;


        if (diff < 0) diff = 0;
        std::string rightSeq = ref.fetchSequence(chr2, variant2Start+diff+1, variant2Start+diff+rightReadSeq.length());

        std::cout << rightSeq.length() <<  " " << rightReadQual.length() << std::endl;
        std::cout << leftSeq.length() <<  " " << leftReadQual.length() << std::endl;

        if (isLeftReverse) leftSeq = reverseComplement(leftSeq);
        if (isRightReverse) rightSeq = reverseComplement(rightSeq);
        if (isLeftFirst) {
            newPair.read1Seq = leftSeq;
            newPair.read1Qual = leftReadQual;
            newPair.read2Seq = rightSeq;
            newPair.read2Qual = rightReadQual;
        }
        else {
            newPair.read2Seq = leftSeq;
            newPair.read2Qual = leftReadQual;
            newPair.read1Seq = rightSeq;
            newPair.read1Qual = rightReadQual;
        }
    }

    else if (overlap_type == "6" ) {

        std::string rightSeq = rightReadSeq;
        int diff = variant1End-leftReadEnd;


        if (diff < 0) diff = 0;
        std::string leftSeq = ref.fetchSequence(chr2, variant2Start+diff+1, variant2Start+diff+leftReadSeq.length());

        std::cout << rightSeq.length() <<  " " << rightReadQual.length() << std::endl;
        std::cout << leftSeq.length() <<  " " << leftReadQual.length() << std::endl;

        if (isLeftReverse) leftSeq = reverseComplement(leftSeq);
        if (isRightReverse) rightSeq = reverseComplement(rightSeq);
        if (isLeftFirst) {
            newPair.read1Seq = leftSeq;
            newPair.read1Qual = leftReadQual;
            newPair.read2Seq = rightSeq;
            newPair.read2Qual = rightReadQual;
        }
        else {
            newPair.read2Seq = leftSeq;
            newPair.read2Qual = leftReadQual;
            newPair.read1Seq = rightSeq;
            newPair.read1Qual = rightReadQual;
        }
    }



    else if (overlap_type == "11" ) {

        int matchBases = variant1Start-leftReadPos;
        int nRightBases = leftReadEnd-variant1Start;

        std::string leftMatch = leftReadSeq.substr(0, matchBases);
        std::string rightSoftClip = ref.fetchSequence(chr2, variant2Start+1,variant2Start+nRightBases);
        std::string leftSeq = leftMatch+rightSoftClip;

        matchBases = variant1Start-rightReadPos;
        nRightBases = rightReadEnd-variant1Start;

        leftMatch = rightReadSeq.substr(0, matchBases);
        rightSoftClip = ref.fetchSequence(chr2, variant2Start+1, variant2Start+nRightBases);
        std::string rightSeq = leftMatch+rightSoftClip;

        if (isLeftReverse) leftSeq = reverseComplement(leftSeq);
        if (isRightReverse) rightSeq = reverseComplement(rightSeq);
        if (isLeftFirst) {
            newPair.read1Seq = leftSeq;
            newPair.read1Qual = leftReadQual;
            newPair.read2Seq = rightSeq;
            newPair.read2Qual = rightReadQual;
        }
        else {
            newPair.read2Seq = leftSeq;
            newPair.read2Qual = leftReadQual;
            newPair.read1Seq = rightSeq;
            newPair.read1Qual = rightReadQual;
        }
    }

    else if (overlap_type == "12" ) {

        int matchBases = leftReadEnd-variant1End;
        int nLeftBases = variant1End-leftReadPos;

        std::string leftMatch = leftReadSeq.substr(nLeftBases);
        std::string leftSoftClip = ref.fetchSequence(chr2, variant2End-nLeftBases+1,variant2End);
        std::string leftSeq = leftSoftClip+leftMatch;

        matchBases = rightReadEnd-variant1End;
        nLeftBases = variant1End-rightReadPos;

        leftMatch = rightReadSeq.substr(nLeftBases);
        leftSoftClip = ref.fetchSequence(chr2, variant2End-nLeftBases+1,variant2End);
        std::string rightSeq = leftSoftClip+leftMatch;

        if (isLeftReverse) leftSeq = reverseComplement(leftSeq);
        if (isRightReverse) rightSeq = reverseComplement(rightSeq);
        if (isLeftFirst) {
            newPair.read1Seq = leftSeq;
            newPair.read1Qual = leftReadQual;
            newPair.read2Seq = rightSeq;
            newPair.read2Qual = rightReadQual;
        }
        else {
            newPair.read2Seq = leftSeq;
            newPair.read2Qual = leftReadQual;
            newPair.read1Seq = rightSeq;
            newPair.read1Qual = rightReadQual;
        }
    }

    else {
        newPair.isSkipped = true;
        return newPair;
    }
    return newPair;
}



ReadPair simulateWholeArmTranslocation(ReadPair& pair, const Variant& variant, const std::string& overlap_type, RefFasta& ref, int adjacency) {
    int64_t leftReadPos;
    int64_t leftReadEnd;
    std::string leftReadSeq;
    std::string leftReadQual;
    bool isLeftReverse;
    bool isLeftFirst;

    int64_t rightReadPos;
    int64_t rightReadEnd;
    std::string rightReadSeq;
    std::string rightReadQual;
    bool isRightReverse;
    bool isRightFirst;

    std::string chr2;
    int64_t variant2Start;
    int64_t variant2End;

    std::string chr1;
    int64_t variant1Start;
    int64_t variant1End;

    if (adjacency == 1) {
        chr1 = variant.chr;
        variant1Start = variant.start;
        variant1End = variant.end;

        chr2 = variant.chr2;
        variant2Start = variant.start2;
        variant2End = variant.end2;
    }
    if (adjacency == 2) {
        chr1 = variant.chr2;
        variant1Start = variant.start2;
        variant1End = variant.end2;

        chr2 = variant.chr;
        variant2Start = variant.start;
        variant2End = variant.end;
    }
    if (pair.read2Pos > pair.read1Pos) {
        leftReadPos = pair.read1Pos;
        leftReadEnd = pair.read1End;
        leftReadSeq = pair.read1Seq;
        leftReadQual = pair.read1Qual;

        isLeftReverse = pair.isRead1Reverse;
        isLeftFirst = true;

        rightReadPos = pair.read2Pos;
        rightReadEnd = pair.read2End;
        rightReadSeq = pair.read2Seq;
        rightReadQual = pair.read2Qual;
        isRightReverse = pair.isRead2Reverse;
        isRightFirst = false;
    }
    else {
        leftReadPos = pair.read2Pos;
        leftReadEnd = pair.read2End;
        leftReadSeq = pair.read2Seq;
        leftReadQual = pair.read2Qual;
        isLeftReverse = pair.isRead2Reverse;
        isLeftFirst = false;

        rightReadPos = pair.read1Pos;
        rightReadEnd = pair.read1End;
        rightReadSeq = pair.read1Seq;
        rightReadQual = pair.read1Qual;

        isRightReverse = pair.isRead1Reverse;
        isRightFirst = true;
    }
    ReadPair newPair;
    newPair.readGroup = "OLAP" + overlap_type + pair.tag  + "ADJ" + std::to_string(adjacency) + pair.readGroup;

    if (overlap_type == "1") {
        std::string leftSeq = leftReadSeq;
        int matchBases = variant1Start-rightReadPos;
        int nRightBases = rightReadEnd-variant1Start;
        std::string leftMatch = rightReadSeq.substr(0,matchBases);
        std::string rightSoftClip = ref.fetchSequence(chr2, variant2Start+1, variant2Start+nRightBases);

        std::string rightSeq = leftMatch+rightSoftClip;
        if (isLeftReverse) leftSeq = reverseComplement(leftSeq);
        if (isRightReverse) rightSeq = reverseComplement(rightSeq);
        if (isLeftFirst) {
            newPair.read1Seq = leftSeq;
            newPair.read1Qual = leftReadQual;
            newPair.read2Seq = rightSeq;
            newPair.read2Qual = rightReadQual;
        }
        else {
            newPair.read2Seq = leftSeq;
            newPair.read2Qual = leftReadQual;
            newPair.read1Seq = rightSeq;
            newPair.read1Qual = rightReadQual;
        } 
    }
    if (overlap_type == "2") {
        int matchBases = abs(variant1Start-leftReadPos);
        int nRightBases = leftReadSeq.length()-matchBases;
        std::string leftMatch = leftReadSeq.substr(0, matchBases);
        std::string rightSoftClip = ref.fetchSequence(chr2, variant2Start+1, variant2Start+nRightBases);
        std::string leftSeq = leftMatch+rightSoftClip;
        
        int diff = abs(rightReadPos-leftReadEnd);
        std::string rightSeq = ref.fetchSequence(chr2, variant2Start+diff+1, variant2Start+diff+rightReadSeq.length());
        if (isLeftReverse) leftSeq = reverseComplement(leftSeq);
        if (isRightReverse) rightSeq = reverseComplement(rightSeq);
        if (isLeftFirst) {
            newPair.read1Seq = leftSeq;
            newPair.read1Qual = leftReadQual;
            newPair.read2Seq = rightSeq;
            newPair.read2Qual = rightReadQual;
        }
        else {
            newPair.read2Seq = leftSeq;
            newPair.read2Qual = leftReadQual;
            newPair.read1Seq = rightSeq;
            newPair.read1Qual = rightReadQual;
        } 

        if (leftSeq.length() > leftReadQual.length()) {
            std::cout << " error: " << variant1Start << " " << leftReadPos << " " << matchBases << " " << nRightBases << std::endl;
        }

        if (rightSeq.length() > rightReadQual.length()) {
            std::cout << " error: " << variant1Start << " " << rightReadPos << " "  << matchBases << " " << nRightBases << std::endl;
        }
    }

    if (overlap_type == "11" || overlap_type == "12") {

        int matchBases = abs(variant1Start-leftReadPos);
        int nRightBases = leftReadSeq.length()-matchBases;

        std::string leftMatch = leftReadSeq.substr(0, matchBases);
        std::string rightSoftClip = ref.fetchSequence(chr2, variant2Start+1, variant2Start+nRightBases);
        std::string leftSeq = leftMatch+rightSoftClip;

        matchBases = abs(variant1Start-rightReadPos);
        nRightBases = rightReadSeq.length()-matchBases;

        leftMatch = rightReadSeq.substr(0, matchBases);
        rightSoftClip = ref.fetchSequence(chr2, variant2Start+1, variant2Start+nRightBases);
        std::string rightSeq = leftMatch+rightSoftClip;
        if (isLeftReverse) leftSeq = reverseComplement(leftSeq);
        if (isRightReverse) rightSeq = reverseComplement(rightSeq);
        if (isLeftFirst) {
            newPair.read1Seq = leftSeq;
            newPair.read1Qual = leftReadQual;
            newPair.read2Seq = rightSeq;
            newPair.read2Qual = rightReadQual;
        }
        else {
            newPair.read2Seq = leftSeq;
            newPair.read2Qual = leftReadQual;
            newPair.read1Seq = rightSeq;
            newPair.read1Qual = rightReadQual;
        }
        if (leftSeq.length() > leftReadQual.length()) {
            std::cout << " error: " << variant1Start << " " << leftReadPos << " " << matchBases << " " << nRightBases << std::endl;
        }

        if (rightSeq.length() > rightReadQual.length()) {
            std::cout << " error: " << variant1Start << " " << rightReadPos << " "  << matchBases << " " << nRightBases << std::endl;
        }


    }
    if (overlap_type == "13") {

        std::string leftSeq = leftReadSeq;
        int diff = abs(rightReadPos-leftReadEnd);
        std::string rightSeq = ref.fetchSequence(chr2, variant2Start+diff+1, variant2Start+diff+rightReadSeq.length());
        if (isLeftReverse) leftSeq = reverseComplement(leftSeq);
        if (isRightReverse) rightSeq = reverseComplement(rightSeq);
        if (isLeftFirst) {
            newPair.read1Seq = leftSeq;
            newPair.read1Qual = leftReadQual;
            newPair.read2Seq = rightSeq;
            newPair.read2Qual = rightReadQual;
        }
        else {
            newPair.read2Seq = leftSeq;
            newPair.read2Qual = leftReadQual;
            newPair.read1Seq = rightSeq;
            newPair.read1Qual = rightReadQual;
        } 
        
    }
    return newPair;

}