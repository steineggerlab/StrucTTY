#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cstdint>
#include <cstdio>

class FoldseekDBReader {
public:
    FoldseekDBReader();
    ~FoldseekDBReader();

    bool open(const std::string& basePath);

    bool prepare(const std::unordered_set<std::string>& target_ids);

    void close();
    bool is_open() const;

    size_t read_entry(const std::string& accession,
                      std::vector<float>& coords,
                      std::string& aa_seq);

private:
    struct IndexEntry { size_t offset; size_t length; };

    std::unordered_map<std::string, unsigned int> lookup_;
    std::unordered_map<unsigned int, IndexEntry>  ca_index_;
    std::unordered_map<unsigned int, IndexEntry>  seq_index_;

    std::string base_path_;
    FILE* ca_data_file_;
    FILE* seq_data_file_;
    bool opened_;

    bool scan_lookup(const std::string& path,
                     const std::unordered_set<std::string>& target_ids);
    bool scan_index(const std::string& path,
                    const std::unordered_set<unsigned int>& target_internal_ids,
                    std::unordered_map<unsigned int, IndexEntry>& index_out);
    bool read_raw(FILE* file, size_t offset, size_t length, std::vector<char>& buf);
    void decode_ca(const char* data, size_t entryLen, size_t chainLen,
                   std::vector<float>& out);
};
