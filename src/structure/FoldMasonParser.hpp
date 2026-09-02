#pragma once
#include <string>
#include <vector>
#include <array>
#include <utility>

struct FoldMasonEntry {
    std::string name;
    std::string aa;
    std::string ss;
    std::vector<std::array<float,3>> ca;
};

class FoldMasonParser {
public:
    bool load_json(const std::string& path);
    bool load_fasta(const std::string& path);

    const std::vector<FoldMasonEntry>& get_entries() const { return entries; }
    const std::vector<float>& get_scores() const { return scores; }
    float get_msa_lddt() const { return msa_lddt; }
    int msa_length() const;
    int entry_count() const { return (int)entries.size(); }

    std::vector<int> build_query_col_map(int ref_idx = 0) const;

    std::vector<std::pair<int,int>> build_aligned_pairs(int ref_idx, int other_idx) const;

    std::vector<float> compute_column_entropy(bool use_ss = false) const;

private:
    std::vector<FoldMasonEntry> entries;
    std::vector<float> scores;
    float msa_lddt = -1.0f;
    std::string format;
};
