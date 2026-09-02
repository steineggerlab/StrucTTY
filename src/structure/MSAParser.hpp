#pragma once

#include <string>
#include <vector>

class MSAParser {
public:
    bool load(const std::string& filepath);

    std::vector<float> compute_conservation();

private:
    std::vector<std::string> sequence_names;
    std::vector<std::string> sequences;
    std::string query_sequence;
    std::vector<int> msa_pos_to_query_idx;
};
