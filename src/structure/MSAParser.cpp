#include "MSAParser.hpp"

#include <fstream>
#include <cmath>
#include <array>
#include <algorithm>

static std::string normalize_a3m(const std::string& seq) {
    std::string result;
    result.reserve(seq.size());
    for (char c : seq) {
        if (c == '-' || (c >= 'A' && c <= 'Z')) {
            result += c;
        }
    }
    return result;
}

bool MSAParser::load(const std::string& filepath) {
    std::ifstream ifs(filepath);
    if (!ifs.is_open()) {
        return false;
    }

    sequence_names.clear();
    sequences.clear();
    query_sequence.clear();
    msa_pos_to_query_idx.clear();

    std::string line;

    while (std::getline(ifs, line)) {
        if (!line.empty()) break;
    }

    if (line.empty()) return false;

    if (line.find("# STOCKHOLM") == 0) {
        return false;
    }

    if (line[0] != '>') {
        return false;
    }

    std::string cur_name = line.substr(1);
    std::string cur_seq;

    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!cur_seq.empty()) {
                sequence_names.push_back(cur_name);
                sequences.push_back(normalize_a3m(cur_seq));
            }
            cur_name = line.substr(1);
            cur_seq.clear();
        } else {
            cur_seq += line;
        }
    }

    if (!cur_seq.empty()) {
        sequence_names.push_back(cur_name);
        sequences.push_back(normalize_a3m(cur_seq));
    }

    if (sequences.empty()) return false;

    query_sequence = sequences[0];
    return true;
}

std::vector<float> MSAParser::compute_conservation() {
    if (sequences.empty()) return {};

    int query_len = (int)query_sequence.size();
    if (query_len == 0) return {};

    std::vector<float> conservation(query_len, -1.0f);

    const float log2_20 = std::log2(20.0f);

    for (int pos = 0; pos < query_len; ++pos) {
        std::array<int, 26> counts{};
        int total = 0;

        for (const std::string& seq : sequences) {
            if (pos >= (int)seq.size()) continue;
            char c = seq[pos];
            if (c >= 'A' && c <= 'Z') {
                counts[static_cast<int>(c - 'A')]++;
                total++;
            }
        }

        if (total == 0) {
            conservation[pos] = -1.0f;
            continue;
        }

        float H = 0.0f;
        for (int aa = 0; aa < 26; ++aa) {
            if (counts[aa] == 0) continue;
            float f = static_cast<float>(counts[aa]) / static_cast<float>(total);
            H -= f * std::log2(f);
        }

        float score = 1.0f - H / log2_20;
        conservation[pos] = std::max(0.0f, std::min(1.0f, score));
    }

    return conservation;
}
