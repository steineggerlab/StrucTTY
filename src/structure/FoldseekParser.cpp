#include "FoldseekParser.hpp"

#include <fstream>
#include <sstream>

static std::vector<std::string> split_tabs(const std::string& line) {
    std::vector<std::string> cols;
    std::istringstream ss(line);
    std::string col;
    while (std::getline(ss, col, '\t')) {
        cols.push_back(col);
    }
    return cols;
}

static std::vector<float> parse_csv_floats(const std::string& s) {
    std::vector<float> result;
    std::istringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        if (tok.empty()) continue;
        try { result.push_back(std::stof(tok)); } catch (...) {}
    }
    return result;
}

static float safe_stof(const std::string& s, float default_val = 0.0f) {
    try { return std::stof(s); } catch (...) { return default_val; }
}

static int safe_stoi(const std::string& s, int default_val = 0) {
    try { return std::stoi(s); } catch (...) { return default_val; }
}

static std::string extract_target_id(const std::string& raw) {
    size_t pos = raw.find(' ');
    if (pos == std::string::npos) return raw;
    return raw.substr(0, pos);
}

bool FoldseekParser::parse_line(const std::vector<std::string>& cols,
                                FoldseekHit& hit, int fmt) {
    if ((int)cols.size() < fmt) return false;

    hit.query  = cols[0];
    hit.target = extract_target_id(cols[1]);
    hit.fident = safe_stof(cols[2]);
    hit.alnlen = safe_stoi(cols[3]);
    hit.qstart = safe_stoi(cols[6], 1);
    hit.qend   = safe_stoi(cols[7], 0);
    hit.tstart = safe_stoi(cols[8], 1);
    hit.tend   = safe_stoi(cols[9], 0);

    if (fmt == 12) {
        hit.evalue = safe_stof(cols[10], 0.0f);
        hit.bits   = safe_stof(cols[11], 0.0f);
        hit.has_transform = false;
        hit.has_aln       = false;
        return true;
    }

    if (fmt == 21) {
        hit.prob    = safe_stof(cols[10], -1.0f);
        hit.evalue  = safe_stof(cols[11], 0.0f);
        hit.bits    = safe_stof(cols[12], 0.0f);
        hit.qlength = safe_stoi(cols[13], 0);
        hit.tlength = safe_stoi(cols[14], 0);
        hit.qaln    = cols[15];
        hit.taln    = cols[16];
        hit.alns    = parse_csv_floats(cols[17]);
        hit.tseq    = cols[18];
        hit.taxid   = cols[19];
        hit.taxname = cols[20];
        hit.has_aln       = (!hit.qaln.empty() && !hit.taln.empty());
        hit.has_transform = false;
        hit.is_alis_format = true;
        return true;
    }

    hit.lddt      = safe_stof(cols[12], -1.0f);
    hit.qtmscore  = safe_stof(cols[13], -1.0f);
    hit.ttmscore  = safe_stof(cols[14], -1.0f);
    hit.evalue    = safe_stof(cols[10], 0.0f);
    hit.bits      = safe_stof(cols[11], 0.0f);

    if (fmt == 17) {
        hit.qaln = cols[15];
        hit.taln = cols[16];
        hit.has_aln       = true;
        hit.has_transform = false;
        return true;
    }

    for (int i = 0; i < 9; ++i) {
        hit.U[i] = safe_stof(cols[15 + i]);
    }
    for (int i = 0; i < 3; ++i) {
        hit.T[i] = safe_stof(cols[24 + i]);
    }
    hit.qaln = cols[27];
    hit.taln = cols[28];
    hit.has_transform = true;
    hit.has_aln       = true;
    return true;
}

bool FoldseekParser::load(const std::string& filepath) {
    hits.clear();

    std::ifstream ifs(filepath);
    if (!ifs.is_open()) {
        return false;
    }

    int fmt = 0;
    std::string line;

    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') continue;

        std::vector<std::string> cols = split_tabs(line);
        int ncols = (int)cols.size();

        if (fmt == 0) {
            if (ncols == 12) {
                fmt = 12;
            } else if (ncols == 17) {
                fmt = 17;
            } else if (ncols == 21) {
                fmt = 21;
            } else if (ncols == 29) {
                fmt = 29;
            } else {
                if (ncols > 21) {
                    fmt = 21;
                    while ((int)cols.size() > 21) {
                        cols[20] += "\t" + cols[21];
                        cols.erase(cols.begin() + 21);
                    }
                    ncols = 21;
                } else {
                    return false;
                }
            }
        }

        if (fmt == 21 && (int)cols.size() > 21) {
            while ((int)cols.size() > 21) {
                cols[20] += "\t" + cols[21];
                cols.erase(cols.begin() + 21);
            }
        }

        FoldseekHit hit;
        if (parse_line(cols, hit, fmt)) {
            hits.push_back(std::move(hit));
        }
    }

    return !hits.empty();
}
