#include "FoldseekParser.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>

// 탭으로 구분된 행을 컬럼 벡터로 분리
static std::vector<std::string> split_tabs(const std::string& line) {
    std::vector<std::string> cols;
    std::istringstream ss(line);
    std::string col;
    while (std::getline(ss, col, '\t')) {
        cols.push_back(col);
    }
    return cols;
}

static float safe_stof(const std::string& s, float default_val = 0.0f) {
    try { return std::stof(s); } catch (...) { return default_val; }
}

static int safe_stoi(const std::string& s, int default_val = 0) {
    try { return std::stoi(s); } catch (...) { return default_val; }
}

bool FoldseekParser::parse_line(const std::vector<std::string>& cols,
                                FoldseekHit& hit, int fmt) {
    // fmt: 12, 17, or 29
    if ((int)cols.size() < fmt) return false;

    hit.query  = cols[0];
    hit.target = cols[1];
    hit.fident = safe_stof(cols[2]);
    hit.alnlen = safe_stoi(cols[3]);
    hit.evalue = safe_stof(cols[10], 0.0f);
    hit.bits   = safe_stof(cols[11], 0.0f);

    if (fmt == 12) {
        hit.has_transform = false;
        hit.has_aln       = false;
        return true;
    }

    // fmt 17 or 29: cols[12]=lddt, [13]=qtmscore, [14]=ttmscore
    hit.lddt      = safe_stof(cols[12], -1.0f);
    hit.qtmscore  = safe_stof(cols[13], -1.0f);
    hit.ttmscore  = safe_stof(cols[14], -1.0f);

    if (fmt == 17) {
        // cols[15]=qaln, [16]=taln
        hit.qaln = cols[15];
        hit.taln = cols[16];
        hit.has_aln       = true;
        hit.has_transform = false;
        return true;
    }

    // fmt == 29: cols[15..23]=U (9 values), cols[24..26]=T (3 values)
    //            cols[27]=qaln, cols[28]=taln
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

    int fmt = 0;  // detected format (12, 17, or 29)
    std::string line;

    while (std::getline(ifs, line)) {
        // skip empty lines
        if (line.empty()) continue;
        // skip comment lines
        if (line[0] == '#') continue;

        std::vector<std::string> cols = split_tabs(line);
        int ncols = (int)cols.size();

        // detect format from the first data line
        if (fmt == 0) {
            if (ncols == 12) {
                fmt = 12;
            } else if (ncols == 17) {
                fmt = 17;
            } else if (ncols == 29) {
                fmt = 29;
            } else {
                // 알 수 없는 형식: 첫 줄의 컬럼 수가 기준
                // 지원하는 형식이 아니면 로드 실패
                return false;
            }
        }

        FoldseekHit hit;
        if (parse_line(cols, hit, fmt)) {
            hits.push_back(std::move(hit));
        }
    }

    return !hits.empty();
}
