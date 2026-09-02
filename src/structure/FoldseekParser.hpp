#pragma once

#include <string>
#include <vector>

struct FoldseekHit {
    std::string query;
    std::string target;
    float fident    = 0.0f;
    int   alnlen    = 0;

    int   qstart    = 1;
    int   qend      = 0;
    int   tstart    = 1;
    int   tend      = 0;

    float prob      = -1.0f;
    float evalue    = 0.0f;
    float bits      = 0.0f;

    int   qlength   = 0;
    int   tlength   = 0;

    float lddt      = -1.0f;
    float qtmscore  = -1.0f;
    float ttmscore  = -1.0f;

    float U[9] = {};
    float T[3] = {};
    bool has_transform = false;

    std::string qaln;
    std::string taln;
    bool has_aln = false;

    std::vector<float> alns;
    std::string tseq;
    std::string taxid;
    std::string taxname;
    bool is_alis_format = false;
};

class FoldseekParser {
public:
    bool load(const std::string& filepath);

    const std::vector<FoldseekHit>& get_hits() const { return hits; }
    int hit_count() const { return (int)hits.size(); }

private:
    std::vector<FoldseekHit> hits;

    bool parse_line(const std::vector<std::string>& cols, FoldseekHit& hit, int fmt);
};
