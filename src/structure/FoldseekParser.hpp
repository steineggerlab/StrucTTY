#pragma once

#include <string>
#include <vector>

// 기능 4: Foldseek 결과 파일(.m8) 파싱
// foldseek easy-search 출력 형식에서 U/T 변환행렬 + alignment string 추출

struct FoldseekHit {
    std::string query;
    std::string target;
    float fident    = 0.0f;
    int   alnlen    = 0;
    float evalue    = 0.0f;
    float bits      = 0.0f;
    float lddt      = -1.0f;
    float qtmscore  = -1.0f;
    float ttmscore  = -1.0f;

    // U: row-major 3x3 rotation matrix (flat, 9 elements)
    // T: translation vector (3 elements)
    // transform: target_transformed = U * target_pos + T  (≈ query_pos)
    float U[9] = {};
    float T[3] = {};
    bool has_transform = false;

    std::string qaln;
    std::string taln;
    bool has_aln = false;
};

class FoldseekParser {
public:
    // filepath: Foldseek easy-search 출력 파일(.m8 형식, 탭 구분, 헤더 없음)
    // 반환: 성공 시 true, 실패 시 false
    bool load(const std::string& filepath);

    const std::vector<FoldseekHit>& get_hits() const { return hits; }
    int hit_count() const { return (int)hits.size(); }

private:
    std::vector<FoldseekHit> hits;

    // 컬럼 수 기반 포맷 자동 감지 후 파싱
    // 지원 포맷:
    //   12 cols: 기본 m8 (query~bits, U/T/qaln/taln 없음)
    //   17 cols: 축소 형식 (lddt, qtmscore, ttmscore, qaln, taln 포함, U/T 없음)
    //   29 cols: 전체 형식 (lddt, qtmscore, ttmscore, U(9), T(3), qaln, taln 포함)
    bool parse_line(const std::vector<std::string>& cols, FoldseekHit& hit, int fmt);
};
