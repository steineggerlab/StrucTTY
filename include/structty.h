#pragma once
#include <string>
#include <vector>

namespace structty {

struct RunOptions {
    std::vector<std::string> input_files;
    std::string mode        = "protein";
    bool show_structure     = false;
    bool no_panel           = false;
    bool benchmark          = false;

    std::string chains_file;
    std::string msa_file;
    std::string foldmason_file;

    // -fst : target 구조 소스. Foldseek DB | 구조 디렉터리 | 구조 파일 | "auto"(다운로드).
    // 종류는 input_probe::probe() 로 판별하므로 호출부는 경로만 채우면 된다.
    std::string foldseek_target;
    // -fsr : Foldseek 결과. m8(12/17/21/29 컬럼) 또는 멀티머 _report(14 컬럼).
    // 컬럼 수로 멀티머 경로 진입이 결정된다(과거 report_format 플래그 대체).
    std::string foldseek_result;

    // -m aligned 에서 정렬됨으로 칠할 CA 쌍의 최대 거리(Å).
    // foldseek backtrace 는 구조가 벌어진 구간도 짝지으므로 이 값으로 잘라낸다.
    float align_cutoff = 5.0f;
};

// Launch the interactive viewer. Blocks until the user presses Q.
void run(const RunOptions& opts);

} // namespace structty
