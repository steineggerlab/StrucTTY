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
    float aspect            = 1.0f;
    std::string theme       = "dark";
    std::string assembly;
    std::string scale_file;

    std::string chains_file;
    std::string msa_file;
    std::string foldmason_file;

    std::string foldseek_target;
    std::string foldseek_result;
};

bool run(const RunOptions& opts);

}
