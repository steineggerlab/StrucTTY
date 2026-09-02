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

    std::string foldseek_target;
    std::string foldseek_result;
};

bool run(const RunOptions& opts);

}
