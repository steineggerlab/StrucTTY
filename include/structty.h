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
    std::string ut_matrix_file;
    std::string msa_file;
    std::string foldseek_file;
    std::string foldseek_db_path;  // --db-path (PDB download directory)
    std::string foldseek_db;       // --db (direct Foldseek DB path)
    std::string foldseek_query_db; // query tmp DB path (read query structures from DB)
    std::string foldmason_file;
    bool report_format     = false; // foldseek_file is a multimer _report (14-col tsv), not .m8 (D6)
};

// Launch the interactive viewer. Blocks until the user presses Q.
void run(const RunOptions& opts);

} // namespace structty
