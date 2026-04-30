#include "structty.h"
#include "Parameters.hpp"

int main(int argc, char* argv[]) {
    Parameters params(argc, argv);
    if (!params.check_arg_okay()) return -1;
    params.print_args();

    structty::RunOptions opts;
    opts.input_files    = params.get_in_file();
    opts.mode           = params.get_mode();
    opts.show_structure = params.get_show_structure();
    opts.no_panel       = params.get_no_panel();
    opts.benchmark      = params.get_benchmark_mode();
    opts.chains_file    = params.get_chainfile();
    opts.ut_matrix_file = params.get_utmatrix();
    opts.msa_file       = params.get_msa_file();
    opts.foldseek_file  = params.get_foldseek_file();
    opts.foldseek_db_path = params.get_db_path();
    opts.foldseek_db    = params.get_foldseek_db();
    opts.foldmason_file = params.get_foldmason_file();

    structty::run(opts);
    return 0;
}
