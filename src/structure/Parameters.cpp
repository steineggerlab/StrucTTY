#include "Parameters.hpp"

void print_help(){
    std::cout << "Usage: StrucTTY <input_files...> [OPTIONS]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -m, --mode <MODE>       Color mode:\n";
    std::cout << "                            protein (default), chain, rainbow,\n";
    std::cout << "                            plddt, interface, conservation, aligned\n";
    std::cout << "  -c, --chains <FILE>     Show only selected chains (see example/chainfile)\n";
    std::cout << "  -s, --structure         Show secondary structure (alpha helix, beta sheet)\n";
    std::cout << "  -ut, --utmatrix <FILE>  Apply rotation/translation matrix (see example/utfile)\n";
    std::cout << "  --msa <FILE>            MSA file for conservation score (FASTA/A3M)\n";
    std::cout << "  -fs, --foldseek <FILE>  Foldseek result file for hit navigation\n";
    std::cout << "  --db-path <DIR>         Target PDB directory for Foldseek hit loading\n";
    std::cout << "  -fm, --foldmason <FILE> FoldMason result (JSON or FASTA MSA)\n";
    std::cout << "  -n, --nopanel           Hide info panel\n";
    std::cout << "  -b, --benchmark         Benchmark mode (measure FPS/latency)\n";
    std::cout << "  --help                  Show this help message\n";
}
Parameters::Parameters(int argc, char* argv[]) {
    arg_okay = true;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--help")) {
            print_help();
            std::exit(0);
        }
    }
    
    if (argc <= 1) {
        std::cerr << "Need input file dir" << std::endl;
        arg_okay = false;
        return;
    }

    for (int i = 1; i < argc; i++) {
        try {
            if (!strcmp(argv[i], "-m") || !strcmp(argv[i], "--mode")) {
                if (i + 1 < argc) {
                    std::string val(argv[i + 1]);
                    std::transform(val.begin(), val.end(), val.begin(), ::tolower); // to lowercase
                    if (val == "chain" || val == "rainbow" || val == "protein" ||
                        val == "plddt" || val == "interface" || val == "conservation" || val == "aligned") {
                        mode = val;
                        i++;
                    } else {
                        throw std::runtime_error("Error: Invalid value for --mode. Use 'protein', 'chain', 'rainbow', 'plddt', 'interface', 'conservation', or 'aligned'.");
                    }
                } else {
                    throw std::runtime_error("Error: Missing value for -m / --mode.");
                }
            } else if (!strcmp(argv[i], "-c") || !strcmp(argv[i], "--chains")) {
                if (i + 1 < argc) {
                    chainfile = argv[++i];
                } else {
                    throw std::runtime_error("Error: Missing value for -c / --chains.");
                }
            }
            else if (!strcmp(argv[i], "-s") || !strcmp(argv[i], "--structure")) {
                show_structure = true;
            } else if (!strcmp(argv[i], "-n") || !strcmp(argv[i], "--nopanel")) {
                no_panel = true;
            } else if (!strcmp(argv[i], "-ut") || !strcmp(argv[i], "--utmatrix")) {
                if (i + 1 < argc) {
                    utmatrix = argv[++i];
                } else {
                    throw std::runtime_error("Error: Missing value for -ut / --utmatrix.");
                }
            } else if (fs::exists(argv[i]) && fs::is_regular_file(argv[i]) && in_file.size() < 9){
                in_file.push_back(argv[i]);
            } else if (!strcmp(argv[i], "--msa")) {
                if (i + 1 < argc) {
                    msa_file = argv[++i];
                } else {
                    throw std::runtime_error("Error: Missing value for --msa.");
                }
            } else if (!strcmp(argv[i], "-fs") || !strcmp(argv[i], "--foldseek")) {
                if (i + 1 < argc) {
                    foldseek_file = argv[++i];
                } else {
                    throw std::runtime_error("Error: Missing value for -fs / --foldseek.");
                }
            } else if (!strcmp(argv[i], "--db-path")) {
                if (i + 1 < argc) {
                    db_path = argv[++i];
                } else {
                    throw std::runtime_error("Error: Missing value for --db-path.");
                }
            } else if (!strcmp(argv[i], "--foldmason") || !strcmp(argv[i], "-fm")) {
                if (i + 1 < argc) {
                    foldmason_file = argv[++i];
                } else {
                    throw std::runtime_error("Error: Missing value for --foldmason / -fm.");
                }
            } else if (!strcmp(argv[i], "-b") || !strcmp(argv[i], "--benchmark")) {
                benchmark_mode = true;
                show_structure = true;
            } else {
                throw std::runtime_error("Error: Unknown parameter: " + std::string(argv[i]));
            }
        }       
        catch (const std::exception& e) {
            std::cerr << "Wrong input parameters: " << e.what() << std::endl;
            std::cerr << "Error at argument: " << argv[i] << std::endl;
            arg_okay = false;
            return;
        }
    }
    if (in_file.size() == 0){
        std::cerr << "Error: Need input file dir" << std::endl;
        arg_okay = false;
        return;
    }
    return;
}

void Parameters::print_args() {
    cout << "Input parameters >> " << endl;
    cout << "  in_file: " << endl;
    for (size_t i = 0; i < in_file.size(); i++) {
        std::cout << "\t" << in_file[i] << '\n';
    }
    cout << "  mode: " << mode << endl;
    cout << "  utmatrix: " << utmatrix << endl;
    cout << "  chainfile: " << chainfile << endl;
    cout << "  show_structure: " << show_structure << endl;
    cout << "  benchmark_mode: " << benchmark_mode << endl;

    cout << "\n";
    return;
}