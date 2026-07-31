#include "Parameters.hpp"

#include "InputProbe.hpp"

void print_help(){
    std::cout << "Usage: StrucTTY <query...> [OPTIONS]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -m, --mode <MODE>       Color mode:\n";
    std::cout << "                            protein (default), chain, rainbow,\n";
    std::cout << "                            plddt, interface, conservation, aligned\n";
    std::cout << "  -c, --chains <FILE>     Show only selected chains (see example/chainfile)\n";
    std::cout << "  -s, --structure         Show secondary structure (alpha helix, beta sheet)\n";
    std::cout << "  --msa <FILE>            MSA file for conservation score (FASTA/A3M)\n";
    std::cout << "  -fst, --foldseek-target <PATH>\n";
    std::cout << "                          Target source for Foldseek hits: Foldseek DB,\n";
    std::cout << "                          structure directory, structure file, or 'auto'\n";
    std::cout << "                          ('auto' downloads hits from public DBs)\n";
    std::cout << "  -fsr, --foldseek-result <FILE>\n";
    std::cout << "                          Foldseek result: m8 (12/17/21/29 columns) or\n";
    std::cout << "                          multimer _report (14 columns)\n";
    std::cout << "                          -fst and -fsr must be given together\n";
    std::cout << "  -fm, --foldmason <FILE> FoldMason result (JSON or FASTA MSA)\n";
    std::cout << "  -n, --nopanel           Hide info panel\n";
    std::cout << "  -b, --benchmark         Benchmark mode (measure FPS/latency)\n";
    std::cout << "  --help                  Show this help message\n";
    std::cout << "\nSupported inputs (4 kinds; detected automatically):\n";
    std::cout << "  1. structure file       .pdb / .cif / .ent (+ .gz)\n";
    std::cout << "  2. structure directory  a directory of those files (-fst only)\n";
    std::cout << "  3. Foldseek DB          base path of a DB built from structures\n";
    std::cout << "                          (needs <db>_ca; sequence-derived DBs have none)\n";
    std::cout << "  4. Foldseek result      m8 (12/17/21/29 columns) or\n";
    std::cout << "                          multimer _report (14 columns)\n";
    std::cout << "  Sequence FASTA is NOT supported -- it carries no 3D coordinates.\n";
    std::cout << "\nRecipes:\n";
    std::cout << "  StrucTTY query.cif\n";
    std::cout << "  StrucTTY query.cif target.cif -m aligned\n";
    std::cout << "  StrucTTY query.cif -fst targetDB   -fsr result.m8   -m aligned\n";
    std::cout << "  StrucTTY query.cif -fst pdb_dir/   -fsr result.m8\n";
    std::cout << "  StrucTTY query.cif -fst auto       -fsr result.m8\n";
    std::cout << "  StrucTTY queryDB   -fst targetDB   -fsr result.m8   (multi-query: ]/[)\n";
    std::cout << "  StrucTTY queryDB   -fst targetDB   -fsr out_report  (multimer)\n";
}

// 좌표가 없는 입력을 렌더 시작 전에 거른다.
// 서열 FASTA: foldseek createdb 의 ProstT5 경로는 3Di(`_ss`)만 예측하고 `_ca` 를 만들지 않는다.
// `_ca` 없는 DB: 서열 유래이거나 `--index-exclude 2` 로 만든 DB.
static bool check_has_coordinates(const std::string& path, const std::string& role) {
    const input_probe::InputKind kind = input_probe::probe(path);

    if (kind == input_probe::InputKind::SequenceFasta) {
        std::cerr << "Error: " << role << " '" << path << "' is a sequence FASTA.\n"
                  << "       StrucTTY needs 3D coordinates: pass PDB/mmCIF files, a directory\n"
                  << "       of them, or a Foldseek DB built from structures.\n"
                  << "       (foldseek createdb --prostt5-model predicts 3Di but writes no _ca)"
                  << std::endl;
        return false;
    }
    if (kind == input_probe::InputKind::FoldseekDB && !input_probe::db_has_ca(path)) {
        std::cerr << "Error: Foldseek DB '" << path << "' has no C-alpha coordinates ("
                  << path << "_ca is missing).\n"
                  << "       It was built from sequences, or with --index-exclude 2."
                  << std::endl;
        return false;
    }
    return true;
}
Parameters::Parameters(int argc, char* argv[]) {
    arg_okay = true;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
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
            } else if (fs::exists(argv[i]) && fs::is_regular_file(argv[i]) && in_file.size() < 9){
                in_file.push_back(argv[i]);
            } else if (!strcmp(argv[i], "--msa")) {
                if (i + 1 < argc) {
                    msa_file = argv[++i];
                } else {
                    throw std::runtime_error("Error: Missing value for --msa.");
                }
            } else if (!strcmp(argv[i], "-fst") || !strcmp(argv[i], "--foldseek-target")) {
                if (i + 1 < argc) {
                    foldseek_target = argv[++i];
                } else {
                    throw std::runtime_error("Error: Missing value for -fst / --foldseek-target.");
                }
            } else if (!strcmp(argv[i], "-fsr") || !strcmp(argv[i], "--foldseek-result")) {
                if (i + 1 < argc) {
                    foldseek_result = argv[++i];
                } else {
                    throw std::runtime_error("Error: Missing value for -fsr / --foldseek-result.");
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
    // 검증 1 — query 위치 인자는 최소 1개
    if (in_file.size() == 0) {
        std::cerr << "Error: Need input file dir" << std::endl;
        arg_okay = false;
        return;
    }

    // 검증 2 — 쌍 규칙: -fst 와 -fsr 은 항상 함께 주어져야 한다.
    // target 소스 없이 결과만 주면 hit 구조를 읽을 수 없고, 결과 없이 target 만 주면
    // 어떤 hit 을 띄울지 알 수 없다.
    if (!foldseek_target.empty() && foldseek_result.empty()) {
        std::cerr << "Error: -fst / --foldseek-target given without -fsr / --foldseek-result.\n"
                  << "       Both must be given together." << std::endl;
        arg_okay = false;
        return;
    }
    if (foldseek_target.empty() && !foldseek_result.empty()) {
        std::cerr << "Error: -fsr / --foldseek-result given without -fst / --foldseek-target.\n"
                  << "       Both must be given together. Use -fst auto to download hit structures."
                  << std::endl;
        arg_okay = false;
        return;
    }

    // 검증 3·4 — query 각각과 target 이 좌표를 가진 입력인지 확인.
    // "auto" 는 예약값이라 경로 판별을 건너뛴다(다운로드 모드).
    for (const std::string& q : in_file) {
        if (!check_has_coordinates(q, "query")) {
            arg_okay = false;
            return;
        }
    }
    if (!foldseek_target.empty() && foldseek_target != "auto") {
        if (!check_has_coordinates(foldseek_target, "-fst target")) {
            arg_okay = false;
            return;
        }
    }

    if (!foldseek_result.empty()) {
        const input_probe::InputKind result_kind = input_probe::probe(foldseek_result);

        // 검증 5 — 결과 파일은 m8(12/17/21/29) 또는 멀티머 _report(14) 여야 한다
        if (result_kind != input_probe::InputKind::ResultM8 &&
            result_kind != input_probe::InputKind::ResultReport) {
            const int ncols = input_probe::tsv_column_count(foldseek_result);
            std::cerr << "Error: -fsr '" << foldseek_result
                      << "' is not a Foldseek result file.\n"
                      << "       Expected tab-separated 12/17/21/29 columns (m8) or 14 columns"
                      << " (multimer _report),\n       ";
            if (ncols == 0) {
                std::cerr << "but the file could not be read or has no data rows.";
            } else {
                std::cerr << "but found " << ncols << " columns.";
            }
            std::cerr << std::endl;
            arg_okay = false;
            return;
        }

        // 검증 6 — 멀티머 _report 는 complex 체인을 query DB 에서 읽으므로 query 가 DB 여야 한다
        if (result_kind == input_probe::InputKind::ResultReport) {
            const input_probe::InputKind query_kind = input_probe::probe(in_file[0]);
            if (query_kind != input_probe::InputKind::FoldseekDB) {
                std::cerr << "Error: -fsr '" << foldseek_result
                          << "' is a multimer _report (14 columns), which needs a Foldseek\n"
                          << "       query DB as the query input (chain entries are read per"
                          << " complex from the DB).\n"
                          << "       Got " << input_probe::kind_name(query_kind) << ": "
                          << in_file[0] << std::endl;
                arg_okay = false;
                return;
            }
        }
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
    cout << "  chainfile: " << chainfile << endl;
    cout << "  show_structure: " << show_structure << endl;
    cout << "  benchmark_mode: " << benchmark_mode << endl;
    if (!foldseek_target.empty() || !foldseek_result.empty()) {
        cout << "  foldseek_target: " << foldseek_target << endl;
        cout << "  foldseek_result: " << foldseek_result << endl;
    }

    cout << "\n";
    return;
}
