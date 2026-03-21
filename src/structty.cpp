#include <iostream>
#include <unistd.h>
#include <clocale>
#include <filesystem>
#include "Curses.hpp"
#include "Common.hpp"
#include "Protein.hpp"
#include "Parameters.hpp"
#include "Screen.hpp"
#include "MSAParser.hpp"
#include "FoldseekParser.hpp"
#include "FoldMasonParser.hpp"
#include "Benchmark.hpp"

int main(int argc, char* argv[]) {
    Parameters params(argc, argv);
    if (!params.check_arg_okay()) return -1; 
    params.print_args();
    
    setlocale(LC_ALL, "");
    initscr();
    cbreak();
    noecho();

    int term_rows, term_cols;
    getmaxyx(stdscr, term_rows, term_cols);

    Screen screen(term_cols, term_rows,
                  params.get_show_structure(),
                  params.get_mode(),
                  params.get_depthcharacter());
    
    Benchmark bm;
    const bool bench = params.get_benchmark_mode();
    using BenchClock = Benchmark::clock;
    BenchClock::time_point t_load0;

    if (bench) {
        fs::path in_file(params.get_in_file(0));
        bm.start(std::filesystem::current_path().string() + "/structty_bench_" + current_timestamp() + "_" + in_file.stem().string() + ".csv");
        std::cout << "CWD = " << std::filesystem::current_path() << std::endl;
        screen.set_benchmark(&bm);
        t_load0 = Benchmark::clock::now();
    }

    screen.set_chainfile(params.get_chainfile(), params.get_in_file().size());
    for (int i = 0; i < params.get_in_file().size(); i++){
        screen.set_protein(params.get_in_file(i), i, params.get_show_structure());
    }
    screen.set_tmatrix();    
    if (params.get_utmatrix() != ""){
        screen.set_utmatrix(params.get_utmatrix(),0);
    }
    screen.normalize_proteins(params.get_utmatrix());
    screen.update_total_len_ca();

    // 기능 1: interface 모드일 때 inter-chain interface 계산 (threshold=8.0Å)
    if (params.get_mode() == "interface") {
        screen.compute_interface_all();
    }

    // 기능 4: aligned 모드 — -ut 만 있는 경우 (nearest-neighbor fallback)
    // -fs 있는 경우는 load_next_hit() 내부에서 U/T + aligned region 계산을 수행하므로 여기서 처리 안 함
    if (params.get_mode() == "aligned" && params.get_foldseek_file().empty()
        && params.get_foldmason_file().empty()) {
        screen.compute_aligned_all();
    }

    // 기능 3: Foldseek hit 탐색 설정 (-fs 파일이 있을 때)
    // load_next_hit() 내부에서 U/T 적용 + aligned region 계산 + align_method 설정 모두 수행
    if (!params.get_foldseek_file().empty()) {
        FoldseekParser fs_nav_parser;
        if (fs_nav_parser.load(params.get_foldseek_file()) && fs_nav_parser.hit_count() > 0) {
            screen.set_foldseek_hits(fs_nav_parser.get_hits());
            screen.set_fs_db_path(params.get_db_path());
            // 첫 번째 hit 자동 로드 (입력 파일이 1개일 때만 — 이미 target이 로드된 경우 제외)
            if ((int)params.get_in_file().size() <= 1) {
                screen.load_next_hit(+1);
            }
        }
    }

    // 기능 5: conservation 모드일 때 MSA 파일 로드 및 conservation score 계산
    if (params.get_mode() == "conservation" && !params.get_msa_file().empty()) {
        MSAParser msa_parser;
        if (msa_parser.load(params.get_msa_file())) {
            std::vector<float> scores = msa_parser.compute_conservation();
            screen.apply_msa_conservation(0, scores);
        } else {
            std::cerr << "Warning: Failed to load MSA file: " << params.get_msa_file() << std::endl;
        }
    }

    // 기능 8: FoldMason MSA 기반 superposition + conservation
    if (!params.get_foldmason_file().empty()) {
        auto fm_parser = std::make_unique<FoldMasonParser>();
        const std::string& fm_path = params.get_foldmason_file();
        bool fm_loaded = false;

        // 확장자 판별
        std::string ext;
        {
            size_t dot = fm_path.rfind('.');
            if (dot != std::string::npos) ext = fm_path.substr(dot + 1);
            for (auto& c : ext) c = (char)std::tolower((unsigned char)c);
        }

        if (ext == "json") {
            fm_loaded = fm_parser->load_json(fm_path);
        } else {
            fm_loaded = fm_parser->load_fasta(fm_path);
        }

        if (!fm_loaded || fm_parser->entry_count() == 0) {
            std::cerr << "Warning: Failed to load FoldMason file: " << fm_path << std::endl;
        } else {
            // entry 매칭: 파일명(확장자 제외) 기준으로 탐색, 실패 시 0/1 순서 가정
            int fm_query_idx = 0;
            int fm_target_idx = 1;
            const auto& entries = fm_parser->get_entries();
            if ((int)params.get_in_file().size() >= 1) {
                fs::path qp(params.get_in_file(0));
                std::string qstem = qp.stem().string();
                for (int ei = 0; ei < fm_parser->entry_count(); ei++) {
                    if (entries[ei].name.find(qstem) != std::string::npos) {
                        fm_query_idx = ei; break;
                    }
                }
            }
            if ((int)params.get_in_file().size() >= 2) {
                fs::path tp(params.get_in_file(1));
                std::string tstem = tp.stem().string();
                for (int ei = 0; ei < fm_parser->entry_count(); ei++) {
                    if (ei == fm_query_idx) continue;
                    if (entries[ei].name.find(tstem) != std::string::npos) {
                        fm_target_idx = ei; break;
                    }
                }
            }

            int total_entries = fm_parser->entry_count();

            // conservation 색상 (단일/다중 구조 모두)
            if (params.get_mode() == "conservation") {
                std::vector<float> entropy = fm_parser->compute_column_entropy(false);
                auto col_map = fm_parser->build_query_col_map(fm_query_idx);
                // col_map: 잔기 인덱스 → MSA 열 인덱스; 잔기 순서대로 entropy 값 매핑
                std::vector<float> mapped_scores(col_map.size(), 0.0f);
                for (int ri = 0; ri < (int)col_map.size(); ri++) {
                    int col = col_map[ri];
                    if (col < (int)entropy.size()) mapped_scores[ri] = entropy[col];
                }
                screen.apply_msa_conservation(0, mapped_scores);
            }

            // 패널 정보 구성
            FoldMasonInfo fm_info;
            fm_info.valid = true;
            fm_info.entry_count = total_entries;

            // 두 구조 superposition (두 번째 PDB 있을 때만)
            bool do_superposition = ((int)params.get_in_file().size() >= 2 && total_entries >= 2);
            if (do_superposition) {
                fm_info.align_method = (params.get_mode() == "aligned") ? "msa-col" : "-";
            } else {
                fm_info.align_method = "-";
            }

            screen.set_foldmason(std::move(fm_parser));
            screen.set_foldmason_panel_info(fm_info);

            if (do_superposition) {
                screen.apply_foldmason_superposition(0, 1, fm_query_idx, fm_target_idx);
            }
        }
    }
    
    if (bench) {
        auto t_load1 = Benchmark::clock::now();
        bm.log("load", -1, Benchmark::ms_since(t_load0, t_load1));

        const std::vector<int> script = {
            'X','Y','Z','A','D','W','S','R','F'
        };

        const int warmup = 200;   
        const int events = 2000; 

        // Warmup run (not measured)
        bool old_enabled = bm.enabled;
        bm.enabled = false;
        for (int i = 0; i < warmup; i++) {
            screen.draw_screen(params.get_no_panel());
            screen.handle_input(script[i % script.size()]);
        }
        bm.enabled = old_enabled;

        // Measured run
        bm.log("bench_begin", -1, 0.0);
        for (int i = 0; i < events; i++) {
            screen.draw_screen(params.get_no_panel());
            screen.handle_input(script[i % script.size()]);
        }
        bm.log("bench_end", -1, 0.0);
    }
    else{
        bool run = true;
        bool needs_redraw = true;
        while(run) {
            if (needs_redraw) {
                screen.draw_screen(params.get_no_panel());
            }
            // KEY_MOUSE 이벤트 시 needs_redraw=false: 패널 부분 갱신만 수행
            // 키보드 이벤트 시 needs_redraw=true: 전체 재렌더링
            run = screen.handle_input(needs_redraw);
        }
    }

    endwin();
    return 0;
}