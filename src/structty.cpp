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

    // 기능 4: aligned 모드
    // -fs 있으면 alignment string 기반, 없으면 nearest-neighbor 기반
    if (params.get_mode() == "aligned") {
        bool handled_by_fs = false;
        if (!params.get_foldseek_file().empty()) {
            FoldseekParser fs_parser;
            if (fs_parser.load(params.get_foldseek_file()) && fs_parser.hit_count() > 0) {
                const FoldseekHit& hit = fs_parser.get_hits()[0];
                // hit의 U/T transform을 protein1(target)에 적용
                if (hit.has_transform) {
                    screen.apply_foldseek_transform(1, hit.U, hit.T);
                }
                if (hit.has_aln) {
                    screen.compute_aligned_from_aln(hit.qaln, hit.taln, 5.0f);
                    screen.set_align_method("aln-string");
                } else {
                    screen.compute_aligned_all();
                    screen.set_align_method("nearest-nbr");
                }
                handled_by_fs = true;
            }
        }
        if (!handled_by_fs) {
            screen.compute_aligned_all();
            // compute_aligned_all() 내부에서 "nearest-nbr" 설정
        }
    }

    // 기능 3: Foldseek hit 탐색 설정 (-fs 파일이 있을 때)
    // aligned 모드의 첫 hit 로드와 별개로, hit 목록 전체를 Screen에 전달
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