#include <iostream>
#include <unistd.h> 
#include <filesystem>
#include "Curses.hpp"
#include "Common.hpp"
#include "Protein.hpp"
#include "Parameters.hpp"
#include "Screen.hpp"
#include "Benchmark.hpp" 

int main(int argc, char* argv[]) {
    Parameters params(argc, argv);
    if (!params.check_arg_okay()) return -1; 
    params.print_args();
    
    initscr();
    cbreak();
    noecho();
    
    Screen screen(params.get_width(), params.get_height(), 
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
        while(run) {
            screen.draw_screen(params.get_no_panel());
            run = screen.handle_input();
        }
    }

    endwin();
    return 0;
}