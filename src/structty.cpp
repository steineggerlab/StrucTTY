#include <iostream>
#include <unistd.h> 
#include "Curses.hpp"
#include "Common.hpp"
#include "Protein.hpp"
#include "Parameters.hpp"
#include "Screen.hpp"
#include "Benchmark.hpp" 

int main(int argc, char* argv[]) {
    Parameters params(argc, argv);

    if (!params.check_arg_okay()) {
        return -1; 
    }
    params.print_args();

    Benchmark bm;
    bm.start("structty_bench_" + current_timestamp() + ".csv");
    initscr();
    cbreak();
    noecho();
    
    Screen screen(params.get_width(), params.get_height(), params.get_show_structure(), params.get_mode(), params.get_depthcharacter()); 
    screen.set_benchmark(&bm);
    auto t_load0 = Benchmark::clock::now();
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

    auto t_load1 = Benchmark::clock::now();
    bm.log("load", -1, Benchmark::ms_since(t_load0, t_load1));

    bool run = true;
    while(run) {
        screen.draw_screen(params.get_no_panel());
        run = screen.handle_input();
    }

    endwin();
    return 0;
}