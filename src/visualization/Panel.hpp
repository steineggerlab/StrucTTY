#pragma once
#include "Palette.hpp"
#include "Curses.hpp"
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <iostream>
#include <cstring>

struct Entry {
    std::string file_name;
    std::map<std::string, int> chain_atom_info;
    std::map<std::string, int> chain_residue_info;
};

class Panel {
public:
    Panel(int width, const std::string& mode);

    void add_panel_info(const std::string& file_name,
                        const std::map<std::string, int>& chain_info,
                        const std::map<std::string, int>& chain_residue_info);

    int get_height() const;
    int get_height_for_width(int max_cols) const; 

    void draw_panel(int start_row, int start_col,
                    int max_rows, int max_cols) const;

private:
    std::vector<Entry> entries;

    int panel_width;
    std::string panel_mode;
};
