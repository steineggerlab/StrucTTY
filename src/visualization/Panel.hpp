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
    Panel(int width, const std::string& mode, bool show_structure = false);

    void add_panel_info(const std::string& file_name,
                        const std::map<std::string, int>& chain_info,
                        const std::map<std::string, int>& chain_residue_info);

    int get_height() const;
    int get_height_for_width(int max_cols) const;

    void draw_panel(int start_row, int start_col,
                    int max_rows, int max_cols) const;

    // 기능 4: 정렬 방식 표시 ("nearest-nbr" or "aln-string")
    void set_align_method(const std::string& method);

    // 기능 6: 마우스 hover — Residue Info 섹션
    // chainID, residue_name(char[4]), residue_number, structure, bfactor, conservation_score
    void set_hover_residue(const std::string& chainID,
                           const char* residue_name,
                           int residue_number,
                           char structure,
                           float bfactor,
                           float conservation_score);
    void clear_hover_residue();

    // Residue Info 섹션의 줄 수 (항상 고정)
    int get_residue_section_height() const;

    // hover 섹션만 부분 갱신 (hover_start_row = "Residue Info" 헤더 행)
    void draw_hover_section(int hover_start_row, int max_cols) const;

private:
    std::vector<Entry> entries;

    int panel_width;
    std::string panel_mode;
    bool panel_show_structure = false;
    std::string align_method;

    // 기능 6: hover 상태
    bool hover_valid = false;
    std::string hover_chain;
    char hover_residue_name[4] = {};
    int hover_residue_number = -1;
    char hover_structure = 0;
    float hover_bfactor = 0.0f;
    float hover_conservation = -1.0f;
};
