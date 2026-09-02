#pragma once
#include "Palette.hpp"
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
    int color_group = -1;
    int chain_color_base = -1;
};

struct FoldMasonInfo {
    bool        valid        = false;
    int         entry_count  = 0;
    std::string align_method;
};

struct FoldseekHitInfo {
    bool     valid        = false;
    int      current_idx  = 0;
    int      total_hits   = 0;
    int      query_idx    = 0;
    int      total_queries = 0;
    std::string target;
    float    evalue       = 0.0f;
    float    prob         = -1.0f;
    float    lddt         = -1.0f;
    float    qtmscore     = -1.0f;
    std::string align_method;
    std::string status_msg;
    bool     multimer     = false;
    float    ttmscore     = -1.0f;
    float    qcov         = -1.0f;
    float    tcov         = -1.0f;
    std::string interface_lddt;
    int      ass_id       = -1;
};

class Panel {
public:
    Panel(int width, const std::string& mode, bool show_structure = false);

    void add_panel_info(const std::string& file_name,
                        const std::map<std::string, int>& chain_info,
                        const std::map<std::string, int>& chain_residue_info,
                        int color_group = -1, int chain_color_base = -1);

    void reset_entries() { entries.clear(); }

    void update_entry(int idx, const std::string& file_name,
                      const std::map<std::string, int>& chain_info,
                      const std::map<std::string, int>& chain_residue_info);

    int get_height() const;
    int get_height_for_width(int max_cols, int compact_level = 0) const;

    void draw_panel(int start_row, int start_col,
                    int max_rows, int max_cols,
                    int compact_level = 0) const;

    void set_align_method(const std::string& method);

    void set_foldseek_hit_info(const FoldseekHitInfo& info);
    void clear_foldseek_hit_info();
    int  get_foldseek_section_height() const;

    void set_hover_residue(const std::string& chainID,
                           const char* residue_name,
                           int residue_number,
                           char structure,
                           float bfactor,
                           float conservation_score);
    void clear_hover_residue();

    int get_residue_section_height() const;

    int get_last_hover_row() const;

    void draw_hover_section(int hover_start_row, int max_cols) const;

    void set_foldmason_info(const FoldMasonInfo& info);
    void clear_foldmason_info();
    int  get_foldmason_section_height() const;

private:
    std::vector<Entry> entries;

    int panel_width;
    std::string panel_mode;
    bool panel_show_structure = false;
    std::string align_method;

    FoldseekHitInfo fs_hit_info;

    FoldMasonInfo fm_info;

    mutable int last_hover_row = -1;
    bool hover_valid = false;
    std::string hover_chain;
    char hover_residue_name[4] = {};
    int hover_residue_number = -1;
    char hover_structure = 0;
    float hover_bfactor = 0.0f;
    float hover_conservation = -1.0f;
};
