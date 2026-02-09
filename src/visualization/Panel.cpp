#include "Panel.hpp"

Panel::Panel(int width, const std::string& mode) : panel_width(width), panel_mode(mode) {}

void Panel::add_panel_info(const std::string& file_name, 
                           const std::map<std::string, int>& chain_info, 
                           const std::map<std::string, int>& chain_residue_info) {
    entries.push_back(Entry{
        file_name,
        chain_info,
        chain_residue_info
    });
}

int Panel::get_height() const {
    int lines = 0;
    lines += 3; 
    for (const auto& entry : entries) {
        lines += 1; 
        int n = (int)entry.chain_atom_info.size();
        int chain_lines = (n == 0) ? 1 : ((n + 2) / 3); // 3 per line
        lines += chain_lines;
        lines += 1; 
    }
    lines += 1; 
    return lines;
}
void Panel::draw_panel(int start_row, int start_col,
                       int max_rows, int max_cols) const {
    const int num_colors = (int)(sizeof(Palettes::UNRAINBOW) / sizeof(int));
    if (max_rows <= 0 || max_cols <= 0) return;

    const int top    = start_row;
    const int left   = start_col;
    const int bottom = start_row + max_rows; // exclusive
    const int right  = start_col + max_cols; // exclusive

    const int right_limit = right - 1;

    auto in_rows = [&](int rr){ return rr >= top && rr < bottom; };
    auto remain_cols = [&](int x){ return right_limit - x; };

    auto clear_line = [&](int rr){
        if (!in_rows(rr)) return;
        move(rr, left);
        clrtoeol();
        move(rr, left);
    };

    auto put_n = [&](int& rr, int& x, const char* s, int n){
        if (!in_rows(rr)) return;
        int rem = remain_cols(x);
        if (rem <= 0 || n <= 0) return;
        int k = std::min(rem, n);
        addnstr(s, k);
        x += k;
    };

    auto put_str = [&](int& rr, int& x, const std::string& s){
        put_n(rr, x, s.c_str(), (int)s.size());
    };

    auto put_cstr = [&](int& rr, int& x, const char* s){
        put_n(rr, x, s, (int)std::strlen(s));
    };

    auto put_indent = [&](int& rr, int& x){
        put_str(rr, x, "  ");
    };

    int r = start_row;

    // Top border
    clear_line(r);
    {
        int x = left;
        put_str(r, x, "*");
        int mid = std::max(0, std::min(panel_width - 2, max_cols - 2));
        put_str(r, x, std::string(mid, '='));
        if (max_cols >= 2) put_str(r, x, "*");
    }
    ++r;
    if (!in_rows(r)) return;

    // Help line
    clear_line(r);
    {
        int x = left;
        put_str(r, x, "W A S D : ^ < v >  ");
        put_str(r, x, "R F : Zoom In/Out  ");
        put_str(r, x, "X Y Z : Rotate X, Y, Z axis  ");
        put_str(r, x, "C : Screenshot  ");
        put_str(r, x, "Q : Quit");
    }
    ++r;
    if (!in_rows(r)) return;

    // Separator
    clear_line(r);
    {
        move(r, left);
        int w = std::min(panel_width, max_cols);
        w = std::min(w, max_cols - 1);
        for (int i = 0; i < w; ++i) addch('-');
    }
    ++r;
    if (!in_rows(r)) return;

    // Body
    int file_idx = 0;
    for (const auto& entry : entries) {
        if (!in_rows(r)) break;

        const std::string& file_name = entry.file_name;
        const auto& chain_info       = entry.chain_atom_info;

        int protein_pair = 0;
        if (panel_mode == "protein" && num_colors > 0) {
            protein_pair = (file_idx % num_colors) + 1;
        }

        // file name line
        clear_line(r);
        {
            int x = left;
            if (protein_pair > 0) attron(COLOR_PAIR(protein_pair));
            put_str(r, x, file_name);
            if (protein_pair > 0) attroff(COLOR_PAIR(protein_pair));
        }
        ++r;
        if (!in_rows(r)) break;

        // chain lines
        clear_line(r);
        move(r, left);
        int x = left;
        put_indent(r, x);

        int count = 0;
        for (const auto& [chainID, length] : chain_info) {
            if (!in_rows(r)) break;

            int residue_cnt = 0;
            auto itC = entry.chain_residue_info.find(chainID);
            if (itC != entry.chain_residue_info.end()) residue_cnt = itC->second;

            char buf[64];
            int token_len = std::snprintf(buf, sizeof(buf), "%s: %d (%d)  ",
                                          chainID.c_str(), residue_cnt, length);
            if (token_len < 0) token_len = 0;
            if (token_len >= (int)sizeof(buf)) token_len = (int)sizeof(buf) - 1;

            if (x + token_len > right_limit) {
                ++r;
                if (!in_rows(r)) break;
                clear_line(r);
                move(r, left);
                x = left;
                put_indent(r, x);
            }

            int chain_pair = 0;
            if (panel_mode == "chain" && num_colors > 0) {
                chain_pair = (file_idx * 10 + (count % num_colors)) + 1;
            }

            int pair_to_use = (panel_mode == "protein") ? protein_pair : chain_pair;

            if (pair_to_use > 0) attron(COLOR_PAIR(pair_to_use));
            put_n(r, x, buf, token_len);
            if (pair_to_use > 0) attroff(COLOR_PAIR(pair_to_use));

            ++count;
        }

        // blank line (2 rows)
        ++r;
        if (!in_rows(r)) break;
        clear_line(r);
        ++r;

        ++file_idx;
    }

    if (!in_rows(r)) return;

    // Bottom border
    clear_line(r);
    {
        int x = left;
        put_str(r, x, "*");
        int mid = std::max(0, std::min(panel_width - 2, max_cols - 2));
        put_str(r, x, std::string(mid, '='));
        if (max_cols >= 2) put_str(r, x, "*");
    }
}

static int count_wrapped_lines_for_chaininfo(
    const std::map<std::string,int>& chain_info,
    const std::map<std::string,int>& chain_residue_info,
    int avail_cols,
    int indent_cols = 2
){
    if (avail_cols <= indent_cols) return 1;

    int lines = 1;
    int x = indent_cols;

    for (const auto& [chainID, length] : chain_info) {
        int residue_cnt = 0;
        auto itC = chain_residue_info.find(chainID);
        if (itC != chain_residue_info.end()) residue_cnt = itC->second;

        char buf[128];
        std::snprintf(buf, sizeof(buf), "%s: %d (%d)  ", chainID.c_str(), residue_cnt, length);

        int token_len = (int)std::strlen(buf);

        if (token_len > avail_cols - indent_cols) {
            if (x > indent_cols) { lines++; x = indent_cols; }
            x = indent_cols + std::min(token_len, avail_cols - indent_cols);
            continue;
        }

        if (x + token_len > avail_cols) {
            lines++;
            x = indent_cols;
        }
        x += token_len;
    }
    return lines;
}
int Panel::get_height_for_width(int max_cols) const {
    int lines = 0;

    lines += 3; // Top border + Help line + Separator

    int avail_cols = max_cols;
    if (avail_cols < 1) avail_cols = 1;

    for (const auto& entry : entries) {
        lines += 1; // file name line

        lines += count_wrapped_lines_for_chaininfo(
            entry.chain_atom_info, entry.chain_residue_info,
            /*avail_cols=*/avail_cols,
            /*indent_cols=*/2
        );

        lines += 2;
    }

    lines += 1; // Bottom border
    return lines;
}

