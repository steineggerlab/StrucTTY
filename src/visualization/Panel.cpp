#include "Panel.hpp"
#include <cstring>  // strncpy
#include <cstdio>   // snprintf

Panel::Panel(int width, const std::string& mode, bool show_structure)
    : panel_width(width), panel_mode(mode), panel_show_structure(show_structure) {}

void Panel::add_panel_info(const std::string& file_name, 
                           const std::map<std::string, int>& chain_info, 
                           const std::map<std::string, int>& chain_residue_info) {
    entries.push_back(Entry{
        file_name,
        chain_info,
        chain_residue_info
    });
}

void Panel::set_align_method(const std::string& method) {
    align_method = method;
}

// 기능 3: entry 갱신 (target protein 교체 시)
void Panel::update_entry(int idx, const std::string& file_name,
                          const std::map<std::string, int>& chain_info,
                          const std::map<std::string, int>& chain_residue_info) {
    if (idx < 0) return;
    if (idx < (int)entries.size()) {
        entries[idx] = Entry{file_name, chain_info, chain_residue_info};
    } else {
        while ((int)entries.size() < idx) {
            entries.push_back(Entry{});
        }
        entries.push_back(Entry{file_name, chain_info, chain_residue_info});
    }
}

// 기능 3: Foldseek hit info 섹션
void Panel::set_foldseek_hit_info(const FoldseekHitInfo& info) {
    fs_hit_info = info;
}

void Panel::clear_foldseek_hit_info() {
    fs_hit_info = FoldseekHitInfo{};
}

int Panel::get_foldseek_section_height() const {
    if (!fs_hit_info.valid && fs_hit_info.total_hits == 0) return 0;
    return 8;  // 고정 8줄
}

// 기능 6: Residue Info hover -----------------------------------------------

void Panel::set_hover_residue(const std::string& chainID,
                               const char* residue_name,
                               int residue_number,
                               char structure,
                               float bfactor,
                               float conservation_score) {
    hover_valid          = true;
    hover_chain          = chainID;
    strncpy(hover_residue_name, residue_name, 3);
    hover_residue_name[3] = '\0';
    hover_residue_number  = residue_number;
    hover_structure       = structure;
    hover_bfactor         = bfactor;
    hover_conservation    = conservation_score;
}

void Panel::clear_hover_residue() {
    hover_valid = false;
}

int Panel::get_residue_section_height() const {
    // 항상 고정: header + chain + res + ss (= 4 기본)
    // + pLDDT 줄 (plddt 모드일 때) + Cons 줄 (conservation 모드일 때)
    int h = 4;
    if (panel_mode == "plddt")        h += 1;
    if (panel_mode == "conservation") h += 1;
    return h;
}

void Panel::draw_hover_section(int hover_start_row, int max_cols) const {
    // hover_start_row: "Residue Info" 헤더 행 (separator는 이미 draw_panel이 그렸으므로 제외)
    // 이 함수는 Residue Info 섹션만 다시 그린다 (bottom border 제외)
    int r = hover_start_row;
    int right_limit = max_cols - 1;

    auto clear_ln = [&](int rr) {
        move(rr, 0);
        clrtoeol();
        move(rr, 0);
    };

    auto put_text = [&](int rr, const char* s) {
        int len = (int)std::strlen(s);
        int k   = std::min(len, right_limit);
        if (k > 0) addnstr(s, k);
    };

    // "Residue Info" 헤더
    clear_ln(r);
    put_text(r, "Residue Info");
    ++r;

    // Chain: X  or  Chain: -
    clear_ln(r);
    if (hover_valid) {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "Chain: %s", hover_chain.empty() ? "-" : hover_chain.c_str());
        put_text(r, buf);
    } else {
        put_text(r, "Chain: -");
    }
    ++r;

    // Res: GLU 42  or  Res: -
    clear_ln(r);
    if (hover_valid && hover_residue_number >= 0) {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "Res:   %s %d",
                      (hover_residue_name[0] ? hover_residue_name : "?"),
                      hover_residue_number);
        put_text(r, buf);
    } else {
        put_text(r, "Res:   -");
    }
    ++r;

    // SS: Helix / Sheet / Coil  or  SS: -
    clear_ln(r);
    if (hover_valid) {
        const char* ss_str = "Coil";
        if      (hover_structure == 'H') ss_str = "Helix";
        else if (hover_structure == 'S') ss_str = "Sheet";
        char buf[32];
        std::snprintf(buf, sizeof(buf), "SS:    %s", ss_str);
        put_text(r, buf);
    } else {
        put_text(r, "SS:    -");
    }
    ++r;

    // pLDDT: 87.3  (plddt 모드일 때만)
    if (panel_mode == "plddt") {
        clear_ln(r);
        if (hover_valid) {
            char buf[32];
            std::snprintf(buf, sizeof(buf), "pLDDT: %.1f", hover_bfactor);
            put_text(r, buf);
        } else {
            put_text(r, "pLDDT: -");
        }
        ++r;
    }

    // Cons: 0.82  (conservation 모드일 때만)
    if (panel_mode == "conservation") {
        clear_ln(r);
        if (hover_valid && hover_conservation >= 0.0f) {
            char buf[32];
            std::snprintf(buf, sizeof(buf), "Cons:  %.2f", hover_conservation);
            put_text(r, buf);
        } else {
            put_text(r, "Cons:  -");
        }
        ++r;
    }
}

// 기능 6 끝 -----------------------------------------------------------------

int Panel::get_height() const {
    int lines = 0;
    lines += 3;
    if (panel_mode == "aligned" && !align_method.empty()) {
        lines += 1;  // "Align: nearest-nbr" or "Align: aln-string"
    }
    for (const auto& entry : entries) {
        lines += 1;
        int n = (int)entry.chain_atom_info.size();
        int chain_lines = (n == 0) ? 1 : ((n + 2) / 3); // 3 per line
        lines += chain_lines;
        lines += 1;
    }
    // 기능 3: Foldseek hit info 섹션
    int fs_h = get_foldseek_section_height();
    if (fs_h > 0) {
        lines += 1;   // separator
        lines += fs_h;
    }
    // 기능 6: separator + Residue Info 섹션
    lines += 1;                             // separator (---)
    lines += get_residue_section_height();  // 고정 줄 수
    lines += 1;                             // bottom border
    return lines;
}
void Panel::draw_panel(int start_row, int start_col,
                       int max_rows, int max_cols) const {
    const int num_protein_colors = 9;
    const int num_chain_colors   = 15;
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

    // 기능 4: aligned 모드일 때 정렬 방식 표시
    if (panel_mode == "aligned" && !align_method.empty()) {
        if (!in_rows(r)) return;
        clear_line(r);
        int x = left;
        put_cstr(r, x, "Align: ");
        put_str(r, x, align_method);
        ++r;
        if (!in_rows(r)) return;
    }

    // Body
    int file_idx = 0;
    for (const auto& entry : entries) {
        if (!in_rows(r)) break;

        const std::string& file_name = entry.file_name;
        const auto& chain_info       = entry.chain_atom_info;

        int protein_pair = 0;
        if (panel_mode == "protein") {
            protein_pair = (file_idx % num_protein_colors) + 1;  // pairs 1-9
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
            if (panel_mode == "chain") {
                chain_pair = 21 + ((file_idx * 10 + count) % num_chain_colors);  // pairs 21-35
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

    // 기능 3: Foldseek hit info 섹션
    {
        int fs_h = get_foldseek_section_height();
        if (fs_h > 0) {
            // separator
            clear_line(r);
            {
                move(r, left);
                int w = std::min(panel_width, max_cols);
                w = std::min(w, max_cols - 1);
                for (int i = 0; i < w; ++i) addch('-');
            }
            ++r;
            if (!in_rows(r)) return;

            const FoldseekHitInfo& fi = fs_hit_info;

            // Line 1: "Foldseek Hits"
            clear_line(r);
            { int x = left; put_cstr(r, x, "Foldseek Hits"); }
            ++r; if (!in_rows(r)) return;

            // Line 2: "[X / Y]" or hit count
            clear_line(r);
            {
                int x = left;
                char buf[32];
                std::snprintf(buf, sizeof(buf), "[%d / %d]", fi.current_idx, fi.total_hits);
                put_cstr(r, x, buf);
            }
            ++r; if (!in_rows(r)) return;

            // Line 3: "Target: ..."
            clear_line(r);
            {
                int x = left;
                put_cstr(r, x, "Target: ");
                if (!fi.target.empty()) put_str(r, x, fi.target);
                else put_cstr(r, x, "-");
            }
            ++r; if (!in_rows(r)) return;

            // Line 4: "E-val:  ..."
            clear_line(r);
            {
                int x = left;
                char buf[32];
                std::snprintf(buf, sizeof(buf), "E-val:  %.2e", fi.evalue);
                put_cstr(r, x, buf);
            }
            ++r; if (!in_rows(r)) return;

            // Line 5: prob or status_msg
            clear_line(r);
            {
                int x = left;
                if (!fi.status_msg.empty()) {
                    put_str(r, x, fi.status_msg);
                } else if (fi.prob >= 0.0f) {
                    char buf[32];
                    std::snprintf(buf, sizeof(buf), "Prob:   %.3f", fi.prob);
                    put_cstr(r, x, buf);
                } else if (fi.qtmscore >= 0.0f) {
                    char buf[32];
                    std::snprintf(buf, sizeof(buf), "TM:     %.3f", fi.qtmscore);
                    put_cstr(r, x, buf);
                } else {
                    put_cstr(r, x, "-");
                }
            }
            ++r; if (!in_rows(r)) return;

            // Line 6: lDDT or second TM score
            clear_line(r);
            {
                int x = left;
                if (fi.lddt >= 0.0f) {
                    char buf[32];
                    std::snprintf(buf, sizeof(buf), "lDDT:   %.3f", fi.lddt);
                    put_cstr(r, x, buf);
                } else if (fi.qtmscore >= 0.0f && fi.prob < 0.0f) {
                    // 29-col 포맷, TM already shown above only if no prob
                    put_cstr(r, x, "-");
                } else {
                    put_cstr(r, x, "-");
                }
            }
            ++r; if (!in_rows(r)) return;

            // Line 7: Align method
            clear_line(r);
            {
                int x = left;
                if (!fi.align_method.empty()) {
                    put_cstr(r, x, "Align:  ");
                    put_str(r, x, fi.align_method);
                } else {
                    put_cstr(r, x, "Align:  -");
                }
            }
            ++r; if (!in_rows(r)) return;

            // Line 8: nav hint
            clear_line(r);
            {
                int x = left;
                put_cstr(r, x, "[N]ext  [P]rev");
            }
            ++r; if (!in_rows(r)) return;
        }
    }

    // 기능 6: Residue Info 섹션 separator
    clear_line(r);
    {
        move(r, left);
        int w = std::min(panel_width, max_cols);
        w = std::min(w, max_cols - 1);
        for (int i = 0; i < w; ++i) addch('-');
    }
    ++r;
    if (!in_rows(r)) return;

    // 기능 6: Residue Info 섹션 (draw_hover_section과 동일한 내용)
    // draw_panel()은 전체 패널을 그리므로 여기서도 Residue Info를 그린다.
    // hover가 없을 때는 모두 "-" 표시.
    draw_hover_section(r, max_cols);
    r += get_residue_section_height();

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

    if (panel_mode == "aligned" && !align_method.empty()) {
        lines += 1;  // "Align: ..." line
    }

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

    // 기능 3: Foldseek hit info 섹션
    int fs_h = get_foldseek_section_height();
    if (fs_h > 0) {
        lines += 1;   // separator
        lines += fs_h;
    }
    // 기능 6: separator + Residue Info 섹션 + bottom border
    lines += 1;                             // separator (---)
    lines += get_residue_section_height();  // 고정 줄 수
    lines += 1;                             // Bottom border
    return lines;
}

