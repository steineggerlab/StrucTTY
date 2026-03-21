#pragma once
#include "Curses.hpp"
#include "Protein.hpp"
#include "Atom.hpp"
#include "RenderPoint.hpp"
#include "Palette.hpp"
#include "Camera.hpp"
#include "Panel.hpp"
#include "Benchmark.hpp"
#include "../structure/FoldseekParser.hpp"
#include "../structure/PDBDownloader.hpp"
#include "../structure/FoldMasonParser.hpp"
#include <vector>
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <cstdlib>
#include <algorithm>   // clamp, max
#include <limits>      // numeric_limits
#include <memory>
#include <string>

class Screen {
public:
    Screen(const int& width, const int& height, const bool& show_structure,
           const std::string& mode, const std::string& depthcharacter);
    ~Screen();

    // 기능 6: 마우스 hover — main loop용 (needs_redraw 반환)
    bool handle_input(bool& needs_redraw);
    // 벤치마크용: needs_redraw 없이 키 직접 처리
    bool handle_input(int key);
    char get_pixel_char_from_depth(float z, float min_z, float max_z);

    void set_protein(const std::string& in_file, int ii, const bool& show_structure);
    void normalize_proteins(const std::string& utmatrix);

    void set_tmatrix();
    void set_utmatrix(const std::string& utmatrix, bool onlyU);
    void set_chainfile(const std::string& chainfile, int filesize);
    void set_zoom_level(float zoom);

    // 기능 1: 로드된 모든 Protein에 대해 inter-chain interface를 계산
    void compute_interface_all(float threshold = 8.0f);

    // 기능 4: 로드된 모든 Protein 쌍에 대해 nearest-neighbor 기반 정렬 잔기를 계산
    void compute_aligned_all(float threshold = 10.0f);

    // 기능 5: 지정 protein에 conservation scores 적용
    void apply_msa_conservation(int protein_idx, const std::vector<float>& scores);

    // 기능 4: -fs 기반 — Foldseek hit의 U/T transform을 지정 protein에 적용
    // U_flat: row-major 3x3 (9 elements), T: translation (3 elements)
    void apply_foldseek_transform(int protein_idx, const float* U_flat, const float* T);

    // 기능 4: -fs 기반 — alignment string으로 aligned 잔기 계산 (protein0 vs protein1)
    void compute_aligned_from_aln(const std::string& qaln, const std::string& taln,
                                  float threshold = 5.0f);

    // 기능 4: 패널에 정렬 방식 표시 설정 ("aln-string" or "nearest-nbr")
    void set_align_method(const std::string& method);

    // 기능 3: Foldseek hit 탐색
    void set_foldseek_hits(const std::vector<FoldseekHit>& hits);
    void set_fs_db_path(const std::string& path);
    void load_next_hit(int delta);  // delta=+1: next, delta=-1: prev, delta=0: first

    // 기능 8: FoldMason MSA 기반 superposition + aligned region 설정
    void set_foldmason(std::unique_ptr<FoldMasonParser> parser);
    void apply_foldmason_superposition(int query_protein_idx, int target_protein_idx,
                                       int fm_query_entry_idx, int fm_target_entry_idx);
    void set_foldmason_panel_info(const FoldMasonInfo& info);

    void draw_screen(bool no_panel);

    // 기능 6: 마우스 hover — 현재 커서 위치의 잔기 정보를 패널에 반영
    void update_hover_info(int mx, int my);
    void init_color_pairs();
    void assign_colors_to_points(std::vector<RenderPoint>& points, int protein_idx);

    void draw_line(std::vector<RenderPoint>& points,
                   int x1, int x2,
                   int y1, int y2,
                   float z1, float z2,
                   const std::string& chainID, char structure,
                   float min_z, float max_z,
                   int max_x = -1, int max_y = -1, int half = 0);
    
    void set_benchmark(Benchmark* b) { bm = b; }
    
    void update_total_len_ca() {
        int64_t sum = 0;
        for (auto* p : data) {
            if (!p) continue;
            sum += (int64_t)p->get_length();  // Cα-only length라고 가정
        }
        total_len_ca = sum;
    }

private:
    int screen_width, screen_height;
    float aspect_ratio;
    bool screen_show_structure;
    bool yesUT = false;
    std::string screen_mode;
    std::string screen_depthcharacter;
    int structNum = -1;

    // Braille sub-pixel rendering
    bool use_braille = true;
    std::vector<RenderPoint> logicalPixels; // 2*width x 4*height logical buffer
    void print_screen_braille(int y_offset);

    float focal_offset = 3.0f;
    float zoom_level;

    std::vector<float> pan_x;
    std::vector<float> pan_y;
    std::vector<std::string> chainVec;
    float** vectorpointer = nullptr;

    std::vector<RenderPoint> screenPixels;
    std::vector<Protein*> data;

    BoundingBox global_bb;
    Camera* camera = nullptr;
    Panel* panel = nullptr;

    bool depth_calibrated = false;
    float depth_base_min_z = 0.0f;
    float depth_base_max_z = 1.0f;

    Benchmark* bm = nullptr;
    bool ttff_logged = false;

    // 기능 6: 마우스 hover — 패널 위치 정보 (hover 갱신 시 필요)
    int last_panel_h = 0;
    int last_panel_start_row = 0;
    int last_panel_cols = 0;
    bool last_no_panel = false;

    // 기능 3: Foldseek hit 탐색
    std::vector<FoldseekHit> foldseek_hits;
    int current_hit_idx = -1;
    std::string fs_db_path;

    // 기능 8: FoldMason MSA
    std::unique_ptr<FoldMasonParser> foldmason_parser;
    float norm_scale = 1.0f;   // normalize_proteins() 에서 저장
    float norm_cx = 0.0f;      // query protein 정규화 전 centroid
    float norm_cy = 0.0f;
    float norm_cz = 0.0f;

    void calibrate_depth_baseline_first_view();

    void project();
    void project(std::vector<RenderPoint>& screenshotPixels, const int proj_width, const int proj_height);
    void clear_screen();
    void print_screen(int panel_lines);

    // 기능 6: handle_input 내부 구현 (두 public 오버로드 공유)
    bool handle_input_impl(int key, bool& needs_redraw);

    int64_t total_len_ca = 0;
};