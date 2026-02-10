#pragma once
#include "Protein.hpp"
#include "Atom.hpp"
#include "RenderPoint.hpp"
#include "Palette.hpp"
#include "Camera.hpp"
#include "Panel.hpp"
#include "Benchmark.hpp"
#include <vector>
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <ncurses.h>
#include <cstdlib>
#include <algorithm>   // clamp, max
#include <limits>      // numeric_limits
#include <string>

class Screen {
public:
    Screen(const int& width, const int& height, const bool& show_structure,
           const std::string& mode, const std::string& depthcharacter);
    ~Screen();

    bool handle_input();
    char get_pixel_char_from_depth(float z, float min_z, float max_z);

    void set_protein(const std::string& in_file, int ii, const bool& show_structure);
    void normalize_proteins(const std::string& utmatrix);

    void set_tmatrix();
    void set_utmatrix(const std::string& utmatrix, bool onlyU);
    void set_chainfile(const std::string& chainfile, int filesize);
    void set_zoom_level(float zoom);

    void draw_screen(bool no_panel);
    void init_color_pairs();
    void assign_colors_to_points(std::vector<RenderPoint>& points, int protein_idx);

    void draw_line(std::vector<RenderPoint>& points,
                   int x1, int x2,
                   int y1, int y2,
                   float z1, float z2,
                   std::string chainID, char structure,
                   float min_z, float max_z);
    
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

    void calibrate_depth_baseline_first_view();

    void project();
    void project(std::vector<RenderPoint>& screenshotPixels, const int proj_width, const int proj_height);
    void clear_screen();
    void print_screen(int panel_lines);

    int64_t total_len_ca = 0;
};