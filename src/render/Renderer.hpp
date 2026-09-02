#pragma once
#include <vector>
#include <string>
#include <chrono>
#include <cstdint>
#include "RenderPoint.hpp"

struct RenderAtom {
    float x, y, z;
    char  structure;
    float bfactor;
    bool  is_interface;
    bool  is_aligned;
    float conservation_score;
    int   residue_number;
    char  residue_name[4];
    int   chain_id = -1;
    int   chain_color_id = -1;
    int   protein_index;
    float pan_x, pan_y;
};

class Renderer {
public:
    Renderer(int width, int height, const std::string& mode, bool show_structure);

    void set_depth_params(float focal_offset, float zoom_level,
                          float depth_min_z, float depth_max_z);
    void render(const std::vector<RenderAtom>& atoms);

    const std::vector<RenderPoint>& get_pixels() const;
    int get_logical_width()  const;
    int get_logical_height() const;

    size_t get_last_point_count() const { return last_point_count_; }

    int64_t get_last_us_clear() const { return last_us_clear_; }
    int64_t get_last_us_fill()  const { return last_us_fill_;  }
    int64_t get_last_us_zbuf()  const { return last_us_zbuf_;  }

private:
    enum class Mode { Protein, Chain, Rainbow, Plddt, Interface, Aligned, Conservation };

    struct PlotStyle {
        const RenderAtom* a;
        int   protein_idx;
        int   chain_color_idx;
        float rainbow_frac;
    };

    int         width_, height_;
    std::string mode_;
    Mode        mode_id_ = Mode::Protein;
    bool        show_structure_;
    float       focal_offset_     = 3.0f;
    float       zoom_level_       = 2.0f;
    float       depth_base_min_z_ = 0.0f;
    float       depth_base_max_z_ = 1.0f;

    std::vector<RenderPoint> logical_pixels_;
    size_t                   last_point_count_ = 0;
    std::chrono::steady_clock::time_point t_enter_;
    int64_t                  last_us_clear_ = 0;
    int64_t                  last_us_fill_  = 0;
    int64_t                  last_us_zbuf_  = 0;

    static constexpr float FOV_ = 90.0f;
    static constexpr float PI_  = 3.14159265359f;

    void clear();
    int  compute_depth_band(float z) const;
    int  color_from_style(const PlotStyle& s, int band) const;
    void plot(int x, int y, float depth, const PlotStyle& s);
    void draw_line_impl(int x1, int x2, int y1, int y2,
                        float z1, float z2,
                        const PlotStyle& s, int half);
    void project_and_fill(const std::vector<RenderAtom>& atoms);
};
