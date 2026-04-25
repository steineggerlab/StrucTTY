#pragma once
#include <vector>
#include <string>
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
    std::string chain_id;
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

private:
    int         width_, height_;
    std::string mode_;
    bool        show_structure_;
    float       focal_offset_     = 3.0f;
    float       zoom_level_       = 2.0f;
    float       depth_base_min_z_ = 0.0f;
    float       depth_base_max_z_ = 1.0f;

    std::vector<RenderPoint> logical_pixels_;

    static constexpr float FOV_ = 90.0f;
    static constexpr float PI_  = 3.14159265359f;

    void clear();
    int  compute_depth_band(float z) const;
    void draw_line_impl(std::vector<RenderPoint>& out,
                        int x1, int x2, int y1, int y2,
                        float z1, float z2,
                        const std::string& chain_id, char structure,
                        int max_x, int max_y, int half) const;
    void project_and_fill(const std::vector<RenderAtom>& atoms,
                          std::vector<RenderPoint>& out) const;
    void assign_colors_impl(std::vector<RenderPoint>& points, int protein_idx) const;
    void zbuffer_resolve(const std::vector<RenderPoint>& points);
};
