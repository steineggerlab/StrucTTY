#include "Renderer.hpp"
#include "Palette.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>

Renderer::Renderer(int width, int height, const std::string& mode, bool show_structure)
    : width_(width), height_(height), mode_(mode), show_structure_(show_structure)
{
    logical_pixels_.resize(width_ * 2 * height_ * 4);
}

void Renderer::set_depth_params(float focal_offset, float zoom_level,
                                float depth_min_z, float depth_max_z) {
    focal_offset_     = focal_offset;
    zoom_level_       = zoom_level;
    depth_base_min_z_ = depth_min_z;
    depth_base_max_z_ = depth_max_z;
}

void Renderer::render(const std::vector<RenderAtom>& atoms) {
    clear();
    std::vector<RenderPoint> final_points;
    final_points.reserve(50000);
    project_and_fill(atoms, final_points);
    zbuffer_resolve(final_points);
}

const std::vector<RenderPoint>& Renderer::get_pixels() const {
    return logical_pixels_;
}

int Renderer::get_logical_width()  const { return width_  * 2; }
int Renderer::get_logical_height() const { return height_ * 4; }

void Renderer::clear() {
    logical_pixels_.assign(width_ * 2 * height_ * 4, RenderPoint());
}

int Renderer::compute_depth_band(float z) const {
    float range = depth_base_max_z_ - depth_base_min_z_;
    if (range <= 0.0f) return 1;
    float t = (z - depth_base_min_z_) / range;
    if      (t < 0.25f) return 0;
    else if (t < 0.60f) return 1;
    else                return 2;
}

void Renderer::draw_line_impl(std::vector<RenderPoint>& out,
                              int x1, int x2, int y1, int y2,
                              float z1, float z2,
                              const std::string& chain_id, char structure,
                              int max_x, int max_y, int half) const {
    int dx = x2 - x1;
    int dy = y2 - y1;
    float dz = z2 - z1;

    int steps = std::max(std::abs(dx), std::abs(dy));
    if (steps == 0) steps = 1;

    float xInc = (float)dx / steps;
    float yInc = (float)dy / steps;
    float zInc = dz / steps;

    float fx = (float)x1;
    float fy = (float)y1;
    float fz = z1;

    for (int i = 0; i <= steps; ++i) {
        int ix = (int)fx;
        int iy = (int)fy;
        int band = compute_depth_band(fz);

        for (int oy = -half; oy <= half; oy++) {
            for (int ox = -half; ox <= half; ox++) {
                if (ox != 0 && oy != 0) continue;
                int nx = ix + ox, ny = iy + oy;
                if (nx >= 0 && nx < max_x && ny >= 0 && ny < max_y) {
                    RenderPoint rp{nx, ny, fz, ' ', 0, chain_id, structure};
                    rp.depth_band = band;
                    out.push_back(rp);
                }
            }
        }

        fx += xInc;
        fy += yInc;
        fz += zInc;
    }
}

void Renderer::project_and_fill(const std::vector<RenderAtom>& atoms,
                                std::vector<RenderPoint>& out) const {
    const float nearPlane = 0.05f;
    float fovRads = 1.0f / std::tan((FOV_ / zoom_level_) * 0.5f / 180.0f * PI_);

    const int logical_w = width_  * 2;
    const int logical_h = height_ * 4;

    std::vector<RenderPoint> chain_points;
    chain_points.reserve(8000);

    size_t n = atoms.size();
    size_t protein_start = 0;

    while (protein_start < n) {
        int cur_protein = atoms[protein_start].protein_index;
        size_t protein_end = protein_start;
        while (protein_end < n && atoms[protein_end].protein_index == cur_protein)
            ++protein_end;

        chain_points.clear();

        size_t chain_start = protein_start;
        while (chain_start < protein_end) {
            const std::string& cur_chain = atoms[chain_start].chain_id;
            size_t chain_end = chain_start;
            while (chain_end < protein_end && atoms[chain_end].chain_id == cur_chain)
                ++chain_end;

            int num_atoms = (int)(chain_end - chain_start);
            int prevScreenX = -1, prevScreenY = -1;
            float prevZ = -1.0f;

            for (int i = 0; i < num_atoms; ++i) {
                const RenderAtom& a = atoms[chain_start + i];
                float x = a.x;
                float y = a.y;
                float z = a.z + focal_offset_;

                if (z < nearPlane) {
                    prevScreenX = prevScreenY = -1;
                    prevZ = -1.0f;
                    continue;
                }

                char structure = a.structure;
                float projectedX = (x / z) * fovRads + a.pan_x;
                float projectedY = (y / z) * fovRads + a.pan_y;

                int screenX = (int)((projectedX + 1.0f) * 0.5f * logical_w);
                int screenY = (int)((1.0f - projectedY) * 0.5f * logical_h);

                if (prevScreenX != -1 && prevScreenY != -1) {
                    size_t before_draw = chain_points.size();
                    if (structure != 'H' && structure != 'S') {
                        const RenderAtom& P0 = atoms[chain_start + std::max(0, i - 2)];
                        const RenderAtom& P1 = atoms[chain_start + (i - 1)];
                        const RenderAtom& P2 = atoms[chain_start + i];
                        const RenderAtom& P3 = atoms[chain_start + std::min(num_atoms - 1, i + 1)];
                        int cr_prevX = prevScreenX, cr_prevY = prevScreenY;
                        float cr_prevZ = prevZ;
                        for (int cr = 1; cr <= 4; ++cr) {
                            float t = cr * 0.25f;
                            float t2 = t * t, t3 = t2 * t;
                            float cx = 0.5f * ((-P0.x + 3*P1.x - 3*P2.x + P3.x)*t3
                                             + ( 2*P0.x - 5*P1.x + 4*P2.x - P3.x)*t2
                                             + (-P0.x + P2.x)*t + 2*P1.x);
                            float cy = 0.5f * ((-P0.y + 3*P1.y - 3*P2.y + P3.y)*t3
                                             + ( 2*P0.y - 5*P1.y + 4*P2.y - P3.y)*t2
                                             + (-P0.y + P2.y)*t + 2*P1.y);
                            float cz = 0.5f * ((-P0.z + 3*P1.z - 3*P2.z + P3.z)*t3
                                             + ( 2*P0.z - 5*P1.z + 4*P2.z - P3.z)*t2
                                             + (-P0.z + P2.z)*t + 2*P1.z) + focal_offset_;
                            if (cz < nearPlane) break;
                            float cpX = (cx / cz) * fovRads + a.pan_x;
                            float cpY = (cy / cz) * fovRads + a.pan_y;
                            int cX = (int)((cpX + 1.0f) * 0.5f * logical_w);
                            int cY = (int)((1.0f - cpY) * 0.5f * logical_h);
                            draw_line_impl(chain_points, cr_prevX, cX, cr_prevY, cY,
                                          cr_prevZ, cz, cur_chain, structure,
                                          logical_w, logical_h, 1);
                            cr_prevX = cX; cr_prevY = cY; cr_prevZ = cz;
                        }
                    } else {
                        draw_line_impl(chain_points, prevScreenX, screenX, prevScreenY, screenY,
                                      prevZ, z, cur_chain, structure,
                                      logical_w, logical_h, 1);
                    }
                    for (size_t k = before_draw; k < chain_points.size(); ++k) {
                        chain_points[k].bfactor            = a.bfactor;
                        chain_points[k].is_interface       = a.is_interface;
                        chain_points[k].is_aligned         = a.is_aligned;
                        chain_points[k].conservation_score = a.conservation_score;
                        chain_points[k].residue_number     = a.residue_number;
                        std::strncpy(chain_points[k].residue_name, a.residue_name, 3);
                        chain_points[k].residue_name[3]    = '\0';
                    }
                }

                if (screenX >= 0 && screenX < logical_w && screenY >= 0 && screenY < logical_h) {
                    int band = compute_depth_band(z);
                    if (structure == 'x') {
                        RenderPoint rp{screenX, screenY, z, ' ', 0, cur_chain, structure};
                        rp.bfactor = a.bfactor; rp.is_interface = a.is_interface;
                        rp.is_aligned = a.is_aligned; rp.conservation_score = a.conservation_score;
                        rp.residue_number = a.residue_number;
                        std::strncpy(rp.residue_name, a.residue_name, 3);
                        rp.residue_name[3] = '\0';
                        rp.depth_band = band;
                        chain_points.push_back(rp);
                        rp.y = screenY - 1;
                        if (rp.y >= 0) chain_points.push_back(rp);
                        rp.y = screenY + 1;
                        if (rp.y < logical_h) chain_points.push_back(rp);
                    } else {
                        for (int oy = -1; oy <= 1; oy++)
                            for (int ox = -1; ox <= 1; ox++) {
                                if (ox != 0 && oy != 0) continue;
                                int nx = screenX + ox, ny = screenY + oy;
                                if (nx >= 0 && nx < logical_w && ny >= 0 && ny < logical_h) {
                                    RenderPoint rp{nx, ny, z, ' ', 0, cur_chain, structure};
                                    rp.bfactor = a.bfactor; rp.is_interface = a.is_interface;
                                    rp.is_aligned = a.is_aligned; rp.conservation_score = a.conservation_score;
                                    rp.residue_number = a.residue_number;
                                    std::strncpy(rp.residue_name, a.residue_name, 3);
                                    rp.residue_name[3] = '\0';
                                    rp.depth_band = band;
                                    chain_points.push_back(rp);
                                }
                            }
                    }
                }

                prevScreenX = screenX;
                prevScreenY = screenY;
                prevZ = z;
            }

            chain_start = chain_end;
        }

        assign_colors_impl(chain_points, cur_protein);
        out.insert(out.end(), chain_points.begin(), chain_points.end());

        protein_start = protein_end;
    }
}

void Renderer::assign_colors_impl(std::vector<RenderPoint>& points, int protein_idx) const {
    if (points.empty()) return;

    if (mode_ == "protein") {
        int idx = protein_idx % 9;
        for (auto& pt : points) {
            if      (pt.depth_band == 0) pt.color_id = idx + 120;
            else if (pt.depth_band == 1) pt.color_id = idx + 1;
            else                         pt.color_id = idx + 200;
        }
    } else if (mode_ == "chain") {
        std::string cur_chain = points[0].chainID;
        int color_idx = 0;
        for (auto& pt : points) {
            if (pt.chainID != cur_chain) { color_idx++; cur_chain = pt.chainID; }
            int ci = (protein_idx * 10 + color_idx) % 15;
            if      (pt.depth_band == 0) pt.color_id = 130 + ci;
            else if (pt.depth_band == 1) pt.color_id = 21  + ci;
            else                         pt.color_id = 145 + ci;
        }
    } else if (mode_ == "rainbow") {
        int num_points = (int)points.size();
        for (int i = 0; i < num_points; i++) {
            int color_idx = (i * 20) / std::max(1, num_points);
            if      (points[i].depth_band == 0) points[i].color_id = color_idx + 160;
            else if (points[i].depth_band == 1) points[i].color_id = color_idx + 51;
            else                                points[i].color_id = color_idx + 180;
        }
    } else if (mode_ == "plddt") {
        for (auto& pt : points) {
            float plddt = pt.bfactor;
            int base;
            if      (plddt >= 90) base = 0;
            else if (plddt >= 70) base = 1;
            else if (plddt >= 50) base = 2;
            else                  base = 3;
            if      (pt.depth_band == 0) pt.color_id = 209 + base;
            else if (pt.depth_band == 1) pt.color_id = 71  + base;
            else                         pt.color_id = 213 + base;
        }
    } else if (mode_ == "interface") {
        for (auto& pt : points) {
            if      (pt.depth_band == 0) pt.color_id = pt.is_interface ? 237 : 238;
            else if (pt.depth_band == 1) pt.color_id = pt.is_interface ? 43  : 44;
            else                         pt.color_id = pt.is_interface ? 239 : 240;
        }
    } else if (mode_ == "aligned") {
        int bright_id = (protein_idx % 9) + 101;
        int near_id   = (protein_idx % 9) + 241;
        for (auto& pt : points) {
            if (pt.is_aligned) {
                if      (pt.depth_band == 0) pt.color_id = near_id;
                else if (pt.depth_band == 1) pt.color_id = bright_id;
                else                         pt.color_id = bright_id;
            } else {
                pt.color_id = (pt.depth_band == 2) ? 250 : 110;
            }
        }
    } else if (mode_ == "conservation") {
        for (auto& pt : points) {
            float score = pt.conservation_score;
            if (score < 0) {
                pt.color_id = 11;
            } else {
                int idx = std::max(0, std::min(9, (int)(score * 9.0f)));
                if      (pt.depth_band == 0) pt.color_id = 217 + idx;
                else if (pt.depth_band == 1) pt.color_id = 75  + idx;
                else                         pt.color_id = 227 + idx;
            }
        }
    }

    if (show_structure_ && mode_ == "protein") {
        int dim_id = (protein_idx % 9) + 11;
        for (auto& pt : points) {
            if      (pt.structure == 'H') pt.color_id = 41;
            else if (pt.structure == 'S') pt.color_id = 42;
            else                          pt.color_id = dim_id;
        }
    }
}

void Renderer::zbuffer_resolve(const std::vector<RenderPoint>& points) {
    const int logical_w = width_  * 2;
    const int logical_h = height_ * 4;
    for (const auto& pt : points) {
        if (pt.x < 0 || pt.x >= logical_w || pt.y < 0 || pt.y >= logical_h) continue;
        int idx = pt.y * logical_w + pt.x;
        if (pt.depth < logical_pixels_[idx].depth) {
            logical_pixels_[idx] = pt;
        }
    }
}
