#include "StructureMaker.hpp"

const float PI = 3.14159265359f;

StructureMaker::StructureMaker() {}
StructureMaker::~StructureMaker() {}

std::vector<std::vector<Atom>> StructureMaker::extract_helix_segments(const Atom* atoms, int num_atoms) {
    std::vector<std::vector<Atom>> helices;
    std::vector<Atom> current;

    for (int i = 0; i < num_atoms; ++i) {
        if (atoms[i].structure == 'H') {
            current.push_back(atoms[i]);
        } else {
            if (current.size() >= 4) {
                helices.push_back(std::move(current));
            }
            current.clear();
        }
    }

    if (current.size() >= 4) {
        helices.push_back(std::move(current));
    }

    return helices;
}


void StructureMaker::compute_helix_axis(const std::vector<Atom>& helix, float (&center)[3], float (&axis)[3]) {
    size_t N = helix.size();
    if (N == 0) return;

    float sum[3] = {0.0f, 0.0f, 0.0f};

    // calculate center coord
    for (const Atom& atom : helix) {
        sum[0] += atom.x;
        sum[1] += atom.y;
        sum[2] += atom.z;
    }

    center[0] = sum[0] / N;
    center[1] = sum[1] / N;
    center[2] = sum[2] / N;

    // calculate covariance matrix
    float cov[3][3] = {0};

    for (const Atom& atom : helix) {
        float dx = atom.x - center[0];
        float dy = atom.y - center[1];
        float dz = atom.z - center[2];

        cov[0][0] += dx * dx;
        cov[0][1] += dx * dy;
        cov[0][2] += dx * dz;

        cov[1][0] += dy * dx;
        cov[1][1] += dy * dy;
        cov[1][2] += dy * dz;

        cov[2][0] += dz * dx;
        cov[2][1] += dz * dy;
        cov[2][2] += dz * dz;
    }

    // main axis change by power iteration
    float vec[3] = {1.0f, 0.0f, 0.0f};  // init vector
    for (int iter = 0; iter < 10; ++iter) {
        float x = cov[0][0]*vec[0] + cov[0][1]*vec[1] + cov[0][2]*vec[2];
        float y = cov[1][0]*vec[0] + cov[1][1]*vec[1] + cov[1][2]*vec[2];
        float z = cov[2][0]*vec[0] + cov[2][1]*vec[1] + cov[2][2]*vec[2];

        float norm = std::sqrt(x*x + y*y + z*z);
        if (norm == 0.0f) break;

        vec[0] = x / norm;
        vec[1] = y / norm;
        vec[2] = z / norm;
    }

    axis[0] = vec[0];
    axis[1] = vec[1];
    axis[2] = vec[2];
}



void StructureMaker::calculate_ss_points(std::map<std::string, std::vector<Atom>>& init_atoms,
                                         std::map<std::string, std::vector<Atom>>& ss_atoms) {
    ss_atoms.clear();

    for (auto& [chainID, atoms] : init_atoms) {
        std::vector<Atom>& output = ss_atoms[chainID];
        size_t i = 0;
        while (i < atoms.size()) {
            char s = atoms[i].structure;

            if (s == 'H') {
                // Find the full helix segment
                size_t start = i;
                while (i < atoms.size() && atoms[i].structure == 'H') ++i;
                size_t end = i;

                if (end - start >= 4) {
                    auto segment = std::vector<Atom>(atoms.begin() + start, atoms.begin() + end);

                    float center[3], axis[3];
                    compute_helix_axis(segment, center, axis);

                    float dx = segment.back().x - segment.front().x;
                    float dy = segment.back().y - segment.front().y;
                    float dz = segment.back().z - segment.front().z;
                    float length = std::sqrt(dx * dx + dy * dy + dz * dz);

                    // Half-residue axial sampling: smooth with 16 stripes,
                    // avoids O(residues * circle_steps) geometry blowup on long helices.
                    const int steps = std::max(8, (int)(end - start) / 2);

                    float up[3] = {0.0f, 0.0f, 1.0f};
                    if (std::abs(axis[2]) > 0.99f) { up[0] = 1.0f; up[2] = 0.0f; }

                    float n1[3] = {
                        axis[1]*up[2] - axis[2]*up[1],
                        axis[2]*up[0] - axis[0]*up[2],
                        axis[0]*up[1] - axis[1]*up[0]
                    };
                    float n1_norm = std::sqrt(n1[0]*n1[0] + n1[1]*n1[1] + n1[2]*n1[2]);
                    for (int k = 0; k < 3; ++k) n1[k] /= n1_norm;

                    float n2[3] = {
                        axis[1]*n1[2] - axis[2]*n1[1],
                        axis[2]*n1[0] - axis[0]*n1[2],
                        axis[0]*n1[1] - axis[1]*n1[0]
                    };

                    // Render as a spiral ribbon that winds around the helix axis.
                    // The radial direction rotates by total_turns full turns as t goes 0→1,
                    // producing a coiled-ribbon appearance.  A few parallel lines give the
                    // ribbon its visible width.
                    float num_residues = (float)(end - start);
                    float total_turns  = num_residues / 3.6f;  // alpha-helix: ~3.6 res/turn
                    const float ribbon_half_width = 1.5f;      // Angstroms
                    const int   ribbon_lines      = 3;         // half-count; 7 lines total
                    const int   spiral_steps      = std::max(100, (int)(end - start) * 8);

                    for (int rs = -ribbon_lines; rs <= ribbon_lines; ++rs) {
                        float rw = ((float)rs / ribbon_lines) * ribbon_half_width;
                        bool forward = ((rs + ribbon_lines) % 2 == 0);

                        for (int si = 0; si <= spiral_steps; ++si) {
                            int   s_idx = forward ? si : (spiral_steps - si);
                            float t     = static_cast<float>(s_idx) / spiral_steps;
                            float angle = total_turns * 2.0f * PI * t;
                            float cos_a = std::cos(angle);
                            float sin_a = std::sin(angle);

                            // Outward radial direction at this angle
                            float rad[3] = {
                                cos_a * n1[0] + sin_a * n2[0],
                                cos_a * n1[1] + sin_a * n2[1],
                                cos_a * n1[2] + sin_a * n2[2]
                            };

                            // Tangential direction on cylinder surface (axis × rad)
                            // – ribbon width lies along this direction
                            float tan_dir[3] = {
                                axis[1]*rad[2] - axis[2]*rad[1],
                                axis[2]*rad[0] - axis[0]*rad[2],
                                axis[0]*rad[1] - axis[1]*rad[0]
                            };

                            float base[3] = {
                                center[0] + axis[0] * (t - 0.5f) * length,
                                center[1] + axis[1] * (t - 0.5f) * length,
                                center[2] + axis[2] * (t - 0.5f) * length,
                            };

                            float px = base[0] + radius * rad[0] + rw * tan_dir[0];
                            float py = base[1] + radius * rad[1] + rw * tan_dir[1];
                            float pz = base[2] + radius * rad[2] + rw * tan_dir[2];
                            output.emplace_back(px, py, pz, 'H');
                        }
                    }
                }
                // i already advanced to end of helix segment by the inner while loop
            }

            else if (s == 'S') {
                // Process the full consecutive sheet segment at once for a uniform ribbon
                size_t seg_start = i;
                while (i < atoms.size() && atoms[i].structure == 'S') ++i;
                size_t seg_end = i;
                int seg_len = (int)(seg_end - seg_start);

                if (seg_len < 2) {
                    output.push_back(atoms[seg_start]);
                    continue;
                }

                // Compute a consistent perpendicular from the overall segment direction
                float dx = atoms[seg_end-1].x - atoms[seg_start].x;
                float dy = atoms[seg_end-1].y - atoms[seg_start].y;
                float dz = atoms[seg_end-1].z - atoms[seg_start].z;
                float len = std::sqrt(dx*dx + dy*dy + dz*dz);
                if (len < 1e-6f) {
                    for (size_t k = seg_start; k < seg_end; ++k) output.push_back(atoms[k]);
                    continue;
                }

                float fwd[3] = {dx/len, dy/len, dz/len};
                float up[3] = {0.0f, 0.0f, 1.0f};
                if (std::abs(fwd[2]) > 0.99f) { up[0] = 1.0f; up[2] = 0.0f; }

                float n1[3] = {
                    fwd[1]*up[2] - fwd[2]*up[1],
                    fwd[2]*up[0] - fwd[0]*up[2],
                    fwd[0]*up[1] - fwd[1]*up[0]
                };
                float n1_norm = std::sqrt(n1[0]*n1[0] + n1[1]*n1[1] + n1[2]*n1[2]);
                if (n1_norm < 1e-6f) {
                    for (size_t k = seg_start; k < seg_end; ++k) output.push_back(atoms[k]);
                    continue;
                }
                for (int k = 0; k < 3; ++k) n1[k] /= n1_norm;

                // Draw ribbon as parallel zigzag stripes along the full segment.
                // Even steps traverse forward; odd steps traverse backward → short
                // cross-connections at the ribbon ends, long clean lines along the strand.
                for (int step = -width; step <= width; ++step) {
                    float offset = step * sheet_step;
                    float ox = n1[0]*offset, oy = n1[1]*offset, oz = n1[2]*offset;
                    bool forward = ((step + width) % 2 == 0);

                    for (int j = 0; j < seg_len - 1; ++j) {
                        int actual = forward ? j : (seg_len - 2 - j);
                        const Atom& pa = atoms[seg_start + actual];
                        const Atom& pb = atoms[seg_start + actual + 1];

                        float pair_dx = pb.x - pa.x;
                        float pair_dy = pb.y - pa.y;
                        float pair_dz = pb.z - pa.z;

                        // Two endpoints per Ca-Ca pair are enough; the renderer connects them
                        // with draw_line, so sub-sampling here just wastes atoms.
                        for (int t = 0; t <= 1; ++t) {
                            float f = static_cast<float>(t);
                            output.emplace_back(
                                pa.x + ox + f * pair_dx,
                                pa.y + oy + f * pair_dy,
                                pa.z + oz + f * pair_dz,
                                'S'
                            );
                        }
                    }
                }
                // i already advanced to seg_end by the inner while loop
            }

            else {
                output.push_back(atoms[i]);
                ++i;
            }
        }
    }
}
