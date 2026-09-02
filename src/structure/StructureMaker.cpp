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

    for (const Atom& atom : helix) {
        sum[0] += atom.x;
        sum[1] += atom.y;
        sum[2] += atom.z;
    }

    center[0] = sum[0] / N;
    center[1] = sum[1] / N;
    center[2] = sum[2] / N;

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

    float vec[3] = {1.0f, 0.0f, 0.0f};
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
                size_t start = i;
                while (i < atoms.size() && atoms[i].structure == 'H') ++i;
                size_t end = i;

                if (end - start < 4) {
                    for (size_t k = start; k < end; ++k) {
                        Atom coil_atom = atoms[k];
                        coil_atom.structure = 'x';
                        output.push_back(coil_atom);
                    }
                } else {
                    Atom junction_start = atoms[start];
                    junction_start.structure = 'x';
                    output.push_back(junction_start);

                    auto segment = std::vector<Atom>(atoms.begin() + start, atoms.begin() + end);

                    float center[3], axis[3];
                    compute_helix_axis(segment, center, axis);

                    float dx = segment.back().x - segment.front().x;
                    float dy = segment.back().y - segment.front().y;
                    float dz = segment.back().z - segment.front().z;
                    float length = std::sqrt(dx * dx + dy * dy + dz * dz);

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

                    float num_residues = (float)(end - start);
                    float total_turns  = num_residues / 3.6f;
                    const float ribbon_half_width = 0.5f;
                    const int   ribbon_lines      = 1;
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

                            float rad[3] = {
                                cos_a * n1[0] + sin_a * n2[0],
                                cos_a * n1[1] + sin_a * n2[1],
                                cos_a * n1[2] + sin_a * n2[2]
                            };

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
                            Atom geom(px, py, pz, 'H');
                            geom.bfactor            = segment.front().bfactor + t * (segment.back().bfactor - segment.front().bfactor);
                            geom.is_interface       = segment.front().is_interface;
                            geom.is_aligned         = segment.front().is_aligned;
                            geom.conservation_score = segment.front().conservation_score + t * (segment.back().conservation_score - segment.front().conservation_score);
                            {
                                int ca_idx = std::min(
                                    (int)std::round(t * (float)(segment.size() - 1)),
                                    (int)segment.size() - 1);
                                geom.residue_number = segment[ca_idx].residue_number;
                                geom.residue_name   = segment[ca_idx].residue_name;
                            }
                            output.push_back(geom);
                        }
                    }

                    Atom junction_end = atoms[end - 1];
                    junction_end.structure = 'x';
                    output.push_back(junction_end);
                }
            }

            else if (s == 'S') {
                size_t seg_start = i;
                while (i < atoms.size() && atoms[i].structure == 'S') ++i;
                size_t seg_end = i;
                int seg_len = (int)(seg_end - seg_start);

                if (seg_len < 2) {
                    Atom coil_atom = atoms[seg_start];
                    coil_atom.structure = 'x';
                    output.push_back(coil_atom);
                    continue;
                }

                {
                    Atom junc = atoms[seg_start];
                    junc.structure = 'x';
                    output.push_back(junc);
                }

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

                        for (int t = 0; t <= 1; ++t) {
                            float f = static_cast<float>(t);
                            Atom geom(pa.x + ox + f * pair_dx,
                                      pa.y + oy + f * pair_dy,
                                      pa.z + oz + f * pair_dz,
                                      'S');
                            geom.bfactor            = pa.bfactor + f * (pb.bfactor - pa.bfactor);
                            geom.is_interface       = pa.is_interface;
                            geom.is_aligned         = pa.is_aligned;
                            geom.conservation_score = pa.conservation_score + f * (pb.conservation_score - pa.conservation_score);
                            geom.residue_number = pa.residue_number;
                            geom.residue_name   = pa.residue_name;
                            output.push_back(geom);
                        }
                    }
                }

                {
                    Atom junc = atoms[seg_end - 1];
                    junc.structure = 'x';
                    output.push_back(junc);
                }
            }

            else {
                output.push_back(atoms[i]);
                ++i;
            }
        }
    }
}
