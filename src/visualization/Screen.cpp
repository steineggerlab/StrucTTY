#include "Screen.hpp"
#include <cstring>  // strncpy, memset

const float FOV = 90.0f;
const float PI  = 3.14159265359f;
const int MAX_STRUCT_NUM = 9;

Screen::Screen(const int& width, const int& height, const bool& show_structure,
               const std::string& mode, const std::string& depthcharacter) {
    screen_width = width;
    screen_height = height;
    screen_show_structure = show_structure;
    screen_mode = mode;
    screen_depthcharacter = depthcharacter;
    aspect_ratio = (float)screen_width / screen_height;
    zoom_level = 2.0f;

    start_color();
    use_default_colors();
    init_color_pairs();

    camera = new Camera(width, height, mode);
    panel  = new Panel(width, mode, show_structure);

    // 기능 6: 마우스 hover 초기화
    // mousemask()는 반환값과 무관하게 항상 escape sequence를 전송한다.
    // VSCode/Windows Terminal 등 terminfo 불일치 터미널 대응.
    keypad(stdscr, TRUE);
    mousemask(ALL_MOUSE_EVENTS | REPORT_MOUSE_POSITION, nullptr);
    mouseinterval(0);
    printf("\033[?1003h");  // xterm-1003: 모든 마우스 이동 추적 활성화
    fflush(stdout);
}

Screen::~Screen() {
    for (Protein* p : data) delete p;
    data.clear();
    delete camera;
    delete panel;

    // 기능 6: 마우스 이동 추적 비활성화
    printf("\033[?1003l");
    fflush(stdout);
}

static float compute_scene_radius_from_render_positions(const std::vector<Protein*>& data) {
    float max_r2 = 0.0f;
    bool any = false;

    for (auto* p : data) {
        for (const auto& [chainID, chain_atoms] : p->get_atoms()) {
            for (const auto& atom : chain_atoms) {
                float* pos = atom.get_position();
                float x = pos[0], y = pos[1], z = pos[2];
                float r2 = x*x + y*y + z*z;
                if (r2 > max_r2) max_r2 = r2;
                any = true;
            }
        }
    }
    if (!any) return 1.0f;
    return std::sqrt(max_r2);
}

void Screen::init_color_pairs() {
    // Protein vivid: pairs 1-9
    for (int i = 0; i < 9; ++i)
        init_pair(i + 1,  Palettes::PROTEIN_COLORS[i],     -1);
    // Protein dim (coil in protein+-s): pairs 11-19
    for (int i = 0; i < 9; ++i)
        init_pair(i + 11, Palettes::PROTEIN_DIM_COLORS[i], -1);
    // Chain colors: pairs 21-35
    for (int i = 0; i < 15; ++i)
        init_pair(i + 21, Palettes::CHAIN_COLORS[i],       -1);
    // Rainbow: pairs 51-70
    for (int i = 0; i < 20; ++i)
        init_pair(i + 51, Palettes::RAINBOW[i],            -1);
    // Secondary structure
    init_pair(41, 226, -1);  // yellow helix
    init_pair(42,  51, -1);  // cyan sheet
    // Interface region: pairs 43-44
    init_pair(43, Palettes::INTERFACE_COLOR,     -1);  // interface residue (bright magenta)
    init_pair(44, Palettes::INTERFACE_DIM_COLOR, -1);  // non-interface dim
    // Aligned region: pairs 45-46
    init_pair(45, Palettes::ALIGNED_COLOR,     -1);  // aligned residue (bright green)
    init_pair(46, Palettes::ALIGNED_DIM_COLOR, -1);  // non-aligned dim
    // pLDDT: pairs 71-74
    for (int i = 0; i < 4; ++i)
        init_pair(i + 71, Palettes::PLDDT_COLORS[i], -1);
    // Conservation gradient: pairs 75-84
    for (int i = 0; i < 10; ++i)
        init_pair(i + 75, Palettes::CONSERVATION_COLORS[i], -1);
}

void Screen::set_protein(const std::string& in_file, int ii, const bool& show_structure) {
    Protein* protein = new Protein(in_file, chainVec.at(ii), show_structure);
    data.push_back(protein);
    pan_x.push_back(0.0f);
    pan_y.push_back(0.0f);

    init_color_pairs();
}

void Screen::set_tmatrix() {
    size_t filenum = data.size();
    vectorpointer = new float*[filenum];
    for (size_t i = 0; i < filenum; i++) {
        vectorpointer[i] = new float[3];
        vectorpointer[i][0] = 0;
        vectorpointer[i][1] = 0;
        vectorpointer[i][2] = 0;
    }
}

void Screen::set_chainfile(const std::string& chainfile, int filesize) {
    for (int i = 0; i < filesize; i++) {
        chainVec.push_back("-");
    }
    if (chainfile.empty()) return;

    std::ifstream file(chainfile);
    if (!file.is_open()) {
        std::cerr << "Failed to open chainfile\n";
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        std::istringstream iss(line);
        int index;
        std::string chainlist;
        iss >> index >> chainlist;
        if (index >= filesize) {
            std::cout << "Index in chain file should be 0..filenum-1\n";
            continue;
        }
        chainVec[index] = chainlist;
    }
    file.close();
}

void Screen::set_utmatrix(const std::string& utmatrix, bool applyUT) {
    yesUT = !utmatrix.empty();

    const size_t filenum = data.size();
    float** matrixpointer = new float*[filenum];
    for (size_t i = 0; i < filenum; i++) {
        matrixpointer[i] = new float[9];
        for (int j = 0; j < 9; j++) {
            matrixpointer[i][j] = (j % 4 == 0) ? 1.f : 0.f; // identity
        }
    }

    if (utmatrix.empty()) {
        for (size_t i = 0; i < filenum; i++) delete[] matrixpointer[i];
        delete[] matrixpointer;
        return;
    }

    std::ifstream file(utmatrix);
    if (!file.is_open()) {
        std::cerr << "Failed to open utmatrix file\n";
        for (size_t i = 0; i < filenum; i++) delete[] matrixpointer[i];
        delete[] matrixpointer;
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        std::istringstream iss(line);
        int index;
        std::string mat9Str;
        std::string mat3Str;

        iss >> index >> mat9Str >> mat3Str;
        if (index < 0 || index >= (int)filenum) continue;

        {
            std::istringstream mss(mat9Str);
            std::string val;
            int count = 0;
            while (std::getline(mss, val, ',') && count < 9) {
                matrixpointer[index][count++] = std::stof(val);
            }
        }

        {
            std::istringstream mss(mat3Str);
            std::string val;
            int count = 0;
            while (std::getline(mss, val, ',') && count < 3) {
                vectorpointer[index][count++] = std::stof(val);
            }
        }
    }

    if (applyUT) {
        for (size_t i = 0; i < filenum; i++) {
            data[i]->do_naive_rotation(matrixpointer[i]); // U → screen_atoms
            data[i]->do_shift(vectorpointer[i]);          // T → screen_atoms
            data[i]->apply_ut_to_init_atoms(matrixpointer[i], vectorpointer[i]); // U,T → init_atoms
        }
    }

    for (size_t i = 0; i < filenum; i++) delete[] matrixpointer[i];
    delete[] matrixpointer;
}

void Screen::normalize_proteins(const std::string& utmatrix) {
    const bool hasUT = !utmatrix.empty();
    for (size_t i = 0; i < data.size(); i++) {
        auto* p = data[i];
        p->load_data(vectorpointer[i], yesUT);
        panel->add_panel_info(p->get_file_name(),
                              p->get_chain_length(),
                              p->get_residue_count());
    }

    if (hasUT) {
        set_utmatrix(utmatrix, true);
    }

    global_bb = BoundingBox();
    for (auto* p : data) {
        p->set_bounding_box();
        global_bb = global_bb + p->get_bounding_box();
    }

    float max_ext = std::max(global_bb.max_x - global_bb.min_x,
                             global_bb.max_y - global_bb.min_y);
    max_ext = std::max(max_ext, global_bb.max_z - global_bb.min_z);
    float scale = (max_ext > 0.f) ? (2.0f / max_ext) : 1.0f;

    if (hasUT) {
        float gx = 0.5f * (global_bb.min_x + global_bb.max_x);
        float gy = 0.5f * (global_bb.min_y + global_bb.max_y);
        float gz = 0.5f * (global_bb.min_z + global_bb.max_z);
        float global_shift[3] = { -gx, -gy, -gz };

        for (auto* p : data) {
            p->set_scale(scale);
            p->do_shift(global_shift);
            p->do_scale(scale);
        }
    } else {
        for (auto* p : data) {
            float center_shift[3] = { -p->cx, -p->cy, -p->cz };
            p->set_scale(scale);
            p->do_shift(center_shift);
            p->do_scale(scale);
        }
    }

    depth_calibrated = false;

    // 기능 3: 정규화 파라미터 저장 (hit 탐색 시 동일 스케일 적용)
    norm_scale = scale;
    if (!data.empty() && data[0]) {
        // cx/cy/cz는 set_bounding_box() 이후, do_shift() 이전 값
        // normalize_proteins() 의 non-UT 경로에서 -p->cx, -p->cy, -p->cz 로 shift
        // 여기서는 shift 이전의 centroid를 저장해야 하므로 do_shift 전에 저장
        // 하지만 이미 do_shift가 호출됐으므로 cx,cy,cz는 변경 전 값 그대로 있음
        norm_cx = data[0]->cx;
        norm_cy = data[0]->cy;
        norm_cz = data[0]->cz;
    }

    float radius = compute_scene_radius_from_render_positions(data);

    focal_offset = std::clamp(2.5f * radius + 1.0f, 2.0f, 8.0f);

    int n = (int)pan_x.size();
    int cols, rows;
    switch (n) {
        case 1:  cols=1; rows=1; break;
        case 2:  cols=2; rows=1; break;
        case 3:  cols=3; rows=1; break;
        case 4:  cols=2; rows=2; break;
        case 5:  cols=3; rows=2; break;
        case 6:  cols=3; rows=2; break;
        case 7:  cols=3; rows=3; break;
        case 8:  cols=3; rows=3; break;
        default: cols=3; rows=(n+2)/3; break;
    }

    if (n > 1 && !hasUT) {
        int max_dim = std::max(cols, rows);
        float step      = (max_dim == 2) ? 0.75f : 0.5f;
        float foc_scale = (max_dim == 2) ? 0.8f  : 0.6f;
        focal_offset *= max_dim * foc_scale;
        for (int i = 0; i < n; i++) {
            int col = i % cols;
            int row = i / cols;
            pan_x[i] =  (col - (cols - 1) / 2.0f) * step;
            pan_y[i] = -((row - (rows - 1) / 2.0f) * step);
        }
    } else {
        for (int i = 0; i < n; i++) {
            pan_x[i] = 0.0f;
            pan_y[i] = 0.0f;
        }
    }
}

char Screen::get_pixel_char_from_depth(float z, float min_z, float max_z) {
    if (screen_depthcharacter.empty()) return '#';
    if (!(max_z > min_z)) return screen_depthcharacter.back();

    float t = (z - min_z) / (max_z - min_z); // 0 near, 1 far
    t = std::clamp(t, 0.0f, 1.0f);

    int L = (int)screen_depthcharacter.size();
    int idx = (int)std::floor(t * (L - 1));
    idx = std::clamp(idx, 0, L - 1);
    return screen_depthcharacter[idx];
}

void Screen::calibrate_depth_baseline_first_view() {
    const float nearPlane = 0.05f;
    float fovRads = 1.0f / std::tan((FOV / zoom_level) * 0.5f / 180.0f * PI);

    float minz = std::numeric_limits<float>::infinity();
    float maxz = -std::numeric_limits<float>::infinity();
    bool any = false;

    for (size_t ii = 0; ii < data.size(); ii++) {
        Protein* target = data[ii];

        for (const auto& [chainID, chain_atoms] : target->get_atoms()) {
            if (chain_atoms.empty()) continue;
            int num_atoms = target->get_chain_length(chainID);

            for (int i = 0; i < num_atoms; ++i) {
                float* position = chain_atoms[i].get_position();
                float x = position[0];
                float y = position[1];
                float z = position[2] + focal_offset;

                if (z < nearPlane) continue;

                float projectedX = (x / z) * fovRads + pan_x[ii];
                float projectedY = (y / z) * fovRads + pan_y[ii];

                int screenX = (int)((projectedX + 1.0f) * 0.5f * screen_width);
                int screenY = (int)((1.0f - projectedY) * 0.5f * screen_height);

                if (screenX < 0 || screenX >= screen_width || screenY < 0 || screenY >= screen_height) continue;

                minz = std::min(minz, z);
                maxz = std::max(maxz, z);
                any = true;
            }
        }
    }

    if (!any) {
        depth_base_min_z = focal_offset;
        depth_base_max_z = focal_offset + 1.0f;
    } else {
        const float minSpan = 0.3f;
        if ((maxz - minz) < minSpan) maxz = minz + minSpan;
        depth_base_min_z = minz;
        depth_base_max_z = maxz;
    }

    depth_calibrated = true;
}

void Screen::draw_line(std::vector<RenderPoint>& points,
                       int x1, int x2,
                       int y1, int y2,
                       float z1, float z2,
                       const std::string& chainID, char structure,
                       float min_z, float max_z,
                       int max_x, int max_y, int half) {
    if (max_x < 0) max_x = screen_width;
    if (max_y < 0) max_y = screen_height;

    int dx = x2 - x1;
    int dy = y2 - y1;
    float dz = z2 - z1;

    int steps = std::max(std::abs(dx), std::abs(dy));
    if (steps == 0) steps = 1;

    float xIncrement = (float)dx / steps;
    float yIncrement = (float)dy / steps;
    float zIncrement = (float)dz / steps;

    float x = (float)x1;
    float y = (float)y1;
    float z = z1;

    for (int i = 0; i <= steps; ++i) {
        int ix = (int)x;
        int iy = (int)y;

        char pix = get_pixel_char_from_depth(z, min_z, max_z);
        for (int oy = -half; oy <= half; oy++) {
            for (int ox = -half; ox <= half; ox++) {
                if (ox != 0 && oy != 0) continue;  // cross (+) pattern, not square
                int nx = ix + ox, ny = iy + oy;
                if (nx >= 0 && nx < max_x && ny >= 0 && ny < max_y)
                    points.push_back({nx, ny, z, pix, 0, chainID, structure});
            }
        }

        x += xIncrement;
        y += yIncrement;
        z += zIncrement;
    }
}

void Screen::assign_colors_to_points(std::vector<RenderPoint>& points, int protein_idx) {
    if (points.empty()) return;

    if (screen_mode == "protein") {
        int idx = protein_idx % 9;
        for (auto& pt : points) pt.color_id = idx + 1;  // pairs 1-9
    } else if (screen_mode == "chain") {
        std::string cur_chain = points[0].chainID;
        int color_idx  = 0;
        for (auto& pt : points) {
            if (pt.chainID != cur_chain) { color_idx++; cur_chain = pt.chainID; }
            pt.color_id = 21 + ((protein_idx * 10 + color_idx) % 15);  // pairs 21-35
        }
    } else if (screen_mode == "rainbow") {
        int num_points = (int)points.size();
        for (int i = 0; i < num_points; i++) {
            int color_idx = (i * 20) / std::max(1, num_points);
            points[i].color_id = color_idx + 51;  // pairs 51-70
        }
    } else if (screen_mode == "plddt") {
        for (auto& pt : points) {
            float plddt = pt.bfactor;
            if      (plddt >= 90) pt.color_id = 71;  // Very high: 파란색
            else if (plddt >= 70) pt.color_id = 72;  // Confident: 청록색
            else if (plddt >= 50) pt.color_id = 73;  // Low: 노란색
            else                  pt.color_id = 74;  // Very low: 주황색
        }
    } else if (screen_mode == "interface") {
        for (auto& pt : points) {
            pt.color_id = pt.is_interface ? 43 : 44;  // 강조(마젠타) or dim
        }
    } else if (screen_mode == "aligned") {
        for (auto& pt : points) {
            pt.color_id = pt.is_aligned ? 45 : 46;  // 정렬됨(초록) or dim
        }
    } else if (screen_mode == "conservation") {
        for (auto& pt : points) {
            float score = pt.conservation_score;
            if (score < 0) {
                pt.color_id = 11;  // 미설정: olive dim
            } else {
                int idx = std::max(0, std::min(9, (int)(score * 9.0f)));
                pt.color_id = 75 + idx;  // conservation gradient pairs 75-84
            }
        }
    } else {
        std::cerr << "Unknown mode: " << screen_mode << std::endl;
    }

    // In protein+-s mode: H=yellow(41), S=cyan(42), coil=dimmed protein color(11-19).
    // Chain and rainbow are never overridden.
    if (screen_show_structure && screen_mode == "protein") {
        int dim_id = (protein_idx % 9) + 11;
        for (auto& pt : points) {
            if      (pt.structure == 'H') pt.color_id = 41;
            else if (pt.structure == 'S') pt.color_id = 42;
            else                          pt.color_id = dim_id;
        }
    }
}

void Screen::project() {
    const float nearPlane = 0.05f;

    float fovRads = 1.0f / std::tan((FOV / zoom_level) * 0.5f / 180.0f * PI);

    if (!depth_calibrated) {
        calibrate_depth_baseline_first_view();
    }

    if (use_braille) {
        // Braille sub-pixel path: render to logicalPixels at 2x width, 4x height.
        // Each terminal cell maps to a 2x4 Braille dot grid giving 8x more detail.
        const int logical_w = screen_width * 2;
        const int logical_h = screen_height * 4;

        std::vector<RenderPoint> finalPoints;
        std::vector<RenderPoint> chainPoints;
        finalPoints.reserve(50000);

        int protein_idx = 0;
        for (size_t ii = 0; ii < data.size(); ii++) {
            Protein* target = data[ii];
            chainPoints.clear();

            for (const auto& [chainID, chain_atoms] : target->get_atoms()) {
                if (chain_atoms.empty()) continue;

                int num_atoms = target->get_chain_length(chainID);
                int prevScreenX = -1, prevScreenY = -1;
                float prevZ = -1.0f;

                for (int i = 0; i < num_atoms; ++i) {
                    float* position = chain_atoms[i].get_position();
                    float x = position[0];
                    float y = position[1];
                    float z = position[2] + focal_offset;

                    if (z < nearPlane) {
                        prevScreenX = prevScreenY = -1;
                        prevZ = -1.0f;
                        continue;
                    }

                    char structure = chain_atoms[i].get_structure();
                    float projectedX = (x / z) * fovRads + pan_x[ii];
                    float projectedY = (y / z) * fovRads + pan_y[ii];

                    int screenX = (int)((projectedX + 1.0f) * 0.5f * logical_w);
                    int screenY = (int)((1.0f - projectedY) * 0.5f * logical_h);

                    const Atom& cur_atom = chain_atoms[i];
                    if (prevScreenX != -1 && prevScreenY != -1) {
                        size_t before_draw = chainPoints.size();
                        if (structure != 'H' && structure != 'S') {
                            // Catmull-Rom for coil/backbone only — SS geometry is already
                            // dense from StructureMaker so applying it there is wasteful.
                            const Atom& P0 = chain_atoms[std::max(0, i-2)];
                            const Atom& P1 = chain_atoms[i-1];
                            const Atom& P2 = chain_atoms[i];
                            const Atom& P3 = chain_atoms[std::min(num_atoms-1, i+1)];
                            int cr_prevX = prevScreenX, cr_prevY = prevScreenY;
                            float cr_prevZ = prevZ;
                            for (int cr = 1; cr <= 4; ++cr) {
                                float t = cr * 0.25f;
                                float t2 = t*t, t3 = t2*t;
                                float cx = 0.5f*((-P0.x+3*P1.x-3*P2.x+P3.x)*t3 + (2*P0.x-5*P1.x+4*P2.x-P3.x)*t2 + (-P0.x+P2.x)*t + 2*P1.x);
                                float cy = 0.5f*((-P0.y+3*P1.y-3*P2.y+P3.y)*t3 + (2*P0.y-5*P1.y+4*P2.y-P3.y)*t2 + (-P0.y+P2.y)*t + 2*P1.y);
                                float cz = 0.5f*((-P0.z+3*P1.z-3*P2.z+P3.z)*t3 + (2*P0.z-5*P1.z+4*P2.z-P3.z)*t2 + (-P0.z+P2.z)*t + 2*P1.z) + focal_offset;
                                if (cz < nearPlane) break;
                                float cpX = (cx/cz)*fovRads + pan_x[ii];
                                float cpY = (cy/cz)*fovRads + pan_y[ii];
                                int cX = (int)((cpX+1.0f)*0.5f*logical_w);
                                int cY = (int)((1.0f-cpY)*0.5f*logical_h);
                                draw_line(chainPoints, cr_prevX, cX, cr_prevY, cY, cr_prevZ, cz, chainID, structure, depth_base_min_z, depth_base_max_z, logical_w, logical_h, 1);
                                cr_prevX = cX; cr_prevY = cY; cr_prevZ = cz;
                            }
                        } else {
                            draw_line(chainPoints, prevScreenX, screenX, prevScreenY, screenY, prevZ, z, chainID, structure, depth_base_min_z, depth_base_max_z, logical_w, logical_h, 1);
                        }
                        // Propagate current atom's metadata to all points added by draw_line
                        for (size_t k = before_draw; k < chainPoints.size(); ++k) {
                            chainPoints[k].bfactor            = cur_atom.bfactor;
                            chainPoints[k].is_interface       = cur_atom.is_interface;
                            chainPoints[k].is_aligned         = cur_atom.is_aligned;
                            chainPoints[k].conservation_score = cur_atom.conservation_score;
                            chainPoints[k].residue_number     = cur_atom.residue_number;
                            strncpy(chainPoints[k].residue_name, cur_atom.residue_name.c_str(), 3);
                            chainPoints[k].residue_name[3]    = '\0';
                        }
                    }

                    if (screenX >= 0 && screenX < logical_w && screenY >= 0 && screenY < logical_h) {
                        char pix = get_pixel_char_from_depth(z, depth_base_min_z, depth_base_max_z);
                        if (structure == 'x') {
                            // Coil: 3-pixel vertical cross node (center + top + bottom)
                            RenderPoint rp{screenX, screenY, z, pix, 0, chainID, structure};
                            rp.bfactor = cur_atom.bfactor; rp.is_interface = cur_atom.is_interface;
                            rp.is_aligned = cur_atom.is_aligned; rp.conservation_score = cur_atom.conservation_score;
                            rp.residue_number = cur_atom.residue_number;
                            strncpy(rp.residue_name, cur_atom.residue_name.c_str(), 3);
                            rp.residue_name[3] = '\0';
                            chainPoints.push_back(rp);
                            rp.y = screenY - 1;
                            if (rp.y >= 0) chainPoints.push_back(rp);
                            rp.y = screenY + 1;
                            if (rp.y < logical_h) chainPoints.push_back(rp);
                        } else {
                            // Helix/Sheet geometry: 5-pixel cross node
                            for (int oy = -1; oy <= 1; oy++)
                                for (int ox = -1; ox <= 1; ox++) {
                                    if (ox != 0 && oy != 0) continue;
                                    int nx = screenX + ox, ny = screenY + oy;
                                    if (nx >= 0 && nx < logical_w && ny >= 0 && ny < logical_h) {
                                        RenderPoint rp{nx, ny, z, pix, 0, chainID, structure};
                                        rp.bfactor = cur_atom.bfactor; rp.is_interface = cur_atom.is_interface;
                                        rp.is_aligned = cur_atom.is_aligned; rp.conservation_score = cur_atom.conservation_score;
                                        rp.residue_number = cur_atom.residue_number;
                                        strncpy(rp.residue_name, cur_atom.residue_name.c_str(), 3);
                                        rp.residue_name[3] = '\0';
                                        chainPoints.push_back(rp);
                                    }
                                }
                        }
                    }

                    prevScreenX = screenX;
                    prevScreenY = screenY;
                    prevZ = z;
                }
            }

            assign_colors_to_points(chainPoints, protein_idx);
            finalPoints.insert(finalPoints.end(), chainPoints.begin(), chainPoints.end());
            protein_idx++;
        }

        // z-buffer resolve into logicalPixels
        for (const auto& pt : finalPoints) {
            if (pt.x < 0 || pt.x >= logical_w || pt.y < 0 || pt.y >= logical_h) continue;
            int idx = pt.y * logical_w + pt.x;
            if (pt.depth < logicalPixels[idx].depth) {
                logicalPixels[idx] = pt;
            }
        }
        return;
    }

    // Standard (non-braille) path
    for (auto& px : screenPixels) {
        px.depth = std::numeric_limits<float>::infinity();
        px.pixel = ' ';
        px.color_id = 0;
    }

    std::vector<RenderPoint> finalPoints;
    std::vector<RenderPoint> chainPoints;
    finalPoints.reserve(200000);

    int protein_idx = 0;
    for (size_t ii = 0; ii < data.size(); ii++) {
        Protein* target = data[ii];
        chainPoints.clear();

        for (const auto& [chainID, chain_atoms] : target->get_atoms()) {
            if (chain_atoms.empty()) continue;

            int num_atoms = target->get_chain_length(chainID);
            int prevScreenX = -1, prevScreenY = -1;
            float prevZ = -1.0f;

            for (int i = 0; i < num_atoms; ++i) {
                float* position = chain_atoms[i].get_position();

                float x = position[0];
                float y = position[1];
                float z = position[2] + focal_offset;

                if (z < nearPlane) {
                    prevScreenX = prevScreenY = -1;
                    prevZ = -1.0f;
                    continue;
                }

                char structure = chain_atoms[i].get_structure();

                float projectedX = (x / z) * fovRads + pan_x[ii];
                float projectedY = (y / z) * fovRads + pan_y[ii];

                int screenX = (int)((projectedX + 1.0f) * 0.5f * screen_width);
                int screenY = (int)((1.0f - projectedY) * 0.5f * screen_height);

                const Atom& cur_atom = chain_atoms[i];
                if (prevScreenX != -1 && prevScreenY != -1) {
                    size_t before_draw = chainPoints.size();
                    if (structure != 'H' && structure != 'S') {
                        const Atom& P0 = chain_atoms[std::max(0, i-2)];
                        const Atom& P1 = chain_atoms[i-1];
                        const Atom& P2 = chain_atoms[i];
                        const Atom& P3 = chain_atoms[std::min(num_atoms-1, i+1)];
                        int cr_prevX = prevScreenX, cr_prevY = prevScreenY;
                        float cr_prevZ = prevZ;
                        for (int cr = 1; cr <= 4; ++cr) {
                            float t = cr * 0.25f;
                            float t2 = t*t, t3 = t2*t;
                            float cx = 0.5f*((-P0.x+3*P1.x-3*P2.x+P3.x)*t3 + (2*P0.x-5*P1.x+4*P2.x-P3.x)*t2 + (-P0.x+P2.x)*t + 2*P1.x);
                            float cy = 0.5f*((-P0.y+3*P1.y-3*P2.y+P3.y)*t3 + (2*P0.y-5*P1.y+4*P2.y-P3.y)*t2 + (-P0.y+P2.y)*t + 2*P1.y);
                            float cz = 0.5f*((-P0.z+3*P1.z-3*P2.z+P3.z)*t3 + (2*P0.z-5*P1.z+4*P2.z-P3.z)*t2 + (-P0.z+P2.z)*t + 2*P1.z) + focal_offset;
                            if (cz < nearPlane) break;
                            float cpX = (cx/cz)*fovRads + pan_x[ii];
                            float cpY = (cy/cz)*fovRads + pan_y[ii];
                            int cX = (int)((cpX+1.0f)*0.5f*screen_width);
                            int cY = (int)((1.0f-cpY)*0.5f*screen_height);
                            draw_line(chainPoints, cr_prevX, cX, cr_prevY, cY, cr_prevZ, cz, chainID, structure, depth_base_min_z, depth_base_max_z);
                            cr_prevX = cX; cr_prevY = cY; cr_prevZ = cz;
                        }
                    } else {
                        draw_line(chainPoints, prevScreenX, screenX, prevScreenY, screenY, prevZ, z, chainID, structure, depth_base_min_z, depth_base_max_z);
                    }
                    // Propagate current atom's metadata to all points added by draw_line
                    for (size_t k = before_draw; k < chainPoints.size(); ++k) {
                        chainPoints[k].bfactor            = cur_atom.bfactor;
                        chainPoints[k].is_interface       = cur_atom.is_interface;
                        chainPoints[k].is_aligned         = cur_atom.is_aligned;
                        chainPoints[k].conservation_score = cur_atom.conservation_score;
                        chainPoints[k].residue_number     = cur_atom.residue_number;
                        strncpy(chainPoints[k].residue_name, cur_atom.residue_name.c_str(), 3);
                        chainPoints[k].residue_name[3]    = '\0';
                    }
                }

                if (screenX >= 0 && screenX < screen_width && screenY >= 0 && screenY < screen_height) {
                    RenderPoint rp{screenX, screenY, z,
                                   get_pixel_char_from_depth(z, depth_base_min_z, depth_base_max_z),
                                   0, chainID, structure};
                    rp.bfactor = cur_atom.bfactor; rp.is_interface = cur_atom.is_interface;
                    rp.is_aligned = cur_atom.is_aligned; rp.conservation_score = cur_atom.conservation_score;
                    rp.residue_number = cur_atom.residue_number;
                    strncpy(rp.residue_name, cur_atom.residue_name.c_str(), 3);
                    rp.residue_name[3] = '\0';
                    chainPoints.push_back(rp);
                }

                prevScreenX = screenX;
                prevScreenY = screenY;
                prevZ = z;
            }
        }

        assign_colors_to_points(chainPoints, protein_idx);
        finalPoints.insert(finalPoints.end(), chainPoints.begin(), chainPoints.end());
        protein_idx++;
    }

    // z-buffer resolve — 전체 복사 필수: residue_number 등 모든 필드 포함
    for (const auto& pt : finalPoints) {
        if (pt.x < 0 || pt.x >= screen_width || pt.y < 0 || pt.y >= screen_height) continue;
        int idx = pt.y * screen_width + pt.x;
        if (pt.depth < screenPixels[idx].depth) {
            screenPixels[idx] = pt;  // RenderPoint 전체 복사 (residue_number 포함)
        }
    }
}

void Screen::project(std::vector<RenderPoint>& projectPixels, const int proj_width, const int proj_height) {
    const float nearPlane = 0.05f;
    float fovRads = 1.0f / std::tan((FOV / zoom_level) * 0.5f / 180.0f * PI);

    // init
    if ((int)projectPixels.size() != proj_width * proj_height) {
        projectPixels.assign(proj_width * proj_height, RenderPoint());
    }
    for (auto& px : projectPixels) {
        px.depth = std::numeric_limits<float>::infinity();
        px.pixel = ' ';
        px.color_id = 0;
    }

    if (!depth_calibrated) {
        calibrate_depth_baseline_first_view();
    }

    std::vector<RenderPoint> finalPoints;
    std::vector<RenderPoint> chainPoints;
    finalPoints.reserve(200000);

    int protein_idx = 0;
    for (size_t ii = 0; ii < data.size(); ii++) {
        Protein* target = data[ii];
        chainPoints.clear();

        for (const auto& [chainID, chain_atoms] : target->get_atoms()) {
            if (chain_atoms.empty()) continue;

            int num_atoms = target->get_chain_length(chainID);
            int prevScreenX = -1, prevScreenY = -1;
            float prevZ = -1.0f;

            for (int i = 0; i < num_atoms; ++i) {
                float* position = chain_atoms[i].get_position();

                float x = position[0];
                float y = position[1];
                float z = position[2] + focal_offset;

                if (z < nearPlane) {
                    prevScreenX = prevScreenY = -1;
                    prevZ = -1.0f;
                    continue;
                }

                char structure = chain_atoms[i].get_structure();

                float projectedX = (x / z) * fovRads + pan_x[ii];
                float projectedY = (y / z) * fovRads + pan_y[ii];

                int screenX = (int)((projectedX + 1.0f) * 0.5f * proj_width);
                int screenY = (int)((1.0f - projectedY) * 0.5f * proj_height);

                const Atom& cur_atom = chain_atoms[i];
                if (prevScreenX != -1 && prevScreenY != -1) {
                    size_t before_draw = chainPoints.size();
                    if (structure != 'H' && structure != 'S') {
                        const Atom& P0 = chain_atoms[std::max(0, i-2)];
                        const Atom& P1 = chain_atoms[i-1];
                        const Atom& P2 = chain_atoms[i];
                        const Atom& P3 = chain_atoms[std::min(num_atoms-1, i+1)];
                        int cr_prevX = prevScreenX, cr_prevY = prevScreenY;
                        float cr_prevZ = prevZ;
                        for (int cr = 1; cr <= 4; ++cr) {
                            float t = cr * 0.25f;
                            float t2 = t*t, t3 = t2*t;
                            float cx = 0.5f*((-P0.x+3*P1.x-3*P2.x+P3.x)*t3 + (2*P0.x-5*P1.x+4*P2.x-P3.x)*t2 + (-P0.x+P2.x)*t + 2*P1.x);
                            float cy = 0.5f*((-P0.y+3*P1.y-3*P2.y+P3.y)*t3 + (2*P0.y-5*P1.y+4*P2.y-P3.y)*t2 + (-P0.y+P2.y)*t + 2*P1.y);
                            float cz = 0.5f*((-P0.z+3*P1.z-3*P2.z+P3.z)*t3 + (2*P0.z-5*P1.z+4*P2.z-P3.z)*t2 + (-P0.z+P2.z)*t + 2*P1.z) + focal_offset;
                            if (cz < nearPlane) break;
                            float cpX = (cx/cz)*fovRads + pan_x[ii];
                            float cpY = (cy/cz)*fovRads + pan_y[ii];
                            int cX = (int)((cpX+1.0f)*0.5f*proj_width);
                            int cY = (int)((1.0f-cpY)*0.5f*proj_height);
                            draw_line(chainPoints, cr_prevX, cX, cr_prevY, cY, cr_prevZ, cz, chainID, structure, depth_base_min_z, depth_base_max_z);
                            cr_prevX = cX; cr_prevY = cY; cr_prevZ = cz;
                        }
                    } else {
                        draw_line(chainPoints, prevScreenX, screenX, prevScreenY, screenY, prevZ, z, chainID, structure, depth_base_min_z, depth_base_max_z);
                    }
                    // Propagate current atom's metadata to all points added by draw_line
                    for (size_t k = before_draw; k < chainPoints.size(); ++k) {
                        chainPoints[k].bfactor            = cur_atom.bfactor;
                        chainPoints[k].is_interface       = cur_atom.is_interface;
                        chainPoints[k].is_aligned         = cur_atom.is_aligned;
                        chainPoints[k].conservation_score = cur_atom.conservation_score;
                        chainPoints[k].residue_number     = cur_atom.residue_number;
                        strncpy(chainPoints[k].residue_name, cur_atom.residue_name.c_str(), 3);
                        chainPoints[k].residue_name[3]    = '\0';
                    }
                }

                if (screenX >= 0 && screenX < proj_width && screenY >= 0 && screenY < proj_height) {
                    RenderPoint rp{screenX, screenY, z,
                                   get_pixel_char_from_depth(z, depth_base_min_z, depth_base_max_z),
                                   0, chainID, structure};
                    rp.bfactor = cur_atom.bfactor; rp.is_interface = cur_atom.is_interface;
                    rp.is_aligned = cur_atom.is_aligned; rp.conservation_score = cur_atom.conservation_score;
                    rp.residue_number = cur_atom.residue_number;
                    strncpy(rp.residue_name, cur_atom.residue_name.c_str(), 3);
                    rp.residue_name[3] = '\0';
                    chainPoints.push_back(rp);
                }

                prevScreenX = screenX;
                prevScreenY = screenY;
                prevZ = z;
            }
        }

        assign_colors_to_points(chainPoints, protein_idx);
        finalPoints.insert(finalPoints.end(), chainPoints.begin(), chainPoints.end());
        protein_idx++;
    }

    for (const auto& pt : finalPoints) {
        if (pt.x < 0 || pt.x >= proj_width || pt.y < 0 || pt.y >= proj_height) continue;
        int idx = pt.y * proj_width + pt.x;

        if (pt.depth < projectPixels[idx].depth) {
            projectPixels[idx].depth = pt.depth;
            projectPixels[idx].pixel = pt.pixel;
            projectPixels[idx].color_id = pt.color_id;
        }
    }
}

void Screen::clear_screen() {
    screenPixels.assign(screen_width * screen_height, RenderPoint());
    if (use_braille) {
        logicalPixels.assign(screen_width * 2 * screen_height * 4, RenderPoint());
    }
}

void Screen::draw_screen(bool no_panel) {
    auto t0 = Benchmark::clock::now();
    clear_screen();
    project();

    int rows, cols;
    getmaxyx(stdscr, rows, cols);
    int panel_cols = std::min(cols, screen_width);

    int panel_h = panel->get_height_for_width(panel_cols);
    if (panel_h > rows) panel_h = rows;

    int offset = 0;
    if (!no_panel) {
        offset += panel_h;
    }
    if (offset > rows) offset = rows;

    clear();
    print_screen(offset);

    int start_row = rows;
    if (!no_panel) { 
        start_row -= panel_h;
    }
    if (start_row < 0) start_row = 0;

    for (int r = start_row; r < rows; ++r) {
        move(r, 0);
        clrtoeol();
    }
    if (!no_panel){
        panel->draw_panel(start_row, 0, panel_h, panel_cols);
    }
    // 기능 6: 패널 위치 저장 (hover 갱신 시 부분 재렌더링에 사용)
    last_panel_h         = no_panel ? 0 : panel_h;
    last_panel_start_row = no_panel ? rows : start_row;
    last_panel_cols      = panel_cols;
    last_no_panel        = no_panel;
    refresh();

    auto t1 = Benchmark::clock::now();
    int64_t render_dt_ms = Benchmark::ms_since(t0, t1);

    if (bm && bm->enabled) {
        if (!ttff_logged) {
            bm->log("ttff", -1, Benchmark::ms_since(bm->t0, t1), total_len_ca, (int64_t)data.size());
            ttff_logged = true;
        }
        bm->mark_frame_end(render_dt_ms, total_len_ca, (int64_t)data.size());
    }
}

void Screen::print_screen_braille(int y_offset) {
    // Render logicalPixels (2*W x 4*H) to terminal using Unicode Braille characters.
    // Each terminal cell covers a 2-wide x 4-tall sub-pixel block:
    //
    //   subcol=0  subcol=1          Braille dot numbering (bit positions):
    //   subrow=0: dot1(0)  dot4(3)
    //   subrow=1: dot2(1)  dot5(4)
    //   subrow=2: dot3(2)  dot6(5)
    //   subrow=3: dot7(6)  dot8(7)
    //
    // Unicode Braille: U+2800 + bitmask
    static const int dot_bits[2][4] = {
        {0, 1, 2, 6},  // left column  (subcol=0)
        {3, 4, 5, 7}   // right column (subcol=1)
    };

    int rows, cols;
    getmaxyx(stdscr, rows, cols);

    const int logical_w = screen_width * 2;
    const int logical_h = screen_height * 4;

    for (int ty = 0; ty < screen_height; ++ty) {
        int row = ty - (y_offset / 2) - 3;
        if (row < 0) continue;
        if (row >= rows) break;

        int max_cols = std::min(screen_width, cols);

        for (int tx = 0; tx < max_cols; ++tx) {
            int bitmask = 0;
            int best_color_id = 0;
            float best_depth = std::numeric_limits<float>::infinity();

            for (int sc = 0; sc < 2; ++sc) {
                for (int sr = 0; sr < 4; ++sr) {
                    int lx = tx * 2 + sc;
                    int ly = ty * 4 + sr;
                    if (lx >= logical_w || ly >= logical_h) continue;

                    const RenderPoint& lp = logicalPixels[ly * logical_w + lx];
                    if (lp.color_id > 0) {
                        bitmask |= (1 << dot_bits[sc][sr]);
                        if (lp.depth < best_depth) {
                            best_depth = lp.depth;
                            best_color_id = lp.color_id;
                        }
                    }
                }
            }

            if (bitmask > 0 && best_color_id > 0) {
                // Encode U+2800+bitmask as UTF-8 (3 bytes) directly — faster than setcchar/mvadd_wch
                char utf8[4];
                utf8[0] = (char)0xE2;
                utf8[1] = (char)(0xA0 | (bitmask >> 6));
                utf8[2] = (char)(0x80 | (bitmask & 0x3F));
                utf8[3] = '\0';
                attron(COLOR_PAIR(best_color_id));
                mvaddstr(row, tx, utf8);
                attroff(COLOR_PAIR(best_color_id));
            }
        }
    }
}

void Screen::print_screen(int y_offset) {
    if (use_braille) {
        print_screen_braille(y_offset);
        return;
    }

    int rows, cols;
    getmaxyx(stdscr, rows, cols);

    for (int i = 0; i < screen_height; ++i) {
        int row = i - (y_offset/2)-3;
        if (row < 0) continue;
        if (row >= rows) break;

        int max_width = std::min(screen_width, cols);

        for (int j = 0; j < max_width; ++j) {
            int idx = i * screen_width + j;
            const RenderPoint& px = screenPixels[idx];

            if (px.color_id > 0) {
                attron(COLOR_PAIR(px.color_id));
                mvaddch(row, j, px.pixel);
                attroff(COLOR_PAIR(px.color_id));
            } else {
                mvaddch(row, j, px.pixel);
            }
        }
    }
}

void Screen::set_zoom_level(float zoom){
    if ((zoom_level + zoom > 0.5)&&(zoom_level + zoom < 50)){
        float f_old = 1.0f / std::tan((FOV / zoom_level) * 0.5f / 180.0f * PI);
        zoom_level += zoom;
        float f_new = 1.0f / std::tan((FOV / zoom_level) * 0.5f / 180.0f * PI);
        float r = f_new / f_old;
        for (size_t i = 0; i < pan_x.size(); i++) {
            pan_x[i] *= r;
            pan_y[i] *= r;
        }
    }
}

// 기능 6: 마우스 커서 위치의 잔기 정보 검색 → 패널 hover 섹션 부분 갱신
void Screen::update_hover_info(int mx, int my) {
    // 터미널 좌표 → 스크린 좌표 변환
    // print_screen() 공식: row = screen_i - (y_offset/2) - 3
    // 역변환: screen_i = terminal_row + (panel_h/2) + 3
    int panel_offset = last_panel_h / 2 + 3;

    const RenderPoint* best = nullptr;
    float best_depth = std::numeric_limits<float>::infinity();

    if (use_braille) {
        int ty_screen = my + panel_offset;  // 브레일 스크린 행
        if (ty_screen >= 0 && ty_screen < screen_height && mx >= 0 && mx < screen_width) {
            const int logical_w = screen_width * 2;
            const int logical_h = screen_height * 4;
            for (int sc = 0; sc < 2; ++sc) {
                for (int sr = 0; sr < 4; ++sr) {
                    int lx = mx * 2 + sc;
                    int ly = ty_screen * 4 + sr;
                    if (lx >= logical_w || ly >= logical_h) continue;
                    const RenderPoint& lp = logicalPixels[ly * logical_w + lx];
                    if (lp.residue_number >= 0 && lp.depth < best_depth) {
                        best_depth = lp.depth;
                        best = &lp;
                    }
                }
            }
        }
    } else {
        int screen_y = my + panel_offset;
        int screen_x = mx;
        if (screen_x >= 0 && screen_x < screen_width && screen_y >= 0 && screen_y < screen_height) {
            const RenderPoint& px = screenPixels[screen_y * screen_width + screen_x];
            if (px.residue_number >= 0) {
                best = &px;
            }
        }
    }

    if (best) {
        panel->set_hover_residue(best->chainID, best->residue_name,
                                 best->residue_number, best->structure,
                                 best->bfactor, best->conservation_score);
    } else {
        panel->clear_hover_residue();
    }

    // 패널이 보이는 경우에만 hover 섹션을 부분 갱신
    if (!last_no_panel && last_panel_h > 0) {
        int rsh = panel->get_residue_section_height();
        // Residue Info 헤더 행 = 패널 하단에서 (rsh + 1) 줄 위 (1=bottom border)
        int hover_row = last_panel_start_row + last_panel_h - 1 - rsh;
        if (hover_row >= last_panel_start_row && hover_row + rsh <= last_panel_start_row + last_panel_h - 1) {
            panel->draw_hover_section(hover_row, last_panel_cols);
            refresh();
        }
    }
}

bool Screen::handle_input(bool& needs_redraw) {
    int key = getch();
    return handle_input_impl(key, needs_redraw);
}

bool Screen::handle_input(int key) {
    bool dummy = true;
    return handle_input_impl(key, dummy);
}

bool Screen::handle_input_impl(int key, bool& needs_redraw) {
    needs_redraw = true;  // 기본: 전체 재렌더링 필요

    // KEY_MOUSE: flushinp() 호출 전에 getmouse()를 먼저 처리
    // (flushinp()이 내부 mouse event 큐를 비워 getmouse()가 실패하는 것을 방지)
    if (key == KEY_MOUSE) {
        MEVENT event;
        if (getmouse(&event) == OK) {
            update_hover_info(event.x, event.y);
        }
        needs_redraw = false;  // 마우스 이동은 전체 재렌더링 불필요
        return true;
    }

    bool keep_show = true;

    auto pan_step_x = 2.0f * 4.0f / screen_width;
    auto pan_step_y = 2.0f * 2.0f / screen_height;

    auto apply_pan = [&](int idx, float dx, float dy){
        if (idx < 0 || idx >= (int)pan_x.size()) return;
        pan_x[idx] += dx;
        pan_y[idx] += dy;
    };

    if (bm && bm->enabled) bm->mark_event(key);
    flushinp();
    switch(key){
        // select protein
        case 48:
        case 49:
        case 50:
        case 51:
        case 52:
        case 53:
        case 54:
            if (key == 48){
                structNum = -1;
            }
            else if (key - 48 <= data.size()) {
                structNum = key - 49;
            }
            break;
        // A, a (minus x-axis)
        case 65:
        case 97:
            if (structNum != -1) apply_pan(structNum, -pan_step_x, 0.0f);
            else for (int i = 0; i < (int)data.size(); i++) apply_pan(i, -pan_step_x, 0.0f);
            break;
        // D, d (plus x-axis)
        case 68:
        case 100:
            if (structNum != -1) apply_pan(structNum, +pan_step_x, 0.0f);
            else for (int i = 0; i < (int)data.size(); i++) apply_pan(i, +pan_step_x, 0.0f);
            break;
        // S, s (minus y-axis)
        case 83:
        case 115:
            if (structNum != -1) apply_pan(structNum, 0.0f, -pan_step_y);
            else for (int i = 0; i < (int)data.size(); i++) apply_pan(i, 0.0f, -pan_step_y);
            break;    
        // W, w (plus y-axis)
        case 87:
        case 119:
            if (structNum != -1) apply_pan(structNum, 0.0f, +pan_step_y);
            else for (int i = 0; i < (int)data.size(); i++) apply_pan(i, 0.0f, +pan_step_y);
            break;

        // X, x (rotate x-centered)
        case 88:
        case 120:
            if (structNum != -1) {
                data[structNum]->set_rotate(1, 0, 0);
            } else if (yesUT) {
                float c = cos(PI / 48.0f), s = sin(PI / 48.0f);
                float m[9] = {1,0,0, 0,c,-s, 0,s,c};
                for (auto* p : data) p->do_naive_rotation(m);
            } else {
                for (int i = 0; i < (int)data.size(); i++) data[i]->set_rotate(1, 0, 0);
            }
            break;
        // Y, y (rotate y-centered)
        case 89:
        case 121:
            if (structNum != -1) {
                data[structNum]->set_rotate(0, 1, 0);
            } else if (yesUT) {
                float c = cos(PI / 48.0f), s = sin(PI / 48.0f);
                float m[9] = {c,0,s, 0,1,0, -s,0,c};
                for (auto* p : data) p->do_naive_rotation(m);
            } else {
                for (int i = 0; i < (int)data.size(); i++) data[i]->set_rotate(0, 1, 0);
            }
            break;
        // Z, z (rotate z-centered)
        case 90:
        case 122:
            if (structNum != -1) {
                data[structNum]->set_rotate(0, 0, 1);
            } else if (yesUT) {
                float c = cos(PI / 48.0f), s = sin(PI / 48.0f);
                float m[9] = {c,-s,0, s,c,0, 0,0,1};
                for (auto* p : data) p->do_naive_rotation(m);
            } else {
                for (int i = 0; i < (int)data.size(); i++) data[i]->set_rotate(0, 0, 1);
            }
            break;

        // F, f (zoom out)
        case 70:
        case 102:
            set_zoom_level(-0.3);
            break;   
        // R, R (zoom in)
        case 82:
        case 114:
            set_zoom_level(0.3);
            break;   

        // C, c (camera)
        case 67:
        case 99:     
        {     
            std::vector<RenderPoint> screenshotPixels;
            screenshotPixels.assign(screen_width * screen_height, RenderPoint());
            project(screenshotPixels, screen_width, screen_height);
            camera->screenshot(screenshotPixels);
            // camera->screenshot(screenPixels);
            break;
        }
        // N, n (next Foldseek hit)
        case 78:
        case 110:
            if (!foldseek_hits.empty()) load_next_hit(+1);
            break;

        // P, p (prev Foldseek hit)
        case 80:
        case 112:
            if (!foldseek_hits.empty()) load_next_hit(-1);
            break;

        // Q, q
        case 81:
        case 113:
            keep_show = false;
            break;

        default:
            break;
    }

    return keep_show;
}

// 기능 1: 로드된 모든 Protein에 대해 inter-chain interface를 계산
void Screen::compute_interface_all(float threshold) {
    for (Protein* p : data) {
        if (p) p->compute_interface(threshold);
    }
}

// 기능 5: 지정 protein에 conservation scores 적용
void Screen::apply_msa_conservation(int protein_idx, const std::vector<float>& scores) {
    if (protein_idx >= 0 && protein_idx < (int)data.size() && data[protein_idx]) {
        data[protein_idx]->apply_conservation_scores(scores);
    }
}

// 기능 4: 로드된 모든 Protein 쌍에 대해 nearest-neighbor 기반 정렬 잔기를 계산
void Screen::compute_aligned_all(float threshold) {
    if (data.size() < 2) {
        // 단일 단백질인 경우: 전체 잔기를 aligned로 표시
        for (Protein* p : data) {
            if (p) p->compute_aligned_regions_nn(*p, threshold);
        }
        return;
    }
    for (size_t i = 0; i < data.size(); ++i) {
        for (size_t j = i + 1; j < data.size(); ++j) {
            if (data[i] && data[j]) {
                data[i]->compute_aligned_regions_nn(*data[j], threshold);
            }
        }
    }
    panel->set_align_method("nearest-nbr");
}

// 기능 4: -fs 기반 — Foldseek hit의 U/T transform을 지정 protein에 적용
// T_norm: 정규화 공간 T (screen_atoms 이동용)
// T_angstrom: Å 공간 T (init_atoms 갱신용). nullptr이면 init_atoms 미갱신.
//   → compute_aligned_regions_from_aln은 init_atoms 기준 거리 비교를 하므로
//     반드시 Å 공간 T를 전달해야 aligned 색상이 올바르게 동작함.
void Screen::apply_foldseek_transform(int protein_idx, const float* U_flat,
                                      const float* T_norm, const float* T_angstrom) {
    if (protein_idx < 0 || protein_idx >= (int)data.size() || !data[protein_idx]) return;
    // screen_atoms에 적용 (시각적 정렬)
    data[protein_idx]->do_naive_rotation(const_cast<float*>(U_flat));
    data[protein_idx]->do_shift(const_cast<float*>(T_norm));
    // init_atoms에 Å 공간 T로 적용 (거리 비교 기준)
    if (T_angstrom) {
        data[protein_idx]->apply_ut_to_init_atoms(U_flat, T_angstrom);
    }
    yesUT = true;
}

// 기능 4: -fs 기반 — alignment string으로 aligned 잔기 계산 (protein0 vs protein1)
void Screen::compute_aligned_from_aln(const std::string& qaln, const std::string& taln,
                                      float threshold) {
    if ((int)data.size() < 2 || !data[0] || !data[1]) return;
    data[0]->compute_aligned_regions_from_aln(*data[1], qaln, taln, threshold);
}

// 기능 4: 패널에 정렬 방식 표시 설정
void Screen::set_align_method(const std::string& method) {
    if (panel) panel->set_align_method(method);
}

// ── 기능 3: Foldseek hit 탐색 ──────────────────────────────────────────────────

// 3×3 대칭 행렬 Jacobi 고유값 분해
// a: 입력 (수정됨), d: 고유값, v: 고유벡터 (열 = 고유벡터)
static void jacobi3_sym(float a[3][3], float d[3], float v[3][3]) {
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) v[i][j] = (i == j) ? 1.f : 0.f;
    for (int iter = 0; iter < 100; ++iter) {
        int p = 0, q = 1;
        float max_val = std::abs(a[0][1]);
        if (std::abs(a[0][2]) > max_val) { max_val = std::abs(a[0][2]); p = 0; q = 2; }
        if (std::abs(a[1][2]) > max_val) { max_val = std::abs(a[1][2]); p = 1; q = 2; }
        if (max_val < 1e-9f) break;
        float diff = a[q][q] - a[p][p];
        float t;
        if (std::abs(diff) < 1e-9f * std::abs(a[p][q])) {
            t = 1.f;
        } else {
            float phi = diff / (2.f * a[p][q]);
            t = (phi > 0) ? (1.f / (phi + std::sqrt(1.f + phi*phi)))
                           : (1.f / (phi - std::sqrt(1.f + phi*phi)));
        }
        float c = 1.f / std::sqrt(1.f + t*t);
        float s = t * c;
        a[p][p] -= t * a[p][q];
        a[q][q] += t * a[p][q];
        a[p][q] = a[q][p] = 0.f;
        int r3 = 3 - p - q;
        float apr = a[p][r3], aqr = a[q][r3];
        a[p][r3] = a[r3][p] =  c * apr - s * aqr;
        a[q][r3] = a[r3][q] =  s * apr + c * aqr;
        for (int i = 0; i < 3; i++) {
            float vp = v[i][p], vq = v[i][q];
            v[i][p] = c * vp - s * vq;
            v[i][q] = s * vp + c * vq;
        }
    }
    for (int i = 0; i < 3; i++) d[i] = a[i][i];
}

static float mat3_det_flat(const float m[9]) {
    return m[0]*(m[4]*m[8]-m[5]*m[7])
          -m[1]*(m[3]*m[8]-m[5]*m[6])
          +m[2]*(m[3]*m[7]-m[4]*m[6]);
}

// Kabsch 알고리즘: U*Q[i]+T ≈ P[i] 최소화
// P: 참조 좌표 (Nx3 flat), Q: 회전 대상 좌표 (Nx3 flat), N: 쌍 수
// 결과: U[9] (row-major 3×3), T[3]
static void kabsch(const std::vector<float>& P, const std::vector<float>& Q,
                   int N, float U[9], float T[3]) {
    // 단위 행렬로 초기화
    for (int i = 0; i < 9; i++) U[i] = (i % 4 == 0) ? 1.f : 0.f;
    T[0] = T[1] = T[2] = 0.f;
    if (N < 3) return;

    // 무게중심
    float Pc[3] = {}, Qc[3] = {};
    for (int i = 0; i < N; i++) {
        Pc[0] += P[i*3]; Pc[1] += P[i*3+1]; Pc[2] += P[i*3+2];
        Qc[0] += Q[i*3]; Qc[1] += Q[i*3+1]; Qc[2] += Q[i*3+2];
    }
    for (int k = 0; k < 3; k++) { Pc[k] /= N; Qc[k] /= N; }

    // H = Σ (Q[i]-Qc)(P[i]-Pc)^T  [3×3]
    float H[3][3] = {};
    for (int i = 0; i < N; i++) {
        float q[3] = { Q[i*3]-Qc[0], Q[i*3+1]-Qc[1], Q[i*3+2]-Qc[2] };
        float p[3] = { P[i*3]-Pc[0], P[i*3+1]-Pc[1], P[i*3+2]-Pc[2] };
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++)
                H[r][c] += q[r] * p[c];
    }

    // H^T H (대칭) → 고유값 분해로 V (오른쪽 특이벡터) 계산
    float HtH[3][3] = {};
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++)
            for (int k = 0; k < 3; k++)
                HtH[r][c] += H[k][r] * H[k][c];

    float d[3];
    float Vmat[3][3];
    jacobi3_sym(HtH, d, Vmat);

    // 고유값 내림차순 정렬
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2-i; j++) {
            if (d[j] < d[j+1]) {
                std::swap(d[j], d[j+1]);
                for (int k = 0; k < 3; k++) std::swap(Vmat[k][j], Vmat[k][j+1]);
            }
        }
    }

    // 특이값
    float sigma[3];
    for (int i = 0; i < 3; i++) sigma[i] = (d[i] > 0.f) ? std::sqrt(d[i]) : 0.f;

    // 왼쪽 특이벡터: Umat[:,i] = H * Vmat[:,i] / sigma[i]
    float Umat[3][3] = {};
    for (int i = 0; i < 3; i++) {
        if (sigma[i] < 1e-7f) continue;
        for (int r = 0; r < 3; r++) {
            for (int k = 0; k < 3; k++) Umat[r][i] += H[r][k] * Vmat[k][i];
            Umat[r][i] /= sigma[i];
        }
    }
    // 退화 처리: 세 번째 열 = 교차곱
    if (sigma[2] < 1e-7f) {
        Umat[0][2] = Umat[1][0]*Umat[2][1] - Umat[2][0]*Umat[1][1];
        Umat[1][2] = Umat[2][0]*Umat[0][1] - Umat[0][0]*Umat[2][1];
        Umat[2][2] = Umat[0][0]*Umat[1][1] - Umat[1][0]*Umat[0][1];
    }

    // det(V*U^T) 반사 보정
    float VUt[9];
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++) {
            VUt[r*3+c] = 0;
            for (int k = 0; k < 3; k++) VUt[r*3+c] += Vmat[r][k] * Umat[c][k];
        }
    float det_sign = (mat3_det_flat(VUt) < 0) ? -1.f : 1.f;

    // R = V * diag(1,1,det_sign) * U^T
    float diag[3] = {1.f, 1.f, det_sign};
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++) {
            U[r*3+c] = 0;
            for (int k = 0; k < 3; k++) U[r*3+c] += Vmat[r][k] * diag[k] * Umat[c][k];
        }

    // T = Pc - R*Qc
    for (int r = 0; r < 3; r++) {
        T[r] = Pc[r];
        for (int c = 0; c < 3; c++) T[r] -= U[r*3+c] * Qc[c];
    }
}

void Screen::set_foldseek_hits(const std::vector<FoldseekHit>& hits) {
    foldseek_hits = hits;
    current_hit_idx = -1;
    // 패널에 총 hit 수 표시
    if (panel && !hits.empty()) {
        FoldseekHitInfo fi;
        fi.valid = true;
        fi.current_idx = 0;
        fi.total_hits = (int)hits.size();
        panel->set_foldseek_hit_info(fi);
    }
}

void Screen::set_fs_db_path(const std::string& path) {
    fs_db_path = path;
}

void Screen::load_next_hit(int delta) {
    if (foldseek_hits.empty()) return;

    int new_idx;
    if (current_hit_idx < 0) {
        // 초기 로드: delta > 0 이면 0번, delta < 0 이면 마지막
        new_idx = (delta >= 0) ? 0 : (int)foldseek_hits.size() - 1;
    } else {
        new_idx = current_hit_idx + delta;
        new_idx = std::max(0, std::min((int)foldseek_hits.size() - 1, new_idx));
        if (new_idx == current_hit_idx) return;
    }
    current_hit_idx = new_idx;

    const FoldseekHit& hit = foldseek_hits[current_hit_idx];

    // 파일 경로 결정 (다운로드 포함)
    std::string status_msg;
    std::string file_path = PDBDownloader::resolve_target_file(
        hit.target, fs_db_path, status_msg);

    // 패널 정보 갱신
    FoldseekHitInfo fs_info;
    fs_info.valid      = true;
    fs_info.current_idx = current_hit_idx + 1;  // 1-based
    fs_info.total_hits  = (int)foldseek_hits.size();
    fs_info.target      = hit.target;
    fs_info.evalue      = hit.evalue;
    fs_info.prob        = hit.prob;
    fs_info.lddt        = hit.lddt;
    fs_info.qtmscore    = hit.qtmscore;
    fs_info.status_msg  = status_msg;
    if (panel) panel->set_foldseek_hit_info(fs_info);

    if (file_path.empty()) return;

    // 기존 target protein (index 1+) 제거
    while ((int)data.size() > 1) {
        delete data.back();
        data.pop_back();
        pan_x.pop_back();
        pan_y.pop_back();
    }
    while ((int)chainVec.size() > 1) chainVec.pop_back();

    // 체인 필터 결정
    DBType db_type = PDBDownloader::detect_db_type(hit.target);
    std::string chain_filter = PDBDownloader::extract_chain(hit.target, db_type);
    chainVec.push_back(chain_filter);

    // 새 target protein 로드
    Protein* target_protein = new Protein(file_path, chain_filter, screen_show_structure);
    data.push_back(target_protein);
    pan_x.push_back(0.0f);
    pan_y.push_back(0.0f);

    float tvec[3] = {0.f, 0.f, 0.f};
    target_protein->load_data(tvec, false);

    // 패널 entry 갱신
    if (panel) {
        panel->update_entry(1, target_protein->get_file_name(),
                            target_protein->get_chain_length(),
                            target_protein->get_residue_count());
    }

    // query와 동일한 norm_scale 로 target 정규화
    target_protein->set_bounding_box();
    target_protein->set_scale(norm_scale);
    float t_cx = target_protein->cx;
    float t_cy = target_protein->cy;
    float t_cz = target_protein->cz;
    float t_shift[3] = { -t_cx, -t_cy, -t_cz };
    target_protein->do_shift(t_shift);
    target_protein->do_scale(norm_scale);

    // U/T transform 계산
    float U[9] = {1,0,0, 0,1,0, 0,0,1};
    float T[3]  = {0,0,0};
    float T_ang[3] = {0,0,0};  // Å 공간 T (init_atoms 갱신용)
    bool computed_transform = false;
    std::string align_method_str;

    if (hit.has_transform) {
        // 29컬럼 포맷: U/T를 정규화 공간으로 변환
        // target_norm = (target_orig - t_centroid) * norm_scale
        // alns_norm   = (alns_orig   - q_centroid) * norm_scale
        // U_norm      = U_fs (회전은 스케일 불변)
        // T_norm      = (U_fs*t_centroid + T_fs - q_centroid) * norm_scale
        const float* Uf = hit.U;
        const float* Tf = hit.T;
        float Utc[3] = {
            Uf[0]*t_cx + Uf[1]*t_cy + Uf[2]*t_cz,
            Uf[3]*t_cx + Uf[4]*t_cy + Uf[5]*t_cz,
            Uf[6]*t_cx + Uf[7]*t_cy + Uf[8]*t_cz
        };
        for (int i = 0; i < 9; i++) U[i] = Uf[i];
        T[0] = (Utc[0] + Tf[0] - norm_cx) * norm_scale;
        T[1] = (Utc[1] + Tf[1] - norm_cy) * norm_scale;
        T[2] = (Utc[2] + Tf[2] - norm_cz) * norm_scale;
        // init_atoms용: Foldseek 원본 Å 공간 T (hit.T)를 그대로 사용
        T_ang[0] = Tf[0]; T_ang[1] = Tf[1]; T_ang[2] = Tf[2];
        computed_transform = true;
        align_method_str = "aln-string";

    } else if (hit.is_alis_format && !hit.alns.empty() && hit.has_aln) {
        // alis 21컬럼 포맷: Kabsch SVD로 U/T 역산
        // alns는 원래 query frame 좌표 → 정규화: (alns - q_centroid) * norm_scale
        // target CA 좌표는 이미 정규화됨 (do_shift + do_scale 후)

        const auto& atoms_map = target_protein->get_atoms();

        // target CA atom 플랫 리스트 (정규화된 좌표)
        std::vector<std::array<float,3>> target_cas;
        for (const auto& [cid, chain] : atoms_map) {
            for (const auto& atom : chain) {
                float* pos = const_cast<Atom&>(atom).get_position();
                target_cas.push_back({pos[0], pos[1], pos[2]});
            }
        }

        // taln 순회: 비-갭 위치마다 쌍 수집
        std::vector<float> P_norm, Q_norm;
        int aln_idx = 0;  // alns 인덱스 (3 float/잔기)
        int t_seq_idx = hit.tstart - 1;  // target CA 0-based 인덱스

        for (size_t ai = 0; ai < hit.taln.size(); ai++) {
            if (hit.taln[ai] == '-') continue;
            if (t_seq_idx < (int)target_cas.size() &&
                aln_idx * 3 + 2 < (int)hit.alns.size()) {
                // 정규화된 alns
                P_norm.push_back((hit.alns[aln_idx*3]   - norm_cx) * norm_scale);
                P_norm.push_back((hit.alns[aln_idx*3+1] - norm_cy) * norm_scale);
                P_norm.push_back((hit.alns[aln_idx*3+2] - norm_cz) * norm_scale);
                // 정규화된 target CA
                Q_norm.push_back(target_cas[t_seq_idx][0]);
                Q_norm.push_back(target_cas[t_seq_idx][1]);
                Q_norm.push_back(target_cas[t_seq_idx][2]);
            }
            aln_idx++;
            t_seq_idx++;
        }

        int N = (int)std::min(P_norm.size(), Q_norm.size()) / 3;
        if (N >= 3) {
            kabsch(P_norm, Q_norm, N, U, T);
            // BUG-A 2단계: kabsch는 정규화 공간에서 계산됨.
            // init_atoms는 Å 공간이므로 T_Å = T_norm/norm_scale + q_centroid - U*t_centroid
            {
                const float q_cen[3] = {norm_cx, norm_cy, norm_cz};
                const float t_cen[3] = {t_cx, t_cy, t_cz};
                for (int r = 0; r < 3; r++) {
                    T_ang[r] = T[r] / norm_scale + q_cen[r];
                    for (int c = 0; c < 3; c++) T_ang[r] -= U[r*3+c] * t_cen[c];
                }
            }
            computed_transform = true;
            align_method_str = "kabsch-alns";
        }
    }

    if (computed_transform) {
        apply_foldseek_transform(1, U, T, T_ang);
    }

    // 정렬 방식 패널 갱신
    fs_info.align_method = align_method_str;
    if (panel) panel->set_foldseek_hit_info(fs_info);

    // aligned 모드일 때 is_aligned 계산
    if (screen_mode == "aligned") {
        if (hit.has_aln) {
            compute_aligned_from_aln(hit.qaln, hit.taln, 5.0f);
            set_align_method("aln-string");
        } else {
            compute_aligned_all();
        }
    }

    depth_calibrated = false;
}

// ── 기능 8: FoldMason MSA 기반 superposition ─────────────────────────────────

void Screen::set_foldmason(std::unique_ptr<FoldMasonParser> parser) {
    foldmason_parser = std::move(parser);
}

void Screen::set_foldmason_panel_info(const FoldMasonInfo& info) {
    if (panel) panel->set_foldmason_info(info);
}

void Screen::apply_foldmason_superposition(int query_protein_idx, int target_protein_idx,
                                           int fm_query_entry_idx, int fm_target_entry_idx) {
    if (!foldmason_parser) return;
    if (query_protein_idx < 0 || query_protein_idx >= (int)data.size() || !data[query_protein_idx]) return;
    if (target_protein_idx < 0 || target_protein_idx >= (int)data.size() || !data[target_protein_idx]) return;

    const auto& entries = foldmason_parser->get_entries();
    if (fm_query_entry_idx >= (int)entries.size() || fm_target_entry_idx >= (int)entries.size()) return;

    auto pairs = foldmason_parser->build_aligned_pairs(fm_query_entry_idx, fm_target_entry_idx);
    if (pairs.empty()) return;

    // query/target protein CA atoms 플랫화 (screen_atoms 순서)
    std::vector<std::array<float,3>> query_cas;
    for (auto& [cid, chain] : data[query_protein_idx]->get_atoms()) {
        for (const auto& a : chain) query_cas.push_back({a.x, a.y, a.z});
    }
    std::vector<std::array<float,3>> target_cas;
    for (auto& [cid, chain] : data[target_protein_idx]->get_atoms()) {
        for (const auto& a : chain) target_cas.push_back({a.x, a.y, a.z});
    }

    const int q_size = (int)query_cas.size();
    const int t_size = (int)target_cas.size();

    // P = query (참조), Q = target (회전 대상)
    std::vector<float> P_flat, Q_flat;
    for (const auto& [ref_res, oth_res] : pairs) {
        if (ref_res >= q_size || oth_res >= t_size) continue;
        P_flat.push_back(query_cas[ref_res][0]);
        P_flat.push_back(query_cas[ref_res][1]);
        P_flat.push_back(query_cas[ref_res][2]);
        Q_flat.push_back(target_cas[oth_res][0]);
        Q_flat.push_back(target_cas[oth_res][1]);
        Q_flat.push_back(target_cas[oth_res][2]);
    }

    int N = (int)std::min(P_flat.size(), Q_flat.size()) / 3;
    if (N < 3) return;

    float U[9], T[3];
    kabsch(P_flat, Q_flat, N, U, T);

    // BUG-A 2단계: kabsch는 정규화 screen_atoms 공간에서 계산됨.
    // init_atoms는 Å 공간이므로 T_Å = T_norm/norm_scale + q_centroid - U*t_centroid
    float T_ang[3];
    {
        const float q_cen[3] = {norm_cx, norm_cy, norm_cz};
        const float t_cen[3] = {data[target_protein_idx]->cx,
                                 data[target_protein_idx]->cy,
                                 data[target_protein_idx]->cz};
        for (int r = 0; r < 3; r++) {
            T_ang[r] = T[r] / norm_scale + q_cen[r];
            for (int c = 0; c < 3; c++) T_ang[r] -= U[r*3+c] * t_cen[c];
        }
    }
    apply_foldseek_transform(target_protein_idx, U, T, T_ang);

    // aligned 모드일 때 is_aligned 잔기 설정
    // MSA aa strings을 qaln/taln으로 사용 (gap 형식 동일)
    if (screen_mode == "aligned") {
        data[query_protein_idx]->compute_aligned_regions_from_aln(
            *data[target_protein_idx],
            entries[fm_query_entry_idx].aa,
            entries[fm_target_entry_idx].aa,
            5.0f);
        set_align_method("msa-col");
    }

    depth_calibrated = false;
}