#pragma once
#include "Protein.hpp"
#include "Atom.hpp"
#include "RenderPoint.hpp"
#include "Palette.hpp"
#include "Camera.hpp"
#include "Panel.hpp"
#include "Benchmark.hpp"
#include "Renderer.hpp"
#include "../structure/FoldseekParser.hpp"
#include "../structure/PDBDownloader.hpp"
#include "../structure/FoldMasonParser.hpp"
#include "../structure/FoldseekDBReader.hpp"
#include "../structure/MultimerReportParser.hpp"
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <limits>
#include <map>
#include <unordered_map>
#include <memory>
#include <string>

class Screen {
public:
    Screen(const int& width, const int& height, const bool& show_structure,
           const std::string& mode);
    ~Screen();

    bool handle_input(bool& needs_redraw);
    bool handle_input(int key);

    void set_protein(const std::string& in_file, int ii, const bool& show_structure);

    std::string apply_accession_chain(int idx, const std::string& accession,
                                     const std::string& file_path);

    static std::string chain_from_accession(const std::string& accession,
                                            const std::string& file_path);

    bool set_query_from_db(const std::string& query_db_path,
                           const std::string& accession,
                           const bool& show_structure);

    void set_query_nav(const std::vector<std::string>& query_ids,
                       const std::map<std::string, std::vector<FoldseekHit>>& hits_by_query,
                       const std::string& query_db_path,
                       const std::string& target_db_path,
                       const bool& show_structure);
    void set_query_nav_from_file(const std::vector<std::string>& query_ids,
                                 const std::map<std::string, std::vector<FoldseekHit>>& hits_by_query,
                                 const std::string& query_source,
                                 const std::string& target_db_path,
                                 const bool& show_structure,
                                 const bool& source_is_directory = false);
    void switch_query(int delta);
    int  query_count() const { return (int)query_ids_.size(); }

    void set_multimer_report(const std::vector<MultimerHit>& hits,
                             const std::string& query_source,
                             bool query_is_db,
                             bool query_is_dir,
                             const std::string& target_db_path,
                             const bool& show_structure);

    void normalize_proteins();

    void set_tmatrix();
    void set_chainfile(const std::string& chainfile, int filesize);
    void set_zoom_level(float zoom);

    void compute_interface_all(float threshold = 8.0f);

    void compute_aligned_all(float threshold = 10.0f);

    void apply_msa_conservation(int protein_idx, const std::vector<float>& scores);

    void apply_foldseek_transform(int protein_idx, const float* U_flat, const float* T_norm,
                                  const float* T_angstrom = nullptr);

    void compute_aligned_from_aln(const std::string& qaln, const std::string& taln,
                                  int q_start = 1, int t_start = 1,
                                  float threshold = 5.0f,
                                  bool skip_distance_check = false);

    void set_align_method(const std::string& method);

    void set_foldseek_hits(const std::vector<FoldseekHit>& hits);
    void set_fs_db_path(const std::string& path);
    void load_next_hit(int delta);

    bool open_foldseek_db(const std::string& db_base_path);
    bool prepare_foldseek_db(const std::vector<FoldseekHit>& hits);

    void apply_hit_transform(int target_protein_idx, const FoldseekHit& hit);

    void set_foldmason(std::unique_ptr<FoldMasonParser> parser);
    void apply_foldmason_superposition(int query_protein_idx, int target_protein_idx,
                                       int fm_query_entry_idx, int fm_target_entry_idx);
    void set_foldmason_panel_info(const FoldMasonInfo& info);

    void draw_screen(bool no_panel);

    void update_hover_info(int mx, int my);

    void set_benchmark(Benchmark* b) { bm = b; }

    void update_total_len_ca() {
        int64_t sum = 0;
        for (auto* p : data) {
            if (!p) continue;
            sum += (int64_t)p->get_length();
        }
        total_len_ca = sum;
    }

private:
    int screen_width, screen_height;
    float aspect_ratio;
    bool screen_show_structure;
    bool yesUT = false;
    std::string screen_mode;
    int structNum = -1;

    Renderer renderer_;
    void print_screen_braille(int y_offset);

    std::vector<std::string> chain_names_;
    std::unordered_map<std::string, int> chain_intern_;
    int intern_chain(const std::string& id);
    const std::string& chain_name(int idx) const;

    float focal_offset = 3.0f;
    float zoom_level;

    std::vector<float> pan_x;
    std::vector<float> pan_y;
    std::vector<std::string> chainVec;
    float** vectorpointer = nullptr;

    std::vector<Protein*> data;

    BoundingBox global_bb;
    Camera* camera = nullptr;
    Panel* panel = nullptr;

    bool depth_calibrated = false;
    float depth_base_min_z = 0.0f;
    float depth_base_max_z = 1.0f;

    Benchmark* bm = nullptr;
    bool ttff_logged = false;

    int last_panel_h = 0;
    int last_panel_start_row = 0;
    int last_panel_cols = 0;
    bool last_no_panel = false;

    std::vector<FoldseekHit> foldseek_hits;
    int current_hit_idx = -1;
    std::string fs_db_path;

    FoldseekDBReader fs_db_reader_;

    FoldseekDBReader query_db_reader_;

    std::vector<std::string> query_ids_;
    std::map<std::string, std::vector<FoldseekHit>> hits_by_query_;
    int  current_query_idx_ = -1;
    std::string query_db_path_;
    std::string query_file_;
    bool query_file_is_dir_ = false;
    std::string target_db_path_;
    bool multi_query_show_structure_ = false;

    void activate_query(int idx);
    bool load_query_into_data0(const std::string& accession, const bool& show_structure);

    bool multimer_mode_ = false;
    std::vector<std::string> mm_query_complexes_;
    std::map<std::string, std::vector<MultimerHit>> mm_hits_by_query_;
    int mm_current_query_idx_ = -1;
    int mm_current_hit_idx_ = -1;
    int mm_query_chain_count_ = 0;
    std::string mm_query_source_;
    bool mm_query_is_db_ = true;
    bool mm_query_is_dir_ = false;
    std::vector<std::string> mm_chain_labels_;
    void activate_multimer_query(int idx);
    void load_multimer_hit(int delta);
    void switch_multimer_query(int delta);
    void normalize_complex();
    bool load_chain_into_data(FoldseekDBReader& reader, const std::string& accession,
                              const bool& show_structure);
    int load_complex_chains(const std::string& complex,
                            const std::vector<std::string>& chains,
                            bool from_db, FoldseekDBReader& reader,
                            const std::string& source_path, bool source_is_dir,
                            const bool& show_structure);
    void transform_target_chain(int idx, const float U[9], const float T[3]);
    void free_tmatrix();
    std::pair<int,int> mm_complex_range(int complex_idx) const;
    size_t vectorpointer_len_ = 0;

    std::unique_ptr<FoldMasonParser> foldmason_parser;
    float norm_scale = 1.0f;
    float norm_cx = 0.0f;
    float norm_cy = 0.0f;
    float norm_cz = 0.0f;
    float rot_pivot_[3] = {0.f, 0.f, 0.f};

    void calibrate_depth_baseline_first_view();

    std::vector<RenderAtom> to_render_atoms();
    void print_screen(int panel_lines);

    bool handle_input_impl(int key, bool& needs_redraw);

    int64_t total_len_ca = 0;
};
