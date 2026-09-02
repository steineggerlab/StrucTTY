#pragma once

#include <string>
#include <map>
#include <vector>
#include <tuple>
#include <cmath>
#include <cctype>
#include <limits>
#include <algorithm>

#include <gemmi/mmread.hpp>
#include <gemmi/model.hpp>
#include <gemmi/metadata.hpp>

#include "Atom.hpp"
#include "StructureMaker.hpp"
#include "SSPredictor.hpp"

struct BoundingBox {
    float min_x = std::numeric_limits<float>::max();
    float min_y = std::numeric_limits<float>::max();
    float min_z = std::numeric_limits<float>::max();
    float max_x = std::numeric_limits<float>::lowest();
    float max_y = std::numeric_limits<float>::lowest();
    float max_z = std::numeric_limits<float>::lowest();

    void expand(float x, float y, float z) {
        min_x = std::min(min_x, x);
        max_x = std::max(max_x, x);
        min_y = std::min(min_y, y);
        max_y = std::max(max_y, y);
        min_z = std::min(min_z, z);
        max_z = std::max(max_z, z);
    }

    BoundingBox operator+(const BoundingBox& other) const {
        BoundingBox result;
        result.min_x = std::min(min_x, other.min_x);
        result.min_y = std::min(min_y, other.min_y);
        result.min_z = std::min(min_z, other.min_z);
        result.max_x = std::max(max_x, other.max_x);
        result.max_y = std::max(max_y, other.max_y);
        result.max_z = std::max(max_z, other.max_z);
        return result;
    }
};

class Protein {
public:
    Protein(const std::string& in_file_, const std::string& target_chains_, const bool& show_structure_);
    Protein(const std::string& name, bool show_structure);
    ~Protein();

    std::map<std::string, std::vector<Atom>>& get_atoms();

    const std::map<std::string, std::vector<Atom>>& get_ca_atoms() const { return init_atoms; }
    std::map<std::string, int> get_residue_count();
    std::map<std::string, int> get_chain_length();
    int get_chain_length(std::string chainID);
    int get_length();
    void set_bounding_box();

    BoundingBox& get_bounding_box();
    void set_scale(float scale_);
    std::string get_file_name() { return in_file; }

    void load_data(float * vectorpointers, bool yesUT);
    void load_from_ca(const std::vector<float>& coords, size_t n_residues,
                      const std::string& aa_seq);
    bool is_ca_only() const { return ca_only_; }

    void compute_interface(float threshold = 8.0f);

    void apply_ut_to_init_atoms(const float* U, const float* T);

    void compute_aligned_regions_from_aln(Protein& other,
                                          const std::string& qaln,
                                          const std::string& taln,
                                          int q_start = 1,
                                          int t_start = 1,
                                          float threshold = 10.0f,
                                          bool skip_distance_check = false);

    void compute_aligned_regions_nn(Protein& other, float threshold = 10.0f);

    void apply_conservation_scores(const std::vector<float>& scores);

    void set_rotate(int x_rotate, int y_rotate, int z_rotate);
    void set_shift(float shift_x, float shift_y, float shift_z);
    void do_naive_rotation(float* rotate_mat);
    void do_rotation(float* rotate_mat);
    void do_shift(float* shift_mat);
    void do_scale(float sclae);
    float cx, cy, cz, scale;

private:
    void count_seqres(const std::string& file);
    bool is_ss_in_file(const std::string& in_file);
    void load_ss_info(const std::string& in_file,
                      const std::string& target_chains,
                      std::vector<std::tuple<std::string, int, std::string, int, char>>& ss_info);
    void load_init_atoms(const std::string& in_file,
                             const std::string& target_chains,
                             const std::vector<std::tuple<std::string, int, std::string, int, char>>& ss_info, float * vectorpointers, bool yesUT);
    void load_init_atoms(const std::string& in_file,
                             const std::string& target_chains, float * vectorpointers, bool yesUT);

    void pred_ss_info(std::map<std::string, std::vector<Atom>>& init_atoms);
    void compute_interface_pair(const std::string& chain_A, const std::string& chain_B, float threshold);
    void sync_interface_to_screen();
    void sync_aligned_to_screen();
    void sync_conservation_to_screen();

    std::map<std::string, std::vector<Atom>> init_atoms;
    std::map<std::string, std::vector<Atom>> screen_atoms;

    std::map<std::string, int> chain_res_count;

    std::string in_file;
    std::string target_chains;
    bool show_structure;
    bool ca_only_ = false;

    BoundingBox bounding_box;

    StructureMaker structureMaker;
    SSPredictor ssPredictor;
};

