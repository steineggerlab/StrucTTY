#pragma once
#include <vector>
#include <map>
#include <cmath>
#include <cstddef>
#include "Atom.hpp"

class SSPredictor {
public:
    float scale = 1.0f;

    float helix_tau_mid   = 89.0f,  helix_tau_tol   = 12.0f;
    float helix_alpha_mid = 50.0f,  helix_alpha_tol = 20.0f;
    float helix_d2_mid    = 5.5f,   helix_d2_tol    = 0.5f;
    float helix_d3_mid    = 5.3f,   helix_d3_tol    = 0.5f;
    float helix_d4_mid    = 6.4f,   helix_d4_tol    = 0.6f;

    float beta_tau_mid    = 124.0f, beta_tau_tol    = 14.0f;
    float beta_alpha_mid  = -170.0f, beta_alpha_tol = 45.0f;
    float beta_d2_mid     = 6.7f,   beta_d2_tol     = 0.6f;
    float beta_d3_mid     = 9.9f,   beta_d3_tol     = 0.9f;
    float beta_d4_mid     = 12.4f,  beta_d4_tol     = 1.1f;

    int   helix_min_len = 5;
    int   beta_min_len  = 4;
    int   short_beta_len = 3;

    float sheet_contact_min = 4.2f;
    float sheet_contact_max = 5.2f;
    int   sheet_contact_count = 5;

    float break_gap = 4.8f;

    void run(std::map<std::string, std::vector<Atom>>& atoms);

    void run_chain(std::vector<Atom>& chain_atoms);

    void set_scale(float scale_) { scale = scale_; }

private:
    struct Geom {
        bool  valid = false;
        float d2 = 0.0f, d3 = 0.0f, d4 = 0.0f;
        float tau = 0.0f, alpha = 0.0f;
    };

    static inline float dist(const Atom& a, const Atom& b) {
        float dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }

    static inline void cross(const float a[3], const float b[3], float out[3]) {
        out[0] = a[1]*b[2] - a[2]*b[1];
        out[1] = a[2]*b[0] - a[0]*b[2];
        out[2] = a[0]*b[1] - a[1]*b[0];
    }

    static inline float dot(const float a[3], const float b[3]) {
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    }

    static inline void sub(const Atom& p, const Atom& q, float out[3]) {
        out[0] = p.x - q.x; out[1] = p.y - q.y; out[2] = p.z - q.z;
    }

    static inline void normalize(float v[3]) {
        float n = std::sqrt(dot(v, v));
        if (n > 1e-8f) { v[0]/=n; v[1]/=n; v[2]/=n; }
    }

    static inline bool within(float value, float mid, float tol) {
        return value >= mid - tol && value <= mid + tol;
    }

    static bool within_angle(float value, float mid, float tol);

    static float angle_deg(const Atom& a, const Atom& b, const Atom& c);

    static float torsion_deg(const Atom& a, const Atom& b, const Atom& c, const Atom& d);

    std::vector<char> compute_breaks(const std::vector<Atom>& A);

    std::vector<Geom> compute_geometry(const std::vector<Atom>& A,
                                       const std::vector<char>& is_break);

    bool helix_site(const Geom& g) const;
    bool beta_site(const Geom& g) const;
    bool helix_extend_site(const Geom& g) const;
    bool beta_extend_site(const Geom& g) const;

    static std::vector<char> mask_consecutive(const std::vector<char>& in, int number);

    bool in_sheet(const std::vector<Atom>& A, size_t first, size_t last) const;

    void assign_helices(const std::vector<Geom>& geom,
                        const std::vector<char>& is_break,
                        std::vector<char>& lab);

    void assign_strands(const std::vector<Atom>& A,
                        const std::vector<Geom>& geom,
                        const std::vector<char>& is_break,
                        std::vector<char>& lab);
};
