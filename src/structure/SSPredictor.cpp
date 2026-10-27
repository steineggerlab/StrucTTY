#include "SSPredictor.hpp"
#include <algorithm>
#include <iostream>

static const float PI_DEG = 180.0f / 3.14159265358979323846f;

bool SSPredictor::within_angle(float value, float mid, float tol) {
    float diff = value - mid;
    while (diff > 180.0f)  diff -= 360.0f;
    while (diff < -180.0f) diff += 360.0f;
    return std::fabs(diff) <= tol;
}

float SSPredictor::angle_deg(const Atom& a, const Atom& b, const Atom& c) {
    float u[3], v[3];
    sub(a, b, u);
    sub(c, b, v);
    normalize(u);
    normalize(v);
    float cosine = dot(u, v);
    if (cosine >  1.0f) cosine =  1.0f;
    if (cosine < -1.0f) cosine = -1.0f;
    return std::acos(cosine) * PI_DEG;
}

float SSPredictor::torsion_deg(const Atom& a, const Atom& b, const Atom& c, const Atom& d) {
    float b1[3], b2[3], b3[3];
    sub(b, a, b1); sub(c, b, b2); sub(d, c, b3);

    float n1[3], n2[3];
    cross(b1, b2, n1);
    cross(b2, b3, n2);

    float b2n[3] = { b2[0], b2[1], b2[2] };
    normalize(b2n);

    float m1[3];
    cross(n1, b2n, m1);

    float x = dot(n1, n2);
    float y = dot(m1, n2);

    return std::atan2(y, x) * PI_DEG;
}

std::vector<char> SSPredictor::compute_breaks(const std::vector<Atom>& A) {
    const size_t n = A.size();
    std::vector<char> br(n, 0);
    if (n < 2) return br;
    for (size_t i = 0; i + 1 < n; ++i) {
        if (dist(A[i], A[i+1]) * scale > break_gap) {
            br[i] = 1;
        }
    }
    return br;
}

std::vector<SSPredictor::Geom>
SSPredictor::compute_geometry(const std::vector<Atom>& A,
                              const std::vector<char>& is_break) {
    const size_t n = A.size();
    std::vector<Geom> geom(n);
    if (n < 5) return geom;

    auto spans = [&](size_t lo, size_t hi) {
        for (size_t k = lo; k < hi; ++k) if (is_break[k]) return false;
        return true;
    };

    for (size_t i = 1; i + 3 < n; ++i) {
        if (!spans(i - 1, i + 3)) continue;
        Geom& g = geom[i];
        g.valid = true;
        g.d2 = dist(A[i-1], A[i+1]) * scale;
        g.d3 = dist(A[i-1], A[i+2]) * scale;
        g.d4 = dist(A[i-1], A[i+3]) * scale;
        g.tau = angle_deg(A[i-1], A[i], A[i+1]);
        g.alpha = torsion_deg(A[i-1], A[i], A[i+1], A[i+2]);
    }
    return geom;
}

bool SSPredictor::helix_site(const Geom& g) const {
    if (!g.valid) return false;
    const bool by_distance = within(g.d3, helix_d3_mid, helix_d3_tol) &&
                             within(g.d4, helix_d4_mid, helix_d4_tol);
    const bool by_angle = within_angle(g.tau, helix_tau_mid, helix_tau_tol) &&
                          within_angle(g.alpha, helix_alpha_mid, helix_alpha_tol);
    return by_distance || by_angle;
}

bool SSPredictor::beta_site(const Geom& g) const {
    if (!g.valid) return false;
    const bool by_distance = within(g.d2, beta_d2_mid, beta_d2_tol) &&
                             within(g.d3, beta_d3_mid, beta_d3_tol) &&
                             within(g.d4, beta_d4_mid, beta_d4_tol);
    const bool by_angle = within_angle(g.tau, beta_tau_mid, beta_tau_tol) &&
                          within_angle(g.alpha, beta_alpha_mid, beta_alpha_tol);
    return by_distance || by_angle;
}

bool SSPredictor::helix_extend_site(const Geom& g) const {
    if (!g.valid) return false;
    return within(g.d3, helix_d3_mid, helix_d3_tol) ||
           within_angle(g.tau, helix_tau_mid, helix_tau_tol);
}

bool SSPredictor::beta_extend_site(const Geom& g) const {
    if (!g.valid) return false;
    return within(g.d3, beta_d3_mid, beta_d3_tol);
}

std::vector<char> SSPredictor::mask_consecutive(const std::vector<char>& in, int number) {
    const size_t n = in.size();
    std::vector<char> out(n, 0);
    if (number <= 0 || n < (size_t)number) return out;
    for (size_t i = 0; i + (size_t)number <= n; ++i) {
        bool all = true;
        for (int k = 0; k < number; ++k) if (!in[i + k]) { all = false; break; }
        if (all) for (int k = 0; k < number; ++k) out[i + k] = 1;
    }
    return out;
}

bool SSPredictor::in_sheet(const std::vector<Atom>& A, size_t first, size_t last) const {
    const size_t n = A.size();
    int contacts = 0;
    for (size_t i = first; i <= last; ++i) {
        for (size_t j = 0; j < n; ++j) {
            const float d = dist(A[i], A[j]) * scale;
            if (d > sheet_contact_min && d <= sheet_contact_max) ++contacts;
        }
    }
    return contacts >= sheet_contact_count;
}

void SSPredictor::assign_helices(const std::vector<Geom>& geom,
                                 const std::vector<char>& is_break,
                                 std::vector<char>& lab) {
    const size_t n = lab.size();
    std::vector<char> strict(n, 0);
    for (size_t i = 0; i < n; ++i) strict[i] = helix_site(geom[i]) ? 1 : 0;

    std::vector<char> mask = mask_consecutive(strict, helix_min_len);

    size_t i = 0;
    while (i < n) {
        if (!mask[i]) { ++i; continue; }
        size_t j = i;
        while (j + 1 < n && mask[j + 1]) ++j;
        if (i > 0 && !is_break[i - 1] && helix_extend_site(geom[i - 1])) mask[i - 1] = 1;
        if (j + 1 < n && !is_break[j] && helix_extend_site(geom[j + 1])) mask[j + 1] = 1;
        i = j + 2;
    }

    for (size_t k = 0; k < n; ++k) if (mask[k]) lab[k] = 'H';
}

void SSPredictor::assign_strands(const std::vector<Atom>& A,
                                 const std::vector<Geom>& geom,
                                 const std::vector<char>& is_break,
                                 std::vector<char>& lab) {
    const size_t n = lab.size();
    std::vector<char> strict(n, 0);
    for (size_t i = 0; i < n; ++i) strict[i] = beta_site(geom[i]) ? 1 : 0;

    std::vector<char> mask = mask_consecutive(strict, beta_min_len);
    std::vector<char> cand = mask_consecutive(strict, short_beta_len);

    size_t i = 0;
    while (i < n) {
        if (!cand[i]) { ++i; continue; }
        size_t j = i;
        while (j + 1 < n && cand[j + 1]) ++j;
        if (in_sheet(A, i, j)) {
            for (size_t k = i; k <= j; ++k) mask[k] = 1;
        }
        i = j + 1;
    }

    i = 0;
    while (i < n) {
        if (!mask[i]) { ++i; continue; }
        size_t j = i;
        while (j + 1 < n && mask[j + 1]) ++j;
        if (i > 0 && !is_break[i - 1] && beta_extend_site(geom[i - 1])) mask[i - 1] = 1;
        if (j + 1 < n && !is_break[j] && beta_extend_site(geom[j + 1])) mask[j + 1] = 1;
        i = j + 2;
    }

    for (size_t k = 0; k < n; ++k) {
        if (mask[k] && lab[k] != 'H') lab[k] = 'S';
    }
}

void SSPredictor::run_chain(std::vector<Atom>& chain_atoms) {
    const size_t n = chain_atoms.size();
    if (n == 0) return;

    std::vector<char> is_break = compute_breaks(chain_atoms);
    std::vector<Geom> geom = compute_geometry(chain_atoms, is_break);
    std::vector<char> lab(n, 'x');

    assign_helices(geom, is_break, lab);
    assign_strands(chain_atoms, geom, is_break, lab);

    for (size_t i = 0; i < n; ++i) {
        chain_atoms[i].set_structure(lab[i]);
    }
}

void SSPredictor::run(std::map<std::string, std::vector<Atom>>& atoms) {
    std::cout << "  predict secondary structure\n";
    for (auto& chain : atoms) {
        run_chain(chain.second);
    }
}
