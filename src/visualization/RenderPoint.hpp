#pragma once
#include <limits>
#include <string>

struct RenderPoint {
    int         x = 0;
    int         y = 0;
    float       depth = std::numeric_limits<float>::infinity();
    char        pixel = ' ';
    int         color_id = 0;
    int         chainID = -1;
    char        structure = 0;

    bool        is_interface = false;

    bool        is_aligned = false;

    float       bfactor = 0.0f;

    float       conservation_score = -1.0f;

    int         residue_number = -1;
    char        residue_name[4] = {};

    int         depth_band = 0;
};
