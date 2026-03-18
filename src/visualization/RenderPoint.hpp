#pragma once
#include <limits>
#include <string>

struct RenderPoint {
    int         x = 0;
    int         y = 0;
    float       depth = std::numeric_limits<float>::infinity();
    char        pixel = ' ';
    int         color_id = 0;
    std::string chainID;
    char        structure = 0;

    // 기능 6: 잔기 정보
    int         residue_number = -1;
    std::string residue_name;

    // 기능 2: pLDDT
    float       bfactor = 0.0f;

    // 기능 5: conservation
    float       conservation_score = -1.0f;
};
