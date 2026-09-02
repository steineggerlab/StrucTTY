#pragma once
#include <array>
#include <cstdint>
#include <string>

struct RGBA { uint8_t r,g,b,a; };
static_assert(sizeof(RGBA) == 4, "RGBA must be 4 bytes");

namespace Palettes {
    inline const std::array<int, 20> RAINBOW = {
        196, 202, 208, 214, 220, 226, 190, 154, 118, 82,
        49,  51,  45,  39,  33,  27,  21,  93,  129, 201
    };

    inline const std::array<int, 9> PROTEIN_COLORS = {
        100,
         80,
         27,
        129,
        213,
        209,
        130,
        214,
        160,
    };

    inline const std::array<int, 9> PROTEIN_DIM_COLORS = {
         58,
         30,
         18,
         54,
        132,
        131,
         94,
        172,
         88,
    };

    inline const std::array<int, 9> PROTEIN_BRIGHT_COLORS = {
        148,
        123,
         69,
        171,
        219,
        216,
        178,
        221,
        203,
    };

    inline constexpr int ALIGNED_NONALIGNED_DIM = 240;

    inline const std::array<int, 9> PROTEIN_NEAR_COLORS = {
        184,
        159,
        105,
        177,
        225,
        223,
        220,
        229,
        210,
    };

    inline const std::array<int, 9> PROTEIN_FAR_COLORS = {
         58,
         30,
         18,
         54,
         96,
        130,
         94,
        136,
         88,
    };

    inline const std::array<int, 15> CHAIN_NEAR_COLORS = {
        184,
        159,
        105,
        177,
        225,
        223,
        220,
        229,
        210,
         87,
        192,
        213,
        228,
        189,
        217,
    };

    inline const std::array<int, 15> CHAIN_FAR_COLORS = {
         58,
         30,
         18,
         54,
         96,
        130,
         94,
        136,
         88,
         23,
         64,
         90,
        136,
         60,
        131,
    };

    inline const std::array<int, 20> RAINBOW_NEAR = {
        210, 216, 222, 228, 229, 230, 192, 157, 157, 120,
         85,  87,  81,  75, 105,  63,  57,  99, 177, 213
    };

    inline const std::array<int, 20> RAINBOW_FAR = {
        124, 130, 136, 172, 178, 142, 106,  70,  34,  28,
         23,  24,  24,  18,  18,  18,  54,  54,  90, 127
    };

    inline const std::array<int, 4> PLDDT_NEAR = {69, 117, 229, 216};
    inline const std::array<int, 4> PLDDT_FAR = {24, 31, 58, 95};

    inline const std::array<int, 10> CONSERVATION_NEAR = {57, 63, 69, 75, 85, 157, 192, 229, 221, 210};
    inline const std::array<int, 10> CONSERVATION_FAR = {18, 18, 19, 24, 28, 30, 64, 58, 94, 88};

    inline constexpr int INTERFACE_NEAR_COLOR     = 213;
    inline constexpr int INTERFACE_DIM_NEAR_COLOR = 243;
    inline constexpr int INTERFACE_FAR_COLOR      =  90;
    inline constexpr int INTERFACE_DIM_FAR_COLOR  = 236;

    inline const std::array<int, 9> ALIGNED_BRIGHT_NEAR = {
        190,
        159,
        111,
        183,
        231,
        224,
        227,
        230,
        211,
    };
    inline constexpr int ALIGNED_NONALIGNED_FAR = 235;

    inline const std::array<int, 9> ALIGNED_BRIGHT_FAR = {
        106,
         73,
         26,
         98,
        176,
        173,
        137,
        172,
        167,
    };

    inline constexpr int ALIGNED_NONALIGNED_NEAR = 245;

    inline constexpr int SS_HELIX_NEAR = 228;
    inline constexpr int SS_HELIX_FAR  = 142;
    inline constexpr int SS_SHEET_NEAR =  87;
    inline constexpr int SS_SHEET_FAR  =  37;

    inline const std::array<int, 9> PROTEIN_COIL_FAR = {
        236,
         23,
         17,
         53,
         96,
         95,
         58,
        130,
         52,
    };

    inline constexpr int CONSERVATION_UNSCORED_NEAR = 249;
    inline constexpr int CONSERVATION_UNSCORED_MID  = 244;
    inline constexpr int CONSERVATION_UNSCORED_FAR  = 238;

    inline constexpr int INTERFACE_COLOR     = 201;
    inline constexpr int INTERFACE_DIM_COLOR = 58;

    inline constexpr int ALIGNED_COLOR     = 46;
    inline constexpr int ALIGNED_DIM_COLOR = 58;

    inline const std::array<int, 4> PLDDT_COLORS = {26, 81, 220, 209};

    inline const std::array<int, 10> CONSERVATION_COLORS = {21, 27, 33, 39, 49, 118, 190, 226, 214, 196};

    inline const std::array<int, 15> CHAIN_COLORS = {
        100,
         80,
         27,
        129,
        213,
        209,
        130,
        214,
        160,
         37,
        118,
        201,
        220,
        183,
        210,
    };

    inline constexpr RGBA ID2RGBA[256] = {
        RGBA{0, 0, 0, 255},
        RGBA{128, 0, 0, 255},
        RGBA{0, 128, 0, 255},
        RGBA{128, 128, 0, 255},
        RGBA{0, 0, 128, 255},
        RGBA{128, 0, 128, 255},
        RGBA{0, 128, 128, 255},
        RGBA{192, 192, 192, 255},
        RGBA{128, 128, 128, 255},
        RGBA{255, 0, 0, 255},
        RGBA{0, 255, 0, 255},
        RGBA{255, 255, 0, 255},
        RGBA{0, 0, 255, 255},
        RGBA{255, 0, 255, 255},
        RGBA{0, 255, 255, 255},
        RGBA{255, 255, 255, 255},
        RGBA{0, 0, 0, 255},
        RGBA{0, 0, 95, 255},
        RGBA{0, 0, 135, 255},
        RGBA{0, 0, 175, 255},
        RGBA{0, 0, 215, 255},
        RGBA{0, 0, 255, 255},
        RGBA{0, 95, 0, 255},
        RGBA{0, 95, 95, 255},
        RGBA{0, 95, 135, 255},
        RGBA{0, 95, 175, 255},
        RGBA{0, 95, 215, 255},
        RGBA{0, 95, 255, 255},
        RGBA{0, 135, 0, 255},
        RGBA{0, 135, 95, 255},
        RGBA{0, 135, 135, 255},
        RGBA{0, 135, 175, 255},
        RGBA{0, 135, 215, 255},
        RGBA{0, 135, 255, 255},
        RGBA{0, 175, 0, 255},
        RGBA{0, 175, 95, 255},
        RGBA{0, 175, 135, 255},
        RGBA{0, 175, 175, 255},
        RGBA{0, 175, 215, 255},
        RGBA{0, 175, 255, 255},
        RGBA{0, 215, 0, 255},
        RGBA{0, 215, 95, 255},
        RGBA{0, 215, 135, 255},
        RGBA{0, 215, 175, 255},
        RGBA{0, 215, 215, 255},
        RGBA{0, 215, 255, 255},
        RGBA{0, 255, 0, 255},
        RGBA{0, 255, 95, 255},
        RGBA{0, 255, 135, 255},
        RGBA{0, 255, 175, 255},
        RGBA{0, 255, 215, 255},
        RGBA{0, 255, 255, 255},
        RGBA{95, 0, 0, 255},
        RGBA{95, 0, 95, 255},
        RGBA{95, 0, 135, 255},
        RGBA{95, 0, 175, 255},
        RGBA{95, 0, 215, 255},
        RGBA{95, 0, 255, 255},
        RGBA{95, 95, 0, 255},
        RGBA{95, 95, 95, 255},
        RGBA{95, 95, 135, 255},
        RGBA{95, 95, 175, 255},
        RGBA{95, 95, 215, 255},
        RGBA{95, 95, 255, 255},
        RGBA{95, 135, 0, 255},
        RGBA{95, 135, 95, 255},
        RGBA{95, 135, 135, 255},
        RGBA{95, 135, 175, 255},
        RGBA{95, 135, 215, 255},
        RGBA{95, 135, 255, 255},
        RGBA{95, 175, 0, 255},
        RGBA{95, 175, 95, 255},
        RGBA{95, 175, 135, 255},
        RGBA{95, 175, 175, 255},
        RGBA{95, 175, 215, 255},
        RGBA{95, 175, 255, 255},
        RGBA{95, 215, 0, 255},
        RGBA{95, 215, 95, 255},
        RGBA{95, 215, 135, 255},
        RGBA{95, 215, 175, 255},
        RGBA{95, 215, 215, 255},
        RGBA{95, 215, 255, 255},
        RGBA{95, 255, 0, 255},
        RGBA{95, 255, 95, 255},
        RGBA{95, 255, 135, 255},
        RGBA{95, 255, 175, 255},
        RGBA{95, 255, 215, 255},
        RGBA{95, 255, 255, 255},
        RGBA{135, 0, 0, 255},
        RGBA{135, 0, 95, 255},
        RGBA{135, 0, 135, 255},
        RGBA{135, 0, 175, 255},
        RGBA{135, 0, 215, 255},
        RGBA{135, 0, 255, 255},
        RGBA{135, 95, 0, 255},
        RGBA{135, 95, 95, 255},
        RGBA{135, 95, 135, 255},
        RGBA{135, 95, 175, 255},
        RGBA{135, 95, 215, 255},
        RGBA{135, 95, 255, 255},
        RGBA{135, 135, 0, 255},
        RGBA{135, 135, 95, 255},
        RGBA{135, 135, 135, 255},
        RGBA{135, 135, 175, 255},
        RGBA{135, 135, 215, 255},
        RGBA{135, 135, 255, 255},
        RGBA{135, 175, 0, 255},
        RGBA{135, 175, 95, 255},
        RGBA{135, 175, 135, 255},
        RGBA{135, 175, 175, 255},
        RGBA{135, 175, 215, 255},
        RGBA{135, 175, 255, 255},
        RGBA{135, 215, 0, 255},
        RGBA{135, 215, 95, 255},
        RGBA{135, 215, 135, 255},
        RGBA{135, 215, 175, 255},
        RGBA{135, 215, 215, 255},
        RGBA{135, 215, 255, 255},
        RGBA{135, 255, 0, 255},
        RGBA{135, 255, 95, 255},
        RGBA{135, 255, 135, 255},
        RGBA{135, 255, 175, 255},
        RGBA{135, 255, 215, 255},
        RGBA{135, 255, 255, 255},
        RGBA{175, 0, 0, 255},
        RGBA{175, 0, 95, 255},
        RGBA{175, 0, 135, 255},
        RGBA{175, 0, 175, 255},
        RGBA{175, 0, 215, 255},
        RGBA{175, 0, 255, 255},
        RGBA{175, 95, 0, 255},
        RGBA{175, 95, 95, 255},
        RGBA{175, 95, 135, 255},
        RGBA{175, 95, 175, 255},
        RGBA{175, 95, 215, 255},
        RGBA{175, 95, 255, 255},
        RGBA{175, 135, 0, 255},
        RGBA{175, 135, 95, 255},
        RGBA{175, 135, 135, 255},
        RGBA{175, 135, 175, 255},
        RGBA{175, 135, 215, 255},
        RGBA{175, 135, 255, 255},
        RGBA{175, 175, 0, 255},
        RGBA{175, 175, 95, 255},
        RGBA{175, 175, 135, 255},
        RGBA{175, 175, 175, 255},
        RGBA{175, 175, 215, 255},
        RGBA{175, 175, 255, 255},
        RGBA{175, 215, 0, 255},
        RGBA{175, 215, 95, 255},
        RGBA{175, 215, 135, 255},
        RGBA{175, 215, 175, 255},
        RGBA{175, 215, 215, 255},
        RGBA{175, 215, 255, 255},
        RGBA{175, 255, 0, 255},
        RGBA{175, 255, 95, 255},
        RGBA{175, 255, 135, 255},
        RGBA{175, 255, 175, 255},
        RGBA{175, 255, 215, 255},
        RGBA{175, 255, 255, 255},
        RGBA{215, 0, 0, 255},
        RGBA{215, 0, 95, 255},
        RGBA{215, 0, 135, 255},
        RGBA{215, 0, 175, 255},
        RGBA{215, 0, 215, 255},
        RGBA{215, 0, 255, 255},
        RGBA{215, 95, 0, 255},
        RGBA{215, 95, 95, 255},
        RGBA{215, 95, 135, 255},
        RGBA{215, 95, 175, 255},
        RGBA{215, 95, 215, 255},
        RGBA{215, 95, 255, 255},
        RGBA{215, 135, 0, 255},
        RGBA{215, 135, 95, 255},
        RGBA{215, 135, 135, 255},
        RGBA{215, 135, 175, 255},
        RGBA{215, 135, 215, 255},
        RGBA{215, 135, 255, 255},
        RGBA{215, 175, 0, 255},
        RGBA{215, 175, 95, 255},
        RGBA{215, 175, 135, 255},
        RGBA{215, 175, 175, 255},
        RGBA{215, 175, 215, 255},
        RGBA{215, 175, 255, 255},
        RGBA{215, 215, 0, 255},
        RGBA{215, 215, 95, 255},
        RGBA{215, 215, 135, 255},
        RGBA{215, 215, 175, 255},
        RGBA{215, 215, 215, 255},
        RGBA{215, 215, 255, 255},
        RGBA{215, 255, 0, 255},
        RGBA{215, 255, 95, 255},
        RGBA{215, 255, 135, 255},
        RGBA{215, 255, 175, 255},
        RGBA{215, 255, 215, 255},
        RGBA{215, 255, 255, 255},
        RGBA{255, 0, 0, 255},
        RGBA{255, 0, 95, 255},
        RGBA{255, 0, 135, 255},
        RGBA{255, 0, 175, 255},
        RGBA{255, 0, 215, 255},
        RGBA{255, 0, 255, 255},
        RGBA{255, 95, 0, 255},
        RGBA{255, 95, 95, 255},
        RGBA{255, 95, 135, 255},
        RGBA{255, 95, 175, 255},
        RGBA{255, 95, 215, 255},
        RGBA{255, 95, 255, 255},
        RGBA{255, 135, 0, 255},
        RGBA{255, 135, 95, 255},
        RGBA{255, 135, 135, 255},
        RGBA{255, 135, 175, 255},
        RGBA{255, 135, 215, 255},
        RGBA{255, 135, 255, 255},
        RGBA{255, 175, 0, 255},
        RGBA{255, 175, 95, 255},
        RGBA{255, 175, 135, 255},
        RGBA{255, 175, 175, 255},
        RGBA{255, 175, 215, 255},
        RGBA{255, 175, 255, 255},
        RGBA{255, 215, 0, 255},
        RGBA{255, 215, 95, 255},
        RGBA{255, 215, 135, 255},
        RGBA{255, 215, 175, 255},
        RGBA{255, 215, 215, 255},
        RGBA{255, 215, 255, 255},
        RGBA{255, 255, 0, 255},
        RGBA{255, 255, 95, 255},
        RGBA{255, 255, 135, 255},
        RGBA{255, 255, 175, 255},
        RGBA{255, 255, 215, 255},
        RGBA{255, 255, 255, 255},
        RGBA{8, 8, 8, 255},
        RGBA{18, 18, 18, 255},
        RGBA{28, 28, 28, 255},
        RGBA{38, 38, 38, 255},
        RGBA{48, 48, 48, 255},
        RGBA{58, 58, 58, 255},
        RGBA{68, 68, 68, 255},
        RGBA{78, 78, 78, 255},
        RGBA{88, 88, 88, 255},
        RGBA{98, 98, 98, 255},
        RGBA{108, 108, 108, 255},
        RGBA{118, 118, 118, 255},
        RGBA{128, 128, 128, 255},
        RGBA{138, 138, 138, 255},
        RGBA{148, 148, 148, 255},
        RGBA{158, 158, 158, 255},
        RGBA{168, 168, 168, 255},
        RGBA{178, 178, 178, 255},
        RGBA{188, 188, 188, 255},
        RGBA{198, 198, 198, 255},
        RGBA{208, 208, 208, 255},
        RGBA{218, 218, 218, 255},
        RGBA{228, 228, 228, 255},
        RGBA{238, 238, 238, 255}
    };

    inline int palette_to_xterm256_fg(int cid) {
        if      (cid == 10)                return CONSERVATION_UNSCORED_NEAR;
        else if (cid == 20)                return CONSERVATION_UNSCORED_FAR;
        else if (cid == 36)                return CONSERVATION_UNSCORED_MID;
        else if (cid == 47)                return SS_HELIX_NEAR;
        else if (cid == 48)                return SS_HELIX_FAR;
        else if (cid == 49)                return SS_SHEET_NEAR;
        else if (cid == 50)                return SS_SHEET_FAR;
        else if (cid >= 85  && cid <= 93)  return PROTEIN_COIL_FAR[cid - 85];
        else if (cid == 111)               return ALIGNED_NONALIGNED_NEAR;
        else if (cid >= 251 && cid <= 259) return ALIGNED_BRIGHT_FAR[cid - 251];
        else if (cid >= 1   && cid <= 9)   return PROTEIN_COLORS[cid - 1];
        else if (cid >= 11  && cid <= 19)  return PROTEIN_DIM_COLORS[cid - 11];
        else if (cid >= 21  && cid <= 35)  return CHAIN_COLORS[cid - 21];
        else if (cid == 41)                return 226;
        else if (cid == 42)                return 51;
        else if (cid == 43)                return INTERFACE_COLOR;
        else if (cid == 44)                return INTERFACE_DIM_COLOR;
        else if (cid == 45)                return ALIGNED_COLOR;
        else if (cid == 46)                return ALIGNED_DIM_COLOR;
        else if (cid >= 51  && cid <= 70)  return RAINBOW[cid - 51];
        else if (cid >= 71  && cid <= 74)  return PLDDT_COLORS[cid - 71];
        else if (cid >= 75  && cid <= 84)  return CONSERVATION_COLORS[cid - 75];
        else if (cid >= 101 && cid <= 109) return PROTEIN_BRIGHT_COLORS[cid - 101];
        else if (cid == 110)               return ALIGNED_NONALIGNED_DIM;
        else if (cid >= 120 && cid <= 128) return PROTEIN_NEAR_COLORS[cid - 120];
        else if (cid >= 130 && cid <= 144) return CHAIN_NEAR_COLORS[cid - 130];
        else if (cid >= 145 && cid <= 159) return CHAIN_FAR_COLORS[cid - 145];
        else if (cid >= 160 && cid <= 179) return RAINBOW_NEAR[cid - 160];
        else if (cid >= 180 && cid <= 199) return RAINBOW_FAR[cid - 180];
        else if (cid >= 200 && cid <= 208) return PROTEIN_FAR_COLORS[cid - 200];
        else if (cid >= 209 && cid <= 212) return PLDDT_NEAR[cid - 209];
        else if (cid >= 213 && cid <= 216) return PLDDT_FAR[cid - 213];
        else if (cid >= 217 && cid <= 226) return CONSERVATION_NEAR[cid - 217];
        else if (cid >= 227 && cid <= 236) return CONSERVATION_FAR[cid - 227];
        else if (cid == 237)               return INTERFACE_NEAR_COLOR;
        else if (cid == 238)               return INTERFACE_DIM_NEAR_COLOR;
        else if (cid == 239)               return INTERFACE_FAR_COLOR;
        else if (cid == 240)               return INTERFACE_DIM_FAR_COLOR;
        else if (cid >= 241 && cid <= 249) return ALIGNED_BRIGHT_NEAR[cid - 241];
        else if (cid == 250)               return ALIGNED_NONALIGNED_FAR;
        else                               return 231;
    }

    inline std::string palette_to_ansi_fg_str(int color_id) {
        return "\033[38;5;" + std::to_string(palette_to_xterm256_fg(color_id)) + "m";
    }

    inline const char* palette_to_ansi_reset() {
        return "\033[0m";
    }
}
