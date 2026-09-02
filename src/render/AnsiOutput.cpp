#include "AnsiOutput.hpp"
#include "Palette.hpp"
#include <clocale>
#include <cstdio>
#include <limits>

static constexpr int DOT_BITS[2][4] = {
    {0, 1, 2, 6},
    {3, 4, 5, 7}
};

std::string AnsiOutput::to_ansi_string(
    const std::vector<RenderPoint>& pixels,
    int logical_width,
    int logical_height)
{
    const int term_w = logical_width  / 2;
    const int term_h = logical_height / 4;

    std::string out;
    out.reserve(static_cast<size_t>(term_w * term_h * 20 + term_h));

    for (int ty = 0; ty < term_h; ++ty) {
        for (int tx = 0; tx < term_w; ++tx) {
            int   bitmask       = 0;
            int   best_color_id = 0;
            float best_depth    = std::numeric_limits<float>::infinity();

            for (int sc = 0; sc < 2; ++sc) {
                for (int sr = 0; sr < 4; ++sr) {
                    int lx = tx * 2 + sc;
                    int ly = ty * 4 + sr;
                    if (lx >= logical_width || ly >= logical_height) continue;

                    const RenderPoint& lp = pixels[ly * logical_width + lx];
                    if (lp.color_id > 0) {
                        bitmask |= (1 << DOT_BITS[sc][sr]);
                        if (lp.depth < best_depth) {
                            best_depth    = lp.depth;
                            best_color_id = lp.color_id;
                        }
                    }
                }
            }

            if (bitmask > 0 && best_color_id > 0) {
                out += Palettes::palette_to_ansi_fg_str(best_color_id);
                out += static_cast<char>(0xE2);
                out += static_cast<char>(0xA0 | (bitmask >> 6));
                out += static_cast<char>(0x80 | (bitmask & 0x3F));
                out += Palettes::palette_to_ansi_reset();
            } else {
                out += ' ';
            }
        }
        out += '\n';
    }

    return out;
}

void AnsiOutput::print_to_stdout(
    const std::vector<RenderPoint>& pixels,
    int logical_width,
    int logical_height)
{
    setlocale(LC_ALL, "");
    const std::string s = to_ansi_string(pixels, logical_width, logical_height);
    std::fwrite(s.data(), 1, s.size(), stdout);
}
