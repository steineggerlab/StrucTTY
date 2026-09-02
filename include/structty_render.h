#pragma once
#include <string>
#include <vector>
#include "Renderer.hpp"
#include "AnsiOutput.hpp"

namespace structty {

inline std::string render_to_ansi(
    const std::vector<RenderAtom>& atoms,
    int width,
    int height,
    const std::string& mode = "protein",
    bool show_structure = false
) {
    Renderer r(width, height, mode, show_structure);
    r.render(atoms);
    return AnsiOutput::to_ansi_string(
        r.get_pixels(), r.get_logical_width(), r.get_logical_height());
}

inline void render_to_stdout(
    const std::vector<RenderAtom>& atoms,
    int width,
    int height,
    const std::string& mode = "protein",
    bool show_structure = false
) {
    Renderer r(width, height, mode, show_structure);
    r.render(atoms);
    AnsiOutput::print_to_stdout(
        r.get_pixels(), r.get_logical_width(), r.get_logical_height());
}

}
