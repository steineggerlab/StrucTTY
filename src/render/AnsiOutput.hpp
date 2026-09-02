#pragma once
#include <string>
#include <vector>
#include "RenderPoint.hpp"

class AnsiOutput {
public:
    static std::string to_ansi_string(
        const std::vector<RenderPoint>& pixels,
        int logical_width,
        int logical_height
    );

    static void print_to_stdout(
        const std::vector<RenderPoint>& pixels,
        int logical_width,
        int logical_height
    );
};
