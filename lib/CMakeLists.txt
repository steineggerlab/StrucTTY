# --------------------------------------------------
# Third-party (submodules)
# --------------------------------------------------

# gemmi: provides CMake targets (e.g., gemmi::gemmi_cpp)
add_subdirectory(gemmi)

# lodepng: simple single-translation-unit library
add_library(lodepng STATIC
    ${CMAKE_CURRENT_SOURCE_DIR}/lodepng/lodepng.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/lodepng/lodepng.h
)

set_target_properties(lodepng PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(lodepng
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/lodepng
)
