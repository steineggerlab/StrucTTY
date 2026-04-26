add_library(structty_render STATIC
    Renderer.cpp
    AnsiOutput.cpp
)

target_include_directories(structty_render
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${PROJECT_SOURCE_DIR}/src/visualization
)

# No ncurses linkage — intentional; consumers verify with:
#   nm -u libstructty_render.a | grep -i ncurses
target_link_libraries(structty_render
    PRIVATE
        gemmi::gemmi_cpp
)
