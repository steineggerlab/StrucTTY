cmake_minimum_required(VERSION 3.15)
project(StrucTTY LANGUAGES C CXX)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set OFF when embedding via add_subdirectory without ncurses (e.g. Foldseek).
option(STRUCTTY_BUILD_APP  "Build StrucTTY terminal application" ON)
# Set ON to also build the render_test integration test binary.
option(BUILD_RENDER_TEST   "Build render_test integration test binary" OFF)

find_package(ZLIB REQUIRED)

# Third-party deps are built in lib/
add_subdirectory(lib)

# Render library (ncurses-free static lib; Renderer.cpp + AnsiOutput.cpp)
add_subdirectory(src/render)

if (STRUCTTY_BUILD_APP)
    # App sources (under src/)
    # src/render/*.cpp is compiled into structty_render — excluded here to avoid duplicates.
    file(GLOB_RECURSE APP_SOURCES CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/src/structure/*.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/visualization/*.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/utils/*.cpp"
    )

    file(GLOB_RECURSE APP_HEADERS CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/src/utils/*.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/structure/*.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/visualization/*.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/render/*.hpp"
    )

    add_executable(StrucTTY
        ${CMAKE_CURRENT_SOURCE_DIR}/src/structty.cpp
        ${APP_SOURCES}
        ${APP_HEADERS}
    )

    target_include_directories(StrucTTY
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/src/utils
            ${CMAKE_CURRENT_SOURCE_DIR}/src/structure
            ${CMAKE_CURRENT_SOURCE_DIR}/src/visualization
            ${CMAKE_CURRENT_SOURCE_DIR}/src/render
    )

    target_link_libraries(StrucTTY
        PRIVATE
            structty_render
            gemmi::gemmi_cpp
            lodepng
            ZLIB::ZLIB
    )
endif()

if (BUILD_RENDER_TEST)
    add_executable(render_test example/render_test.cpp)
    target_link_libraries(render_test PRIVATE structty_render gemmi::gemmi_cpp)
endif()
