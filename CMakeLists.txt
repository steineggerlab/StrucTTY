cmake_minimum_required(VERSION 3.15)
project(StrucTTY LANGUAGES C CXX)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CURSES_NEED_WIDE TRUE)
find_package(Curses REQUIRED)

# Third-party deps are built in lib/
add_subdirectory(lib)

# App sources (under src/)
file(GLOB_RECURSE APP_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/structure/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/visualization/*.cpp"
)

file(GLOB_RECURSE APP_HEADERS CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Benchmark.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/utils/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/structure/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/visualization/*.hpp"
)

add_executable(StrucTTY
    ${CMAKE_CURRENT_SOURCE_DIR}/src/structty.cpp
    ${APP_SOURCES}
    ${APP_HEADERS}
)

target_include_directories(StrucTTY
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/src/utils
        ${CMAKE_CURRENT_SOURCE_DIR}/src/structure
        ${CMAKE_CURRENT_SOURCE_DIR}/src/visualization
)

target_link_libraries(StrucTTY
    PRIVATE
        gemmi::gemmi_cpp
        lodepng
        ${CURSES_LIBRARIES}
)
