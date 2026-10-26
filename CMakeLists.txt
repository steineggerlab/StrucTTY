cmake_minimum_required(VERSION 3.15)
project(StrucTTY LANGUAGES C CXX)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(STRUCTTY_BUILD_APP  "Build StrucTTY terminal application" ON)
option(BUILD_RENDER_TEST   "Build render_test integration test binary" OFF)

find_package(ZLIB)

add_subdirectory(lib)

add_subdirectory(src/render)

file(GLOB_RECURSE STRUCTTY_LIB_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/structty.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/structure/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/visualization/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/utils/*.cpp"
)

file(GLOB_RECURSE STRUCTTY_LIB_HEADERS CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/include/structty.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/utils/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/structure/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/visualization/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/render/*.hpp"
)

add_library(structty STATIC
    ${STRUCTTY_LIB_SOURCES}
    ${STRUCTTY_LIB_HEADERS}
)

target_include_directories(structty
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/src/utils
        ${CMAKE_CURRENT_SOURCE_DIR}/src/structure
        ${CMAKE_CURRENT_SOURCE_DIR}/src/visualization
        ${CMAKE_CURRENT_SOURCE_DIR}/src/render
)

target_link_libraries(structty
    PUBLIC
        structty_render
        gemmi::gemmi_cpp
        stb
)
if(ZLIB_FOUND)
    target_link_libraries(structty PUBLIC ZLIB::ZLIB)
endif()

if (STRUCTTY_BUILD_APP)
    add_executable(structty_app ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)
    set_target_properties(structty_app PROPERTIES OUTPUT_NAME StrucTTY)
    target_link_libraries(structty_app PRIVATE structty)
endif()

if (BUILD_RENDER_TEST)
    add_executable(render_test example/render_test.cpp)
    target_link_libraries(render_test PRIVATE structty_render gemmi::gemmi_cpp)
endif()
