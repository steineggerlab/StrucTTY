
if(NOT TARGET gemmi_cpp)
    add_library(gemmi_cpp INTERFACE)
    target_include_directories(gemmi_cpp SYSTEM INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/gemmi/include)
    target_compile_features(gemmi_cpp INTERFACE cxx_std_17)
    target_compile_options(gemmi_cpp INTERFACE -include cstdint)
    if(ZLIB_FOUND)
        target_link_libraries(gemmi_cpp INTERFACE ZLIB::ZLIB)
    else()
        target_include_directories(gemmi_cpp SYSTEM INTERFACE
            ${CMAKE_CURRENT_SOURCE_DIR}/gemmi/third_party/zlib)
    endif()
    add_library(gemmi::gemmi_cpp ALIAS gemmi_cpp)
endif()

add_library(stb INTERFACE)

target_include_directories(stb
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/stb
)
