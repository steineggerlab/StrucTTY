add_library(structty_render STATIC
    Renderer.cpp
    AnsiOutput.cpp
)

target_include_directories(structty_render
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/visualization>
        $<INSTALL_INTERFACE:include>
)

include(GNUInstallDirs)

install(TARGETS structty_render
    EXPORT  structty_renderTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(FILES
    ${PROJECT_SOURCE_DIR}/include/structty_render.h
    ${CMAKE_CURRENT_SOURCE_DIR}/Renderer.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/AnsiOutput.hpp
    ${PROJECT_SOURCE_DIR}/src/visualization/RenderPoint.hpp
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT structty_renderTargets
    FILE      structty_renderTargets.cmake
    NAMESPACE structty_render::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/structty_render
)

install(FILES
    ${PROJECT_SOURCE_DIR}/cmake/structty_renderConfig.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/structty_render
)
