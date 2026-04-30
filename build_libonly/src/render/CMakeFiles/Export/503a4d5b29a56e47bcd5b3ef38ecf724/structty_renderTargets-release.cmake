#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "structty_render::structty_render" for configuration "Release"
set_property(TARGET structty_render::structty_render APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(structty_render::structty_render PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libstructty_render.a"
  )

list(APPEND _cmake_import_check_targets structty_render::structty_render )
list(APPEND _cmake_import_check_files_for_structty_render::structty_render "${_IMPORT_PREFIX}/lib/libstructty_render.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
