#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "sleef::sleef" for configuration "Release"
set_property(TARGET sleef::sleef APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(sleef::sleef PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C;CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libsleef.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS sleef::sleef )
list(APPEND _IMPORT_CHECK_FILES_FOR_sleef::sleef "${_IMPORT_PREFIX}/lib/libsleef.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
