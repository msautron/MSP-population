#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "SEVNLIBS::sevn_lib_static" for configuration ""
set_property(TARGET SEVNLIBS::sevn_lib_static APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(SEVNLIBS::sevn_lib_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "CXX"
  IMPORTED_LINK_INTERFACE_LIBRARIES_NOCONFIG "/usr/lib/x86_64-linux-gnu/libgsl.so;/usr/lib/x86_64-linux-gnu/libgslcblas.so"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/sevn/libsevn_lib_static.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS SEVNLIBS::sevn_lib_static )
list(APPEND _IMPORT_CHECK_FILES_FOR_SEVNLIBS::sevn_lib_static "${_IMPORT_PREFIX}/lib/sevn/libsevn_lib_static.a" )

# Import target "SEVNLIBS::sevn_lib_shared" for configuration ""
set_property(TARGET SEVNLIBS::sevn_lib_shared APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(SEVNLIBS::sevn_lib_shared PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/sevn/libsevn_lib_shared.so"
  IMPORTED_SONAME_NOCONFIG "libsevn_lib_shared.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS SEVNLIBS::sevn_lib_shared )
list(APPEND _IMPORT_CHECK_FILES_FOR_SEVNLIBS::sevn_lib_shared "${_IMPORT_PREFIX}/lib/sevn/libsevn_lib_shared.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
