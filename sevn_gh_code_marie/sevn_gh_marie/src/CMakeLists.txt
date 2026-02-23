
#TODO at certain point give a look to cpack to distribute the libraries
#include(InstallRequiredSystemLibraries)
#set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
#set(CPACK_PACKAGE_VERSION_MAJOR "${Tutorial_VERSION_MAJOR}")
#set(CPACK_PACKAGE_VERSION_MINOR "${Tutorial_VERSION_MINOR}")
#set(CPACK_SOURCE_GENERATOR "TGZ")
#include(CPack)

file(GLOB SEVN_BASE "star/*.cpp" "star/lambdas/*.cpp" "star/procs/*.cpp" "star/procs/supernova/*.cpp" "star/procs/pairinstability/*.cpp" "star/procs/kicks/*.cpp" "star/procs/neutrinomassloss/*.cpp" "general/*.cpp" "general/utils/*.cpp" "binstar/*.cpp" "binstar/procs/*.cpp")


#if(h5)
#    file(GLOB H5OUT_BASE "h5out/*.cpp")
#    list(APPEND SEVN_BASE ${H5OUT_BASE})
#endif(h5)


message(${SEVN_BASE})

if(NOT staticlib AND NOT sharedlib)
    message( FATAL_ERROR "You disabled both the static and the shared library. At least one of them needs to be enabled" )
endif()


file(GLOB_RECURSE SEVN_HEADER "${PROJECT_SOURCE_DIR}/include/*h" "${PROJECT_SOURCE_DIR}/src/*h")

include(GNUInstallDirs)


if(staticlib)
    add_library(sevn_lib_static STATIC ${SEVN_BASE})
    set_target_properties(sevn_lib_static PROPERTIES PUBLIC_HEADER  "${SEVN_HEADER}")
    install(TARGETS sevn_lib_static
            EXPORT SEVNlibs
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/sevn
            PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/sevn")
    target_link_libraries(sevn_lib_static ${GSL_LIBRARIES_GLOBAL})
endif(staticlib)

if(sharedlib)
    # use, i.e. don't skip the full RPATH for the build tree
    if (APPLE)
          set(CMAKE_MACOSX_RPATH ON)
    endif (APPLE)


    add_library(sevn_lib_shared SHARED ${SEVN_BASE})
    set_target_properties(sevn_lib_shared PROPERTIES PUBLIC_HEADER  "${SEVN_HEADER}")
    install(TARGETS sevn_lib_shared
            EXPORT SEVNlibs
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/sevn
            PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/sevn)
endif(sharedlib)

if (linklib STREQUAL "auto" AND staticlib)
    set(SEVNLIBL sevn_lib_static)
elseif (linklib STREQUAL "auto" AND sharedlib)
    set(SEVNLIBL sevn_lib_shared)
elseif(linklib STREQUAL "shared" AND sharedlib)
    set(SEVNLIBL sevn_lib_shared)
elseif(linklib STREQUAL "shared" AND NOT sharedlib)
    message( FATAL_ERROR "You choose to make the exectuable linking the shared SEVN libraries, but these are not enabled" )
elseif(linklib STREQUAL "static" AND staticlib)
    set(SEVNLIBL sevn_lib_static)
elseif(linklib STREQUAL "static" AND NOT staticlib)
    message( FATAL_ERROR "You choose to make the exectuable linking the static SEVN libraries, but these are not enabled" )
endif()


#Cmake Package staff
set(INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_INCLUDEDIR}/sevn)
set(LIB_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR}/sevn)

#Libraries Cmake
install(
        EXPORT SEVNlibs
        FILE SEVNlibs.cmake
        NAMESPACE SEVNLIBS::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/sevn
)


#Package Cmake
include(CMakePackageConfigHelpers)
# generate the config file that is includes the exports
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/SEVNConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/SEVNConfig.cmake
        INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/sevn"
        PATH_VARS INCLUDE_INSTALL_DIR LIB_INSTALL_DIR)
write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/SEVNConfigVersion.cmake"
        VERSION "2.1.1"
        COMPATIBILITY AnyNewerVersion
)

#Install SEVN cmake to be load by other projects using Cmake
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/SEVNConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/SEVNConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/sevn)


#Alternative way to simplify the SEVN lib usage without Cmake
#Define a bashscript that sets the SEVN env variables
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/sevnenv.sh
        "#!/usr/bin/env bash
export SEVNFLAGS=\"-I${CMAKE_INSTALL_FULL_INCLUDEDIR}/sevn -L${CMAKE_INSTALL_FULL_LIBDIR}/sevn -l${SEVNLIBL} ${CMAKE_CXX_FLAGS}\"
export SEVNDIR=\"${PROJECT_SOURCE_DIR}\"
export SEVNCXX=\"${CMAKE_CXX_COMPILER}\"
alias sevncxx=\"${CMAKE_CXX_COMPILER}\"
        ")
#Move the bash script to the bin folder so that it can be directly sourced by the users
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/sevnenv.sh DESTINATION bin RENAME sevnenv)

