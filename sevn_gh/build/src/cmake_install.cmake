# Install script for directory: /home/matteo.sautron/sevn_gh/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/sevn" TYPE STATIC_LIBRARY FILES "/home/matteo.sautron/sevn_gh/build/src/libsevn_lib_static.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/sevn" TYPE FILE FILES
    "/home/matteo.sautron/sevn_gh/include/sevn.h"
    "/home/matteo.sautron/sevn_gh/include/sevn_version.h"
    "/home/matteo.sautron/sevn_gh/include/static_main.h"
    "/home/matteo.sautron/sevn_gh/src/binstar/BinaryProperty.h"
    "/home/matteo.sautron/sevn_gh/src/binstar/MTstability.h"
    "/home/matteo.sautron/sevn_gh/src/binstar/Orbit.h"
    "/home/matteo.sautron/sevn_gh/src/binstar/Processes.h"
    "/home/matteo.sautron/sevn_gh/src/binstar/binstar.h"
    "/home/matteo.sautron/sevn_gh/src/binstar/procs/Circularisation.h"
    "/home/matteo.sautron/sevn_gh/src/binstar/procs/Collision.h"
    "/home/matteo.sautron/sevn_gh/src/binstar/procs/Hardening.h"
    "/home/matteo.sautron/sevn_gh/src/binstar/procs/Windaccretion.h"
    "/home/matteo.sautron/sevn_gh/src/general/IO.h"
    "/home/matteo.sautron/sevn_gh/src/general/lookup_and_phases.h"
    "/home/matteo.sautron/sevn_gh/src/general/params.h"
    "/home/matteo.sautron/sevn_gh/src/general/starparameter.h"
    "/home/matteo.sautron/sevn_gh/src/general/types.h"
    "/home/matteo.sautron/sevn_gh/src/general/utils/errhand.h"
    "/home/matteo.sautron/sevn_gh/src/general/utils/evolve.h"
    "/home/matteo.sautron/sevn_gh/src/general/utils/sevnlog.h"
    "/home/matteo.sautron/sevn_gh/src/general/utils/utilities.h"
    "/home/matteo.sautron/sevn_gh/src/star/BSEintegrator.h"
    "/home/matteo.sautron/sevn_gh/src/star/lambdas/lambda_base.h"
    "/home/matteo.sautron/sevn_gh/src/star/lambdas/lambda_klencki21.h"
    "/home/matteo.sautron/sevn_gh/src/star/lambdas/lambda_nanjing.h"
    "/home/matteo.sautron/sevn_gh/src/star/procs/kicks/hobbs.h"
    "/home/matteo.sautron/sevn_gh/src/star/procs/kicks/kicks.h"
    "/home/matteo.sautron/sevn_gh/src/star/procs/kicks/specialkicks.h"
    "/home/matteo.sautron/sevn_gh/src/star/procs/kicks/unified.h"
    "/home/matteo.sautron/sevn_gh/src/star/procs/kicks/zeros.h"
    "/home/matteo.sautron/sevn_gh/src/star/procs/neutrinomassloss/neutrinomassloss.h"
    "/home/matteo.sautron/sevn_gh/src/star/procs/neutrinomassloss/nmllattimer89.h"
    "/home/matteo.sautron/sevn_gh/src/star/procs/pairinstability/pairinstability.h"
    "/home/matteo.sautron/sevn_gh/src/star/procs/pairinstability/pifarmer19.h"
    "/home/matteo.sautron/sevn_gh/src/star/procs/pairinstability/piiorio22.h"
    "/home/matteo.sautron/sevn_gh/src/star/procs/pairinstability/pimapelli20.h"
    "/home/matteo.sautron/sevn_gh/src/star/procs/supernova/compactness.h"
    "/home/matteo.sautron/sevn_gh/src/star/procs/supernova/deathmatrix.h"
    "/home/matteo.sautron/sevn_gh/src/star/procs/supernova/delayed.h"
    "/home/matteo.sautron/sevn_gh/src/star/procs/supernova/directcollapse.h"
    "/home/matteo.sautron/sevn_gh/src/star/procs/supernova/rapid.h"
    "/home/matteo.sautron/sevn_gh/src/star/procs/supernova/supernova.h"
    "/home/matteo.sautron/sevn_gh/src/star/property.h"
    "/home/matteo.sautron/sevn_gh/src/star/remnant.h"
    "/home/matteo.sautron/sevn_gh/src/star/star.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/sevn/libsevn_lib_shared.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/sevn/libsevn_lib_shared.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/sevn/libsevn_lib_shared.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/sevn" TYPE SHARED_LIBRARY FILES "/home/matteo.sautron/sevn_gh/build/src/libsevn_lib_shared.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/sevn/libsevn_lib_shared.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/sevn/libsevn_lib_shared.so")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/sevn/libsevn_lib_shared.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/sevn" TYPE FILE FILES
    "/home/matteo.sautron/sevn_gh/include/sevn.h"
    "/home/matteo.sautron/sevn_gh/include/sevn_version.h"
    "/home/matteo.sautron/sevn_gh/include/static_main.h"
    "/home/matteo.sautron/sevn_gh/src/binstar/BinaryProperty.h"
    "/home/matteo.sautron/sevn_gh/src/binstar/MTstability.h"
    "/home/matteo.sautron/sevn_gh/src/binstar/Orbit.h"
    "/home/matteo.sautron/sevn_gh/src/binstar/Processes.h"
    "/home/matteo.sautron/sevn_gh/src/binstar/binstar.h"
    "/home/matteo.sautron/sevn_gh/src/binstar/procs/Circularisation.h"
    "/home/matteo.sautron/sevn_gh/src/binstar/procs/Collision.h"
    "/home/matteo.sautron/sevn_gh/src/binstar/procs/Hardening.h"
    "/home/matteo.sautron/sevn_gh/src/binstar/procs/Windaccretion.h"
    "/home/matteo.sautron/sevn_gh/src/general/IO.h"
    "/home/matteo.sautron/sevn_gh/src/general/lookup_and_phases.h"
    "/home/matteo.sautron/sevn_gh/src/general/params.h"
    "/home/matteo.sautron/sevn_gh/src/general/starparameter.h"
    "/home/matteo.sautron/sevn_gh/src/general/types.h"
    "/home/matteo.sautron/sevn_gh/src/general/utils/errhand.h"
    "/home/matteo.sautron/sevn_gh/src/general/utils/evolve.h"
    "/home/matteo.sautron/sevn_gh/src/general/utils/sevnlog.h"
    "/home/matteo.sautron/sevn_gh/src/general/utils/utilities.h"
    "/home/matteo.sautron/sevn_gh/src/star/BSEintegrator.h"
    "/home/matteo.sautron/sevn_gh/src/star/lambdas/lambda_base.h"
    "/home/matteo.sautron/sevn_gh/src/star/lambdas/lambda_klencki21.h"
    "/home/matteo.sautron/sevn_gh/src/star/lambdas/lambda_nanjing.h"
    "/home/matteo.sautron/sevn_gh/src/star/procs/kicks/hobbs.h"
    "/home/matteo.sautron/sevn_gh/src/star/procs/kicks/kicks.h"
    "/home/matteo.sautron/sevn_gh/src/star/procs/kicks/specialkicks.h"
    "/home/matteo.sautron/sevn_gh/src/star/procs/kicks/unified.h"
    "/home/matteo.sautron/sevn_gh/src/star/procs/kicks/zeros.h"
    "/home/matteo.sautron/sevn_gh/src/star/procs/neutrinomassloss/neutrinomassloss.h"
    "/home/matteo.sautron/sevn_gh/src/star/procs/neutrinomassloss/nmllattimer89.h"
    "/home/matteo.sautron/sevn_gh/src/star/procs/pairinstability/pairinstability.h"
    "/home/matteo.sautron/sevn_gh/src/star/procs/pairinstability/pifarmer19.h"
    "/home/matteo.sautron/sevn_gh/src/star/procs/pairinstability/piiorio22.h"
    "/home/matteo.sautron/sevn_gh/src/star/procs/pairinstability/pimapelli20.h"
    "/home/matteo.sautron/sevn_gh/src/star/procs/supernova/compactness.h"
    "/home/matteo.sautron/sevn_gh/src/star/procs/supernova/deathmatrix.h"
    "/home/matteo.sautron/sevn_gh/src/star/procs/supernova/delayed.h"
    "/home/matteo.sautron/sevn_gh/src/star/procs/supernova/directcollapse.h"
    "/home/matteo.sautron/sevn_gh/src/star/procs/supernova/rapid.h"
    "/home/matteo.sautron/sevn_gh/src/star/procs/supernova/supernova.h"
    "/home/matteo.sautron/sevn_gh/src/star/property.h"
    "/home/matteo.sautron/sevn_gh/src/star/remnant.h"
    "/home/matteo.sautron/sevn_gh/src/star/star.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/sevn/SEVNlibs.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/sevn/SEVNlibs.cmake"
         "/home/matteo.sautron/sevn_gh/build/src/CMakeFiles/Export/lib/cmake/sevn/SEVNlibs.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/sevn/SEVNlibs-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/sevn/SEVNlibs.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/sevn" TYPE FILE FILES "/home/matteo.sautron/sevn_gh/build/src/CMakeFiles/Export/lib/cmake/sevn/SEVNlibs.cmake")
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^()$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/sevn" TYPE FILE FILES "/home/matteo.sautron/sevn_gh/build/src/CMakeFiles/Export/lib/cmake/sevn/SEVNlibs-noconfig.cmake")
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/sevn" TYPE FILE FILES
    "/home/matteo.sautron/sevn_gh/build/src/SEVNConfig.cmake"
    "/home/matteo.sautron/sevn_gh/build/src/SEVNConfigVersion.cmake"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE PROGRAM RENAME "sevnenv" FILES "/home/matteo.sautron/sevn_gh/build/src/sevnenv.sh")
endif()

