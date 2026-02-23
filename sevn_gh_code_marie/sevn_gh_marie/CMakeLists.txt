cmake_minimum_required(VERSION 2.8)

###VERSIONING
set(SEVN_DEVELOPMENT_LINE "zelaous_redgiant")
set(SEVN_MAJOR_VERSION "2")
set(SEVN_MINOR_VERSION "10")
set(SEVN_PATCH_VERSION "0")
set(SEVN_VERSION_STRING "${SEVN_MAJOR_VERSION}.${SEVN_MINOR_VERSION}.${SEVN_PATCH_VERSION}")

if (CMAKE_VERSION VERSION_LESS 3.0)
    project(SEVN CXX)
    set(PROJECT_VERSION_MAJOR "${SEVN_MAJOR_VERSION}")
    set(PROJECT_VERSION_MINOR "${SEVN_MINOR_VERSION}")
    set(PROJECT_VERSION_PATCH "${SEVN_PATCH_VERSION}")
    set(PROJECT_VERSION "${SEVN_VERSION_STRING}")
else()
    cmake_policy(SET CMP0048 NEW)
    PROJECT(SEVN VERSION "${SEVN_VERSION_STRING}" LANGUAGES CXX)
endif()

find_package(Git)
find_file(GITDIR NAMES .git PATHS ${CMAKE_SOURCE_DIR} NO_DEFAULT_PATH)

# Trouver la bibliothèque GSL
find_package(GSL REQUIRED)
message(STATUS "GSL Libraries: ${GSL_LIBRARIES}")

if(Git_FOUND AND GITDIR)
    message("Using a git cloned repository")

    #Magic to make the dependence to the git head so that is is updated when we push a commit
    #It used a fake targed that depends on the git HEAD that is updated each commit
    add_custom_target(
            fake_target
            DEPENDS ${GITDIR}/logs/HEAD
    )

    #Find branch name
    execute_process(
            COMMAND  git rev-parse --abbrev-ref HEAD
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            OUTPUT_VARIABLE  GIT_BRANCH
            OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    ##LOCAL###
    #Find commit hash
    execute_process(
            COMMAND  git rev-parse HEAD
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            OUTPUT_VARIABLE  GIT_SHA_LOCAL
            OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    #Find commit date
    execute_process(
            COMMAND  git log -1 --format=%ci
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            OUTPUT_VARIABLE  GIT_SHATIME_LOCAL
            OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    #Find commit timestamp
    execute_process(
            COMMAND  git log -1 --format=%ct
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            OUTPUT_VARIABLE  GIT_SHATIMESTAMP_LOCAL
            OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    #Find commit counter
    execute_process(
            COMMAND  git rev-list HEAD --count
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            OUTPUT_VARIABLE  GIT_COUNTER_LOCAL
            OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    ##REMOTE###
    #Find commit hash
    execute_process(
            COMMAND  git rev-parse "origin/${GIT_BRANCH}"
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            OUTPUT_VARIABLE  GIT_SHA_REMOTE
            OUTPUT_STRIP_TRAILING_WHITESPACE
            RESULT_VARIABLE STATUS
    )
    if(STATUS  GREATER 0)
        set(GIT_SHA_REMOTE "NOT AVAILABLE")
    endif()
    #Find commit date
    execute_process(
            COMMAND  git log "origin/${GIT_BRANCH}" -1 --format=%ci
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            OUTPUT_VARIABLE  GIT_SHATIME_REMOTE
            OUTPUT_STRIP_TRAILING_WHITESPACE
            RESULT_VARIABLE STATUS
    )
    if(STATUS  GREATER 0)
        set(GIT_SHATIME_REMOTE "NOT AVAILABLE")
    endif()
    #Find commit timestamp
    execute_process(
            COMMAND  git log "origin/${GIT_BRANCH}" -1 --format=%ct
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            OUTPUT_VARIABLE  GIT_SHATIMESTAMP_REMOTE
            OUTPUT_STRIP_TRAILING_WHITESPACE
            RESULT_VARIABLE STATUS
    )
    if(STATUS  GREATER 0)
        set(GIT_SHATIMESTAMP_REMOTE "0")
    endif()
    #Find commit counter
    execute_process(
            COMMAND  git rev-list "origin/${GIT_BRANCH}" --count
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            OUTPUT_VARIABLE  GIT_COUNTER_REMOTE
            OUTPUT_STRIP_TRAILING_WHITESPACE
            RESULT_VARIABLE STATUS
    )
    if(STATUS  GREATER 0)
        set(GIT_COUNTER_REMOTE "0")
    endif()

    message("Current branch name ${GIT_BRANCH}")
    message("Current Commit hascode ${GIT_SHA}")
    message("Date of last commit  ${GIT_SHATIME} (Timestamp: ${GIT_SHATIMESTAMP})")
else()
    message("Not Using a git cloned repository")
    set(GIT_BRANCH "NOT AVAILABLE")
    set(GIT_SHA_LOCAL "NOT AVAILABLE")
    set(GIT_SHATIME_LOCAL "NOT AVAILABLE")
    set(GIT_SHATIMESTAMP_LOCAL "0")
    set(GIT_COUNTER_LOCAL "0")
    set(GIT_SHA_REMOTE "NOT AVAILABLE")
    set(GIT_BRANCH_REMOTE "NOT AVAILABLE")
    set(GIT_SHATIME_REMOTE "NOT AVAILABLE")
    set(GIT_SHATIMESTAMP_REMOTE "0")
    set(GIT_COUNTER_REMOTE "0")
endif()
#@TODO there is a missing case, i.e. when it is a git repo but It has not a remote (i.e. a new local branch)
#in this case we will have problems in compiling because the sevn.h variable will be not set

configure_file("${PROJECT_SOURCE_DIR}/include/sevn_version.h.in" "${PROJECT_SOURCE_DIR}/include/sevn_version.h")



set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)



#check flags
set(REQUIRED_CXX_FLAGS "")

include(cmake/check_flags.cmake)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MY_FLAGS} -W -Wall  -O3")


if(APPLE)
    set (CMAKE_EXE_LINKER_FLAGS "-undefined dynamic_lookup")
    set (CMAKE_SHARED_LINKER_FLAGS "-undefined dynamic_lookup")
endif(APPLE)


#OPTIONS AND SETS
option(h5 "Enables HDF5 output files" OFF)
option(doc "Build documentation" OFF)
option(debug "Enable debug output" OFF)
option(coverage "Enable option for GNU GCOV utility" OFF)
option(test "Enable units test" OFF)
option(openmp "Enable OpenMP if available" ON)
option(staticlib "Compile the static library" ON)
option(sharedlib "Compile the shared library" ON)
set(linklib "auto"  CACHE STRING "What library use for the exectuable")
option(exe "Compile the executables" ON)
option(bin "Compile sevnB" ON)
option(sin "Compile sevn" ON)
option(verbose "Enable verbose compilation" OFF)

if(verbose)
    set(CMAKE_VERBOSE_MAKEFILE ON)
    message("-- Verbose compilation enabled")
endif(verbose)


#INCLUDE CMAKE FILES
include(cmake/check_opnemp.cmake)
include(cmake/check_debug.cmake)
include(cmake/check_coverage.cmake)
include(cmake/generate_doc.cmake)
#include(cmake/check_stdfunctions.cmake)

#include header files
include_directories(include src)
include_directories(${GSL_INCLUDE_DIRS})
include_directories(src/binstar src/binstar/procs)
include_directories(src/general src/general/utils)
include_directories(src/star src/star/lambdas src/star/procs src/star/procs/supernova src/star/procs/pairinstability  src/star/procs/kicks src/star/procs/neutrinomassloss)

set(GSL_LIBRARIES_GLOBAL ${GSL_LIBRARIES} CACHE INTERNAL "Global GSL libraries")

#Subdirectory
add_subdirectory(src)
if(exe)
    add_subdirectory(main)
endif(exe)

if(${CMAKE_VERSION} VERSION_LESS "3.8.0")
    message("Please consider to switch to CMake 3.8.0")
endif()



if(test AND ${CMAKE_VERSION} VERSION_LESS "3.2.0")
    message(WARNING "You enable the test compilation, but this is possible only with a Cmake version >=3.2, current version is ${CMAKE_VERSION}. The unit test will be disabled.")
    message("-- UNITS TEST compilation disabled")
elseif(test)
    message("-- Units test compilation enabled")
    add_subdirectory(extern/catch)
    add_subdirectory(test)
else()
    message("-- UNITS TEST compilation disabled")
endif()
