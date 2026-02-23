
####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was SEVNConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

#Include Directories
set_and_check(SEVN_INCLUDE_DIR "${PACKAGE_PREFIX_DIR}/include/sevn")
include_directories(${SEVN_INCLUDE_DIR})

#Get the library cmake
include("${CMAKE_CURRENT_LIST_DIR}/SEVNlibs.cmake")

#Include required flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -fopenmp")

set(SEVN_STATICLIB ON)
set(SEVN_SHAREDCLIB ON)

check_required_components(SEVN)
check_required_components(SEVNlibs)
