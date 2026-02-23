cmake_minimum_required(VERSION 2.8)
project(Binevolv)

find_package(SEVN)

if(SEVN_FOUND)
    message("SEVN found!")
else()
    message( FATAL_ERROR "SEVN not found" )
endif()

add_executable(binary_evol.exe binary_evol.cpp)

#Check if static lib are present, otherwise link to shared
if(SEVN_STATICLIB)
    message("Executable will be linked to the SEVN static lib")
    target_link_libraries(binary_evol.exe PRIVATE SEVNLIBS::sevn_lib_static)
elseif(SEVN_SHAREDLIB)
    message("Executable will be linked to the SEVN shared lib")
    target_link_libraries(binary_evol.exe PRIVATE SEVNLIBS::sevn_lib_shared)
else()
    message( FATAL_ERROR "Neither the static nor the shared SEVN libs are available" )
endif()