#Check compiler flags
include(CheckCXXCompilerFlag)

set(MY_FLAGS " ")

check_cxx_compiler_flag(-std=c++14 HAVE_FLAG_STD_CXX14)
check_cxx_compiler_flag(-std=c++1y HAVE_FLAG_STD_CXX1Y)
if(HAVE_FLAG_STD_CXX14)
	set(MY_FLAGS "${MY_FLAGS} -std=c++14")
elseif(HAVE_FLAG_STD_CXX1Y)
	set(MY_FLAGS "${MY_FLAGS}-std=c++1y")
else()
	message(FATAL_ERROR "The compiler have to accept either the flag std=c++14 or std=c++1y")
endif()

check_cxx_compiler_flag(-fdiagnostics-color=auto HAVE_FDIAGNOSTIC_COLOR_FLAG)
if(HAVE_FLAG_STD_CXX14)
	set(MY_FLAGS "${MY_FLAGS} -fdiagnostics-color=auto")
endif()