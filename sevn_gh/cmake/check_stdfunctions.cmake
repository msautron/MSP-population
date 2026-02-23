include(CheckCXXCompilerFlag)

# Work around the fact that cmake does not propagate the language standard flag into
# the CHECK_CXX_SOURCE_COMPILES function. See CMake issue #16456.
# Ensure we do this after the FIND_PACKAGE calls which use C, and will error on a C++
# standards flag.
if(NOT POLICY CMP0067)
   list(APPEND CMAKE_REQUIRED_FLAGS "${CMAKE_CXX${CMAKE_CXX_STANDARD}_EXTENSION_COMPILE_OPTION}")
endif()


#Check if std::make_unique can be used
check_cxx_source_compiles("
#include <memory>
int main () {
 std::unique_ptr<int> foo = std::make_unique<int>();
}
"
HAVE_STD__MAKE_UNIQUE
)

message("CHECK compile ${HAVE_STD__MAKE_UNIQUE}")

if (HAVE_STD__MAKE_UNIQUE)
   add_definitions( -DHAVE_STD__MAKE_UNIQUE)
   message("STD library contains make_unique")
else()
   message("STD library does not contain make_unique")
endif()
