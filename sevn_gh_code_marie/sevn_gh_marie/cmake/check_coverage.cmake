#Option to enable the code coverage through GNU gcov https://gcc.gnu.org/onlinedocs/gcc/Gcov.html
if (coverage)
    message("-- Coverage: enabled")
    SET(GCC_COVERAGE_LINK_FLAGS    "-coverage -lgcov -fprofile-arcs -ftest-coverage -pg")
    SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
else()
    message("-- Coverage: Disabled")
endif(coverage)