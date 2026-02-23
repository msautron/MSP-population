# first we can indicate the documentation build as an option and set it to OFF by default
option(doc "Build documentation" OFF)
if(doc)
    # check if Doxygen is installed
    find_package(Doxygen)
    if (DOXYGEN_FOUND)
        # set input and output files
        set(DOXYGEN_IN ${CMAKE_SOURCE_DIR}/doc/config_file)



        # note the option ALL which allows to build the docs together with the application
        add_custom_target( doc_doxygen ALL
                COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_IN}
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/doc/
                COMMENT "Generating API documentation with Doxygen"
                VERBATIM )
        message("-- Generate Documentation enabled (using Doxygen)")
    else (DOXYGEN_FOUND)
        message("-- Doxygen need to be installed to generate the doxygen documentation")
    endif (DOXYGEN_FOUND)

endif(doc)
