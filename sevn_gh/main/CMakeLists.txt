project(SEVN_app)

if (NOT staticlib AND NOT sharedlib)
    message( FATAL_ERROR "Cannot link the executable, neither staticlib or sharedlib are enabled" )
elseif (linklib STREQUAL "auto" AND staticlib)
    set(SEVNLIB sevn_lib_static)
elseif (linklib STREQUAL "auto" AND sharedlib)
    set(SEVNLIB sevn_lib_shared)
elseif(linklib STREQUAL "shared" AND sharedlib)
    set(SEVNLIB sevn_lib_shared)
elseif(linklib STREQUAL "shared" AND NOT sharedlib)
    message( FATAL_ERROR "You choose to make the exectuable linking the shared SEVN libraries, but these are not enabled" )
elseif(linklib STREQUAL "static" AND staticlib)
    set(SEVNLIB sevn_lib_static)
elseif(linklib STREQUAL "static" AND NOT staticlib)
    message( FATAL_ERROR "You choose to make the exectuable linking the static SEVN libraries, but these are not enabled" )
endif()


if (sin)
    #MAIN program SEVN single evolution
    add_executable(sevn.x sevn.cpp)

    #target_link_libraries(sevn.x ${GSL_LIBRARIES_GLOBAL} sevn_lib_static)

    if(h5)
        target_link_libraries(sevn.x ${SEVNLIB} ${HDF5_LIBRARIES})
	target_link_libraries(sevn.x ${SEVNLIB} ${GSL_LIBRARIES_GLOBAL})
    else()
        target_link_libraries(sevn.x ${SEVNLIB})
    endif(h5)

    add_custom_command(TARGET sevn.x POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E make_directory
            ${CMAKE_SOURCE_DIR}/build/exe/)

    add_custom_command(TARGET sevn.x POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:sevn.x>
            ${CMAKE_SOURCE_DIR}/build/exe/)
endif(sin)



#MAIN program SEVN binary evolution
if (bin)
    add_executable(sevnB.x sevnB.cpp)

    #target_link_libraries(sevnB.x ${SEVNLIB} ${GSL_LIBRARIES_GLOBAL})

    if(h5)
        target_link_libraries(sevnB.x ${SEVNLIB} ${HDF5_LIBRARIES})
	target_link_libraries(sevnB.x ${SEVNLIB} ${GSL_LIBRARIES_GLOBAL})
    else()
        target_link_libraries(sevnB.x ${SEVNLIB})
    endif(h5)

    add_custom_command(TARGET sevnB.x POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E make_directory
            ${CMAKE_SOURCE_DIR}/build/exe/)

    add_custom_command(TARGET sevnB.x POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:sevnB.x>
            ${CMAKE_SOURCE_DIR}/build/exe/)
endif(bin)


