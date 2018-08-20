# here, we add all the files belonging to the networkVMC library to the sources variable


#define macro for adding source files
macro (add_sources)
    set(_relPath "${CMAKE_CURRENT_SOURCE_DIR}")
    # go through each file passed as argument
    foreach (_src ${ARGN})
    # append the file path relative to project root
        if (_relPath)
            list (APPEND SOURCES "${_relPath}/${_src}")
        else()
            list (APPEND SOURCES "${_src}")
        endif()
    endforeach()
    if (_relPath)
        # propagate SOURCES to parent directory
        set (SOURCES ${SOURCES} PARENT_SCOPE)
    endif()
endmacro()

add_subdirectory(CostFunctions)
add_subdirectory(Hamiltonian)
add_subdirectory(HilbertSpace)
add_subdirectory(math)
add_subdirectory(Network)
add_subdirectory(Samplers)
add_subdirectory(Solvers)
add_subdirectory(utilities)
add_sources(Trainer.cxx)
