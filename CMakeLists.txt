cmake_minimum_required (VERSION 3.1)

project(networkVMC CXX)

# set the default build type (RELEASE) if none is specified
set(default_build_type "RELEASE")
# check if the build type is empty, if yes, set it to release
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
       message(STATUS "No build type specified, defaulting to RELEASE")
       set(CMAKE_BUILD_TYPE "${default_build_type}")
endif()

# which compiler?
message(STATUS "Using ${CMAKE_CXX_COMPILER_ID} compiler")
# define the flags for release/debug/default
include(cmake/${CMAKE_CXX_COMPILER_ID}_flags.cmake)

#enable openmp
find_package(OpenMP)
if(OpenMP_FOUND)
	message(STATUS "Using OpenMP ${OpenMP_VERSION}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
	message(STATUS "WARNING: Did not find OpenMP, disabled shared memory paralleization")
endif()

add_subdirectory(src)

set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

#include eigen
include_directories("${PROJECT_SOURCE_DIR}/lib/eigen")

# get the boost installation
find_package(Boost 1.40 COMPONENTS regex REQUIRED)

#include boost
include_directories(${Boost_INCLUDE_DIR})

# define the library
add_library(networkVMC SHARED ${SOURCES})
# link with boost
target_link_libraries(networkVMC LINK_PUBLIC ${Boost_LIBRARIES})

# specify that we always link to networkVMC
#link_libraries(networkVMC)

add_subdirectory(test)
