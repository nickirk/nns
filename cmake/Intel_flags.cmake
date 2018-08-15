# compiler flags for the intel compiler

set(CMAKE_CXX_FLAGS "-std=c++11 -g -Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-std=c++11 -g -debug all -O0 -Wall -traceback -check-uninit")
set(CMAKE_CXX_FLAGS_RELEASE "-std=c++11 -O3 -no-prec-div")
