cmake_minimum_required(VERSION 3.8)
project(EvoAlgoNSS)

set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /usr/local/lib)
set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} /usr/local/include)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_FLAGS "-Wall -std=c++17")
set(CMAKE_EXE_LINKER_FLAGS  "-stdlib=libc++ -lc++abi")
set(CMAKE_C_COMPILER=clang)
set(CMAKE_CXX_COMPILER=clang++)
set(SOURCE_FILES
        src/bond.h
        src/bondhelper.h
        src/differentialevo.h
        src/ealgorithm_base.h
        src/geneticalgo.h
        src/irr.h
        src/local_best_pso.h
        tests/main.cpp
        src/svensson.h
        src/utilities.h
        src/yield_curve_fitting.h)
set(Boost_INCLUDE_DIR /home/test/Documents/boost_1_65_0)
set(Date_INCLUDE_DIR /home/test/Documents/date)
include_directories(${Boost_INCLUDE_DIR} ${Date_INCLUDE_DIR})
add_executable(EvoAlgoNSS ${SOURCE_FILES})