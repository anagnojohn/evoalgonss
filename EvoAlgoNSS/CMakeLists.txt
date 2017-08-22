cmake_minimum_required(VERSION 3.8)
project(EvoAlgoNSS)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS -Wall)
set(SOURCE_FILES
        bond.h
        bondhelper.h
        differentialevo.h
        ealgorithm_base.h
        geneticalgo.h
        irr.h
        local_best_pso.h
        main.cpp
        svensson.h
        utilities.h
        yield_curve_fitting.h)
set(Boost_INCLUDE_DIR /home/test/Documents/boost_1_65_0)
set(Date_INCLUDE_DIR /home/test/Documents/date)
include_directories(${Boost_INCLUDE_DIR} ${Date_INCLUDE_DIR})
add_executable(EvoAlgoNSS ${SOURCE_FILES})
