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
        helper_functions.h
        irr.h
        local_best_pso.h
        main.cpp
        svensson.h
        yield_curve_fitting.h)
set(Boost_INCLUDE_DIR c:/mingw/include/)
set(Date_INCLUDE_DIR c:/users/john/development/repositories/remote/date)
set(Boost_LIBRARY_DIR c:/mingw/lib/)
find_package(Boost COMPONENTS system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIR} ${Date_INCLUDE_DIR})
add_executable(EvoAlgoNSS ${SOURCE_FILES})