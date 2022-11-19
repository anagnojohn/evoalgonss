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
        src/bond/bond.h
        src/bond/bondhelper.h
        src/solver/differentialevo.h
        src/solver/ealgorithm_base.h
        src/solver/geneticalgo.h
        src/model/irr.h
        src/solver/lbestpso.h
        src/solver/pso_sub_swarm.h
        tests/main.cpp
        src/model/svensson.h
        src/utilities.h
        src/model/yield_curve_fitting.h)
set(CMAKE_TOOLCHAIN_FILE "C:/vcpkg/scripts/buildsystems/vcpkg.cmake")
include_directories(${Boost_INCLUDE_DIR} ${Date_INCLUDE_DIR})
add_executable(EvoAlgoNSS ${SOURCE_FILES})