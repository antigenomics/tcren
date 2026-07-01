cmake_minimum_required(VERSION 3.18)
project(tcren LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# pybind11's modern FindPython path (silences CMP0148 deprecation warnings).
set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

# Fitting-alignment discrimination for MHC pseudosequence (MPS) matching. See src/_align/align.cpp.
pybind11_add_module(_align src/_align/align.cpp)
target_compile_options(_align PRIVATE -O3)

# Potential-guided rigid-body peptide refinement (statistical potential + clash MC). See src/_refine/refine.cpp.
pybind11_add_module(_refine src/_refine/refine.cpp)
target_compile_options(_refine PRIVATE -O3)

# CCD loop-closure kernel for open-source anchor-restrained peptide modelling. See src/_fold/fold.cpp.
pybind11_add_module(_fold src/_fold/fold.cpp)
target_compile_options(_fold PRIVATE -O3)

install(TARGETS _align _refine _fold DESTINATION tcren)
