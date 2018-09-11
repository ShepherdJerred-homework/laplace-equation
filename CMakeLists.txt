cmake_minimum_required(VERSION 3.12)
project(laplace-equation C)

set(CMAKE_C_STANDARD 90)

find_package(MPI REQUIRED)
include_directories(/usr/local/Cellar/open-mpi/3.1.1/include)

SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)

set(SOURCE_FILES main.c)
add_executable(laplace-equation src/main.c)