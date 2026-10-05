cmake_minimum_required(VERSION 3.25...3.29)
project(demo_Wannier90 LANGUAGES Fortran)

find_package(Wannier90 REQUIRED)
find_package(MPI REQUIRED COMPONENTS Fortran)

add_executable(demo demo.F90)
target_link_libraries(demo PRIVATE Wannier90::wannier90 MPI::MPI_Fortran)

enable_testing()
find_package(Python 3.8 REQUIRED)
