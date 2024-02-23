cmake_minimum_required(VERSION 3.25...3.29)
project(demo_Wannier90 LANGUAGES Fortran)

find_package(Wannier90 REQUIRED)

add_executable(demo demo.F90)
if (Wannier90_MPI)
	target_compile_definitions(demo PRIVATE Wannier90_MPI)
endif ()
target_link_libraries(demo PRIVATE Wannier90::wannier90)
