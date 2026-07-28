cmake_minimum_required(VERSION 3.25...3.29)
# Fortran stays enabled: the executable embeds the fortran library, so the
# fortran compiler drives the link and supplies its runtime
project(demo_Wannier90 LANGUAGES CXX Fortran)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Wannier90 REQUIRED)

add_executable(demo wannier.cc)
target_link_libraries(demo PRIVATE Wannier90::wannier90)
if (Wannier90_MPI)
  # under MPI the demo initialises MPI and hands mpi_comm_world to the library
  find_package(MPI REQUIRED COMPONENTS CXX)
  target_link_libraries(demo PRIVATE MPI::MPI_CXX)
  target_compile_definitions(demo PRIVATE MPI)
endif ()
set_target_properties(demo PROPERTIES
        LINKER_LANGUAGE Fortran
        OUTPUT_NAME wannier_c
        SUFFIX .x
)
if (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
  # the main program is C++, not fortran: stop ifx supplying its own main
  target_link_options(demo PRIVATE -nofor-main)
endif ()

enable_testing()
# no ctest entry: the demo needs a .win file plus A/M/eigenvalue matrices as
# input, and is run by hand in the test-suite/tests directories (see readme.md)
