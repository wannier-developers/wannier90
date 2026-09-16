cmake_minimum_required(VERSION 3.25...3.29)
project(demo_Wannier90 LANGUAGES Fortran)

find_package(Wannier90 REQUIRED)
find_package(MPI REQUIRED COMPONENTS Fortran)

add_executable(demo demo.F90)
target_link_libraries(demo PRIVATE Wannier90::wannier90 MPI::MPI_Fortran)

enable_testing()
find_package(Python 3.8 REQUIRED)
file(COPY ref/results_ref.dat
        DESTINATION ref/results_ref.dat
)
# TODO: Re-enable these tests
#add_test(NAME run_demo
#        COMMAND $<TARGET_FILE:demo> > results.dat
#)
#add_test(NAME compare_results
#        COMMAND Python::Interpreter ${CMAKE_CURRENT_SOURCE_DIR}/compare_results.py
#)
