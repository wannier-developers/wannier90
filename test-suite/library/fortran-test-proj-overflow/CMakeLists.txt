cmake_minimum_required(VERSION 3.25...3.29)
project(library_mode_test_proj_overflow LANGUAGES Fortran)

find_package(Wannier90 REQUIRED)
find_package(MPI REQUIRED COMPONENTS Fortran)

add_executable(proj_overflow proj_overflow.f90)
target_link_libraries(proj_overflow PRIVATE Wannier90::wannier90 MPI::MPI_Fortran)

enable_testing()
add_test(NAME run_proj_overflow
        COMMAND $<TARGET_FILE:proj_overflow>
)
