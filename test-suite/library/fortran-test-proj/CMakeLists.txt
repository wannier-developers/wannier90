cmake_minimum_required(VERSION 3.25...3.29)
project(w90_library_set_projections LANGUAGES Fortran)

find_package(Wannier90 REQUIRED)
find_package(MPI REQUIRED COMPONENTS Fortran)

add_executable(set_projections set_projections.f90)
target_link_libraries(set_projections PRIVATE Wannier90::wannier90 MPI::MPI_Fortran)

enable_testing()
add_test(NAME run_set_projections
         COMMAND $<TARGET_FILE:set_projections>
)
