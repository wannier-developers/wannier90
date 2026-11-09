cmake_minimum_required(VERSION 3.25...3.29)
project(kblocks_Wannier90 LANGUAGES Fortran)

find_package(Wannier90 REQUIRED)
find_package(MPI REQUIRED COMPONENTS Fortran)

add_executable(kblocks kblocks.F90)
target_link_libraries(kblocks PRIVATE Wannier90::wannier90 MPI::MPI_Fortran)

enable_testing()
# The distributions are built by the program itself, so the only thing that varies here is
# the rank count: one rank for the trivial map, two for the smallest non-trivial one, four
# so that a mesh size which does not divide the rank count is covered as well.
foreach (np IN ITEMS 1 2 4)
    add_test(NAME kblocks_np${np}
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${MPIEXEC_PREFLAGS} $<TARGET_FILE:kblocks> ${MPIEXEC_POSTFLAGS}
    )
endforeach ()
