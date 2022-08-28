function(get_mpi_vendor)

    execute_process(COMMAND mpirun --version OUTPUT_VARIABLE MPIRUN_OUTPUT)

    string(FIND "${MPIRUN_OUTPUT}" "Open MPI" OMPI_POS)
    string(FIND "${MPIRUN_OUTPUT}" "MPICH" MPICH_POS)
    string(FIND "${MPIRUN_OUTPUT}" "Intel(R) MPI" IMPI_POS)

    if(NOT OMPI_POS STREQUAL "-1")
      set(MPI_VENDOR "OpenMPI" PARENT_SCOPE)
    elseif(NOT MPICH_POS STREQUAL "-1")
      set(MPI_VENDOR "MPICH" PARENT_SCOPE)
    elseif(NOT IMPI_POS STREQUAL "-1")
      set(MPI_VENDOR "IntelMPI" PARENT_SCOPE)
    endif()

endfunction()
