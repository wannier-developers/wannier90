target_sources(Wannier90_post PRIVATE
		postw90.F90
		)
target_sources(Wannier90_lib PRIVATE
		berry.F90
		boltzwann.F90
		dos.F90
		geninterp.F90
		get_oper.F90
		gyrotropic.F90
		kpath.F90
		kslice.F90
		postw90.F90
		postw90_common.F90
		postw90_readwrite.F90
		postw90_types.F90
		spin.F90
		wan_ham.F90
		)
target_link_libraries(Wannier90_post PRIVATE
		Wannier90_lib BLAS::BLAS LAPACK::LAPACK)
if (WANNIER90_MPI)
	target_compile_definitions(Wannier90_post PRIVATE
			# TODO: Use mpi_f08 instead. Tests are not adapted yet
			MPI ${_wannier90_mpi_interface}
	)
	if (WANNIER90_MPIH AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
		target_compile_options(Wannier90_post PRIVATE -fallow-argument-mismatch)
	endif ()
endif ()

## Installs
if (WANNIER90_INSTALL)
	install(TARGETS Wannier90_post
			EXPORT Wannier90Targets
	)
endif ()
