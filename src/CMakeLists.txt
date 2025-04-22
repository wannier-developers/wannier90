#!--------------------------------------------------------------!
#! This file is distributed as part of the Wannier90 code and   !
#! under the terms of the                                       !
#!                                                              !
#!        GNU Lesser General Public License, version 2.1        !
#!                                                              !
#! see file `LICENSE' in the root directory of the Wannier90    !
#! distribution, or https://www.gnu.org/licenses/lgpl-2.1.html  !
#!                                                              !
#! The webpage of the Wannier90 code is www.wannier.org         !
#!                                                              !
#! The Wannier90 code is hosted on GitHub                       !
#! https://github.com/wannier-developers/wannier90              !
#!--------------------------------------------------------------!

## Source files definitions
target_sources(Wannier90_lib PRIVATE
		c_interface.F90
		comms.F90
		constants.F90
		disentangle.F90
		error.F90
		error_base.F90
		hamiltonian.F90
		io.F90
		kmesh.F90
		library_extra.F90
		library_interface.F90
		overlap.F90
		plot.F90
		readwrite.F90
		sitesym.F90
		transport.F90
		types.F90
		utility.F90
		wannier90_readwrite.F90
		wannier90_types.F90
		wannierise.F90
		ws_distance.F90
		)
target_sources(Wannier90_exe PRIVATE
		wannier_prog.F90
		)

## Linking targets
target_link_libraries(Wannier90_exe PRIVATE
		Wannier90_lib
)
target_link_libraries(Wannier90_lib PRIVATE
		BLAS::BLAS LAPACK::LAPACK
)
target_include_directories(Wannier90_lib PUBLIC
		"$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
		"$<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>"
		"$<INSTALL_INTERFACE:${CMAKE_INSTALL_MODULEDIR}/Wannier90>"
)
target_sources(Wannier90_lib
		INTERFACE FILE_SET wannier90_C_Headers TYPE HEADERS FILES
		wannier90.h
		INTERFACE FILE_SET wannier90_Fortran_Headers TYPE HEADERS BASE_DIRS ${CMAKE_Fortran_MODULE_DIRECTORY} FILES
		${CMAKE_Fortran_MODULE_DIRECTORY}/w90_library.mod
)
if (WANNIER90_MPI)
	target_link_libraries(Wannier90_lib PUBLIC MPI::MPI_Fortran)
	target_compile_definitions(Wannier90_lib PRIVATE
			# TODO: Use mpi_f08 instead. Tests are not adapted yet
			MPI90 MPI
	)
	target_compile_definitions(Wannier90_exe PRIVATE
			# TODO: Use mpi_f08 instead. Tests are not adapted yet
			MPI90 MPI
	)
endif ()

add_subdirectory(postw90)

## Installs
if (WANNIER90_INSTALL)
	install(TARGETS Wannier90_lib Wannier90_exe
			EXPORT Wannier90Targets
			FILE_SET wannier90_C_Headers
			FILE_SET wannier90_Fortran_Headers DESTINATION ${CMAKE_INSTALL_MODULEDIR}/Wannier90
			)
endif ()
