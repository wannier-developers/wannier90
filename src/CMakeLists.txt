#------------------------------------------------------------#
# Copyright (C) 2026 Wannier Developer Group                 #
#                                                            #
# This library is free software; you can redistribute it     #
# and/or modify it under the terms of the GNU Lesser General #
# Public License as published by the Free Software           #
# Foundation; either version 2.1 of the License, or (at your #
# option) any later version.                                 #
#                                                            #
# This library is distributed in the hope that it will be    #
# useful,but WITHOUT ANY WARRANTY; without even the implied  #
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR    #
# PURPOSE.  See the GNU Lesser General Public License for    #
# more details.                                              #
#                                                            #
# You should have received a copy of the GNU Lesser General  #
# Public License along with this library; if not, see        #
# <https://www.gnu.org/licenses/>.                           #
#                                                            #
# The webpage of the Wannier90 code is                       #
# <https://www.wannier.org>.                                 #
#                                                            #
# The Wannier90 code is hosted on GitHub                     #
# <https://github.com/wannier-developers/wannier90>          #
#------------------------------------------------------------#

## Source files definitions
target_sources(Wannier90_lib PRIVATE
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
if (WANNIER90_WITH_C)
	target_sources(Wannier90_lib PRIVATE
			c_interface.F90
			wannier90.c
			)
endif ()
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
if (WANNIER90_WITH_C)
	target_sources(Wannier90_lib
			INTERFACE FILE_SET wannier90_C_Headers TYPE HEADERS FILES
			wannier90.h
	)
endif ()
target_sources(Wannier90_lib
		INTERFACE FILE_SET wannier90_Fortran_Headers TYPE HEADERS BASE_DIRS ${CMAKE_Fortran_MODULE_DIRECTORY} FILES
		${CMAKE_Fortran_MODULE_DIRECTORY}/w90_library.mod
)
if (WANNIER90_MPI)
	set(_wannier90_mpi_interface MPI90)
	if (WANNIER90_MPIH)
		set(_wannier90_mpi_interface MPIH)
	endif ()

	target_link_libraries(Wannier90_lib PUBLIC MPI::MPI_Fortran)
	target_compile_definitions(Wannier90_lib PRIVATE
			# TODO: Use mpi_f08 instead. Tests are not adapted yet
			MPI ${_wannier90_mpi_interface}
	)
	target_compile_definitions(Wannier90_exe PRIVATE
			# TODO: Use mpi_f08 instead. Tests are not adapted yet
			MPI ${_wannier90_mpi_interface}
	)
	if (WANNIER90_MPIH AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
		target_compile_options(Wannier90_lib PRIVATE -fallow-argument-mismatch)
		target_compile_options(Wannier90_exe PRIVATE -fallow-argument-mismatch)
	endif ()
endif ()

add_subdirectory(postw90)

## Installs
if (WANNIER90_INSTALL)
	if (WANNIER90_WITH_C)
		install(TARGETS Wannier90_lib Wannier90_exe
				EXPORT Wannier90Targets
				FILE_SET wannier90_C_Headers
				FILE_SET wannier90_Fortran_Headers DESTINATION ${CMAKE_INSTALL_MODULEDIR}/Wannier90
				)
	else ()
		install(TARGETS Wannier90_lib Wannier90_exe
				EXPORT Wannier90Targets
				FILE_SET wannier90_Fortran_Headers DESTINATION ${CMAKE_INSTALL_MODULEDIR}/Wannier90
				)
	endif ()
endif ()
