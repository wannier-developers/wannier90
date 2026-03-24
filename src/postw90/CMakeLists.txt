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
			MPI90 MPI
	)
endif ()

## Installs
if (WANNIER90_INSTALL)
	install(TARGETS Wannier90_post
			EXPORT Wannier90Targets
	)
endif ()
