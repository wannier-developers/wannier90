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
		overlap.F90
		plot.F90
		readwrite.F90
		sitesym.F90
		transport.F90
		types.F90
		utility.F90
		wannier90_readwrite.F90
		wannier90_types.F90
		library_interface.F90
		wannierise.F90
		ws_distance.F90
		)
target_sources(Wannier90_exe PRIVATE
		wannier_prog.F90
		)

## Linking targets
target_link_libraries(Wannier90_exe PRIVATE
		Wannier90_lib)
target_link_libraries(Wannier90_lib PRIVATE
		BLAS::BLAS LAPACK::LAPACK)
target_include_directories(Wannier90_lib PUBLIC
		"$<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>")
if (WANNIER90_WITH_MPI)
	target_link_libraries(Wannier90_lib PUBLIC MPI::MPI_Fortran)
endif ()

add_subdirectory(postw90)

## Installs
if (WANNIER90_INSTALL)
	target_include_directories(Wannier90_lib PUBLIC
			"$<INSTALL_INTERFACE:${CMAKE_INSTALL_MODULEDIR}>")
	install(TARGETS Wannier90_lib
			EXPORT Wannier90Targets
			LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Wannier90_Runtime NAMELINK_COMPONENT Wannier90_Development
			ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Wannier90_Development
			)
	install(FILES ${CMAKE_Fortran_MODULE_DIRECTORY}/w90_helper_types.mod
			DESTINATION ${CMAKE_INSTALL_MODULEDIR}
			COMPONENT Wannier90_Development
	)
	install(TARGETS Wannier90_exe
			EXPORT Wannier90Targets
			RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Wannier90_Runtime
			)
	install(TARGETS Wannier90_post
			EXPORT Wannier90Targets
			RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Wannier90_Runtime
			)
endif ()
