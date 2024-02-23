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

## Installs
if (WANNIER90_INSTALL)
	install(TARGETS Wannier90_post
			EXPORT Wannier90Targets
			LIBRARY COMPONENT Wannier90_Runtime NAMELINK_COMPONENT Wannier90_Development
			ARCHIVE COMPONENT Wannier90_Development
			PUBLIC_HEADER COMPONENT Wannier90_Development
			RUNTIME COMPONENT Wannier90_Runtime
	)
endif ()
