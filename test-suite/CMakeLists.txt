set(test_names "")
list(APPEND test_names
		partestw90_mpierr
		testpostw90_boltzwann
		testpostw90_example04_dos
		testpostw90_example04_pdos
		testpostw90_fe_ahc
		testpostw90_fe_ahc_adaptandfermi
		testpostw90_fe_dos_spin
		testpostw90_fe_kpathcurv
		testpostw90_fe_kpathmorbcurv
		testpostw90_fe_kpathmorbcurv_ws
		testpostw90_fe_kslicecurv
		testpostw90_fe_kslicemorb
		testpostw90_fe_kubo_Axy
		testpostw90_fe_kubo_Szz
		testpostw90_fe_kubo_jdos
		testpostw90_fe_morb
		testpostw90_fe_morbandahc
		testpostw90_fe_spin
		testpostw90_gaas_kdotp
		testpostw90_gaas_sc_eta_corr
		testpostw90_gaas_sc_xyz
		testpostw90_gaas_sc_xyz_scphase2
		testpostw90_gaas_sc_xyz_scphase2_ws
		testpostw90_gaas_sc_xyz_ws
		testpostw90_gaas_shc
		testpostw90_pt_kpathbandsshc
		testpostw90_pt_kpathshc
		testpostw90_pt_ksliceshc
		testpostw90_pt_shc
		testpostw90_pt_shc_ryoo
		testpostw90_si_geninterp
		testpostw90_si_geninterp_wsdistance
		testpostw90_te_gyrotropic
		testpostw90_te_gyrotropic_C
		testpostw90_te_gyrotropic_D0
		testpostw90_te_gyrotropic_Dw
		testpostw90_te_gyrotropic_K
		testpostw90_te_gyrotropic_NOA
		testpostw90_te_gyrotropic_dos
		testw90_basic1
		testw90_basic2
		testw90_benzene_gamma_val
		testw90_benzene_gamma_val_hexcell
		testw90_benzene_gamma_valcond
		testw90_bvec
		testw90_cube_format
		testw90_disentanglement_sawfs
		testw90_example01
		testw90_example02
		testw90_example02_restart
		testw90_example03
		testw90_example03_labelinfo
		testw90_example03_optmem
		testw90_example04
		testw90_example05
		testw90_example07
		testw90_example11_1
		testw90_example11_2
		testw90_example21_As_sp
		testw90_example26
		testw90_gaas_disentanglement_issue192
		testw90_lavo3_dissphere
		testw90_na_chain_gamma
		testw90_nnkpt1
		testw90_nnkpt2
		testw90_nnkpt3
		testw90_nnkpt4
		testw90_nnkpt5
		testw90_precond_1
		testw90_precond_2
		testw90_write_u_matrices
		)
set(extra_files "")
list(APPEND extra_files
		checkpoints
		run_tests
		testcode
		tests
		tools
		)

## Prepare test environment
foreach (extra IN LISTS extra_files)
	file(COPY ${extra} DESTINATION test-suite)
endforeach ()
target_link_libraries(w90chk2chk.x PRIVATE Wannier90_lib)
target_sources(w90chk2chk.x PRIVATE
		${PROJECT_SOURCE_DIR}/src/w90chk2chk.F90
		)
target_link_libraries(w90spn2spn.x PRIVATE Wannier90_lib)
target_sources(w90spn2spn.x PRIVATE
		${PROJECT_SOURCE_DIR}/src/w90spn2spn.F90
		)

## Define tests
foreach (test IN LISTS test_names)
	add_test(NAME ${test}
			COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/run_tests -c ${test}
			WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
			)
endforeach ()
