## List tests
set(test_names
                checkpoint0_write
                checkpoint1_read
#               partestw90_mpierr  # Disabled because test intended to fail
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
                testpostw90_fe_kubo_jdos
                testpostw90_fe_kubo_Szz
                testpostw90_fe_morb
                testpostw90_fe_morbandahc
                testpostw90_fe_morb_transl_inv
                testpostw90_fe_morb_transl_inv_higher
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
                testpostw90_pt_shc_higher
                testpostw90_pt_shc_ryoo
                testpostw90_pt_shc_ryoo_transl_inv
                testpostw90_pt_tetra_shcfermi
                testpostw90_pt_tetra_shcfreq
                testpostw90_si_geninterp
                testpostw90_si_geninterp_wsdistance
                testpostw90_te_gyrotropic
                testpostw90_te_gyrotropic_C
                testpostw90_te_gyrotropic_D0
                testpostw90_te_gyrotropic_dos
                testpostw90_te_gyrotropic_Dw
                testpostw90_te_gyrotropic_K
                testpostw90_te_gyrotropic_NOA
                testw90_basic1
                testw90_basic2
                testw90_benzene_gamma_val
                testw90_benzene_gamma_valcond
                testw90_benzene_gamma_val_hexcell
                testw90_bvec
                testw90_cube_format
                testw90_disentanglement_sawfs
                testw90_example01
                testw90_example02
                testw90_example02_restart
                testw90_example03
                testw90_example03_bands_plot
                testw90_example03_explicit_kpath_band_kpt
                testw90_example03_explicit_kpath_bands_plot
                testw90_example03_explicit_kpath_labelinfo
                testw90_example03_labelinfo
                testw90_example03_optmem
                testw90_example04
                testw90_example05
                testw90_example07
                testw90_example11_1
#                testw90_example11_2  # Dropped because CI's old/unfortunate intel/flexiblas versions give too large numerical fuzz
                testw90_example21_As_sp
                testw90_example26
                testw90_example36
                testw90_gaas_disentanglement_issue192
                testw90_graphene_dis_proj_auto
                testw90_graphene_pdwf
                testw90_guidingcentre_selectproj
                testw90_knbo3_higher
                testw90_knbo3_higher_stengel_spaldin
                testw90_lavo3_dissphere
#		testw90_na_chain_gamma  # Dropped because flaky
                testw90_nnkpt1
                testw90_nnkpt2
                testw90_nnkpt3
                testw90_nnkpt6
                testw90_nnkpt7
                testw90_precond_1
                testw90_precond_2
                testw90_proj1
                testw90_proj2
                testw90_proj3
                testw90_rmn
                testw90_write_u_matrices
                testw90_write_u_matrices_disent
                testw90_readin_b_vec
                testw90_nnkpt4 # expects to fail MUST BE AT END OF LIST (MPI subsystem left defunct)
                testw90_nnkpt5 # expects to fail MUST BE AT END OF LIST (MPI subsystem left defunct)
)
if (Wannier90_MPI)
#	list(APPEND test_names
#			partestw90_mpierr # dropped because intended to fail
#	)
endif ()

## Define tests
foreach (test IN LISTS test_names)
	set(extra_args)
	if (Wannier90_MPI)
		list(APPEND extra_args -n 2)
		set(num_procs 2)
	else ()
		set(num_procs 1)
	endif ()
	add_test(NAME ${test}
			COMMAND Python::Interpreter ${CMAKE_CURRENT_BINARY_DIR}/../run_tests -c ${test} ${extra_args}
			WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..
	)
	set_tests_properties(${test} PROPERTIES
			PROCESSORS ${num_procs}
	)
endforeach ()
