set_property(DIRECTORY APPEND
		PROPERTY LABELS functional library
)

set(test_names
		lib-gas-opt-pllel
		libc-gas-mix-serial
		libc-gas-mix-serial-clash
)
if (Wannier90_MPI)
	list(APPEND test_names
			libc-gas-opt-pllel-random
	)
endif ()

foreach (test IN LISTS test_names)
	Wannier90_add_test(${test}
			TEST_COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${test}/demo
	)
endforeach ()
