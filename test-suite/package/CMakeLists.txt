set_property(DIRECTORY APPEND
		PROPERTY LABELS package
)

set(test_list
		FetchContent
		find_package
)
if (NOT Wannier90_IS_TOP_LEVEL)
	list(APPEND test_list pkg-config)
endif ()

foreach (test IN LISTS test_list)
	Wannier90_add_test(${test} TEST_NAME test-${test})
endforeach ()
