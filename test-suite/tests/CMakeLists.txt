## Regression tests
#
# Tests are discovered from tests/*/test.yaml -- there is no hand-maintained list to keep
# in sync. Adding a test means creating a directory with a test.yaml in it, nothing else.
#
# Each ctest case runs pytest on one test directory, out of tree, with its own work
# directory under the build tree.

file(GLOB test_yaml_files CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*/test.yaml")

if (NOT test_yaml_files)
	message(WARNING "No tests/*/test.yaml found; no regression tests will be registered")
endif ()

set(pytest_extra_args)
if (Wannier90_MPI)
	# Matches the previous behaviour: with an MPI build, ctest exercises the parallel path.
	# Tests marked `serial_only` still run serially; the harness enforces that.
	list(APPEND pytest_extra_args --nprocs=2)
	set(num_procs 2)
else ()
	set(num_procs 1)
endif ()

foreach (test_yaml IN LISTS test_yaml_files)
	get_filename_component(test_dir ${test_yaml} DIRECTORY)
	get_filename_component(test ${test_dir} NAME)

	# Labels mirror the pytest markers, so `ctest -L postw90` and `pytest -m postw90`
	# select the same set.
	set(test_labels "")
	if (test MATCHES "^partest")
		set(test_labels parallel)
	elseif (test MATCHES "^testpostw90_")
		set(test_labels postw90)
	elseif (test MATCHES "^testw90_")
		set(test_labels wannier90)
	elseif (test MATCHES "^checkpoint")
		set(test_labels checkpoint)
	endif ()

	add_test(NAME ${test}
			COMMAND ${Python_EXECUTABLE} -m pytest ${test_dir}
			--wannier90-exe=$<TARGET_FILE:Wannier90_exe>
			--postw90-exe=$<TARGET_FILE:Wannier90_post>
			--w90chk2chk-exe=$<TARGET_FILE:w90chk2chk.x>
			--workdir=${CMAKE_CURRENT_BINARY_DIR}/work/${test}
			${pytest_extra_args}
			-p no:cacheprovider
			WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/..
	)
	set_tests_properties(${test} PROPERTIES PROCESSORS ${num_procs})
	if (test_labels)
		set_tests_properties(${test} PROPERTIES LABELS "${test_labels}")
	endif ()

	# `status: skip` in test.yaml disables the test in ctest too, so the two agree.
	file(STRINGS ${test_yaml} status_line REGEX "^status:[ \t]*skip")
	if (status_line)
		set_tests_properties(${test} PROPERTIES DISABLED TRUE)
	endif ()
endforeach ()

list(LENGTH test_yaml_files num_tests)
message(STATUS "Registered ${num_tests} Wannier90 regression tests")
