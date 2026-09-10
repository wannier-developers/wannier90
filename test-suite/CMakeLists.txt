cmake_minimum_required(VERSION 3.25...3.29)

#[=============================================================================[
#                           Basic project definition                           #
]=============================================================================]

list(APPEND CMAKE_MESSAGE_CONTEXT Test)
project(Wannier90_test
		LANGUAGES CXX Fortran
)

enable_testing()

#[=============================================================================[
#                                   Options                                   #
]=============================================================================]

#[=============================================================================[
#                            Project configuration                            #
]=============================================================================]

#[=============================================================================[
#                              External packages                              #
]=============================================================================]

# 3.10 is the floor: it gives every syntax feature the harness uses and is available on
# RHEL/Rocky 9 via appstream, Debian bookworm and Ubuntu 22.04. Raising it further is a
# one-line change.
find_package(Python 3.10 REQUIRED COMPONENTS Interpreter)

# The regression suite is driven by pytest and reads YAML. Rather than failing the whole
# configure for someone who only wants the binaries, warn and skip registering the tests.
execute_process(
		COMMAND ${Python_EXECUTABLE} -c "import pytest, yaml"
		RESULT_VARIABLE _w90_pytest_missing
		OUTPUT_QUIET ERROR_QUIET
)
if (_w90_pytest_missing)
	message(WARNING
			"Python packages 'pytest' and/or 'PyYAML' not found for ${Python_EXECUTABLE}.\n"
			"The regression tests will NOT be registered. To enable them:\n"
			"    pip install -r ${CMAKE_CURRENT_SOURCE_DIR}/requirements.txt\n"
			"then re-run cmake."
	)
endif ()

#[=============================================================================[
#                               Main definition                               #
]=============================================================================]

function(Wannier90_add_test test)
	#[===[.md
	# Wannier90_add_test

	Internal helper for adding functional tests specific for Wannier90

	## Synopsis
	```cmake
	Wannier90_add_test(<name>
			[TEST_NAME <test_name>]
			[LABELS <label1> <label2>]
			[TEST_COMMAND <arg1> <arg2>]
			)
	```

	## Options

	: `<name>`
	  Path to the CMake project to be executed relative to `${CMAKE_CURRENT_SOURCE_DIR}`

	: `TEST_NAME` [Default: `<name>`]
	  Name for the test to be used as the ctest name

	: `LABELS`
	  Additional labels to be added

	: `TEST_COMMAND` [Default: ctest]
	  Test command to use instead of ctest

	]===]

	list(APPEND CMAKE_MESSAGE_CONTEXT "Wannier90_add_test")

	set(ARGS_Options)
	set(ARGS_OneValue
			TEST_NAME
	)
	set(ARGS_MultiValue
			LABELS
			TEST_COMMAND
	)
	cmake_parse_arguments(PARSE_ARGV 1 ARGS "${ARGS_Options}" "${ARGS_OneValue}" "${ARGS_MultiValue}")
	# Check required/optional arguments
	if (ARGC LESS 1)
		message(FATAL_ERROR "Missing test name")
	endif ()
	if (NOT DEFINED ARGS_TEST_NAME)
		set(ARGS_TEST_NAME ${test})
	endif ()
	if (NOT DEFINED ARGS_TEST_COMMAND)
		set(ARGS_TEST_COMMAND
				${CMAKE_CTEST_COMMAND} --test-dir ${CMAKE_CURRENT_BINARY_DIR}/${test} --output-on-failure
		)
	endif ()

	set(extra_args)
	if (Wannier90_IS_TOP_LEVEL)
		list(APPEND extra_args
				-DFETCHCONTENT_TRY_FIND_PACKAGE_MODE=ALWAYS
				# Generated Config file point to binary targets until it is installed
				-DWannier90_ROOT=${Wannier90_BINARY_DIR}
				-DFETCHCONTENT_SOURCE_DIR_WANNIER90=${Wannier90_SOURCE_DIR}
		)
	endif ()
	set(compiler_args)
	foreach (lang IN ITEMS CXX Fortran)
		if (DEFINED CMAKE_${lang}_COMPILER)
			list(APPEND compiler_args -DCMAKE_${lang}_COMPILER=${CMAKE_${lang}_COMPILER})
		endif ()
	endforeach ()

	add_test(NAME ${ARGS_TEST_NAME}
			COMMAND ${CMAKE_CTEST_COMMAND} --build-and-test ${CMAKE_CURRENT_SOURCE_DIR}/${test}
			${CMAKE_CURRENT_BINARY_DIR}/${test}
			# Use the same build environment as the current runner
			--build-generator "${CMAKE_GENERATOR}"
			--build-options
			${compiler_args} ${extra_args}
			--test-command ${ARGS_TEST_COMMAND}
	)
	set_tests_properties(${ARGS_TEST_NAME} PROPERTIES
			LABELS "${ARGS_LABELS}"
	)
endfunction()

## Add regression tests
# Nothing is copied into the build tree: pytest reads the source tree directly and runs
# each test in its own work directory under the build tree. The previous file(COPY ...)
# block duplicated ~240 MB on every configure.
if (NOT _w90_pytest_missing)
	add_subdirectory(tests)
endif ()

## Functional tests
set(functional_tests
)
if (WANNIER90_WITH_C)
	# These test requires that the C interface is built
        # (which requires a --std=f2018 compliant compiler
	list(APPEND functional_tests
			library-mode-test-C-interface
        )
endif ()
if (WANNIER90_MPI)
	# These tests require MPI interface
	list(APPEND functional_tests
			library-mode-test-proj-overflow
			library-mode-test
			library-mode-test-projections
	)
endif ()
foreach (test IN LISTS functional_tests)
	Wannier90_add_test(${test})
endforeach ()
