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

find_package(Python 3.8 REQUIRED)

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

## Prepare test environment
foreach (extra_file IN ITEMS
		checkpoints
		run_tests
		testcode
		tests
		tools
)
	# TODO: Use `PYTHONPATH` to avoid needing to copy and setting workdir
	file(COPY ${extra_file} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endforeach ()

## Add unit tests
add_subdirectory(tests)
## Add functional tests
add_subdirectory(libv2-demo)
