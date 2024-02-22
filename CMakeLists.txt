cmake_minimum_required(VERSION 3.25...3.29)

#[=============================================================================[
#                           Basic project definition                           #
]=============================================================================]

list(APPEND CMAKE_MESSAGE_CONTEXT Wannier90)
project(Wannier90
		VERSION 4.0.0
		DESCRIPTION "Compute maximally-localised Wannier functions."
		HOMEPAGE_URL https://www.wannier90.org
		LANGUAGES Fortran C
		)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

if (NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif ()

#[=============================================================================[
#                                   Options                                   #
]=============================================================================]

option(WANNIER90_MPI "Wannier90: Build with MPI support" OFF)
option(WANNIER90_SHARED_LIBS "Wannier90: Build library as a shared library" ${PROJECT_IS_TOP_LEVEL})
option(WANNIER90_INSTALL "Wannier90: Install files" ${PROJECT_IS_TOP_LEVEL})
option(WANNIER90_TEST "Wannier90: Build test-suite" ${PROJECT_IS_TOP_LEVEL})

#[=============================================================================[
#                            Project configuration                            #
]=============================================================================]

if (NOT CMAKE_Fortran_MODULE_DIRECTORY)
	set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fortran_mods)
endif ()
set(BUILD_SHARED_LIBS ${WANNIER90_SHARED_LIBS})

if (WANNIER90_INSTALL)
	include(GNUInstallDirs)
	include(CMakePackageConfigHelpers)

	# CMake does not properly support fortran module installation paths.
	# https://gitlab.kitware.com/cmake/cmake/-/issues/19608
	# https://discourse.cmake.org/t/api-design-c-modules-source-listings-and-interface-properties/5389/14
	# Note: This does not follow fortran-stdlib format, instead it follows Fedora guidelines
	# https://pagure.io/packaging-committee/issue/1334
	cmake_path(APPEND CMAKE_INSTALL_LIBDIR "${CMAKE_Fortran_COMPILER_ID}-${CMAKE_Fortran_COMPILER_VERSION}" modules
			OUTPUT_VARIABLE _DEFAULT_CMAKE_INSTALL_MODULEDIR
	)
	set(CMAKE_INSTALL_MODULEDIR ${_DEFAULT_CMAKE_INSTALL_MODULEDIR}
			CACHE STRING
			"Fortran module installation path (Not a cmake native variable)"
	)
	cmake_path(IS_ABSOLUTE CMAKE_INSTALL_MODULEDIR _is_absolute)
	if (_is_absolute)
		set(CMAKE_INSTALL_FULL_MODULEDIR ${CMAKE_INSTALL_MODULEDIR})
	else ()
		cmake_path(APPEND CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_MODULEDIR}
				OUTPUT_VARIABLE CMAKE_INSTALL_FULL_MODULEDIR
		)
	endif ()
endif ()

## Populate variables for FetchContent and sub-projects
# All variables have to be consistent with CMakeExtraUtilsConfig.cmake
set(Wannier90_MPI ${WANNIER90_WITH_MPI})

#[=============================================================================[
#                              External packages                              #
]=============================================================================]

## Third-party libraries
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
if (WANNIER90_MPI)
	find_package(MPI REQUIRED)
endif ()

#[=============================================================================[
#                               Main definition                               #
]=============================================================================]

## Main targets
add_executable(Wannier90_exe)
set_target_properties(Wannier90_exe PROPERTIES
		OUTPUT_NAME wannier90.x
		EXPORT_NAME exe
)
add_executable(Wannier90::exe ALIAS Wannier90_exe)
add_executable(Wannier90_post)
set_target_properties(Wannier90_post PROPERTIES
		OUTPUT_NAME postw90.x
		EXPORT_NAME post
)
add_executable(Wannier90::post ALIAS Wannier90_post)
add_library(Wannier90_lib)
set_target_properties(Wannier90_lib PROPERTIES
		OUTPUT_NAME wannier90
		EXPORT_NAME wannier90
		VERSION ${PROJECT_VERSION}
		SOVERSION ${PROJECT_VERSION_MAJOR}
)
add_library(Wannier90::wannier90 ALIAS Wannier90_lib)

## Main implementations
add_subdirectory(src)

## Testing
if (WANNIER90_TEST)
	enable_testing()
	## Define compiled support programs
	# TODO: These should be moved outside of src and hard coded path?
	add_executable(w90chk2chk.x src/w90chk2chk.F90)
	add_executable(w90spn2spn.x src/w90spn2spn.F90)
	target_link_libraries(w90chk2chk.x PRIVATE Wannier90_lib)
	target_link_libraries(w90spn2spn.x PRIVATE Wannier90_lib)

	add_subdirectory(test-suite)
endif ()

#[=============================================================================[
#                              Install or Export                              #
]=============================================================================]

## Installs
if (WANNIER90_INSTALL)
	configure_file(cmake/wannier90.pc.in wannier90.pc)
	install(FILES ${CMAKE_CURRENT_BINARY_DIR}/wannier90.pc
			DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
			COMPONENT Wannier90_Development
	)
	write_basic_package_version_file(
			${CMAKE_CURRENT_BINARY_DIR}/Wannier90ConfigVersion.cmake
			VERSION ${PROJECT_VERSION}
			COMPATIBILITY SameMajorVersion
	)
	configure_package_config_file(
			cmake/Wannier90Config.cmake.in
			Wannier90Config.cmake
			INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Wannier90
	)
	install(FILES ${CMAKE_CURRENT_BINARY_DIR}/Wannier90ConfigVersion.cmake
			${CMAKE_CURRENT_BINARY_DIR}/Wannier90Config.cmake
			DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Wannier90
			COMPONENT Spglib_Development
	)
	install(EXPORT Wannier90Targets
			FILE Wannier90Targets.cmake
			NAMESPACE Wannier90::
			DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Wannier90
			COMPONENT Wannier90_Development
	)
	export(EXPORT Wannier90Targets
			FILE Wannier90Targets.cmake
			NAMESPACE Wannier90::
	)
endif ()

## Make project available to FetchContent
if (NOT PROJECT_IS_TOP_LEVEL)
	# Propagate variables for FetchContent
	return(PROPAGATE
			Wannier90_VERSION
			Wannier90_VERSION_MAJOR
			Wannier90_VERSION_MINOR
			Wannier90_VERSION_PATCH
			Wannier90_VERSION_TWEAK
			Wannier90_MPI
			)
endif ()
