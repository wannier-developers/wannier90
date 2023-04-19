cmake_minimum_required(VERSION 3.18)

project(Wannier90
		VERSION 4.0.0
		DESCRIPTION "Compute maximally-localised Wannier functions."
		HOMEPAGE_URL https://www.wannier90.org
		LANGUAGES Fortran
		)

## Options
option(WANNIER90_WITH_MPI "Wannier90: Build library as a shared library" OFF)
option(WANNIER90_SHARED_LIBS "Wannier90: Build library as a shared library" ${PROJECT_IS_TOP_LEVEL})
option(WANNIER90_INSTALL "Wannier90: Install files" ${PROJECT_IS_TOP_LEVEL})
option(WANNIER90_TEST "Wannier90: Build unit tests" ${PROJECT_IS_TOP_LEVEL})

## Variable setup
if (NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif ()
if (WANNIER90_INSTALL)
	include(GNUInstallDirs)
	include(CMakePackageConfigHelpers)
endif ()
if (NOT CMAKE_Fortran_MODULE_DIRECTORY)
	set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fortran_mods)
endif ()
set(BUILD_SHARED_LIBS ${WANNIER90_SHARED_LIBS})

## Third-party libraries
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
if (WANNIER90_WITH_MPI)
	find_package(MPI REQUIRED)
endif ()

## Main targets
add_executable(Wannier90_exe)
set_target_properties(Wannier90_exe PROPERTIES
		OUTPUT_NAME wannier90.x
		EXPORT_NAME exe)
#add_library(Wannier90::exe ALIAS Wannier90_exe)
add_executable(Wannier90_post)
set_target_properties(Wannier90_post PROPERTIES
		OUTPUT_NAME postw90.x
		EXPORT_NAME post)
#add_library(Wannier90::post ALIAS Wannier90_post)
add_library(Wannier90_lib)
set_target_properties(Wannier90_lib PROPERTIES
		OUTPUT_NAME wannier90
		EXPORT_NAME wannier90
		VERSION ${PROJECT_VERSION}
		SOVERSION ${PROJECT_VERSION_MAJOR})
add_library(Wannier90::wannier90 ALIAS Wannier90_lib)

## Main implementations
add_subdirectory(src)

## Testing
if (WANNIER90_TEST)
	enable_testing()
	add_executable(w90chk2chk.x)
	add_executable(w90spn2spn.x)
	add_subdirectory(test-suite)
endif ()

## Installs
if (WANNIER90_INSTALL)
	configure_file(cmake/wannier90.pc.in wannier90.pc)
	install(FILES ${CMAKE_CURRENT_BINARY_DIR}/wannier90.pc
			DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
			COMPONENT Wannier90_Development)
	write_basic_package_version_file(
			${CMAKE_CURRENT_BINARY_DIR}/Wannier90ConfigVersion.cmake
			VERSION ${PROJECT_VERSION}
			COMPATIBILITY SameMajorVersion)
	configure_package_config_file(
			cmake/Wannier90Config.cmake.in
			Wannier90Config.cmake
			INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Wannier90)
	install(FILES ${CMAKE_CURRENT_BINARY_DIR}/Wannier90ConfigVersion.cmake
			${CMAKE_CURRENT_BINARY_DIR}/Wannier90Config.cmake
			DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Wannier90
			COMPONENT Spglib_Development)
	install(EXPORT Wannier90Targets
			FILE Wannier90Targets.cmake
			NAMESPACE ${ARGS_PROJECT}::
			DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Wannier90
			COMPONENT Wannier90_Development)
	export(EXPORT Wannier90Targets
			FILE Wannier90Targets.cmake
			NAMESPACE Wannier90::)
endif ()

## Make project available to FetchContent
if (NOT PROJECT_IS_TOP_LEVEL)
	# Propagate variables for FetchContent
	# All variables have to be consistent with CMakeExtraUtilsConfig.cmake
	set(Wannier90_MPI ${WANNIER90_WITH_MPI})
	if (CMAKE_VERSION VERSION_LESS 3.25)
		# TODO: Remove when cmake 3.25 is commonly distributed
		set(Wannier90_VERSION ${Wannier90_VERSION} PARENT_SCOPE)
		set(Wannier90_VERSION_MAJOR ${Wannier90_VERSION_MAJOR} PARENT_SCOPE)
		set(Wannier90_VERSION_MINOR ${Wannier90_VERSION_MINOR} PARENT_SCOPE)
		set(Wannier90_VERSION_PATCH ${Wannier90_VERSION_PATCH} PARENT_SCOPE)
		set(Wannier90_VERSION_TWEAK ${Wannier90_VERSION_TWEAK} PARENT_SCOPE)
		set(Wannier90_MPI ${Wannier90_MPI} PARENT_SCOPE)
	endif ()
	return(PROPAGATE
			Wannier90_VERSION
			Wannier90_VERSION_MAJOR
			Wannier90_VERSION_MINOR
			Wannier90_VERSION_PATCH
			Wannier90_VERSION_TWEAK
			Wannier90_MPI
			)
endif ()
