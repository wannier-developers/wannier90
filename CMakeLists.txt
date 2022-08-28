# using 3.18 to have nice lapack imported target
cmake_minimum_required(VERSION 3.18)

project(wannier90
  LANGUAGES C CXX Fortran
  VERSION 3.1.0
  HOMEPAGE_URL http://www.wannier.org/
  DESCRIPTION "Computational chemistry package for computing maximally-localised Wannier functions.")

#
# Build options
#

# Do we want want shared or static libraries ?
option(WANNIER_SHARED_LIB "Build shared library ? (Default: OFF, i.e. static libraries)" OFF)

# build doc ?
option(BUILD_DOC "Build doc ?" OFF)

# Travis CI enabled ?
option(TRAVISCI_ENABLED "Use custom flags for Travis CI ? (Default: OFF)" OFF)

#
# Export compile command to json (for editors like emacs, clion, vim, etc...).
# It allows nice code editing features provided by LSP (Language Server Protocol)
#
set(CMAKE_EXPORT_COMPILE_COMMANDS on)

#
# default local cmake macro repository (located in wannier90 top level sources)
#
list(INSERT CMAKE_MODULE_PATH 0 "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

#
# Enforce user to use a build directory outside of source tree
#
include(prevent_in_source_build)
prevent_in_source_build()

#
# Set default compile optimization flag
#
set(wannier90_BUILD_TYPE "Release" CACHE STRING
  "Optimization flags: set to Debug, Release, RelWithDebInfo, or MinSizeRel")

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "${wannier90_BUILD_TYPE}" CACHE INTERNAL "" FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif ()

#
# Mixing Fortran and C
#
include(FortranCInterface)
FortranCInterface_VERIFY()
FortranCInterface_HEADER(
  fc_mangle.h
  MACRO_NAMESPACE "FC_"
  )

string(TOLOWER ${CMAKE_BUILD_TYPE} build_type)
if (build_type STREQUAL debug)
  set(DEBUG_MODE 1)
endif()

#------------------------------------------------------------------------------#
#
# Additional compiler flags, compiler specific
#
#------------------------------------------------------------------------------#
include(SetupCompiler)


# BLAS / LAPACK setup
#
# this is controlled by cmake var BLA_VENDOR
# available values slightly depends on cmake version
# currently in cmake 3.18, valid values are
#
#    Goto
#    OpenBLAS
#    FLAME
#    ATLAS PhiPACK
#    CXML
#    DXML
#    SunPerf
#    SCSL
#    SGIMATH
#    IBMESSL
#    Intel10_32 (intel mkl v10 32 bit)
#    Intel10_64lp (intel mkl v10+ 64 bit, threaded code, lp64 model)
#    Intel10_64lp_seq (intel mkl v10+ 64 bit, sequential code, lp64 model)
#    Intel10_64ilp (intel mkl v10+ 64 bit, threaded code, ilp64 model)
#    Intel10_64ilp_seq (intel mkl v10+ 64 bit, sequential code, ilp64 model)
#    Intel10_64_dyn (intel mkl v10+ 64 bit, single dynamic library)
#    Intel (obsolete versions of mkl 32 and 64 bit)
#    ACML
#    ACML_MP
#    ACML_GPU
#    Apple
#    NAS
#    Arm
#    Arm_mp
#    Arm_ilp64
#    Arm_ilp64_mp
#    Generic
#
# NOTE:
# - if variable BLA_VENDOR is not set, all vendors are considered.
# - if multiple implementation are available on your system, the first to be detected,
#   will be considered
#
# Recommended value for MKL : Intel10_32
#
# Example
# - cmake configure step (from source top-level)
#   cmake -DBLA_VENDOR=OpenBLAS -S . -B build/openblas
# - cmake build project
#   cmake --build build/openblas
#
# See doc https://cmake.org/cmake/help/latest/manual/cmake.1.html
#
find_package(BLAS)
find_package(LAPACK)
if (NOT BLAS_FOUND OR NOT LAPACK_FOUND)
  message(fatal_error "BLAS or LAPACK were not found, please adjust variable BLA_VENDOR to select blas/lapack implementation and/or install blas/lapack implementation.")
endif()

#
# MPI
#
set(MPI_DETERMINE_LIBRARY_VERSION TRUE)
find_package(MPI COMPONENTS Fortran)

#
# check MPI fortran interface
#
if(MPI_FOUND)
  add_compile_definitions(MPI)

  if (MPI_Fortran_HAVE_F08_MODULE)
    add_compile_definitions(MPI08)
    message("MPI fortran interface f08 supported, MPI_Fortran_MODULE_DIR=${MPI_Fortran_MODULE_DIR}")
  elseif(MPI_Fortran_HAVE_F90_MODULE)
    add_compile_definitions(MPIF90)
  elseif(MPI_Fortran_HAVE_F77_HEADER)
    add_compile_definitions(MPIH)
  else()
    # default to f90 style "use mpi"
    add_compile_definitions(MPIF90)
  endif()

endif(MPI_FOUND)

# determine mpi vendor (MPI_VENDOR)
# valid values are : OpenMPI, MPICH or IntelMPI
include(get_mpi_vendor)
get_mpi_vendor()

if (NOT DEFINED MPI_VERSION AND DEFINED MPI_Fortran_VERSION)
  set(MPI_VERSION ${MPI_Fortran_VERSION})
  set(MPI_VERSION_MAJOR ${MPI_Fortran_VERSION_MAJOR})
  set(MPI_VERSION_MINOR ${MPI_Fortran_VERSION_MINOR})
endif()

#
# OpenMP
#
find_package(OpenMP)


##################### PRINT CONFIGURE STATUS ######################
message("//===================================================")
message("// ${PROJECT_NAME} build configuration:")
message("// ${PROJECT_NAME} version : ${PROJECT_VERSION}")
message("//===================================================")
message("  CMake version          : ${CMAKE_VERSION}")
if (NOT CMAKE_BUILD_TYPE)
  message("  CMake build type       : NOT SET !")
else()
  message("  CMake build type       : ${CMAKE_BUILD_TYPE}")
endif()
message("  CMake install prefix   : ${CMAKE_INSTALL_PREFIX}")
message("  CMake system processor : ${CMAKE_SYSTEM_PROCESSOR}")
message("  CMake system name (OS) : ${CMAKE_SYSTEM_NAME}")
message("")

message("  Fortran compiler Id      : ${CMAKE_Fortran_COMPILER_ID}")
message("  Fortran compiler version : ${CMAKE_Fortran_COMPILER_VERSION}")
message("  Fortran compiler exe     : ${CMAKE_Fortran_COMPILER}")
message("  Fortran flags            : ${CMAKE_Fortran_FLAGS}")
message("  Fortran compiler wrapper : ${CMAKE_Fortran_COMPILER_WRAPPER}")
message("")

message(STATUS "MPI config:")
message("    MPI found            : ${MPI_FOUND}")
message("    MPI standard version : ${MPI_Fortran_VERSION}")
# note : make sure you correctly your C, C++ and Fortran compiler
# through variable CC, CXX and FC
message("    MPI_VENDOR           : ${MPI_VENDOR}")
message("    MPI library version  : ${MPI_Fortran_LIBRARY_VERSION_STRING}")
message("    MPI fortran compiler : ${MPI_Fortran_COMPILER}")
message("    MPI headers          : ${MPI_Fortran_INCLUDE_DIRS}")
message("    MPI libraries        : ${MPI_Fortran_LIBRARIES}")
#message("    MPI link flags       : ${MPI_Fortran_LINK_FLAGS}")
message("    MPI module dir       : ${MPI_Fortran_MODULE_DIR}")
message("")

message(STATUS "OpenMP config :")
message("    OpenMP found         : ${OpenMP_FOUND}")
message("    OpenMP version       : ${OpenMP_Fortran_VERSION}")
message("")

message(STATUS "BLAS config:")
message("  BLA_VENDOR : ${BLA_VENDOR}")
if (BLAS_FOUND)
  message("  BLAS_LINKER_FLAGS : ${BLAS_LINKER_FLAGS}")
  message("  BLAS_LIBRARIES    : ${BLAS_LIBRARIES}")
else()
  message("  BLAS not found ! Please recheck your environment variables.")
endif()
message("")

message(STATUS "LAPACK config:")
message("  BLA_VENDOR : ${BLA_VENDOR}")
if (LAPACK_FOUND)
  message("  LAPACK_LINKER_FLAGS : ${LAPACK_LINKER_FLAGS}")
  message("  LAPACK_LIBRARIES    : ${LAPACK_LIBRARIES}")
else()
  message("  LAPACK not found ! Please recheck your environment variables.")
endif()
message("")

if (DEFINED ENV{PKG_CONFIG_PATH})
  message(STATUS "ENV{PKG_CONFIG_PATH} was :")
  message("    $ENV{PKG_CONFIG_PATH}")
  message("")
else()
  message(STATUS "ENV{PKG_CONFIG_PATH} was empty")
  message("")
endif()

if (DEFINED ENV{CMAKE_PREFIX_PATH})
  message(STATUS "ENV{CMAKE_PREFIX_PATH} was :")
  message("    $ENV{CMAKE_PREFIX_PATH}")
  message("")
else()
  message(STATUS "ENV{CMAKE_PREFIX_PATH} was empty")
  message("")
endif()

# --- make git ignore build directory
if(NOT EXISTS ${PROJECT_BINARY_DIR}/.gitignore)
  file(WRITE ${PROJECT_BINARY_DIR}/.gitignore "*")
endif()

add_subdirectory(src)
add_subdirectory(utility)

if (BUILD_DOC)
  find_package(LATEX COMPONENTS PDFLATEX BIBTEX)

  if (LATEX_FOUND)
    add_subdirectory(doc)
    add_custom_target(doc DEPENDS doc_tutorial doc_user_guide)
    add_custom_target(doc_clean DEPENDS doc_tutorial_clean doc_user_guide_clean)
    add_custom_target(doc_veryclean DEPENDS doc_tutorial_veryclean doc_user_guide_veryclean)
  else()
    message(FATAL_ERROR "pdflatex not found; can't build documentation. Turn OFF BUILD_DOC or install pdflatex.")
  endif()

endif(BUILD_DOC)

#
# dist targets
#
find_package(Git QUIET)

if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
  set(tarfile ${PROJECT_NAME}-${PROJECT_VERSION}.tar.gz)
  add_custom_command(
    OUTPUT ${tarfile}
    COMMAND ${GIT_EXECUTABLE} archive HEAD --prefix=${PROJECT_NAME}-${PROJECT_VERSION}/ -o ${CMAKE_BINARY_DIR}/${tarfile}
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    )
  add_custom_target(dist DEPENDS "${PROJECT_NAME}-${PROJECT_VERSION}.tar.gz")
endif()

# gnu compatibility,
# see https://cmake.org/cmake/help/latest/module/GNUInstallDirs.html
include(GNUInstallDirs)

# used to install libwannier
set(INSTALL_TARGETS wannier)

################################# EXPORT CONFIG #################################
include(CMakePackageConfigHelpers)

# setup some variables
set(version_config ${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake)
set(project_config_src ${PROJECT_SOURCE_DIR}/${PROJECT_NAME}-config.cmake.in)
set(project_config_dst ${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config.cmake)
set(targets_export_name ${PROJECT_NAME}-targets)

# important variables
set(INSTALL_BINDIR ${CMAKE_INSTALL_BINDIR} CACHE STRING
  "Installation directory for executables, relative to ${CMAKE_INSTALL_PREFIX}.")

set(INSTALL_LIBDIR ${CMAKE_INSTALL_LIBDIR} CACHE STRING
  "Installation directory for libraries, relative to ${CMAKE_INSTALL_PREFIX}.")

set(INSTALL_INCLUDEDIR ${CMAKE_INSTALL_INCLUDEDIR} CACHE STRING
  "Installation directory for include files, relative to ${CMAKE_INSTALL_PREFIX}.")

set(INSTALL_PKGCONFIG_DIR ${CMAKE_INSTALL_LIBDIR}/pkgconfig CACHE PATH
  "Installation directory for pkgconfig (.pc) files, relative to ${CMAKE_INSTALL_PREFIX}.")

set(INSTALL_CMAKE_DIR ${CMAKE_INSTALL_LIBDIR}/cmake CACHE STRING
  "Installation directory for cmake files, relative to ${CMAKE_INSTALL_PREFIX}.")

# Generate the version, config and target files into the build directory.
write_basic_package_version_file(
  ${version_config}
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion)

# Generate cmake my_package-config.cmake file
configure_package_config_file(
  ${project_config_src}
  ${project_config_dst}
  INSTALL_DESTINATION ${INSTALL_CMAKE_DIR})
  # Use a namespace because CMake provides better diagnostics
  # for namespaced imported targets.
export(
  TARGETS ${INSTALL_TARGETS} NAMESPACE wannier::
  FILE ${PROJECT_BINARY_DIR}/${targets_export_name}.cmake)

# macro helper to generate pkg-config file
include(generate_pkgconfig)

# generate the pkg-config file wannier90.pc
generate_pkgconfig(wannier90)

################################# INSTALL LIBRARY #################################

# install executables
install(
  TARGETS wannier90.x postw90.x w90chk2chk.x w90spn2spn.x
  RUNTIME DESTINATION ${INSTALL_BINDIR}
  COMPONENT bin
  )

# install libwannier
install(
  TARGETS ${INSTALL_TARGETS}
  EXPORT ${targets_export_name}
  ARCHIVE DESTINATION ${INSTALL_LIBDIR}
  LIBRARY DESTINATION ${INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

# install cmake config and targets
install(
  FILES ${project_config_dst} ${version_config}
  DESTINATION ${INSTALL_CMAKE_DIR})

install(
  EXPORT ${targets_export_name}
  DESTINATION ${INSTALL_CMAKE_DIR}
  NAMESPACE wannier::)

# install pkgconfig
install(
  FILES ${CMAKE_BINARY_DIR}/${PROJECT_NAME}.pc
  DESTINATION "${INSTALL_PKGCONFIG_DIR}")
