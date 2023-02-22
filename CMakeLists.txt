cmake_minimum_required(VERSION 3.12)

project(wannier90 VERSION 4.0.0 LANGUAGES Fortran) # welcome to Wannier90!

include(GNUInstallDirs)

# -Dwith_mpi=off for serial build (default MPI)
# otherwise:
# -Dwith_mpif08=off to disable use of f08 MPI module, to use f90 module instead
# -Dwith_mpih=on to force use of mpif.h include file

option(with_mpi "build MPI version" ON)
if (with_mpi)
  option(with_mpif08 "build MPI version using f08 module" ON)
  option(with_mpih "build MPI version using Fortran 77 mpif.h header" OFF)
  if(with_mpih)
      set(with_mpif08 OFF)
  endif()
endif()

if(with_mpi)
  find_package(MPI REQUIRED)
  include_directories(${MPI_Fortran_INCLUDE_PATH})
  message(STATUS "MPI_Fortran_INCLUDE_PATH ${MPI_Fortran_INCLUDE_PATH}")
  set(LIBS ${LIBS} ${MPI_Fortran_LINK_FLAGS} ${MPI_Fortran_LIBRARIES})
  if(with_mpif08)
      # default to using f08 if available, falling back to f90
      if(MPI_Fortran_HAVE_F08_MODULE)
          add_compile_definitions("MPI08")
      elseif(MPI_Fortran_HAVE_F90_MODULE)
          add_compile_definitions("MPI90")
          message(WARNING "MPI f08 module not found, falling back to f90")
      else()
          message(FATAL_ERROR "MPI f90 module not found and MPI f08 module not found")
      endif()
  elseif(with_mpih)
      # use f77 interface if explicitly requested
      if(MPI_Fortran_HAVE_F77_HEADER)
          add_compile_definitions("MPIH")
      else()
          message(FATAL_ERROR "Fortran 77 mpif.h header requested but not found")
      endif()
  else()
      # if f08 is not requested, then use f90 interface
      if(MPI_Fortran_HAVE_F90_MODULE)
          add_compile_definitions("MPI90")
      else()
          message(FATAL_ERROR "MPI f90 module not found")
      endif()
  endif()
  add_compile_definitions("MPI")
else()
  message(STATUS "(serial build)")
endif()

find_package(BLAS REQUIRED) # binaries require blas
if (BLAS_FOUND)
    set(LIBS ${LIBS} ${BLAS_LIBRARIES})
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# common objects
set(src ${PROJECT_SOURCE_DIR}/src/comms.F90
    ${PROJECT_SOURCE_DIR}/src/constants.F90
    ${PROJECT_SOURCE_DIR}/src/error_base.F90
    ${PROJECT_SOURCE_DIR}/src/error.F90
    ${PROJECT_SOURCE_DIR}/src/io.F90
    ${PROJECT_SOURCE_DIR}/src/kmesh.F90
    ${PROJECT_SOURCE_DIR}/src/readwrite.F90
    ${PROJECT_SOURCE_DIR}/src/types.F90
    ${PROJECT_SOURCE_DIR}/src/utility.F90
    ${PROJECT_SOURCE_DIR}/src/ws_distance.F90)

# wannier90.x specific
set(srcw90 ${PROJECT_SOURCE_DIR}/src/disentangle.F90
    ${PROJECT_SOURCE_DIR}/src/hamiltonian.F90
    ${PROJECT_SOURCE_DIR}/src/overlap.F90
    ${PROJECT_SOURCE_DIR}/src/plot.F90
    ${PROJECT_SOURCE_DIR}/src/sitesym.F90
    ${PROJECT_SOURCE_DIR}/src/transport.F90
    ${PROJECT_SOURCE_DIR}/src/wannier90_readwrite.F90
    ${PROJECT_SOURCE_DIR}/src/wannier90_types.F90
    ${PROJECT_SOURCE_DIR}/src/wannierise.F90)

# postw90.x specific
set(srcp90 ${PROJECT_SOURCE_DIR}/src/postw90/berry.F90
    ${PROJECT_SOURCE_DIR}/src/postw90/boltzwann.F90
    ${PROJECT_SOURCE_DIR}/src/postw90/dos.F90
    ${PROJECT_SOURCE_DIR}/src/postw90/geninterp.F90
    ${PROJECT_SOURCE_DIR}/src/postw90/get_oper.F90
    ${PROJECT_SOURCE_DIR}/src/postw90/gyrotropic.F90
    ${PROJECT_SOURCE_DIR}/src/postw90/kpath.F90
    ${PROJECT_SOURCE_DIR}/src/postw90/kslice.F90
    ${PROJECT_SOURCE_DIR}/src/postw90/postw90_common.F90
    ${PROJECT_SOURCE_DIR}/src/postw90/postw90_readwrite.F90
    ${PROJECT_SOURCE_DIR}/src/postw90/postw90_types.F90
    ${PROJECT_SOURCE_DIR}/src/postw90/spin.F90
    ${PROJECT_SOURCE_DIR}/src/postw90/wan_ham.F90)

# library specific
set(srclib ${PROJECT_SOURCE_DIR}/src/library_interface.F90)

# library objects to avoid recompilation of common objects
add_library(obj_common OBJECT ${src})
add_library(obj_wann90 OBJECT ${srcw90} $<TARGET_OBJECTS:obj_common>)
add_library(obj_post90 OBJECT ${srcp90} $<TARGET_OBJECTS:obj_common>)
add_library(obj_lib OBJECT ${srclib} $<TARGET_OBJECTS:obj_common> $<TARGET_OBJECTS:obj_wann90> )

# shared library
add_library(dynamic_lib SHARED $<TARGET_OBJECTS:obj_lib> $<TARGET_OBJECTS:obj_common> $<TARGET_OBJECTS:obj_wann90>)
set_property(TARGET dynamic_lib PROPERTY VERSION ${PROJECT_VERSION})
#set(soserial "2") # version 2 of library, but wannier90 version 4.0.0, to messy
#set_property(TARGET dynamic_lib PROPERTY SOVERSION ${soserial})
set_property(TARGET dynamic_lib PROPERTY OUTPUT_NAME wannier90)
install(TARGETS dynamic_lib LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

# static library
add_library(static_lib STATIC $<TARGET_OBJECTS:obj_lib> $<TARGET_OBJECTS:obj_common> $<TARGET_OBJECTS:obj_wann90>)
set_property(TARGET static_lib PROPERTY VERSION ${PROJECT_VERSION})
set_property(TARGET static_lib PROPERTY OUTPUT_NAME wannier90)
install(TARGETS static_lib ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

# wannier90.x executable
add_executable(wannier90.x ${PROJECT_SOURCE_DIR}/src/wannier_prog.F90 $<TARGET_OBJECTS:obj_common> $<TARGET_OBJECTS:obj_wann90>)
target_link_libraries(wannier90.x dynamic_lib ${LIBS})
install(TARGETS wannier90.x RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# postw90.x executable
add_executable(postw90.x ${PROJECT_SOURCE_DIR}/src/postw90/postw90.F90 $<TARGET_OBJECTS:obj_common> $<TARGET_OBJECTS:obj_post90>)
target_link_libraries(postw90.x ${LIBS})
install(TARGETS postw90.x RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# include files
install(FILES ${PROJECT_SOURCE_DIR}/src/library_interface.F90 DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/w90_helper_types.mod DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# pkg-config
configure_file(wannier90.pc.cmake ${CMAKE_CURRENT_BINARY_DIR}/wannier90.pc)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/wannier90.pc DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/pkgconfig COMPONENT pkgconfig)

# test system
include(CTest)

# currently, test runner expects these utilities to be in the top level build directory
add_executable(w90chk2chk.x ${PROJECT_SOURCE_DIR}/src/w90chk2chk.F90 $<TARGET_OBJECTS:obj_common>)
target_link_libraries(w90chk2chk.x ${LIBS})

add_executable(w90spn2spn.x ${PROJECT_SOURCE_DIR}/src/w90spn2spn.F90 $<TARGET_OBJECTS:obj_common>)
target_link_libraries(w90spn2spn.x ${LIBS})

# copy entire test bumpf to build directory
#  warning, currently the run_tests script invokes gnumake to unzip some checkpoint files
#  fixme, translate those to explicit targets here (see Makefiles in specific test dirs)
file(COPY ${PROJECT_SOURCE_DIR}/test-suite/checkpoints DESTINATION test-suite)
file(COPY ${PROJECT_SOURCE_DIR}/test-suite/run_tests DESTINATION test-suite)
file(COPY ${PROJECT_SOURCE_DIR}/test-suite/testcode DESTINATION test-suite)
file(COPY ${PROJECT_SOURCE_DIR}/test-suite/tests DESTINATION test-suite)
file(COPY ${PROJECT_SOURCE_DIR}/test-suite/tools DESTINATION test-suite)

# but the run_tests python must be executed in the test-suite directory
add_subdirectory(test-suite)
