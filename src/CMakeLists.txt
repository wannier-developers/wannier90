#
# common lib
#
add_library(wannier_common OBJECT
  constants.F90
  io.F90
  utility.F90
  types.F90
  hamiltonian.F90
  overlap.F90
  kmesh.F90
  disentangle.F90
  ws_distance.F90
  wannierise.F90
  plot.F90
  transport.F90
  sitesym.F90
  comms.F90
  wannier90_types.F90
  wannier90_readwrite.F90
  readwrite.F90
  error.F90
  error_base.F90)

if (MPI_FOUND)
  target_link_libraries(wannier_common PUBLIC
    MPI::MPI_Fortran)
endif()

#
# post lib
#
add_library(wannier_post OBJECT
  ws_distance.F90
  types.F90
  kmesh.F90
  io.F90
  comms.F90
  utility.F90
  postw90/get_oper.F90
  constants.F90
  postw90/postw90_common.F90
  postw90/wan_ham.F90
  postw90/spin.F90
  postw90/dos.F90
  postw90/berry.F90
  postw90/gyrotropic.F90
  postw90/kpath.F90
  postw90/kslice.F90
  postw90/boltzwann.F90
  postw90/geninterp.F90
  postw90/postw90_types.F90
  postw90/postw90_readwrite.F90
  readwrite.F90
  error.F90
  error_base.F90)
if (MPI_FOUND)
  target_link_libraries(wannier_post PUBLIC
    MPI::MPI_Fortran)
endif()


#
# lib wan2
#
add_library(wan2
  $<TARGET_OBJECTS:wannier_common>)
if (MPI_FOUND)
  target_link_libraries(wan2 PUBLIC
    MPI::MPI_Fortran)
endif()
add_library(wannier::wan2 ALIAS wan2)

#
# libwannier
#
add_library(wannier
  $<TARGET_OBJECTS:wannier_common>
  wannier_lib.F90)
if (MPI_FOUND)
  target_link_libraries(wannier PUBLIC
    MPI::MPI_Fortran)
endif()
add_library(wannier::wannier ALIAS wannier)

#
# move libwannier and libwan2 to top-level bin dir
#
set_target_properties(wannier wan2
  PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

#
# Executables
#

# w90chk2chk.x
add_executable(w90chk2chk.x
  w90chk2chk.F90)

target_link_libraries(w90chk2chk.x PUBLIC
  wannier::wan2
  )

# w90spn2spn.x
add_executable(w90spn2spn.x
  w90spn2spn.F90)

target_link_libraries(w90spn2spn.x PUBLIC
  wannier::wan2
  )

# wannier90.x
add_executable(wannier90.x
  wannier_prog.F90)

target_link_libraries(wannier90.x PUBLIC
  wannier::wan2
  BLAS::BLAS
  LAPACK::LAPACK)

# postw90.x
add_executable(postw90.x
  $<TARGET_OBJECTS:wannier_post>
  postw90/postw90.F90)
if (MPI_FOUND)
  target_link_libraries(postw90.x PUBLIC
    MPI::MPI_Fortran
    BLAS::BLAS
    LAPACK::LAPACK)
endif()

#
# move wannier90.x and postw90.x to top-level bin dir
#
set_target_properties(wannier90.x postw90.x
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
