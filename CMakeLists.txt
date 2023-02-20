
cmake_minimum_required(VERSION 2.8.12) # 2.8.8 adds object library

project(wannier90 Fortran)

#set(CMAKE_C_FLAGS_RELEASE "-O2")
#set(CMAKE_C_FLAGS_DEBUG "-g -Wall")
#if(NOT CMAKE_BUILD_TYPE)
#  set(CMAKE_BUILD_TYPE Release)
#endif(NOT CMAKE_BUILD_TYPE)
#message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

include(GNUInstallDirs)

find_package(BLAS)
if (BLAS_FOUND)
    set(LIBS ${LIBS} ${BLAS_LIBRARIES})
else()
    message(WARNING "no beef")
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_MESSAGE_LOG_LEVEL "WARNING")

# Version numbers
file(READ ${PROJECT_SOURCE_DIR}/src/version version_file)
string(REGEX MATCH "WANNIER90_MAJOR_VERSION ([0-9]+)" wannier90_major_version ${version_file})
set(wannier90_major_version ${CMAKE_MATCH_1})
string(REGEX MATCH "WANNIER90_MINOR_VERSION ([0-9]+)" wannier90_minor_version ${version_file})
set(wannier90_minor_version ${CMAKE_MATCH_1})
string(REGEX MATCH "WANNIER90_MICRO_VERSION ([0-9]+)" wannier90_micro_version ${version_file})
set(wannier90_micro_version ${CMAKE_MATCH_1})
set(serial "${wannier90_major_version}.${wannier90_minor_version}.${wannier90_micro_version}")
set(soserial "1")

# common code
set(SRC ${PROJECT_SOURCE_DIR}/src/comms.F90
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
set(SRCW90 ${PROJECT_SOURCE_DIR}/src/disentangle.F90
    ${PROJECT_SOURCE_DIR}/src/hamiltonian.F90
    ${PROJECT_SOURCE_DIR}/src/library_interface.F90
    ${PROJECT_SOURCE_DIR}/src/overlap.F90
    ${PROJECT_SOURCE_DIR}/src/plot.F90
    ${PROJECT_SOURCE_DIR}/src/sitesym.F90
    ${PROJECT_SOURCE_DIR}/src/transport.F90
    ${PROJECT_SOURCE_DIR}/src/wannier90_readwrite.F90
    ${PROJECT_SOURCE_DIR}/src/wannier90_types.F90
    ${PROJECT_SOURCE_DIR}/src/wannierise.F90)

# postw90.x specific
set(SRCP90 ${PROJECT_SOURCE_DIR}/src/postw90/berry.F90
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

# library objects to avoid recompilation of common objects
add_library(obj_common OBJECT ${SRC})
add_library(obj_wann90 OBJECT ${SRCW90} $<TARGET_OBJECTS:obj_common>)
add_library(obj_post90 OBJECT ${SRCP90} $<TARGET_OBJECTS:obj_common>)

# shared library
add_library(dynamic_lib SHARED $<TARGET_OBJECTS:obj_common> $<TARGET_OBJECTS:obj_wann90>)
set_property(TARGET dynamic_lib PROPERTY VERSION ${serial})
set_property(TARGET dynamic_lib PROPERTY SOVERSION ${soserial})
set_property(TARGET dynamic_lib PROPERTY OUTPUT_NAME wannier90)
install(TARGETS dynamic_lib LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

# static library
add_library(static_lib STATIC $<TARGET_OBJECTS:obj_common> $<TARGET_OBJECTS:obj_wann90>)
set_property(TARGET static_lib PROPERTY VERSION ${serial})
set_property(TARGET static_lib PROPERTY SOVERSION ${soserial})
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
