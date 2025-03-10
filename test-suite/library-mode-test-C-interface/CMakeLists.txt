cmake_minimum_required(VERSION 3.25...3.29)
# User projects still need to have Fortran language enabled for the linker
project(demo_Wannier90 LANGUAGES CXX Fortran)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Wannier90 REQUIRED)

add_executable(demo wannier.cc)
target_link_libraries(demo PRIVATE Wannier90::wannier90)
# TODO: user should not need to set this
set_property(TARGET demo PROPERTY LINKER_LANGUAGE Fortran)

enable_testing()
# TODO: The C interface is still buggy and would segfault
#add_test(NAME run_demo
#        COMMAND $<TARGET_FILE:demo>
#)
