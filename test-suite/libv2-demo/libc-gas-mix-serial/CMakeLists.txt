cmake_minimum_required(VERSION 3.25...3.29)
project(demo_Wannier90 LANGUAGES CXX)

find_package(Wannier90 REQUIRED)

add_executable(demo demo.cc wannier.cc)
target_link_libraries(demo PRIVATE Wannier90::wannier90)
