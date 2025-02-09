cmake_minimum_required(VERSION 3.25...3.29)
project(demo_Wannier90 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Wannier90 REQUIRED)

add_executable(demo wannier.cc)
target_link_libraries(demo PRIVATE Wannier90::wannier90)

enable_testing()
add_test(NAME run_demo
        COMMAND $<TARGET_FILE:demo>
)
