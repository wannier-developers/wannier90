cmake_minimum_required(VERSION 3.25...3.29)
project(example_Wannier90 LANGUAGES Fortran)

find_package(PkgConfig REQUIRED)
pkg_check_modules(Wannier90 REQUIRED IMPORTED_TARGET wannier90)

add_executable(main src/main.F90)
target_link_libraries(main PRIVATE PkgConfig::Wannier90)

enable_testing()

add_test(NAME smoke_test_lib
		COMMAND $<TARGET_FILE:main>
)
