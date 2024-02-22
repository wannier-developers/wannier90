cmake_minimum_required(VERSION 3.25...3.29)
project(example_Wannier90 LANGUAGES Fortran)

find_package(Wannier90 REQUIRED)

add_executable(main src/main.F90)
target_link_libraries(main PRIVATE Wannier90::wannier90)

enable_testing()

add_test(NAME smoke_test_lib
		COMMAND $<TARGET_FILE:main>
)
add_test(NAME smoke_test_exe
		COMMAND $<TARGET_FILE:Wannier90::exe> --version
)
add_test(NAME smoke_test_post
		COMMAND $<TARGET_FILE:Wannier90::post> --version
)
