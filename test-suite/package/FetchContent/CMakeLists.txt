cmake_minimum_required(VERSION 3.25...3.29)
project(example_Wannier90 LANGUAGES Fortran)


include(FetchContent)

if (DEFINED ENV{PACKIT_SOURCE_URL})
	set(git_repo $ENV{PACKIT_SOURCE_URL})
else ()
	set(git_repo https://github.com/wannier-developers/wannier90)
endif ()
if (DEFINED ENV{PACKIT_SOURCE_SHA})
	set(git_tag $ENV{PACKIT_SOURCE_SHA})
elseif (DEFINED ENV{PACKIT_COMMIT_SHA})
	set(git_tag $ENV{PACKIT_COMMIT_SHA})
else ()
	set(git_tag develop)
endif ()

FetchContent_Declare(Wannier90
		GIT_REPOSITORY ${git_repo}
		GIT_TAG ${git_tag}
)
FetchContent_MakeAvailable(Wannier90)

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
