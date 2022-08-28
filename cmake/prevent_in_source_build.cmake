function(prevent_in_source_build)

  # make sure the user doesn't play dirty with symlinks
  get_filename_component(srcdir "${CMAKE_SOURCE_DIR}" REALPATH)
  get_filename_component(bindir "${CMAKE_BINARY_DIR}" REALPATH)

  # disallow in-source builds
  if("${srcdir}" STREQUAL "${bindir}")
    message("######################################################")
    message("# ${PROJECT_NAME} should not be configured and built in the source directory")
    message("# You must run cmake in a build directory.")
    message("# For example:")
    message("# mkdir build ; cd build")
    message("# run cmake from the build directory.")
    message("######################################################")
    message(FATAL_ERROR "Quitting configuration")
  endif()

endfunction(prevent_in_source_build)
