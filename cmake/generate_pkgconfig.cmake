macro(generate_pkgconfig name)

  # about pkg-config
  # set variables as required by autotools
  set (VERSION ${PROJECT_VERSION})
  set (prefix ${CMAKE_INSTALL_PREFIX})
  set (exec_prefix ${CMAKE_INSTALL_FULL_BINDIR})
  set (libdir ${CMAKE_INSTALL_FULL_LIBDIR})
  set (includedir ${CMAKE_INSTALL_FULL_INCLUDEDIR})

  # actually produce the pkg-config file
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/${name}.pc.in
    ${CMAKE_CURRENT_BINARY_DIR}/${name}.pc
    @ONLY
    )

  unset (VERSION)
  unset (prefix)
  unset (exec_prefix)
  unset (libdir)
  unset (includedir)

endmacro(generate_pkgconfig)
