if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")

  if (TRAVISCI_ENABLED)
    add_compile_options(
      "$<$<COMPILE_LANGUAGE:Fortran>:-fprofile-arcs;-ftest-coverage;-fstrict-aliasing;-fno-omit-frame-pointer;-fno-realloc-lhs;-fcheck=bounds,do,recursion,pointer;-ffree-form;-Wall;-Waliasing;-Wsurprising;-Wline-truncation;-Wno-tabs;-Wno-uninitialized;-Wno-unused-dummy-argument;-Wno-unused;-Wno-character-truncation;-O1;-g;-fbacktrace>")

    add_link_options(
      "$<$<COMPILE_LANGUAGE:Fortran>:-fprofile-arcs;-ftest-coverage;-fstrict-aliasing;-fno-omit-frame-pointer;-fno-realloc-lhs;-fcheck=bounds,do,recursion,pointer;-ffree-form;-Wall;-Waliasing;-Wsurprising;-Wline-truncation;-Wno-tabs;-Wno-uninitialized;-Wno-unused-dummy-argument;-Wno-unused;-Wno-character-truncation;-O1;-g;-fbacktrace>")
  endif()

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "NAG")

  add_compile_definitions(NAG)

endif()
