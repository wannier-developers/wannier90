add_custom_command(
  OUTPUT tutorial.pdf
  COMMAND ${CMAKE_COMMAND} -E echo "Building tutorial.pdf"
  COMMAND ${PDFLATEX_COMPILER} tutorial.tex
  COMMAND ${BIBTEX_COMPILER}   tutorial
  COMMAND ${PDFLATEX_COMPILER} tutorial.tex
  COMMAND ${BIBTEX_COMPILER}   tutorial
  COMMAND ${PDFLATEX_COMPILER} tutorial.tex
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  DEPENDS tutorial.tex
  COMMENT "Generating tutorial.pdf with pdflatex")

add_custom_target(doc_tutorial DEPENDS tutorial.pdf)

#
# add clean and veryclean target
#
file(GLOB tutorial_clean_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} CONFIGURE_DEPENDS "*.aux" "*.log" "*.bbl" "*.blg" "*.lot" "*.lof" "*.toc" "*.dvi" "*.bak" "*~" "*.ps")
message("[wannier90/doc/tutorial] files to clean : ${tutorial_clean_files}")

add_custom_command(
  OUTPUT clean_tutorial
  COMMAND rm -f ${tutorial_clean_files}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_custom_command(
  OUTPUT veryclean_tutorial
  COMMAND rm -f tutorial.pdf
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_custom_target(doc_tutorial_clean DEPENDS clean_tutorial)
add_custom_target(doc_tutorial_veryclean DEPENDS clean_tutorial veryclean_tutorial)
