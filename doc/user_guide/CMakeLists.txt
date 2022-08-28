add_custom_command(
  OUTPUT user_guide.pdf
  COMMAND ${CMAKE_COMMAND} -E echo "Building user_guide.pdf"
  COMMAND ${PDFLATEX_COMPILER} user_guide.tex
  COMMAND ${BIBTEX_COMPILER}   user_guide
  COMMAND ${PDFLATEX_COMPILER} user_guide.tex
  COMMAND ${BIBTEX_COMPILER}   user_guide
  COMMAND ${PDFLATEX_COMPILER} user_guide.tex
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  DEPENDS user_guide.tex
  COMMENT "Generating user_guide.pdf with pdflatex")

add_custom_target(doc_user_guide DEPENDS user_guide.pdf)

#
# add clean and veryclean target
#
file(GLOB user_guide_clean_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} CONFIGURE_DEPENDS "*.aux" "*.log" "*.bbl" "*.blg" "*.lot" "*.lof" "*.toc" "*.dvi" "*.bak" "*~" "*.ps")
message("[wannier90/doc/user_guide] files to clean : ${user_guide_clean_files}")

add_custom_command(
  OUTPUT clean_user_guide
  COMMAND rm -f ${user_guide_clean_files}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_custom_command(
  OUTPUT veryclean_user_guide
  COMMAND rm -f user_guide.pdf
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_custom_target(doc_user_guide_clean DEPENDS clean_user_guide)
add_custom_target(doc_user_guide_veryclean DEPENDS clean_user_guide veryclean_user_guide)
