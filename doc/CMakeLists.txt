IF(NOT(PDFLATEX_COMPILER STREQUAL "PDFLATEX_COMPILER-NOTFOUND"))

  ADD_LATEX_DOCUMENT(discrover-manual.tex
    CONFIGURE discrover-manual.tex
    BIBFILES references.bib
    # IMAGE_DIRS images
    DEPENDS ${CMAKE_BINARY_DIR}/plasma-cli-help.txt ${CMAKE_BINARY_DIR}/discrover-cli-help.txt
    DEFAULT_PDF)

  INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/discrover-manual.pdf
    DESTINATION "${CMAKE_INSTALL_PREFIX}/share/doc/${PROJECT_NAME}")

ENDIF()
