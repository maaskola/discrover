IF(NOT(PDFLATEX_COMPILER STREQUAL "PDFLATEX_COMPILER-NOTFOUND"))

  include(${CMAKE_BINARY_DIR}/discrover-targets.cmake)

  ADD_CUSTOM_COMMAND(
    OUTPUT plasma-cli-help.txt
    COMMENT "Rebuilding plasma CLI help"
    COMMAND plasma -h > plasma-cli-help.txt && echo Success || echo Failure
    DEPENDS plasma)

  ADD_CUSTOM_COMMAND(
    OUTPUT discrover-cli-help.txt
    COMMENT "Rebuilding discrover CLI help"
    COMMAND discrover-bin -hv > discrover-cli-help.txt && echo Success || echo Failure
    DEPENDS discrover-bin)

  ADD_LATEX_DOCUMENT(discrover-manual.tex
    CONFIGURE discrover-manual.tex
    BIBFILES references.bib
    # IMAGE_DIRS images
    DEPENDS plasma-cli-help.txt discrover-cli-help.txt
    DEFAULT_PDF)

  INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/discrover-manual.pdf
    DESTINATION "${CMAKE_INSTALL_PREFIX}/share/doc/${PROJECT_NAME}")

ENDIF()
