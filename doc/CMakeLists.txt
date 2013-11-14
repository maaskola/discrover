INCLUDE(UseLATEX)

ADD_LATEX_DOCUMENT(documentation.tex
  CONFIGURE documentation.tex
  BIBFILES documentation.bib
  # IMAGE_DIRS images
  DEFAULT_PDF
)
ADD_CUSTOM_COMMAND(TARGET pdf POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy
  ${CMAKE_CURRENT_BINARY_DIR}/documentation.pdf
  ${CMAKE_BINARY_DIR}/doc/documentation.pdf)
