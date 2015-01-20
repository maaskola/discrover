
SET(MANUAL_LOCATION
  "${CMAKE_INSTALL_PREFIX}/share/doc/discrover/discrover-manual.pdf")

CONFIGURE_FILE(discrover.1.in discrover.1)
CONFIGURE_FILE(discrover-logo.1.in discrover-logo.1)
CONFIGURE_FILE(discrover-shuffle.1.in discrover-shuffle.1)
CONFIGURE_FILE(plasma.1.in plasma.1)

INSTALL(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/discrover.1
  ${CMAKE_CURRENT_BINARY_DIR}/discrover-logo.1
  ${CMAKE_CURRENT_BINARY_DIR}/discrover-shuffle.1
  ${CMAKE_CURRENT_BINARY_DIR}/plasma.1
DESTINATION "${CMAKE_INSTALL_PREFIX}/share/man/man1/")
