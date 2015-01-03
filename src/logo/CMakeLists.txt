# build a shared library
ADD_LIBRARY(discroverlogo SHARED logo.cpp options.cpp cli.cpp ../aux.cpp "${CMAKE_CURRENT_BINARY_DIR}/../GitSHA1.cpp")

ADD_EXECUTABLE(discrover-logo main.cpp ../sha1.cpp ../terminal.cpp
  ../plasma/code.cpp ../executioninformation.cpp)

TARGET_LINK_LIBRARIES(discroverlogo
  ${CAIRO_LIBRARIES}
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_IOSTREAMS_LIBRARY})

TARGET_LINK_LIBRARIES(discrover-logo discroverlogo ${Boost_PROGRAM_OPTIONS_LIBRARY})

INSTALL(TARGETS discrover-logo discroverlogo
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION ${LIBDIR}
  ARCHIVE DESTINATION ${LIBDIR}/static
  DESTINATION "${CMAKE_INSTALL_PREFIX}"
  )
