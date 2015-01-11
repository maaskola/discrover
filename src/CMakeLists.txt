# add definitions, compiler switches, etc.
ADD_DEFINITIONS(-Wall -funroll-loops -march=native -Wno-unused-local-typedefs)

ADD_SUBDIRECTORY(shuffle)
ADD_SUBDIRECTORY(plasma)
ADD_SUBDIRECTORY(hmm)

# The directory mcmc contains a demo application that can be built by
# uncommenting the next line
# ADD_SUBDIRECTORY(mcmc)

ADD_LIBRARY(discrover-common OBJECT aux.cpp executioninformation.cpp matrix.cpp
  random_seed.cpp sha1.cpp terminal.cpp timer.cpp mcmc/montecarlo.cpp)

LIST(APPEND DISCROVER_OBJECTS
  $<TARGET_OBJECTS:discrover-common>
  $<TARGET_OBJECTS:discrover-hmm>
  $<TARGET_OBJECTS:discrover-plasma>
  $<TARGET_OBJECTS:discrover-shuffle>)

IF(CAIRO_FOUND)
  ADD_SUBDIRECTORY(logo)
  LIST(APPEND DISCROVER_OBJECTS $<TARGET_OBJECTS:discrover-logo>)
ENDIF()

IF(NOT(LIBR_FOUND AND LIBR_MATHLIB_LIBRARY))
  ADD_SUBDIRECTORY(stats)
  LIST(APPEND DISCROVER_OBJECTS $<TARGET_OBJECTS:discrover-stats>)
ENDIF()

ADD_LIBRARY(discrover SHARED ${DISCROVER_OBJECTS})

SET_TARGET_PROPERTIES(discrover
  PROPERTIES VERSION ${DISCROVER_VERSION}
  SONAME ${DISCROVER_VERSION})

TARGET_LINK_LIBRARIES(discrover
  m
  -L${Boost_LIBRARY_DIRS}
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_REGEX_LIBRARY}
  ${Boost_IOSTREAMS_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${Boost_DATE_TIME_LIBRARY})

if(CAIRO_FOUND)
  TARGET_LINK_LIBRARIES(discrover ${CAIRO_LIBRARIES})
ENDIF()

IF(LIBR_FOUND AND LIBR_MATHLIB_LIBRARY)
  TARGET_LINK_LIBRARIES(discrover ${LIBR_LIBRARIES})
ENDIF()

CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/GitSHA1.cpp.in"
  "${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp" @ONLY)

CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/stats_config.hpp.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/stats_config.hpp" @ONLY)

CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/logo_config.hpp.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/logo_config.hpp" @ONLY)

CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/dreme_config.hpp.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/dreme_config.hpp" @ONLY)

INSTALL(TARGETS discrover
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION ${LIBDIR}
    ARCHIVE DESTINATION ${LIBDIR}/static
    DESTINATION "${CMAKE_INSTALL_PREFIX}"
    )
