# add definitions, compiler switches, etc.
ADD_DEFINITIONS(-Wall -funroll-loops -march=native -Wno-unused-local-typedefs)

# Allow to #include <discrover_config.hpp>, if doing out-of-source build.
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

ADD_SUBDIRECTORY(shuffle)
ADD_SUBDIRECTORY(plasma)
ADD_SUBDIRECTORY(hmm)

ADD_LIBRARY(discrover-common OBJECT aux.cpp executioninformation.cpp matrix.cpp
  random_distributions.cpp random_seed.cpp sha1.cpp terminal.cpp timer.cpp
  mcmc/montecarlo.cpp)

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

CONFIGURE_FILE(discrover_config.hpp.in discrover_config.hpp)
CONFIGURE_FILE(git_config.hpp.in git_config.hpp)

INSTALL(TARGETS discrover DESTINATION ${LIBDIR})
