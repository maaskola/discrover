# add definitions, compiler switches, etc.
ADD_DEFINITIONS(-Wall -funroll-loops -march=native -Wno-unused-local-typedefs)
# ADD_DEFINITIONS(-Wall -funroll-loops -march=native -flto -Wno-unused-local-typedefs)

ADD_SUBDIRECTORY(shuffle)
ADD_SUBDIRECTORY(plasma)
ADD_SUBDIRECTORY(hmm)
if(CAIRO_FOUND)
  ADD_SUBDIRECTORY(logo)
ENDIF()

# The directory mcmc contains a demo application that can be built by
# uncommenting the next line
# ADD_SUBDIRECTORY(mcmc)

CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/GitSHA1.cpp.in"
  "${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp" @ONLY)

CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/stats_config.hpp.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/stats_config.hpp" @ONLY)

CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/dreme_config.hpp.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/dreme_config.hpp" @ONLY)
