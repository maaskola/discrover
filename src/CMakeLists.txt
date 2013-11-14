# add definitions, compiler switches, etc.
ADD_DEFINITIONS(-Wall -funroll-loops -march=native -Wno-unused-local-typedefs)
# ADD_DEFINITIONS(-Wall -funroll-loops -march=native -flto -Wno-unused-local-typedefs)

ADD_SUBDIRECTORY(plasma)
ADD_SUBDIRECTORY(hmm)

CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/GitSHA1.cpp.in"
  "${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp" @ONLY)

CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/stats_config.hpp.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/stats_config.hpp" @ONLY)

CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/dreme_config.hpp.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/dreme_config.hpp" @ONLY)

ADD_EXECUTABLE(rundreme dreme.cpp "${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp")
