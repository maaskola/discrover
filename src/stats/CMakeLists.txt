
# CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/config.h.cmake ${CMAKE_BINARY_DIR}/src/config.h)

#add definitions, compiler switches, etc.
ADD_DEFINITIONS(-Wall -I ${Boost_INCLUDE_DIR} -funroll-loops -mtune=native)

#build a shared library
ADD_LIBRARY(discroverstats SHARED pgamma.cpp chisq.cpp)

TARGET_LINK_LIBRARIES(discroverstats m)

INSTALL(TARGETS discroverstats
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib/static
    DESTINATION "${CMAKE_INSTALL_PREFIX}"
    )
