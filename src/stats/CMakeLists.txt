#build a shared library
ADD_LIBRARY(discroverstats SHARED pgamma.cpp chisq.cpp)

TARGET_LINK_LIBRARIES(discroverstats m)

INSTALL(TARGETS discroverstats
    LIBRARY DESTINATION ${LIBDIR}
    ARCHIVE DESTINATION ${LIBDIR}/static
    DESTINATION "${CMAKE_INSTALL_PREFIX}"
    )
