ADD_LIBRARY(discrover-shuffle OBJECT dinucleotide_shuffle.cpp)

ADD_EXECUTABLE(discrover-shuffle-bin main.cpp)
SET_TARGET_PROPERTIES(discrover-shuffle-bin PROPERTIES OUTPUT_NAME discrover-shuffle)
TARGET_LINK_LIBRARIES(discrover-shuffle-bin discrover)

INSTALL(TARGETS discrover-shuffle-bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION ${LIBDIR}
  ARCHIVE DESTINATION ${LIBDIR}/static
  DESTINATION "${CMAKE_INSTALL_PREFIX}"
  )
