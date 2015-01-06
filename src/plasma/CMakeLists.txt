ADD_LIBRARY(discrover-plasma OBJECT options.cpp measure.cpp align.cpp fasta.cpp
  specification.cpp plasma.cpp mask.cpp count.cpp results.cpp cli.cpp
  score.cpp code.cpp data.cpp motif.cpp correction.cpp
  plasma_stats.cpp randomization_test.cpp dreme/dreme.cpp)

ADD_SUBDIRECTORY(dreme)

ADD_EXECUTABLE(plasma main.cpp "${CMAKE_CURRENT_BINARY_DIR}/../GitSHA1.cpp")
TARGET_LINK_LIBRARIES(plasma discrover)

IF(GOOGLEPERFTOOLS_FOUND)
  TARGET_LINK_LIBRARIES(plasma ${PROFILER_LIBRARY} ${TCMALLOC_LIBRARY})
ENDIF()

INSTALL(TARGETS plasma
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION ${LIBDIR}
  ARCHIVE DESTINATION ${LIBDIR}/static
  DESTINATION "${CMAKE_INSTALL_PREFIX}"
  )

