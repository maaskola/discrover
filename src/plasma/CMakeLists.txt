ADD_LIBRARY(discrover-plasma OBJECT options.cpp measure.cpp align.cpp fasta.cpp
  specification.cpp plasma.cpp mask.cpp count.cpp results.cpp cli.cpp
  score.cpp code.cpp data.cpp motif.cpp correction.cpp
  plasma_stats.cpp dreme/dreme.cpp)

SET_TARGET_PROPERTIES(discrover-plasma PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

# un-comment to build a test program for the DREME driver code
# ADD_SUBDIRECTORY(dreme)

ADD_EXECUTABLE(plasma main.cpp)
SET_TARGET_PROPERTIES(plasma PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
TARGET_LINK_LIBRARIES(plasma discrover)

IF(GOOGLEPERFTOOLS_FOUND)
  TARGET_LINK_LIBRARIES(plasma ${PROFILER_LIBRARY} ${TCMALLOC_LIBRARY})
ENDIF()

INSTALL(TARGETS plasma DESTINATION bin)
