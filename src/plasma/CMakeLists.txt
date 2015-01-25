ADD_LIBRARY(discrover-plasma OBJECT align.cpp cli.cpp code.cpp correction.cpp
  count.cpp data.cpp fasta.cpp harmonization.cpp mask.cpp measure.cpp
  motif.cpp options.cpp plasma.cpp plasma_stats.cpp results.cpp score.cpp
  specification.cpp dreme/dreme.cpp)

# un-comment to build a test program for the DREME driver code
# ADD_SUBDIRECTORY(dreme)

ADD_EXECUTABLE(plasma main.cpp)
TARGET_LINK_LIBRARIES(plasma discrover)

IF(GOOGLEPERFTOOLS_FOUND)
  TARGET_LINK_LIBRARIES(plasma ${PROFILER_LIBRARY} ${TCMALLOC_LIBRARY})
ENDIF()

INSTALL(TARGETS plasma DESTINATION bin)
