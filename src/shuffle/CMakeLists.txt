ADD_LIBRARY(discrover-shuffle OBJECT dinucleotide_shuffle.cpp)

SET_TARGET_PROPERTIES(discrover-shuffle PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

ADD_EXECUTABLE(discrover-shuffle-bin main.cpp)
SET_TARGET_PROPERTIES(discrover-shuffle-bin
  PROPERTIES OUTPUT_NAME discrover-shuffle)
SET_TARGET_PROPERTIES(discrover-shuffle-bin
  PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
TARGET_LINK_LIBRARIES(discrover-shuffle-bin discrover)

INSTALL(TARGETS discrover-shuffle-bin DESTINATION bin)
