ADD_LIBRARY(discrover-shuffle OBJECT dinucleotide_shuffle.cpp)

ADD_EXECUTABLE(discrover-shuffle-bin main.cpp)
SET_TARGET_PROPERTIES(discrover-shuffle-bin
  PROPERTIES OUTPUT_NAME discrover-shuffle)

IF(COMPILER_SUPPORTS_PIC)
  SET_TARGET_PROPERTIES(discrover-shuffle PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
ENDIF()
IF(COMPILER_SUPPORTS_PIE)
  SET_TARGET_PROPERTIES(discrover-shuffle-bin
    PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
ENDIF()

TARGET_LINK_LIBRARIES(discrover-shuffle-bin discrover)

INSTALL(TARGETS discrover-shuffle-bin DESTINATION bin)
