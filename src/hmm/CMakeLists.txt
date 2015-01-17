ADD_LIBRARY(discrover-hmm OBJECT association.cpp basedefs.cpp
  analysis.cpp bitmask.cpp
  conditional_mutual_information.cpp
  conditional_decoder.cpp
  hmm.cpp hmm_core.cpp hmm_init.cpp hmm_learn.cpp hmm_linesearch.cpp
  hmm_aux.cpp hmm_gradient.cpp hmm_mcmc.cpp hmm_score.cpp hmm_options.cpp
  polyfit.cpp registration.cpp report.cpp results.cpp
  sequence.cpp subhmm.cpp trainingmode.cpp )

ADD_EXECUTABLE(discrover-bin main.cpp)
SET_TARGET_PROPERTIES(discrover-bin PROPERTIES OUTPUT_NAME discrover)
TARGET_LINK_LIBRARIES(discrover-bin discrover)

# ADD_EXECUTABLE(mcmc mcmc/montecarlo.cpp) # this is a test program for the Gibbs sampling code
# ADD_EXECUTABLE(polyfit polyfittest.cpp polyfit.cpp "${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp")

IF(GOOGLEPERFTOOLS_FOUND)
  TARGET_LINK_LIBRARIES(discrover-bin ${PROFILER_LIBRARY} ${TCMALLOC_LIBRARY})
ENDIF()

INSTALL(TARGETS discrover-bin DESTINATION bin)
