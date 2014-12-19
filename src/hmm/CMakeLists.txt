# build a shared library
ADD_LIBRARY(discrover SHARED association.cpp ../aux.cpp basedefs.cpp
  analysis.cpp bitmask.cpp
  conditional_mutual_information.cpp
  conditional_decoder.cpp
  ../executioninformation.cpp ../random_seed.cpp
  hmm.cpp hmm_core.cpp hmm_init.cpp hmm_learn.cpp hmm_linesearch.cpp
  hmm_aux.cpp hmm_gradient.cpp hmm_mcmc.cpp hmm_score.cpp hmm_options.cpp 
  ../mcmc/montecarlo.cpp polyfit.cpp registration.cpp report.cpp results.cpp
  sequence.cpp ../sha1.cpp subhmm.cpp trainingmode.cpp ../terminal.cpp
  ../timer.cpp)

ADD_EXECUTABLE(discrover-bin main.cpp
  "${CMAKE_CURRENT_BINARY_DIR}/../GitSHA1.cpp")

# ADD_EXECUTABLE(mcmc mcmc/montecarlo.cpp) # this is a test program for the Gibbs sampling code
# ADD_EXECUTABLE(polyfit polyfittest.cpp polyfit.cpp "${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp")

SET_TARGET_PROPERTIES(discrover-bin
  PROPERTIES OUTPUT_NAME discrover)

TARGET_LINK_LIBRARIES(discrover
  discroverplasma
  -L${Boost_LIBRARY_DIRS}
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_REGEX_LIBRARY}
  ${Boost_IOSTREAMS_LIBRARY})

IF(LIBR_FOUND AND LIBR_MATHLIB_LIBRARY)
  TARGET_LINK_LIBRARIES(discrover m ${LIBR_LIBRARIES})
ELSE()
  SUBDIRS(../stats)
  TARGET_LINK_LIBRARIES(discrover discroverstats)
ENDIF()

TARGET_LINK_LIBRARIES(discrover-bin
  discrover
  -L${Boost_LIBRARY_DIRS}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${Boost_DATE_TIME_LIBRARY})

IF(GOOGLEPERFTOOLS_FOUND)
  TARGET_LINK_LIBRARIES(discrover-bin ${PROFILER_LIBRARY} ${TCMALLOC_LIBRARY})
ENDIF()

# TARGET_LINK_LIBRARIES(polyfit -L${Boost_LIBRARY_DIRS} ${Boost_PROGRAM_OPTIONS_LIBRARY})

INSTALL(TARGETS discrover-bin discrover
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib/static
    DESTINATION "${CMAKE_INSTALL_PREFIX}"
    )
