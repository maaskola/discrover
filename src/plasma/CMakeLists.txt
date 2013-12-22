
# CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/config.h.cmake ${CMAKE_BINARY_DIR}/src/config.h)

# build a shared library
ADD_LIBRARY(plasma SHARED options.cpp measure.cpp align.cpp fasta.cpp specification.cpp find.cpp
  mask.cpp count.cpp results.cpp plasma_cli.cpp score.cpp code.cpp data.cpp motif.cpp
  ../aux.cpp correction.cpp plasma.cpp plasma_stats.cpp ../random_seed.cpp randomization_test.cpp
  dreme/dreme.cpp ../sha1.cpp ../timer.cpp ../shuffle/dinucleotide_shuffle.cpp
  ../executioninformation.cpp ../mcmc/montecarlo.cpp
  "${CMAKE_CURRENT_BINARY_DIR}/../GitSHA1.cpp")

ADD_SUBDIRECTORY(dreme)

ADD_EXECUTABLE(plasma-bin main.cpp)
SET_TARGET_PROPERTIES(plasma-bin
  PROPERTIES OUTPUT_NAME plasma)

TARGET_LINK_LIBRARIES(plasma
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_IOSTREAMS_LIBRARY})
IF(LIBR_FOUND AND LIBR_MATHLIB_LIBRARY)
  TARGET_LINK_LIBRARIES(plasma m ${LIBR_LIBRARIES})
ELSE()
  TARGET_LINK_LIBRARIES(plasma stats)
ENDIF()

TARGET_LINK_LIBRARIES(plasma-bin plasma ${Boost_PROGRAM_OPTIONS_LIBRARY})
IF(GOOGLEPERFTOOLS_FOUND)
  TARGET_LINK_LIBRARIES(plasma-bin ${PROFILER_LIBRARY} ${TCMALLOC_LIBRARY})
ENDIF()

INSTALL(TARGETS plasma-bin plasma
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static
  DESTINATION "${CMAKE_INSTALL_PREFIX}"
  )

