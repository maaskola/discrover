# build a shared library
ADD_LIBRARY(discrovershuffle SHARED dinucleotide_shuffle.cpp "${CMAKE_CURRENT_BINARY_DIR}/../GitSHA1.cpp")

ADD_EXECUTABLE(dinucleotide_shuffle main.cpp ../plasma/fasta.cpp)

SET_TARGET_PROPERTIES(dinucleotide_shuffle
  PROPERTIES OUTPUT_NAME discrover-shuffle)

TARGET_LINK_LIBRARIES(dinucleotide_shuffle
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_IOSTREAMS_LIBRARY})

TARGET_LINK_LIBRARIES(dinucleotide_shuffle discrovershuffle ${Boost_PROGRAM_OPTIONS_LIBRARY})

INSTALL(TARGETS dinucleotide_shuffle discrovershuffle
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION ${LIBDIR}
  ARCHIVE DESTINATION ${LIBDIR}/static
  DESTINATION "${CMAKE_INSTALL_PREFIX}"
  )
