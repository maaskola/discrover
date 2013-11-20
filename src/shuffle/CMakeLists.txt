# build a shared library
ADD_LIBRARY(shuffle SHARED dinucleotide_shuffle.cpp "${CMAKE_CURRENT_BINARY_DIR}/../GitSHA1.cpp")

ADD_EXECUTABLE(dinucleotide_shuffle main.cpp ../plasma/fasta.cpp)

TARGET_LINK_LIBRARIES(dinucleotide_shuffle
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_IOSTREAMS_LIBRARY})

TARGET_LINK_LIBRARIES(dinucleotide_shuffle shuffle ${Boost_PROGRAM_OPTIONS_LIBRARY})

INSTALL(TARGETS dinucleotide_shuffle shuffle
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static
  DESTINATION "${CMAKE_INSTALL_PREFIX}"
  )
