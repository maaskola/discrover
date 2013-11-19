

ADD_EXECUTABLE(test_dreme test_dreme.cpp dreme.cpp
  "${CMAKE_CURRENT_BINARY_DIR}/../../GitSHA1.cpp")

TARGET_LINK_LIBRARIES(test_dreme
  ${Boost_FILESYSTEM_LIBRARY})

