

cmake_minimum_required(VERSION 2.4.8)

IF(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
ENDIF(COMMAND cmake_policy)

PROJECT(dlhmm)

# if you want the full compiler output, uncomment the following line
# SET(CMAKE_VERBOSE_MAKEFILE ON)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")




# Change build configuration
#
# SET(CMAKE_BUILD_TYPE Debug)
SET(CMAKE_BUILD_TYPE Release)

# Set installation destination and prefix path to find libraries installed on
# this sytem. This is the only required manual configuration step.
#
# NOTE: Please modify this to suit your system
# ============================================
#
SET(LOCAL_PREFIX "~/local")
SET(CMAKE_INSTALL_PREFIX ${LOCAL_PREFIX})
SET(CMAKE_PREFIX_PATH ${LOCAL_PREFIX})




SET(CMAKE_LINKER ld.gold)

INCLUDE(InstallRequiredSystemLibraries)

SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Binding site pattern discovery from nucleic acid sequences by discriminative learning of hidden Markov models")
SET(CPACK_PACKAGE_VENDOR "Jonas Maaskola")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
SET(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README")
SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "0")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "CMake ${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")
SET(CPACK_SOURCE_GENERATOR TGZ;TBZ2;ZIP)
SET(CPACK_GENERATOR TGZ;TBZ2;ZIP)
SET(CPACK_PACKAGE_CONTACT "Jonas Maaskola")
# SET(CPACK_STRIP_FILES "bin/MyExecutable")
# SET(CPACK_SOURCE_STRIP_FILES "")
# SET(CPACK_PACKAGE_EXECUTABLES "MyExecutable" "My Executable")
INCLUDE(CPack)

SET(DLHMM_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")


# Check for C++11 features of the compiler
#
INCLUDE(CheckCXX11Features)

# Boost configuration
# You may try uncommenting some of the commented options below if you have
# trouble automatically finding Boost on your system.

# Toggle verbose output of the boost configuration
#
SET(Boost_DEBUG OFF)

# If you have your Boost installation in a non-standard location, that is not
# part below ${LOCAL_PREFIX} and the automatic finding fails, you may set the
# Boost root directory here
#
# SET(BOOST_ROOT "~/local/boost-1.49.0")

# Uncomment if linking should not use system Boost installation but instead
# should prefer one in a non-standard location.
#
# SET(Boost_NO_SYSTEM_PATHS TRUE)

# SET(Boost_USE_STATIC_LIBS  OFF)
# SET(Boost_USE_STATIC_RUNTIME ON)
# SET(Boost_ADDITIONAL_VERSIONS "1.47" "1.47.0" "1.46" "1.46.0" "1.42.0" "1.40.0")
FIND_PACKAGE(Boost COMPONENTS iostreams regex system filesystem
  unit_test_framework program_options serialization date_time REQUIRED)
IF(Boost_FOUND)
  ADD_DEFINITIONS( "-DHAS_BOOST" )
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
  LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})
ENDIF()

FIND_PACKAGE(OpenMP)

IF(NOT(DEFINED WITH_TCMALLOC) OR WITH_TCMALLOC)
  FIND_PACKAGE(LibGooglePerfTools)
  IF(NOT(GOOGLEPERFTOOLS_FOUND))
    IF(WITH_TCMALLOC)
      MESSAGE(FATAL_ERROR " TCMALLOC required but not found.")
    ELSE()
      MESSAGE(INFO " TCMALLOC not found.")
    ENDIF()
  ENDIF()
ENDIF()

IF(NOT(DEFINED WITH_RMATHLIB) OR WITH_RMATHLIB)
  FIND_PACKAGE(LibR)
  IF(LIBR_FOUND)
    INCLUDE_DIRECTORIES(${LIBR_INCLUDE_DIRS})
    SET(HAVE_LIBR 1)
  ELSE()
    IF(WITH_RMATHLIB)
      MESSAGE(FATAL_ERROR " standalone Rmathlib required but not found.")
    ELSE()
      SET(HAVE_LIBR 0)
    ENDIF()
  ENDIF()
ENDIF()

IF(NOT(DEFINED WITH_DREME) OR WITH_DREME)
  FIND_PACKAGE(Dreme)
  IF(NOT DREME_FOUND)
    IF(WITH_DREME)
      MESSAGE(FATAL_ERROR " Dreme required but not found")
    ELSE()
      MESSAGE(INFO " Dreme was not found")
    ENDIF()
  ENDIF()
ENDIF()

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -std=c++11")

INCLUDE(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
git_describe(GIT_DESCRIPTION)
git_branch(GIT_BRANCH)

ADD_CUSTOM_TARGET(ctags ALL
    COMMAND ctags -R --c++-kinds=+p --fields=+iaS --extra=+q src/
    COMMENT "Creating Ctags."
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} )

ADD_SUBDIRECTORY(src build)
ADD_SUBDIRECTORY(scripts)

IF(NOT(DEFINED WITH_DOC) OR WITH_DOC)
  ADD_SUBDIRECTORY(doc build/doc)
ENDIF()

INSTALL(FILES ${CMAKE_CURRENT_SOURCE_DIR}/README
  ${CMAKE_CURRENT_SOURCE_DIR}/COPYING
  ${CMAKE_CURRENT_SOURCE_DIR}/FAQ
  ${CMAKE_CURRENT_SOURCE_DIR}/TODO
  ${CMAKE_CURRENT_SOURCE_DIR}/AUTHORS
  DESTINATION "${CMAKE_INSTALL_PREFIX}/share/doc/${PROJECT_NAME}"
  )
