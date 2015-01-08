

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.8)

IF(COMMAND cmake_policy)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND cmake_policy)

PROJECT(discrover)

# if you want the full compiler output, uncomment the following line
# SET(CMAKE_VERBOSE_MAKEFILE ON)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")




# Change build configuration
#
# SET(CMAKE_BUILD_TYPE Debug)
SET(CMAKE_BUILD_TYPE Release)

# Set installation destination and prefix path to find libraries installed on
# this sytem. If this is not specified, /usr/local/ will be used.
#
# SET(LOCAL_PREFIX "/usr/local")
# SET(CMAKE_INSTALL_PREFIX ${LOCAL_PREFIX})
# SET(CMAKE_PREFIX_PATH ${LOCAL_PREFIX})




SET(CMAKE_LINKER ld.gold)

INCLUDE(InstallRequiredSystemLibraries)

SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Binding site pattern discovery from nucleic acid sequences by discriminative learning of hidden Markov models")
SET(CPACK_PACKAGE_VENDOR "Jonas Maaskola")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
SET(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "4")
SET(CPACK_PACKAGE_VERSION_PATCH "1")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "CMake ${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")
SET(CPACK_SOURCE_GENERATOR TGZ;TBZ2;ZIP)
SET(CPACK_GENERATOR TGZ;TBZ2;ZIP)
SET(CPACK_PACKAGE_CONTACT "Jonas Maaskola")
# SET(CPACK_STRIP_FILES "bin/MyExecutable")
# SET(CPACK_SOURCE_STRIP_FILES "")
# SET(CPACK_PACKAGE_EXECUTABLES "MyExecutable" "My Executable")
INCLUDE(CPack)

SET(DISCROVER_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")


# Check for C++11 features of the compiler
INCLUDE(CheckCXX11Features)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX11_FLAGS}")

# Check whether the compiler supports link-time optimization

INCLUDE(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG(-flto COMPILER_SUPPORTS_LTO)
IF(COMPILER_SUPPORTS_LTO)
  MESSAGE(STATUS "Enabled: link-time optimization")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")
ELSEIF(COMPILER_SUPPORTS_CXX0X)
  MESSAGE(STATUS "Disabled: link-time optimization - no compiler support")
ENDIF()


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
FIND_PACKAGE(Ctags)
FIND_PACKAGE(Ruby)

IF(NOT(DEFINED WITH_CAIRO) OR WITH_CAIRO)
  FIND_PACKAGE(Cairo)
  IF(CAIRO_FOUND)
    INCLUDE_DIRECTORIES(${CAIRO_INCLUDE_DIRS})
    MESSAGE(STATUS "Enabled: logo plotting capabilities using Cairo")
  ELSE()
    IF(WITH_CAIRO)
      MESSAGE(FATAL_ERROR "Cairo required but not found")
    ELSE()
      MESSAGE(STATUS "Disabled: logo plotting capabilities - Cairo not found")
    ENDIF()
  ENDIF()
ELSE()
  MESSAGE(STATUS "Disabled: logo plotting capabilities (use -DWITH_CAIRO=ON to enable)")
ENDIF()

INCLUDE(UseLATEX)
IF(NOT(PDFLATEX_COMPILER STREQUAL "PDFLATEX_COMPILER-NOTFOUND"))
  MESSAGE(STATUS "Found pdflatex: ${PDFLATEX_COMPILER}")
ENDIF()

IF(NOT(DEFINED WITH_TCMALLOC) OR WITH_TCMALLOC)
  FIND_PACKAGE(LibGooglePerfTools)
  IF(GOOGLEPERFTOOLS_FOUND)
    MESSAGE(STATUS "Enabled: usage of TCMALLOC")
  ELSE()
    IF(WITH_TCMALLOC)
      MESSAGE(FATAL_ERROR "TCMALLOC required but not found")
    ELSE()
      MESSAGE(STATUS "Disabled: usage of TCMALLOC - TCMALLOC not found")
    ENDIF()
  ENDIF()
ELSE()
  MESSAGE(STATUS "Disabled: usage of TCMALLOC (use -DWITH_TCMALLOC=ON to enable)")
ENDIF()

SET(HAVE_LIBR 0)
IF(NOT(DEFINED WITH_RMATHLIB) OR WITH_RMATHLIB)
  FIND_PACKAGE(LibR)
  IF(LIBR_FOUND)
    MESSAGE(STATUS "Enabled: linking to standalone Rmathlib")
    INCLUDE_DIRECTORIES(${LIBR_INCLUDE_DIRS})
    SET(HAVE_LIBR 1)
  ELSE()
    IF(WITH_RMATHLIB)
      MESSAGE(FATAL_ERROR "Standalone Rmathlib required but not found")
    ELSE()
      MESSAGE(STATUS "Disabled: linking to standalone Rmathlib - Rmathlib not found")
    ENDIF()
  ENDIF()
ELSE()
  MESSAGE(STATUS "Disabled: linking to standalone Rmathlib (use -DWITH_RMATHLIB=ON to enable)")
ENDIF()

IF(NOT(DEFINED WITH_DREME) OR WITH_DREME)
  FIND_PACKAGE(Dreme)
  IF(DREME_FOUND)
    MESSAGE(STATUS "Enabled: DREME seeding method")
  ELSE()
    IF(WITH_DREME)
      MESSAGE(FATAL_ERROR "DREME required but not found")
    ELSE()
      MESSAGE(STATUS "Disabled: DREME seeding method - DREME not found")
    ENDIF()
  ENDIF()
ELSE()
  MESSAGE(STATUS "Disabled: DREME seeding method (use -DWITH_DREME=ON to enable)")
ENDIF()

SET(INSTALL_LOGO_SCRIPT OFF)
IF(WITH_TIKZLOGO)
  IF(NOT(PDFLATEX_COMPILER STREQUAL "PDFLATEX_COMPILER-NOTFOUND"))
    IF(RUBY_FOUND)
      MESSAGE(STATUS "Enabled: sequence logo drawing script")
      SET(INSTALL_LOGO_SCRIPT ON)
    ELSE()
      IF(WITH_TIKZLOGO)
        MESSAGE(FATAL_ERROR "Ruby not found! Sequence logo drawing script needs Ruby!")
      ELSE()
        MESSAGE(STATUS "Disabled: sequence logo drawing script - Ruby not found")
      ENDIF()
    ENDIF()
  ELSE()
    IF(WITH_TIKZLOGO)
      MESSAGE(FATAL_ERROR "PDFlatex not found! Sequence logo drawing script needs PDFlatex!")
    ELSE()
      MESSAGE(STATUS "Disabled: sequence logo drawing script - PDFlatex not found")
    ENDIF()
  ENDIF()
ELSE()
  MESSAGE(STATUS "Disabled: sequence logo drawing script (use -DWITH_TIKZLOGO=ON to enable)")
ENDIF()

SET(INSTALL_MISC_SCRIPTS OFF)
IF(WITH_MISC_SCRIPTS)
  IF(NOT RUBY_FOUND)
    MESSAGE(STATUS "Enabled: miscellaneous scripts")
    SET(INSTALL_MISC_SCRIPTS ON)
  ELSE()
    MESSAGE(FATAL_ERROR "Ruby not found! Miscellaneous scripts disabled!")
  ENDIF()
ELSE()
  MESSAGE(STATUS "Disabled: miscellaneous scripts (use -DWITH_MISC_SCRIPTS=ON to enable)")
ENDIF()

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

INCLUDE(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
IF(GIT_SHA1 STREQUAL "GITDIR-NOTFOUND")
  MESSAGE(STATUS "This is not a build from a git clone.")
  SET(GIT_SHA1 "Not compiled from a git repository.")
  SET(GIT_DESCRIPTION
    "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
  SET(GIT_BRANCH "release")
ELSE()
  SET(GIT_BUILD 1)
  git_describe(GIT_DESCRIPTION)
  git_branch(GIT_BRANCH)
ENDIF()

IF(CTAGS_FOUND AND NOT(APPLE))
  ADD_CUSTOM_TARGET(ctags ALL
    COMMAND ctags -R --c++-kinds=+p --fields=+iaS --extra=+q src/
    COMMENT "Creating Ctags."
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} )
ENDIF()

# Set LIB_SUFFIX if you want to install libraries into a specific directory.
# By default, lib is used.
# This is useful if your distribution prefers 64bit to go into lib64.
# Debian and Ubuntu: there is only /usr/lib.
# Gentoo: /usr/lib is a symlink to /usr/lib64.
# Fedora:
#   architecture-independent libraries go into /usr/lib
#   64-bit libraries go into /usrlib64. In this case set LIB_SUFFIX to "64".
IF(DEFINED LIB_SUFFIX)
  SET(LIBDIR "lib${LIB_SUFFIX}")
ELSE()
  SET(LIBDIR "lib")
ENDIF()
MESSAGE(STATUS "Libraries will be installed into ${CMAKE_INSTALL_PREFIX}/${LIBDIR}")

ADD_SUBDIRECTORY(src build)
ADD_SUBDIRECTORY(scripts)

IF(NOT(DEFINED WITH_DOC) OR WITH_DOC)
  IF(NOT(PDFLATEX_COMPILER STREQUAL "PDFLATEX_COMPILER-NOTFOUND"))
    MESSAGE(STATUS "Enabled: building user manual")
    ADD_SUBDIRECTORY(doc build/doc)
  ELSE()
    IF(WITH_DOC)
      MESSAGE(FATAL_ERROR "PDFlatex not found! Building manual requires PDFlatex")
    ELSE()
      MESSAGE(WARNING "Disabled: building user manual - PDFlatex not found")
    ENDIF()
  ENDIF()
ELSE()
  MESSAGE(STATUS "Disabled: building user manual (use -DWITH_DOC=ON to enable)")
ENDIF()

INSTALL(FILES ${CMAKE_CURRENT_SOURCE_DIR}/README.md
  ${CMAKE_CURRENT_SOURCE_DIR}/COPYING
  ${CMAKE_CURRENT_SOURCE_DIR}/FAQ
  ${CMAKE_CURRENT_SOURCE_DIR}/TODO
  ${CMAKE_CURRENT_SOURCE_DIR}/AUTHORS
  DESTINATION "${CMAKE_INSTALL_PREFIX}/share/doc/${PROJECT_NAME}"
  )
