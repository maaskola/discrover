# - Try to find LibLBFGS
# Once done this will define
#  LIBLBFGS_FOUND - System has LibLBFGS
#  LIBLBFGS_INCLUDE_DIRS - The LibLBFGS include directories
#  LIBLBFGS_LIBRARIES - The libraries needed to use LibLBFGS

find_path(LIBLBFGS_INCLUDE_DIR NAMES lbfgs.h)
# find_path(LIBSLBFGS_INCLUDE_DIR NAMES lbfgs.h PATHS ${CMAKE_INCLUDE_PATH})

find_library(LIBLBFGS_LIBRARY NAMES lbfgs liblbfgs)
#  xml2 libxml2
#             HINTS ${PC_LIBXML_LIBDIR} ${PC_LIBXML_LIBRARY_DIRS} )

set(LIBLBFGS_LIBRARIES ${LIBLBFGS_LIBRARY} )
set(LIBLBFGS_INCLUDE_DIRS ${LIBLBFGS_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBLBFGS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(LibLBFGS DEFAULT_MSG
                                  LIBLBFGS_LIBRARY LIBLBFGS_INCLUDE_DIR)


