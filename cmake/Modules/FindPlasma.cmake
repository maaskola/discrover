# - Try to find LibPlasma
# Once done this will define
#  LIBPLASMA_FOUND - System has LibPlasma
#  LIBPLASMA_INCLUDE_DIRS - The LibPlasma include directories
#  LIBPLASMA_LIBRARIES - The libraries needed to use LibPlasma

find_path(LIBPLASMA_INCLUDE_DIR 
    NAMES plasma/plasma.hpp
    )

find_library(LIBPLASMA_LIBRARY NAMES plasma libplasma)
#  xml2 libxml2
#             HINTS ${PC_LIBXML_LIBDIR} ${PC_LIBXML_LIBRARY_DIRS} )

set(LIBPLASMA_LIBRARIES ${LIBPLASMA_LIBRARY} )
set(LIBPLASMA_INCLUDE_DIRS ${LIBPLASMA_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBPLASMA_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(LibPlasma DEFAULT_MSG
                                  LIBPLASMA_LIBRARY LIBPLASMA_INCLUDE_DIR)


