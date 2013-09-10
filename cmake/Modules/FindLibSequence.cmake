# - Try to find LibSequence
# Once done this will define
#  LIBSEQUENCE_FOUND - System has LibSequence
#  LIBSEQUENCE_INCLUDE_DIRS - The LibSequence include directories
#  LIBSEQUENCE_LIBRARIES - The libraries needed to use LibSequence

find_path(LIBSEQUENCE_INCLUDE_DIR 
    NAMES Sequence/Fasta.hpp
    )

find_library(LIBSEQUENCE_LIBRARY NAMES sequence libsequence)
#  xml2 libxml2
#             HINTS ${PC_LIBXML_LIBDIR} ${PC_LIBXML_LIBRARY_DIRS} )

set(LIBSEQUENCE_LIBRARIES ${LIBSEQUENCE_LIBRARY} )
set(LIBSEQUENCE_INCLUDE_DIRS ${LIBSEQUENCE_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBSEQUENCE_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(LibSequence DEFAULT_MSG
                                  LIBSEQUENCE_LIBRARY LIBSEQUENCE_INCLUDE_DIR)


