# - Try to find DREME from the MEME suite
# Once done this will define
#  DREME_FOUND - System has Dreme
#  DREME_DIR - The directory with the dreme binary

find_path(DREME_DIR NAMES dreme meme-dreme)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBGOOGLEPERFTOOLS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(DREME DEFAULT_MSG DREME_DIR)

