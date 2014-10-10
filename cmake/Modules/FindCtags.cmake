# - Try to find CTAGS
# Once done this will define
#  CTAGS_FOUND - System has ctags
#  CTAGS_DIR - The directory with the dreme binary

find_path(CTAGS_DIR NAMES ctags)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set CTAGS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(CTAGS DEFAULT_MSG CTAGS_DIR)

