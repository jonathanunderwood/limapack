# Try to find libconfig
#
# This will define:
#    libconfig_FOUND
#    libconfig_INCLUDE_DIRS
#    libconfig_LIBRARIES
#

# Use pkg-config to take a guess at header and library locations
find_package(PkgConfig QUIET)
pkg_check_modules(PC_LIBCONFIG QUIET libconfig)

find_path(LIBCONFIG_INCLUDE_DIR
    NAMES libconfig.h
    HINTS ${PC_LIBCONFIG_INCLUDEDIR} ${PC_LIBCONFIG_INCLUDE_DIRS}
)

find_library(LIBCONFIG_LIBRARIES
    NAMES config libconfig
    HINTS  ${PC_LIBCONFIG_LIBDIR} ${PC_LIBCONFIG_LIBRARY_DIRS}
)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(libconfig
                                  REQUIRED_VARS LIBCONFIG_LIBRARIES LIBCONFIG_INCLUDE_DIR
				  )


