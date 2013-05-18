cmake_minimum_required (VERSION 2.8)

### Project information
project (limapack C)
set (LIMAPACK_MAJOR_VERSION "0")
set (LIMAPACK_MINOR_VERSION "1")
set (LIMAPACK_PATCH_VERSION "0")
set (LIMAPACK_VERSION 
     "${LIMAPACK_MAJOR_VERSION}.${LIMAPACK_MINOR_VERSION}.${LIMAPACK_PATCH_VERSION}")

### Print out compiler and linker commands
set (CMAKE_VERBOSE_MAKEFILE "TRUE")

### Add project cmake module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

#### Find needed packages
find_package(libconfig REQUIRED)
find_package(GSL REQUIRED)
find_package(MPI)
find_package(HDF5 REQUIRED)
find_package(Threads)
find_library (MATH_LIBRARIES m)

## Compiler flags
if(CMAKE_COMPILER_IS_GNUCC)
#    set(CMAKE_C_FLAGS "-Wall -Werror -g -O2")
    set(CMAKE_C_FLAGS "-Wall -Werror -g -O")
#    set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary
endif()

## Turn on tests
#include(CTest)
#ENABLE_TESTING()

add_subdirectory (src)
#add_subdirectory (tests)

### Source tarball generation
# See http://agateau.com/2009/08/09/cmake-and-make-dist-the-simple-version/
set (TARBALL_NAME "${CMAKE_PROJECT_NAME}-${LIMAPACK_VERSION}")

# Use CPACK to create a make dist target to produce a tarball
set (CPACK_SOURCE_GENERATOR "TBZ2")
set (CPACK_SOURCE_PACKAGE_FILE_NAME "${TARBALL_NAME}")
set (CPACK_SOURCE_IGNORE_FILES "${CMAKE_CURRENT_BINARY_DIR}/;.git/;CMakeFiles/;Makefile;.gitignore;~$;${CPACK_SOURCE_IGNORE_FILES}")
include (CPack)
add_custom_target (dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)

# Alternative - create a make gitdist tarket which uses git to create the tarball
add_custom_target(gitdist
    COMMAND git archive --prefix=${TARBALL_NAME}/ HEAD
        | bzip2 > ${CMAKE_BINARY_DIR}/${TARBALL_NAME}.tar.bz2
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
