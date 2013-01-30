cmake_minimum_required (VERSION 2.8)
project (LIMAPACK)

## Print out compiler and linker commands
set (CMAKE_VERBOSE_MAKEFILE "TRUE")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

find_package(libconfig)

#set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake) 
#FIND_PACKAGE(GSL REQUIRED)

# Ascertain if we have MPI available
find_package(MPI)

# Find HDF5 package, required for data output
find_package(HDF5 REQUIRED)

## Compiler flags
if(CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_C_FLAGS "-Wall -Werror -g -O2")
#    set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary
endif()

## Turn on tests
#include(CTest)
#ENABLE_TESTING()

add_subdirectory (src)
#add_subdirectory (tests)
