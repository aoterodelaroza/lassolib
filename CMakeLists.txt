## project and version
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

## some project variables
set(lassolib_MAJOR_VERSION 1)
set(lassolib_MINOR_VERSION 0)
set(lassolib_URL "https://github.com/aoterodelaroza/lassolib")
set(lassolib_EMAIL "aoterodelaroza@gmail.com")
set(lassolib_LICENSE "GNU/GPL version 3")
set(lassolib_DESCRIPTION "A Fortran library for LASSO fitting")
string(TIMESTAMP DATE "%a %d %b %Y %H:%M:%S")

## project definition
project(lassolib VERSION "${lassolib_MAJOR_VERSION}.${lassolib_MINOR_VERSION}" LANGUAGES Fortran)

## global flags for cmake
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_VERBOSE_MAKEFILE ON)
include(GNUInstallDirs)

## default build type, and type options
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
endif()
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release)

## compiler flags, fortran, and check versions
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
##  set(FFLAGS_DEBUG "-g -O0 -C -traceback -debug all -fpp -check all")
##  set(FFLAGS_RELEASE "-O3")
##  if ("${CMAKE_Fortran_COMPILER_VERSION}" VERSION_EQUAL 20.2.1.20201112)
##    set(FFLAGS_DEBUG "-g -O0 -traceback -debug all -fpp")
##  elseif ("${CMAKE_Fortran_COMPILER_VERSION}" VERSION_LESS 19.0)
##    set(_message "\nA version of intel fortran >=19.5 is required to compile the development version of lassolib. "
##      "Please, visit: https://aoterodelaroza.github.io/lassolib/installation/#whichcompilerswork"
##      "for more information.\n")
##    message(FATAL_ERROR ${_message})
##  endif()
elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  set(FFLAGS_DEBUG "-g -g -Og -fcheck=all -fbounds-check -Wall -Wunused-parameter -Wno-maybe-uninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace -fdump-core -cpp")
  set(FFLAGS_RELEASE "-O3")
endif()

## set the compiler flags
string(TOUPPER ${CMAKE_BUILD_TYPE} BT)
set(CMAKE_Fortran_FLAGS_${BT} "${FFLAGS_${BT}}")

## lapack
find_package(LAPACK)
if (NOT LAPACK_FOUND)
   message(FATAL_ERROR "LAPACK was not found and it is required to compile this program.")
endif()

## process subdirectories
add_subdirectory(src)

## summary
set(_message
"\n\n### ${CMAKE_PROJECT_NAME}, ${CMAKE_PROJECT_VERSION}###\n"
"Date: ${DATE}\n"
"Build type: ${CMAKE_BUILD_TYPE}\n"
"Fortran compiler: ${CMAKE_Fortran_COMPILER} (${CMAKE_Fortran_COMPILER_ID}, ${CMAKE_Fortran_COMPILER_VERSION})\n"
"Fortran flags: ${CMAKE_Fortran_FLAGS_${BT}}\n"
"Install directory: ${CMAKE_INSTALL_PREFIX}\n"
"${_message}Lapack: lib=${LAPACK_LIBRARIES}\n"
)
set(_message "${_message}################################\n")
message(STATUS ${_message})
