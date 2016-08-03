################################################################################
# This file is part of Shadowfax
# Copyright (C) 2015 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
#               2016 Bert Vandenbroucke
#
# Shadowfax is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Shadowfax is distributed in the hope that it will be useful,
# but WITOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with Shadowfax. If not, see <http://www.gnu.org/licenses/>.
################################################################################

cmake_minimum_required(VERSION 2.6)
project(shadowfax)

# The version number.
set(SHADOWFAX_VERSION_MAJOR 1)
set(SHADOWFAX_VERSION_MINOR 2)

# Variables used throughout the build process
set(COOLINGTABLE_LOCATION "${PROJECT_BINARY_DIR}/rundir/coolingtables/")
set(COOLINGTABLE_URL
    "http://users.ugent.be/~sdrijcke/CoolingHeatingTables.tar.gz")

enable_testing()

# suppress warnings about duplicate build directories (caused by the foreach
# over dimensions)
# these warnings are unnecessary in this case and displaying them may confuse
# the enduser
cmake_policy(SET CMP0013 OLD)

add_definitions(-std=c++0x)

set(CMAKE_CXX_FLAGS_PROFILE "-pg -g")
set(CMAKE_C_FLAGS_PROFILE "-pg -g")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wignored-qualifiers -Werror")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs")

if(NO_SANITIZER)
    set(SANITIZER False)
else(NO_SANITIZER)
    set(SANITIZER True)
endif(NO_SANITIZER)

# currently the only way to link in the new GNU address sanitizer
if(SANITIZER)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        if(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.8")
            set(CMAKE_CXX_FLAGS_DEBUG
                "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address")
            set(CMAKE_CXX_FLAGS_DEBUG
                "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer")
            set(CMAKE_EXE_LINKER_FLAGS_DEBUG
                "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -fsanitize=address")
            add_definitions(-DSANITIZER)
        endif(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.8")
    endif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        set(CMAKE_CXX_FLAGS_DEBUG
            "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address")
        set(CMAKE_CXX_FLAGS_DEBUG
            "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer")
        set(CMAKE_EXE_LINKER_FLAGS_DEBUG
            "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -fsanitize=address")
        add_definitions(-DSANITIZER)
    endif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
endif(SANITIZER)

# fix for strange errors on the UGent HPC
if(NO_MPI_SEEK)
    add_definitions(-DMPICH_IGNORE_CXX_SEEK)
endif(NO_MPI_SEEK)

# HDF5 somewhere in between versions 1.6 and 1.8 changed its API
# old systems are not aware of this and need the 1.6 API. Luckily, the 1.8
# API provides macros that convert 1.6 API function calls to the correct
# 1.8 version. By specifying this flag, we indicate that we will use 1.6 version
# function calls everywhere. Old systems will completely ignore this flag and
# everything will work. New systems will be compiled using the compatibility
# macros and will hence also work. Everybody happy.
# If at some point in the future, HDF5 stops supporting the old API, some work
# will have to be done to replace old function calls by new ones and this flag
# may be removed.
#
# more info: http://www.hdfgroup.org/HDF5/doc/RM/APICompatMacros.html
add_definitions(-DH5_USE_16_API)

# the include directories
include_directories(${PROJECT_BINARY_DIR}/inc)
include_directories(${PROJECT_SOURCE_DIR}/src)

find_package(MPI REQUIRED)
message(STATUS "Found MPI_EXEC: " ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG}
        ${MPIEXEC_PREFLAGS} ${MPIEXEC_POSTFLAGS})
include_directories(${MPI_INCLUDE_PATH})

find_package(Boost 1.53.0 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})
set(ALL_Boost_LIBRARIES ${Boost_LIBRARIES})

find_package(HDF5 REQUIRED)
message(STATUS "Found HDF5_DIFF: " ${HDF5_DIFF_EXECUTABLE})
include_directories(${HDF5_INCLUDE_DIR})
link_directories(${HDF5_LIBRARY_DIRS})

if(PYTHON_MODULES)
    find_package(Boost 1.53.0 COMPONENTS "python" REQUIRED)
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
    set(Boost_Python_LIBRARIES ${Boost_LIBRARIES})

    find_package(PythonLibs 2.7 REQUIRED)
    include_directories(${PYTHON_INCLUDE_PATH})
endif(PYTHON_MODULES)

# select the programs to compile
if(ONLY_SIMULATION_PROGRAM)
    set(DO_SHADOWFAX True)
    set(DO_ICMAKER False)
    set(DO_VTKMAKER False)
else(ONLY_SIMULATION_PROGRAM)
    set(DO_SHADOWFAX True)
    set(DO_ICMAKER True)
    set(DO_VTKMAKER True)
endif(ONLY_SIMULATION_PROGRAM)

if(GIZMO)
    set(DO_GIZMO True)
else(GIZMO)
    set(DO_GIZMO False)
endif(GIZMO)

if(TENET)
    set(DO_TENET True)
else(TENET)
    set(DO_TENET False)
endif(TENET)

if(LLOYD)
    set(DO_LLOYD True)
else(LLOYD)
    set(DO_LLOYD False)
endif(LLOYD)

# put executable in the rundirs
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/rundir)

# now do the source dir, creates the executable
# we perform a loop to make separate executables for 2 and 3 dimensions
if(ONLY_3D)
    set(NDIMS 3)
elseif(ONLY_2D)
    set(NDIMS 2)
else(ONLY_3D)
    set(NDIMS 2 3)
endif(ONLY_3D)

foreach(NDIM ${NDIMS})
    add_subdirectory(src)
endforeach(NDIM)

add_subdirectory(test)

# copy the testsuite
message(STATUS "Copying testsuite...")
file(GLOB_RECURSE TESTSUITE_FILES
     "${PROJECT_SOURCE_DIR}/rundir/testsuite/*")
foreach(TESTSUITE_FILE ${TESTSUITE_FILES})
    file(RELATIVE_PATH TESTSUITE_FILE_ONLYNAME ${PROJECT_SOURCE_DIR}
         ${TESTSUITE_FILE})
    # configure_file only works for text files if COPYONLY is not specified
    get_filename_component(TESTSUITE_FILE_EXT ${TESTSUITE_FILE} EXT)
    if(${TESTSUITE_FILE_EXT} STREQUAL ".hdf5")
        configure_file(${TESTSUITE_FILE}
                       ${PROJECT_BINARY_DIR}/${TESTSUITE_FILE_ONLYNAME}
                       COPYONLY)
    else(${TESTSUITE_FILE_EXT} STREQUAL ".hdf5")
        configure_file(${TESTSUITE_FILE}
                       ${PROJECT_BINARY_DIR}/${TESTSUITE_FILE_ONLYNAME} @ONLY)
    endif(${TESTSUITE_FILE_EXT} STREQUAL ".hdf5")
endforeach(TESTSUITE_FILE ${TESTSUITE_FILES})

# doxygen
macro(add_doc_debug)
    set(oneValueArgs FOLDER EXTENSION)
    cmake_parse_arguments(DOC_DEBUG "${options}" "${oneValueArgs}"
                                    "${multiValueArgs}" ${ARGN})
    file(GLOB DOC_SRCS "${DOC_DEBUG_FOLDER}/*.${DOC_DEBUG_EXTENSION}")
    foreach(DOC_SRC ${DOC_SRCS})
        get_filename_component(DOC_SRC_NAME_UNPARSED ${DOC_SRC} NAME_WE)
        string(REGEX REPLACE "_" "__" DOC_SRC_NAME ${DOC_SRC_NAME_UNPARSED})
        set(DOC_HTML_NAME "${DOC_SRC_NAME}_8${DOC_DEBUG_EXTENSION}.html")
        set(DOC_OUTPUT_NAME ${PROJECT_BINARY_DIR}/doc/html/${DOC_HTML_NAME})
        if(${DOC_DEBUG_EXTENSION} STREQUAL "cpp")
            file(GLOB_RECURSE DOC_SRC_EXTRA
                 "${DOC_DEBUG_FOLDER}/${DOC_SRC_NAME_UNPARSED}.hpp")
        else(${DOC_DEBUG_EXTENSION} STREQUAL "cpp")
            file(GLOB_RECURSE DOC_SRC_EXTRA
                 "${DOC_DEBUG_FOLDER}/${DOC_SRC_NAME_UNPARSED}.cpp")
        endif(${DOC_DEBUG_EXTENSION} STREQUAL "cpp")
        set(DOC_SRC "${DOC_SRC};${DOC_SRC_EXTRA}")
        set(DOC_FILTER
            "${DOC_SRC_NAME_UNPARSED}.cpp ${DOC_SRC_NAME_UNPARSED}.hpp")
        set(DOC_COMMAND "doxygen ${PROJECT_BINARY_DIR}/Doxyfile.debug")
        add_custom_command(
            OUTPUT ${DOC_OUTPUT_NAME}
            COMMAND eval "DOXYGEN_INPUT_FILE=\"${DOC_FILTER}\" ${DOC_COMMAND}"
            DEPENDS ${DOC_SRC}
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            VERBATIM)
        set(DOC_NAMES ${DOC_NAMES} ${DOC_OUTPUT_NAME})
    endforeach(DOC_SRC ${DOC_SRCS})
endmacro(add_doc_debug)

macro(add_doc_debug_folder)
    set(oneValueArgs FOLDER)
    cmake_parse_arguments(DOC_DEBUG_F "${options}" "${oneValueArgs}"
                                      "${multiValueArgs}" ${ARGN})
    add_doc_debug(FOLDER ${DOC_DEBUG_F_FOLDER} EXTENSION "cpp")
    add_doc_debug(FOLDER ${DOC_DEBUG_F_FOLDER} EXTENSION "hpp")
endmacro(add_doc_debug_folder)

find_package(Doxygen)
if(DOXYGEN_FOUND)
configure_file(${PROJECT_SOURCE_DIR}/Doxyfile.in
               ${PROJECT_BINARY_DIR}/Doxyfile @ONLY)
add_custom_command(OUTPUT ${PROJECT_BINARY_DIR}/doc/html/index.html
                   COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
                   WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                   COMMENT "Generating API documentation with Doxygen" VERBATIM)
add_custom_target(doc DEPENDS ${PROJECT_BINARY_DIR}/doc/html/index.html)

configure_file(${PROJECT_SOURCE_DIR}/Doxyfile.debug.in
               ${PROJECT_BINARY_DIR}/Doxyfile.debug @ONLY)
add_doc_debug_folder(FOLDER "${PROJECT_SOURCE_DIR}/src")
add_doc_debug_folder(FOLDER "${PROJECT_SOURCE_DIR}/src/exact_arithmetics")
add_doc_debug_folder(FOLDER "${PROJECT_SOURCE_DIR}/src/fortune")
add_doc_debug_folder(FOLDER "${PROJECT_SOURCE_DIR}/src/io")
add_doc_debug_folder(FOLDER "${PROJECT_SOURCE_DIR}/src/python")
add_doc_debug_folder(FOLDER "${PROJECT_SOURCE_DIR}/src/riemann")
add_doc_debug_folder(FOLDER "${PROJECT_SOURCE_DIR}/src/utilities")
add_doc_debug_folder(FOLDER "${PROJECT_SOURCE_DIR}/test")
add_custom_target(doc_debug DEPENDS ${DOC_NAMES}
                                    ${PROJECT_BINARY_DIR}/doc/html/index.html)
endif(DOXYGEN_FOUND)

if(EXISTS ${COOLINGTABLE_LOCATION}RadLoss_0.00_0.00_0.00_1.00e+00.rates)
    message(STATUS "Cooling tables already present")
else(EXISTS ${COOLINGTABLE_LOCATION}RadLoss_0.00_0.00_0.00_1.00e+00.rates)
    message(STATUS "Downloading and unpacking coolingtables...")
    file(DOWNLOAD ${COOLINGTABLE_URL} ${COOLINGTABLE_LOCATION}cool.tar.gz)
    execute_process(
    COMMAND ${CMAKE_COMMAND} -E tar xzf ${COOLINGTABLE_LOCATION}cool.tar.gz
    WORKING_DIRECTORY ${COOLINGTABLE_LOCATION})
    file(GLOB COOLINGTABLES ${COOLINGTABLE_LOCATION}RadLoss/*.rates)
    file(COPY ${COOLINGTABLES} DESTINATION ${COOLINGTABLE_LOCATION})
    file(REMOVE_RECURSE ${COOLINGTABLE_LOCATION}RadLoss)
    file(REMOVE ${COOLINGTABLE_LOCATION}README)
    file(REMOVE ${COOLINGTABLE_LOCATION}cool.tar.gz)
    configure_file(${PROJECT_SOURCE_DIR}/rundir/coolingtables/filecopy.py
                   ${COOLINGTABLE_LOCATION}filecopy.py COPYONLY)
    find_program(PYTHON python)
    execute_process(COMMAND ${PYTHON} ${COOLINGTABLE_LOCATION}filecopy.py
                    WORKING_DIRECTORY ${COOLINGTABLE_LOCATION})
endif(EXISTS ${COOLINGTABLE_LOCATION}RadLoss_0.00_0.00_0.00_1.00e+00.rates)
