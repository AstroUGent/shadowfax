################################################################################
# This file is part of Shadowfax
# Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

# add a new unit test
# a new target with the test sources is constructed, and a CTest test with the
# same name is created
# the new test is also added to the global list of tests contained in the check
# target
macro(add_unit_test)
    set(oneValueArgs NAME NDIM)
    set(multiValueArgs SOURCES LIBS)
    cmake_parse_arguments(TEST "${options}" "${oneValueArgs}"
                               "${multiValueArgs}" ${ARGN})
    message(STATUS "generating " ${TEST_NAME})
    add_executable(${TEST_NAME} EXCLUDE_FROM_ALL ${TEST_SOURCES})
    set_target_properties(${TEST_NAME} PROPERTIES
                          COMPILE_FLAGS -Dndim_=${TEST_NDIM})
    set_target_properties(${TEST_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                          ${PROJECT_BINARY_DIR}/rundir/test)
    target_link_libraries(${TEST_NAME} ${TEST_LIBS})

    add_test(NAME ${TEST_NAME}
             WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/rundir/test
             COMMAND ${TEST_NAME})

    set(TESTNAMES ${TESTNAMES} ${TEST_NAME})
endmacro(add_unit_test)

# RiemannSolver test
set(TESTRIEMANNSOLVER_SOURCES
    myAssert.hpp
    testRiemannSolver.cpp
    ../src/ApproximateSolver.cpp
    ../src/ExactRiemannSolver.cpp
    ../src/MPIGlobal.cpp
    ../src/ProgramLog.cpp
    )

add_unit_test(NAME testRiemannSolver
              SOURCES ${TESTRIEMANNSOLVER_SOURCES}
              NDIM 3
              LIBS ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES} utilities3d)

# ExArith test
set(TESTEXARITH_SOURCES
    testExArith.cpp
    ../src/ExArith.cpp
    ../src/MPIGlobal.cpp
    )

add_unit_test(NAME testExArith2d
              SOURCES ${TESTEXARITH_SOURCES}
              NDIM 2
              LIBS ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES} utilities2d)

add_unit_test(NAME testExArith3d
              SOURCES ${TESTEXARITH_SOURCES}
              NDIM 3
              LIBS ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES} utilities3d)

# DelTess test
set(TESTDELTESS_SOURCES
    testDelTess.cpp
    ../src/DelTess.cpp
    ../src/CubicBox.cpp
    ../src/Line.cpp
    ../src/MPIGlobal.cpp
    ../src/Plane.cpp
    ../src/ProgramLog.cpp
    ../src/RectangularBox.cpp
    ../src/Simplex.cpp
    ../src/VorGen.cpp
    ../src/ExArith.cpp
    )

add_unit_test(NAME testDelTess2d
              SOURCES ${TESTDELTESS_SOURCES}
              NDIM 2
              LIBS ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES} ${HDF5_LIBRARIES}
                   utilities2d)

add_unit_test(NAME testDelTess3d
              SOURCES ${TESTDELTESS_SOURCES}
              NDIM 3
              LIBS ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES} ${HDF5_LIBRARIES}
                   utilities3d)

# Add custom 'make check' target that builds all tests and executes them
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS ${TESTNAMES}
                                                               version_number)
