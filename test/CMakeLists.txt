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

macro(add_unit_test_bash)
    set(oneValueArgs NAME NDIM SCRIPT)
    set(multiValueArgs SOURCES LIBS FILES)
    find_program(BASH_PROGRAM bash)
    cmake_parse_arguments(TEST "${options}" "${oneValueArgs}"
                               "${multiValueArgs}" ${ARGN})
    message(STATUS "generating " ${TEST_NAME})
    add_executable(${TEST_NAME} EXCLUDE_FROM_ALL ${TEST_SOURCES})
    set_target_properties(${TEST_NAME} PROPERTIES
                          COMPILE_FLAGS -Dndim_=${TEST_NDIM})
    set_target_properties(${TEST_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                          ${PROJECT_BINARY_DIR}/rundir/test)
    target_link_libraries(${TEST_NAME} ${TEST_LIBS})
    foreach(TEST_FILE ${TEST_FILES})
        message(STATUS ${TEST_FILE})
        configure_file(${TEST_FILE}
                       ${PROJECT_BINARY_DIR}/rundir/test/${TEST_FILE}
                       COPYONLY)
    endforeach(TEST_FILE ${TEST_FILES})
    add_test(NAME ${TEST_NAME}
             WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/rundir/test
             COMMAND ${BASH_PROGRAM} ${TEST_SCRIPT})

    set(TESTNAMES ${TESTNAMES} ${TEST_NAME})
endmacro(add_unit_test_bash)

# RiemannSolver test
set(TESTRIEMANNSOLVER_SOURCES
    myAssert.hpp
    testRiemannSolver.cpp
    ../src/MPIGlobal.cpp
    ../src/ProgramLog.cpp
    )

add_unit_test(NAME testRiemannSolver
              SOURCES ${TESTRIEMANNSOLVER_SOURCES}
              NDIM 3
              LIBS ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES} riemann3d
              utilities3d)

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

# Unit test (for the Unit class)
set(TESTUNIT_SOURCES
    testUnit.cpp

    ../src/io/Unit.hpp
)

add_unit_test(NAME testUnit
              SOURCES ${TESTUNIT_SOURCES}
              NDIM 2
              LIBS ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES})

#NDTable test
set(TESTNDTABLE_SOURCES
    testNDTable.cpp
    ../src/NDTable.cpp
    ../src/MPIGlobal.cpp
    ../src/io/UnitSet.cpp
    ../src/utilities/Timer.cpp
    )

set(TESTNDTABLE_FILES
    testND.sh
    gentabs.py
    )

add_unit_test_bash(NAME testNDTable
              SCRIPT testND.sh
              SOURCES ${TESTNDTABLE_SOURCES}
              NDIM 2
              FILES ${TESTNDTABLE_FILES}
              LIBS ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES})

#NDIrregTable test
set(TESTNDIRREGTABLE_SOURCES
    testNDIrregTable.cpp
    ../src/NDIrregTable.cpp
    ../src/MPIGlobal.cpp
    ../src/io/UnitSet.cpp
    ../src/utilities/Timer.cpp
    )

set(TESTNDIRREGTABLE_FILES
    testNDIrreg.sh
    gentabs2.py
    )

#GasCooling test
set(TESTGASCOOLING_SOURCES
    testGasCooling.cpp

    ../src/GasCooling.cpp
    ../src/CoolingTable.cpp
    ../src/NDIrregTable.cpp
    ../src/NDTable.cpp
    ../src/ParameterFile.cpp
    ../src/utilities/GasParticle.cpp
    ../src/utilities/Particle.cpp
    ../src/utilities/Hilbert.cpp
    ../src/utilities/HelperFunctions.cpp
    ../src/MPIGlobal.cpp
    ../src/Physics.cpp
    ../src/VorGen.cpp
    ../src/io/UnitSet.cpp
    ../src/utilities/Timer.cpp
    )

set(TESTGASCOOLING_FILES
    testGasCooling.sh
    difcool.py
    )

add_unit_test_bash(NAME testGasCooling
              SCRIPT testGasCooling.sh
              SOURCES ${TESTGASCOOLING_SOURCES}
              NDIM 3
              FILES ${TESTGASCOOLING_FILES}
              LIBS ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES})

# GSL test
set(TESTGSL_SOURCES
    testGSL.cpp

    ../src/GSL.cpp
    ../src/GSL.hpp
    ../src/MPIGlobal.cpp
)

add_unit_test(NAME testGSL
              SOURCES ${TESTGSL_SOURCES}
              NDIM 2
              LIBS ${MPI_CXX_LIBRARIES} utilities2d)

# Cosmology test
set(TESTCOSMOLOGY_SOURCES
    testCosmology.cpp

    ../src/Cosmology.cpp
    ../src/Cosmology.hpp
    ../src/GSL.cpp
    ../src/MPIGlobal.cpp
)

add_unit_test(NAME testCosmology
              SOURCES ${TESTCOSMOLOGY_SOURCES}
              NDIM 2
              LIBS ${MPI_CXX_LIBRARIES} utilities2d)

#SpecificICGenerator test
set(TESTSPECIFICICGENERATOR_SOURCES
    testSpecificICGenerator.cpp

    ../src/ColorMap.cpp
    ../src/DelTess.cpp
    ../src/DelTess_mpi.cpp
    ../src/ExArith.cpp
    ../src/Line.cpp
    ../src/MPIGlobal.cpp
    ../src/NgbSearch.cpp
    ../src/Plane.cpp
    ../src/RectangularBox.cpp
    ../src/Simplex.cpp
    ../src/SpecificICGenerator.cpp
    ../src/VorCell.cpp
    ../src/VorFace.cpp
    ../src/VorGen.cpp
    ../src/VorTess.cpp
    ../src/VorTessManager.cpp
    )

add_unit_test(NAME testSpecificICGenerator
              SOURCES ${TESTSPECIFICICGENERATOR_SOURCES}
              NDIM 3
              LIBS ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES} utilities3d io3d)

# Legendre test
set(TESTLEGENDRE_SOURCES
    testLegendre.cpp
    ../src/Legendre.cpp
    ../src/utilities/Timer.cpp

    ../src/Legendre.hpp
)

add_unit_test(NAME testLegendre
              SOURCES ${TESTLEGENDRE_SOURCES}
              NDIM 3
              LIBS ${MPI_CXX_LIBRARIES})

# Runge-Kutta test
set(TESTRK_SOURCES
    testRungeKutta.cpp

    ../src/RungeKutta.hpp
    ../src/RungeKuttaFlux.hpp
)

add_unit_test(NAME testRungeKutta
              SOURCES ${TESTRK_SOURCES}
              NDIM 2)

# Add custom 'make check' target that builds all tests and executes them
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS ${TESTNAMES}
                                                               version_number)
