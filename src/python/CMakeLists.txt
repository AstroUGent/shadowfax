################################################################################
# This file is part of Shadowfax
# Copyright (C) 2015 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

set(PYTHON_RIEMANNNAME python_riemann${NDIM}d)
set(PYTHON_IONAME python_io${NDIM}d)
set(PYTHON_GRIDNAME python_grid${NDIM}d)
set(PYTHON_LLOYDNAME python_lloyd${NDIM}d)
set(SHAPYNAME shapy${NDIM}d)

set(PYTHON_RIEMANNSOURCES
    RiemannModule.cpp
    ../riemann/ApproximateSolver.cpp
    ../riemann/ExactRiemannSolver.cpp
    ../riemann/RiemannSolver.cpp
    ../MPIGlobal.cpp
    ../ProgramLog.cpp
    ../utilities/HelperFunctions.cpp
    ../utilities/Timer.cpp
    )

add_library(${PYTHON_RIEMANNNAME} SHARED ${PYTHON_RIEMANNSOURCES})

set_target_properties(${PYTHON_RIEMANNNAME} PROPERTIES
                      COMPILE_FLAGS "-DPYTHON_MODULE -Dndim_=${NDIM}")
set_target_properties(${PYTHON_RIEMANNNAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY
                      ${PROJECT_BINARY_DIR}/${SHAPYNAME})

target_link_libraries(${PYTHON_RIEMANNNAME} ${MPI_CXX_LIBRARIES})
target_link_libraries(${PYTHON_RIEMANNNAME} ${MPI_C_LIBRARIES})
target_link_libraries(${PYTHON_RIEMANNNAME} ${ALL_Boost_LIBRARIES})
target_link_libraries(${PYTHON_RIEMANNNAME} ${Boost_Python_LIBRARIES})

set(PYTHON_IO_SOURCES
IOModule.cpp
../MPIGlobal.cpp
../io/AsciiInput.cpp
../io/AsciiOutput.cpp
../io/Block.cpp
../io/FileInput.cpp
../io/FileOutput.cpp
../io/Header.cpp
../io/UnitSet.cpp
../utilities/Timer.cpp
)

add_library(${PYTHON_IONAME} SHARED ${PYTHON_IO_SOURCES})
set_target_properties(${PYTHON_IONAME} PROPERTIES
                      COMPILE_FLAGS "-DPYTHON_MODULE -Dndim_=${NDIM}")
set_target_properties(${PYTHON_IONAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY
                      ${PROJECT_BINARY_DIR}/${SHAPYNAME})

target_link_libraries(${PYTHON_IONAME} ${MPI_CXX_LIBRARIES})
target_link_libraries(${PYTHON_IONAME} ${MPI_C_LIBRARIES})
target_link_libraries(${PYTHON_IONAME} ${ALL_Boost_LIBRARIES})
target_link_libraries(${PYTHON_IONAME} ${Boost_Python_LIBRARIES})
target_link_libraries(${PYTHON_IONAME} ${HDF5_LIBRARIES})

set(PYTHON_GRIDSOURCES
GridModule.cpp
../ColorMap.cpp
../CubicBox.cpp
../DelTess.cpp
../DelTess_mpi.cpp
../ExArith.cpp
../Line.cpp
../MPIGlobal.cpp
../NgbSearch.cpp
../Plane.cpp
../ProgramLog.cpp
../RectangularBox.cpp
../Simplex.cpp
../VorCell.cpp
../VorFace.cpp
../VorGen.cpp
../VorTess.cpp
../io/Header.cpp
../utilities/DMParticle.cpp
../utilities/EwaldTable.cpp
../utilities/GasParticle.cpp
../utilities/HelperFunctions.cpp
../utilities/Hilbert.cpp
../utilities/Particle.cpp
../utilities/ParticleVector.cpp
../utilities/StarParticle.cpp
../utilities/Timer.cpp
../utilities/Tree.cpp
)

add_library(${PYTHON_GRIDNAME} SHARED ${PYTHON_GRIDSOURCES})
set_target_properties(${PYTHON_GRIDNAME} PROPERTIES
                      COMPILE_FLAGS "-DPYTHON_MODULE -Dndim_=${NDIM}")
set_target_properties(${PYTHON_GRIDNAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY
                      ${PROJECT_BINARY_DIR}/${SHAPYNAME})

target_link_libraries(${PYTHON_GRIDNAME} ${MPI_CXX_LIBRARIES})
target_link_libraries(${PYTHON_GRIDNAME} ${MPI_C_LIBRARIES})
target_link_libraries(${PYTHON_GRIDNAME} ${ALL_Boost_LIBRARIES})
target_link_libraries(${PYTHON_GRIDNAME} ${Boost_Python_LIBRARIES})
target_link_libraries(${PYTHON_GRIDNAME} ${HDF5_LIBRARIES})

set(PYTHON_LLOYDSOURCES
LloydModule.cpp
../ColorMap.cpp
../CubicBox.cpp
../DelTess.cpp
../DelTess_mpi.cpp
../ExArith.cpp
../GadgetSnapshotReader.cpp
../Line.cpp
../Lloyd.cpp
../MPIGlobal.cpp
../NgbSearch.cpp
../Plane.cpp
../ProgramLog.cpp
../RectangularBox.cpp
../ShadowfaxSnapshotReader.cpp
../Simplex.cpp
../SnapshotHandler.cpp
../VorCell.cpp
../VorFace.cpp
../VorGen.cpp
../VorTess.cpp
../io/Block.cpp
../io/FileInput.cpp
../io/FileOutput.cpp
../io/Header.cpp
../io/UnitSet.cpp
../utilities/DMParticle.cpp
../utilities/EwaldTable.cpp
../utilities/GasParticle.cpp
../utilities/HelperFunctions.cpp
../utilities/Hilbert.cpp
../utilities/Particle.cpp
../utilities/ParticleVector.cpp
../utilities/StarParticle.cpp
../utilities/Timer.cpp
../utilities/Tree.cpp
)

add_library(${PYTHON_LLOYDNAME} SHARED ${PYTHON_LLOYDSOURCES})

set_target_properties(${PYTHON_LLOYDNAME} PROPERTIES
                      COMPILE_FLAGS "-DPYTHON_MODULE -Dndim_=${NDIM}")
set_target_properties(${PYTHON_LLOYDNAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY
                      ${PROJECT_BINARY_DIR}/${SHAPYNAME})

target_link_libraries(${PYTHON_LLOYDNAME} ${MPI_CXX_LIBRARIES})
target_link_libraries(${PYTHON_LLOYDNAME} ${MPI_C_LIBRARIES})
target_link_libraries(${PYTHON_LLOYDNAME} ${HDF5_LIBRARIES})
target_link_libraries(${PYTHON_LLOYDNAME} ${ALL_Boost_LIBRARIES})
target_link_libraries(${PYTHON_LLOYDNAME} ${Boost_Python_LIBRARIES})

file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/${SHAPYNAME})
configure_file(${PROJECT_SOURCE_DIR}/${SHAPYNAME}/__init__.py
               ${PROJECT_BINARY_DIR}/${SHAPYNAME} COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/${SHAPYNAME}/si_units.py
               ${PROJECT_BINARY_DIR}/${SHAPYNAME} COPYONLY)
