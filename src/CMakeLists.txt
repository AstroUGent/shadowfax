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

set(IONAME io${NDIM}d)
set(RIEMANNNAME riemann${NDIM}d)
set(UTILNAME utilities${NDIM}d)

add_subdirectory(io)
add_subdirectory(riemann)
add_subdirectory(utilities)

if(DO_SHADOWFAX)
    set(SHADOWFAXNAME shadowfax${NDIM}d)

    set(SHADOWFAXSOURCES
    CircularBox.cpp
    CoolingTable.cpp
    CubicBox.cpp
    DelTess.cpp
    DelTess_mpi.cpp
    DelTessCheck.cpp
    ExArith.cpp
    FixedGrid.cpp
    GadgetSnapshotReader.cpp
    GadgetSnapshotWriter.cpp
    GasCooling.cpp
    GravityWalker.cpp
    Image.cpp
    Line.cpp
    LogFiles.cpp
    MPIGlobal.cpp
    NDIrregTable.cpp
    NDTable.cpp
    NgbSearch.cpp
    ParameterFile.cpp
    Physics.hpp
    Plane.cpp
    Plotter.cpp
    ProgramLog.cpp
    RectangularBox.cpp
    RestartFile.cpp
    Shadowfax.cpp
    ShadowfaxSnapshotReader.cpp
    ShadowfaxSnapshotWriter.cpp
    SidePrograms.cpp
    Simplex.cpp
    Simulation.cpp
    SnapshotHandler.cpp
    StarFormationParticleConverter.cpp
    StellarFeedback.cpp
    TimeLine.cpp
    TimeStepWalker.cpp
    VorCell.cpp
    VorFace.cpp
    VorGen.cpp
    VorTess.cpp
    VorTessManager.cpp

    CoolingLocation.hpp
    CoolingTable.hpp
    DelCont.hpp
    DelTess.hpp
    Error.hpp
    ExArith.hpp
    FixedGrid.hpp
    GadgetSnapshotReader.hpp
    GadgetSnapshotWriter.hpp
    GasCooling.hpp
    GravityWalker.hpp
    Image.hpp
    Line.hpp
    LogFiles.hpp
    MPIGlobal.hpp
    MPIMethods.hpp
    NDIrregTable.hpp
    NDTable.hpp
    NgbSearch.hpp
    ParameterFile.hpp
    Physics.cpp
    Plane.hpp
    Plotter.hpp
    ProgramLog.hpp
    RestartFile.hpp
    SidePrograms.hpp
    Simplex.hpp
    Simulation.hpp
    SnapshotHandler.hpp
    SnapshotReaderFactory.hpp
    SnapshotWriterFactory.hpp
    StarFormationParticleConverter.hpp
    StateVector.hpp
    StellarFeedback.hpp
    TimeLine.hpp
    TimeStepWalker.hpp
    Vec.hpp
    VorCell.hpp
    VorFace.hpp
    VorGen.hpp
    VorTess.hpp
    VorTessManager.hpp
    )

    #add the location of the cooling tables
    configure_file(${PROJECT_SOURCE_DIR}/src/CoolingLocation.hpp.in
                   ${PROJECT_BINARY_DIR}/inc/CoolingLocation.hpp)

    # create a new version file if one of the source files changes
    # this way, we will always have the correct compilation time in our
    # restart file
    add_custom_command(OUTPUT ${PROJECT_BINARY_DIR}/inc/ShadowfaxVersion.hpp
                       COMMAND
                       cmake -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                       -DPROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}
                       -DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
                       -P ${PROJECT_SOURCE_DIR}/write_version.cmake
                       WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                       DEPENDS ${SHADOWFAXSOURCES} ${IONAME} ${UTILNAME})

    add_executable(${SHADOWFAXNAME} ${SHADOWFAXSOURCES})
    set_target_properties(${SHADOWFAXNAME} PROPERTIES
                          COMPILE_FLAGS -Dndim_=${NDIM})

    target_link_libraries(${SHADOWFAXNAME} ${IONAME})
    target_link_libraries(${SHADOWFAXNAME} ${RIEMANNNAME})
    target_link_libraries(${SHADOWFAXNAME} ${UTILNAME})
    target_link_libraries(${SHADOWFAXNAME} ${MPI_CXX_LIBRARIES})
    target_link_libraries(${SHADOWFAXNAME} ${MPI_C_LIBRARIES})
    target_link_libraries(${SHADOWFAXNAME} ${HDF5_LIBRARIES})
    target_link_libraries(${SHADOWFAXNAME} ${ALL_Boost_LIBRARIES})

    # make sure version_number is called, but only once
    if(${NDIM} EQUAL 2 OR ONLY_3D)
        add_custom_target(version_number
                          DEPENDS
                          ${PROJECT_BINARY_DIR}/inc/ShadowfaxVersion.hpp)
        add_dependencies(${SHADOWFAXNAME} version_number)
    endif(${NDIM} EQUAL 2 OR ONLY_3D)
endif(DO_SHADOWFAX)

if(DO_ICMAKER)
    set(ICMAKERNAME icmaker${NDIM}d)

    set(ICMAKERSOURCES
    BlockICGenerator.cpp
    CircularBox.cpp
    CubicBox.cpp
    DelTess.cpp
    DelTess_mpi.cpp
    DelTessCheck.cpp
    ExArith.cpp
    GadgetSnapshotWriter.cpp
    ICMaker.cpp
    ICRegion.cpp
    Line.cpp
    MPIGlobal.cpp
    NgbSearch.cpp
    Plane.cpp
    ProgramLog.cpp
    RectangularBox.cpp
    ShadowfaxSnapshotReader.cpp
    ShadowfaxSnapshotWriter.cpp
    Simplex.cpp
    SnapshotHandler.cpp
    SpecificICGenerator.cpp
    TimeLine.cpp
    TimeStepWalker.cpp
    VorCell.cpp
    VorFace.cpp
    VorGen.cpp
    VorTess.cpp
    VorTessManager.cpp
    
    DelCont.hpp
    DelTess.hpp
    ExArith.hpp
    ICGenerator.hpp
    ICMaker.hpp
    ICRegion.hpp
    Line.hpp
    MPIGlobal.hpp
    MPIMethods.hpp
    Plane.hpp
    ProgramLog.hpp
    Simplex.hpp
    ShadowfaxSnapshotReader.hpp
    ShadowfaxSnapshotWriter.hpp
    SnapshotHandler.hpp
    StateVector.hpp
    SymbolicFunction.hpp
    Vec.hpp
    VorCell.hpp
    VorFace.hpp
    VorGen.hpp
    VorTess.hpp
    )

    add_executable(${ICMAKERNAME} ${ICMAKERSOURCES})

    set_target_properties(${ICMAKERNAME} PROPERTIES
                          COMPILE_FLAGS "-Dndim_=${NDIM} -DICMAKER")

    target_link_libraries(${ICMAKERNAME} ${IONAME})
    target_link_libraries(${ICMAKERNAME} ${UTILNAME})
    target_link_libraries(${ICMAKERNAME} ${MPI_CXX_LIBRARIES})
    target_link_libraries(${ICMAKERNAME} ${MPI_C_LIBRARIES})
    target_link_libraries(${ICMAKERNAME} ${HDF5_LIBRARIES})
    target_link_libraries(${ICMAKERNAME} ${ALL_Boost_LIBRARIES})
endif(DO_ICMAKER)

if(DO_VTKMAKER)
    set(VTKMAKERNAME vtkmaker${NDIM}d)
    add_executable(${VTKMAKERNAME}
    CubicBox.cpp
    DelTess.cpp
    DelTess_mpi.cpp
    ExArith.cpp
    GadgetSnapshotReader.cpp
    GadgetSnapshotWriter.cpp
    Image.cpp
    Line.cpp
    MPIGlobal.cpp
    NgbSearch.cpp
    Plane.cpp
    ProgramLog.cpp
    RectangularBox.cpp
    Simplex.cpp
    ShadowfaxSnapshotReader.cpp
    ShadowfaxSnapshotWriter.cpp
    SnapshotHandler.cpp
    TimeLine.cpp
    TimeStepWalker.cpp
    VorCell.cpp
    VorFace.cpp
    VorGen.cpp
    VorTess.cpp
    VTKMaker.cpp
    
    MPIGlobal.hpp
    MPIMethods.hpp
    ProgramLog.hpp
    VTKMaker.hpp
    )
    set_target_properties(${VTKMAKERNAME} PROPERTIES
                          COMPILE_FLAGS "-Dndim_=${NDIM}")

    target_link_libraries(${VTKMAKERNAME} ${IONAME})
    target_link_libraries(${VTKMAKERNAME} ${UTILNAME})
    target_link_libraries(${VTKMAKERNAME} ${MPI_CXX_LIBRARIES})
    target_link_libraries(${VTKMAKERNAME} ${MPI_C_LIBRARIES})
    target_link_libraries(${VTKMAKERNAME} ${HDF5_LIBRARIES})
    target_link_libraries(${VTKMAKERNAME} ${ALL_Boost_LIBRARIES})
endif(DO_VTKMAKER)

if(${NDIM} EQUAL 3)
if(DO_GIZMO)
    set(GIZMONAME myGIZMO${NDIM}d)

    add_executable(${GIZMONAME}
    CubicBox.cpp
    DelTess.cpp
    DelTess_mpi.cpp
    ExArith.cpp
    GadgetSnapshotReader.cpp
    GadgetSnapshotWriter.cpp
    GIZMO.cpp
    Image.cpp
    Line.cpp
    MPIGlobal.cpp
    NgbSearch.cpp
    Plane.cpp
    ProgramLog.cpp
    RectangularBox.cpp
    ShadowfaxSnapshotReader.cpp
    ShadowfaxSnapshotWriter.cpp
    Simplex.cpp
    SnapshotHandler.cpp
    TimeLine.cpp
    TimeStepWalker.cpp
    VorCell.cpp
    VorFace.cpp
    VorGen.cpp
    VorTess.cpp

    GIZMO.hpp
    )
    set_target_properties(${GIZMONAME} PROPERTIES
                          COMPILE_FLAGS "-Dndim_=${NDIM}")

    target_link_libraries(${GIZMONAME} ${IONAME})
    target_link_libraries(${GIZMONAME} ${UTILNAME})
    target_link_libraries(${GIZMONAME} ${RIEMANNNAME})
    target_link_libraries(${GIZMONAME} ${MPI_CXX_LIBRARIES})
    target_link_libraries(${GIZMONAME} ${MPI_C_LIBRARIES})
    target_link_libraries(${GIZMONAME} ${HDF5_LIBRARIES})
    target_link_libraries(${GIZMONAME} ${ALL_Boost_LIBRARIES})

endif(DO_GIZMO)
endif(${NDIM} EQUAL 3)

if(DO_TENET)
    set(TENETNAME myTenet${NDIM}d)

    add_executable(${TENETNAME}
    Legendre.cpp
    MPIGlobal.cpp
    ProgramLog.cpp
    Tenet.cpp
    TenetCell.cpp
    utilities/HelperFunctions.cpp
    utilities/Timer.cpp

    Legendre.hpp
    RungeKutta.hpp
    RungeKuttaFlux.hpp
    Tenet.hpp
    TenetBasisFunctions.hpp
    TenetCell.hpp
    TenetCellWeights.hpp
    TenetGrid.hpp
    TenetMinModSlopeLimiter.hpp
    TenetRungeKuttaFlux.hpp
    TenetSlopeLimiter.hpp
    TenetSystemState.hpp
    )
    set_target_properties(${TENETNAME} PROPERTIES
                          COMPILE_FLAGS "-Dndim_=${NDIM}")

    target_link_libraries(${TENETNAME} ${RIEMANNNAME})
    target_link_libraries(${TENETNAME} ${MPI_CXX_LIBRARIES})
    target_link_libraries(${TENETNAME} ${MPI_C_LIBRARIES})

endif(DO_TENET)
