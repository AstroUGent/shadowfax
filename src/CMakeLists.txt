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

set(IONAME io${NDIM}d)
set(UTILNAME utilities${NDIM}d)

add_subdirectory(io)
add_subdirectory(utilities)

if(DO_SHADOWFAX)
    set(SHADOWFAXNAME shadowfax${NDIM}d)

    set(SHADOWFAXSOURCES
    ParameterFile.cpp
    Plane.cpp
    CircularBox.cpp
    Plotter.cpp
    CubicBox.cpp
    RectangularBox.cpp
    DelTessCheck.cpp
    ExactRiemannSolver.cpp
    DelTess.cpp
    SidePrograms.cpp
    Simplex.cpp
    Simulation.cpp
    SnapshotHandler.cpp
    ShadowfaxSnapshotReader.cpp
    ShadowfaxSnapshotWriter.cpp
    TimeLine.cpp
    VorCell.cpp
    Image.cpp
    VorFace.cpp
    Line.cpp
    VorGen.cpp
    Shadowfax.cpp
    VorTess.cpp
    NgbSearch.cpp
    TimeStepWalker.cpp
    GravityWalker.cpp
    ../inc/GravityWalker.hpp
    ../inc/DelCont.hpp
    ../inc/DelTess.hpp
    ../inc/ExArith.h
    ../inc/Image.hpp
    ../inc/Line.hpp
    ../inc/MPIGlobal.hpp
    MPIGlobal.cpp
    ../inc/MPIMethods.hpp
    ../inc/MPIMethods.hpp
    ../inc/Error.hpp
    ../inc/ParameterFile.hpp
    ../inc/Plane.hpp
    ../inc/Plotter.hpp
    ../inc/RiemannSolver.hpp
    ../inc/ExactRiemannSolver.hpp
    ../inc/SidePrograms.hpp
    ../inc/Simplex.hpp
    ../inc/Simulation.hpp
    ../inc/SnapshotHandler.hpp
    ../inc/StateVector.hpp
    ../inc/TimeLine.hpp
    ../inc/Vec.hpp
    ../inc/VorCell.hpp
    ../inc/VorFace.hpp
    ../inc/VorGen.hpp
    ../inc/VorTess.hpp
    VorTessManager.cpp
    ../inc/VorTessManager.hpp
    ApproximateSolver.cpp
    ../inc/ApproximateSolver.hpp
    ../inc/RiemannSolverFactory.hpp
    GadgetSnapshotReader.cpp
    ../inc/GadgetSnapshotReader.hpp
    GadgetSnapshotWriter.cpp
    ../inc/GadgetSnapshotWriter.hpp
    FixedGrid.cpp
    ../inc/FixedGrid.hpp
    ../inc/TimeStepWalker.hpp
    ../inc/SnapshotWriterFactory.hpp
    ../inc/SnapshotReaderFactory.hpp
    RestartFile.cpp
    ../inc/RestartFile.hpp
    LogFiles.cpp
    ../inc/LogFiles.hpp
    ExArith.cpp
    Physics.cpp
    ../inc/Physics.hpp
    ProgramLog.cpp
    ../inc/ProgramLog.hpp
    )

    # create a new version file if one of the source files changes
    # this way, we will always have the correct compilation time in our
    # restart file
    add_custom_command(OUTPUT ${PROJECT_SOURCE_DIR}/inc/ShadowfaxVersion.hpp
                       COMMAND cmake -P ${PROJECT_SOURCE_DIR}/write_version.cmake
                       WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                       DEPENDS ${SHADOWFAXSOURCES} ${IONAME} ${UTILNAME})

    add_executable(${SHADOWFAXNAME} ${SHADOWFAXSOURCES})
    set_target_properties(${SHADOWFAXNAME} PROPERTIES COMPILE_FLAGS -Dndim_=${NDIM})

    target_link_libraries(${SHADOWFAXNAME} ${IONAME})
    target_link_libraries(${SHADOWFAXNAME} ${UTILNAME})
    target_link_libraries(${SHADOWFAXNAME} ${MPI_CXX_LIBRARIES})
    target_link_libraries(${SHADOWFAXNAME} ${MPI_C_LIBRARIES})
    target_link_libraries(${SHADOWFAXNAME} ${HDF5_LIBRARIES})
    target_link_libraries(${SHADOWFAXNAME} ${ALL_Boost_LIBRARIES})

    # make sure version_number is called, but only once
    if(${NDIM} EQUAL 2 OR ONLY_3D)
        add_custom_target(version_number
                          DEPENDS ${PROJECT_SOURCE_DIR}/inc/ShadowfaxVersion.hpp)
        add_dependencies(${SHADOWFAXNAME} version_number)
    endif(${NDIM} EQUAL 2 OR ONLY_3D)
endif(DO_SHADOWFAX)

if(DO_ICMAKER)
    set(ICMAKERNAME icmaker${NDIM}d)

    set(ICMAKERSOURCES
    BlockICGenerator.cpp
    Plane.cpp
    CircularBox.cpp
    CubicBox.cpp
    RectangularBox.cpp
    DelTessCheck.cpp
    DelTess.cpp
    Simplex.cpp
    SnapshotHandler.cpp
    ShadowfaxSnapshotReader.cpp
    ShadowfaxSnapshotWriter.cpp
    GadgetSnapshotWriter.cpp
    SpecificICGenerator.cpp
    ICMaker.cpp
    ICRegion.cpp
    VorCell.cpp
    VorFace.cpp
    Line.cpp
    VorGen.cpp
    VorTess.cpp
    VorTessManager.cpp
    TimeLine.cpp
    TimeStepWalker.cpp
    NgbSearch.cpp
    ../inc/DelCont.hpp
    ../inc/DelTess.hpp
    ../inc/ExArith.h
    ../inc/ICGenerator.hpp
    ../inc/ICMaker.hpp
    ../inc/ICRegion.hpp
    ../inc/Line.hpp
    ../inc/MPIGlobal.hpp
    MPIGlobal.cpp
    ../inc/MPIMethods.hpp
    ../inc/Plane.hpp
    ../inc/Simplex.hpp
    ../inc/SnapshotHandler.hpp
    ../inc/ShadowfaxSnapshotWriter.hpp
    ../inc/ShadowfaxSnapshotReader.hpp
    ../inc/StateVector.hpp
    ../inc/Vec.hpp
    ../inc/VorCell.hpp
    ../inc/VorFace.hpp
    ../inc/VorGen.hpp
    ../inc/VorTess.hpp
    ../inc/SymbolicFunction.hpp
    ExArith.cpp
    ProgramLog.cpp
    ../inc/ProgramLog.hpp
    )

    add_executable(${ICMAKERNAME} ${ICMAKERSOURCES})

    set_target_properties(${ICMAKERNAME} PROPERTIES COMPILE_FLAGS "-Dndim_=${NDIM} -DICMAKER")

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
    Plane.cpp
    CubicBox.cpp
    RectangularBox.cpp
    DelTess.cpp
    Simplex.cpp
    SnapshotHandler.cpp
    ShadowfaxSnapshotReader.cpp
    GadgetSnapshotReader.cpp
    VorCell.cpp
    VorFace.cpp
    Line.cpp
    ../inc/MPIGlobal.hpp
    MPIGlobal.cpp
    ../inc/MPIMethods.hpp
    VorGen.cpp
    VorTess.cpp
    NgbSearch.cpp
    VTKMaker.cpp
    ../inc/VTKMaker.hpp
    Image.cpp
    TimeLine.cpp
    ShadowfaxSnapshotWriter.cpp
    GadgetSnapshotWriter.cpp
    TimeStepWalker.cpp
    ExArith.cpp
    ProgramLog.cpp
    ../inc/ProgramLog.hpp
    )
    set_target_properties(${VTKMAKERNAME} PROPERTIES COMPILE_FLAGS "-Dndim_=${NDIM}")

    target_link_libraries(${VTKMAKERNAME} ${IONAME})
    target_link_libraries(${VTKMAKERNAME} ${UTILNAME})
    target_link_libraries(${VTKMAKERNAME} ${MPI_CXX_LIBRARIES})
    target_link_libraries(${VTKMAKERNAME} ${MPI_C_LIBRARIES})
    target_link_libraries(${VTKMAKERNAME} ${HDF5_LIBRARIES})
    target_link_libraries(${VTKMAKERNAME} ${ALL_Boost_LIBRARIES})
endif(DO_VTKMAKER)
