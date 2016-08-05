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

add_library(${RIEMANNNAME}
    AdvectionRiemannSolver.cpp
    ApproximateSolver.cpp
    CollelaGlazRiemannSolver.cpp
    ExactRiemannSolver.cpp
    GottliebGrothRiemannSolver.cpp
    HLLCRiemannSolver.cpp
    NoPressureRiemannSolver.cpp
    RiemannSolver.cpp

    AdvectionRiemannSolver.hpp
    ApproximateSolver.hpp
    CollelaGlazRiemannSolver.hpp
    ExactRiemannSolver.hpp
    GottliebGrothRiemannSolver.hpp
    HLLCRiemannSolver.hpp
    NoPressureRiemannSolver.hpp
    RiemannSolver.hpp
    RiemannSolverFactory.hpp
)

set_target_properties(${RIEMANNNAME} PROPERTIES COMPILE_FLAGS -Dndim_=${NDIM})
