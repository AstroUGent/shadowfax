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

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/data)
set(EWALDTABLE_LOCATION ${CMAKE_BINARY_DIR}/data/)
configure_file(EwaldTableLocation.hpp.in
               ${CMAKE_BINARY_DIR}/inc/EwaldTableLocation.hpp @ONLY)

add_library(${UTILNAME}
DMParticle.cpp
EwaldTable.cpp
GasParticle.cpp
HelperFunctions.cpp
Hilbert.cpp
Particle.cpp
ParticleVector.cpp
StarParticle.cpp
Timer.cpp
Tree.cpp

Box.hpp
Cuboid.hpp
DMParticle.hpp
EwaldTable.hpp
EwaldTableLocation.hpp.in
GasParticle.hpp
HelperFunctions.hpp
Hilbert.hpp
ParallelSorter.hpp
Particle.hpp
ParticleConverter.hpp
ParticleFactory.hpp
ParticleTypes.hpp
ParticleVector.hpp
StarParticle.hpp
Timer.hpp
Tree.hpp
TreeRoute.hpp
TreeWalker.hpp
)
set_target_properties(${UTILNAME} PROPERTIES COMPILE_FLAGS -Dndim_=${NDIM})

include_directories("${PROJECT_SOURCE_DIR}/src")
