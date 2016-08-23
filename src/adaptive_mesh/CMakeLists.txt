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

if(${NDIM} EQUAL 2)

add_library(${ADAPTIVE_MESHNAME}
AdaptiveFace2d.cpp
AdaptiveMeshUtils.cpp
AdaptiveVorCell2d.cpp
AdaptiveVorFace2d.cpp
AdaptiveVorTess2d.cpp
predicates.cpp

AdaptiveCellList.hpp
AdaptiveCellParallelSorter.hpp
AdaptiveFace2d.hpp
AdaptiveMeshException.hpp
AdaptiveMeshUtils.hpp
AdaptiveVorCell2d.hpp
AdaptiveVorFace2d.hpp
AdaptiveVorTess2d.hpp
predicates.hpp
)

else(${NDIM} EQUAL 2)

add_library(${ADAPTIVE_MESHNAME}
AdaptiveFace3d.cpp
AdaptiveMeshUtils.cpp
AdaptiveVorCell3d.cpp
AdaptiveVorFace3d.cpp
AdaptiveVorTess3d.cpp
predicates.cpp

AdaptiveCellList.hpp
AdaptiveFace3d.hpp
AdaptiveMeshException.hpp
AdaptiveMeshUtils.hpp
AdaptiveVorCell3d.hpp
AdaptiveVorFace3d.hpp
AdaptiveVorTess3d.hpp
predicates.hpp
)

endif(${NDIM} EQUAL 2)

set_target_properties(${ADAPTIVE_MESHNAME} PROPERTIES
                      COMPILE_FLAGS -Dndim_=${NDIM})

include_directories("${PROJECT_SOURCE_DIR}/src")
