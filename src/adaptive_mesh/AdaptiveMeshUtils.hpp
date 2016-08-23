/*******************************************************************************
 * This file is part of Shadowfax
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *
 * Shadowfax is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Shadowfax is distributed in the hope that it will be useful,
 * but WITOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with Shadowfax. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

/**
 * @file AdaptiveMeshUtils.hpp
 *
 * @brief Utility functions for the mesh evolution algorithm: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef HEAD_ADAPTIVEMESHUTILS
#define HEAD_ADAPTIVEMESHUTILS

class Cuboid;

/**
 * @brief Utility functions to work with box walls in the evolution algorithm
 */
namespace AdaptiveMeshUtils {
void get_wall_position(double* pos, double* newpos, int wall, Cuboid& box);
int get_wall(double* pos, double* opos);
void get_periodic_position(double* pos, int wall, Cuboid& box);
int get_wallpos(int left, int right);
int trim_wallpos(double* pos, Cuboid& box);
int swap_wallpos(int oldvalue, int wallpos);
int replace_wallpos(int oldvalue, int wallpos);
unsigned long get_key(double* pos, Cuboid& box);
}

#endif
