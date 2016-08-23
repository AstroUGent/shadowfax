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
 * @file AdaptiveMeshUtils.cpp
 *
 * @brief Utility functions for the mesh evolution algorithm: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "AdaptiveMeshUtils.hpp"
#include "Error.hpp"              // for my_exit
#include "Vec.hpp"                // for Vec
#include "utilities/Cuboid.hpp"   // for Cuboid
#include "utilities/Hilbert.hpp"  // for get_key
#include <algorithm>              // for max
#include <iostream>               // for operator<<, etc
using namespace std;

#if ndim_ == 3
/*! @brief Table used to correct wall positions in a cell for the movement of
 *  the cell through a given wall. This table was constructed using the python
 *  script construct_3D_tables.py */
static const int swap_pos_table[27][27] = {
        {0,   -14, -14, -16, -17, -17, -16, -17, -17, -22, -23, -23, -25, -26,
         -25, -26, -26, -22, -23, -23, -25, -26, -26, -25, -26, -26, -26},
        {-13, 0,   -14, -15, -16, -17, -15, -16, -17, -21, -22, -23, -24, -26,
         -24, -25, -26, -21, -22, -23, -24, -25, -26, -24, -25, -26, -25},
        {-13, -13, 0,   -15, -15, -16, -15, -15, -16, -21, -21, -22, -24, -25,
         -24, -24, -25, -21, -21, -22, -24, -24, -25, -24, -24, -25, -24},
        {-11, -12, -12, 0,   -14, -14, -16, -17, -17, -19, -20, -20, -22, -23,
         -25, -26, -26, -19, -20, -20, -22, -23, -23, -25, -26, -26, -23},
        {-10, -11, -12, -13, 0,   -14, -15, -16, -17, -18, -19, -20, -21, -23,
         -24, -25, -26, -18, -19, -20, -21, -22, -23, -24, -25, -26, -22},
        {-10, -10, -11, -13, -13, 0,   -15, -15, -16, -18, -18, -19, -21, -22,
         -24, -24, -25, -18, -18, -19, -21, -21, -22, -24, -24, -25, -21},
        {-11, -12, -12, -11, -12, -12, 0,   -14, -14, -19, -20, -20, -19, -20,
         -22, -23, -23, -19, -20, -20, -19, -20, -20, -22, -23, -23, -20},
        {-10, -11, -12, -10, -11, -12, -13, 0,   -14, -18, -19, -20, -18, -20,
         -21, -22, -23, -18, -19, -20, -18, -19, -20, -21, -22, -23, -19},
        {-10, -10, -11, -10, -10, -11, -13, -13, 0,   -18, -18, -19, -18, -19,
         -21, -21, -22, -18, -18, -19, -18, -18, -19, -21, -21, -22, -18},
        {-5,  -6,  -6,  -8,  -9,  -9,  -8,  -9,  -9,  0,   -14, -14, -16, -17,
         -16, -17, -17, -22, -23, -23, -25, -26, -26, -25, -26, -26, -17},
        {-4,  -5,  -6,  -7,  -8,  -9,  -7,  -8,  -9,  -13, 0,   -14, -15, -17,
         -15, -16, -17, -21, -22, -23, -24, -25, -26, -24, -25, -26, -16},
        {-4,  -4,  -5,  -7,  -7,  -8,  -7,  -7,  -8,  -13, -13, 0,   -15, -16,
         -15, -15, -16, -21, -21, -22, -24, -24, -25, -24, -24, -25, -15},
        {-2,  -3,  -3,  -5,  -6,  -6,  -8,  -9,  -9,  -11, -12, -12, 0,  -14,
         -16, -17, -17, -19, -20, -20, -22, -23, -23, -25, -26, -26, -14},
        {-1,  -1,  -2,  -4,  -4,  -5,  -7,  -7,  -8,  -10, -10, -11, -13, 0,
         -15, -15, -16, -18, -18, -19, -21, -21, -22, -24, -24, -25, -13},
        {-2, -3,  -3,  -2,  -3,  -3,  -5,  -6,  -6,  -11, -12, -12, -11, -12,
         0,  -14, -14, -19, -20, -20, -19, -20, -20, -22, -23, -23, -12},
        {-1,  -2, -3,  -1,  -2,  -3,  -4,  -5,  -6,  -10, -11, -12, -10, -12,
         -13, 0,  -14, -18, -19, -20, -18, -19, -20, -21, -22, -23, -11},
        {-1,  -1,  -2, -1,  -1,  -2,  -4,  -4,  -5,  -10, -10, -11, -10, -11,
         -13, -13, 0,  -18, -18, -19, -18, -18, -19, -21, -21, -22, -10},
        {-5, -6, -6, -8, -9,  -9,  -8,  -9,  -9,  -5,  -6,  -6,  -8, -9,
         -8, -9, -9, 0,  -14, -14, -16, -17, -17, -16, -17, -17, -9},
        {-4, -5, -6, -7,  -8, -9,  -7,  -8,  -9,  -4,  -5,  -6,  -7, -9,
         -7, -8, -9, -13, 0,  -14, -15, -16, -17, -15, -16, -17, -8},
        {-4, -4, -5, -7,  -7,  -8, -7,  -7,  -8,  -4,  -4,  -5,  -7, -8,
         -7, -7, -8, -13, -13, 0,  -15, -15, -16, -15, -15, -16, -7},
        {-2, -3, -3, -5,  -6,  -6,  -8, -9,  -9,  -2,  -3,  -3,  -5, -6,
         -8, -9, -9, -11, -12, -12, 0,  -14, -14, -16, -17, -17, -6},
        {-1, -2, -3, -4,  -5,  -6,  -7,  -8, -9,  -1,  -2,  -3,  -4, -6,
         -7, -8, -9, -10, -11, -12, -13, 0,  -14, -15, -16, -17, -5},
        {-1, -1, -2, -4,  -4,  -5,  -7,  -7,  -8, -1,  -1,  -2,  -4, -5,
         -7, -7, -8, -10, -10, -11, -13, -13, 0,  -15, -15, -16, -4},
        {-2, -3, -3, -2,  -3,  -3,  -5,  -6,  -6,  -2, -3,  -3,  -2, -3,
         -5, -6, -6, -11, -12, -12, -11, -12, -12, 0,  -14, -14, -3},
        {-1, -2, -3, -1,  -2,  -3,  -4,  -5,  -6,  -1,  -2, -3,  -1, -3,
         -4, -5, -6, -10, -11, -12, -10, -11, -12, -13, 0,  -14, -2},
        {-1, -1, -2, -1,  -1,  -2,  -4,  -4,  -5,  -1,  -1,  -2, -1, -2,
         -4, -4, -5, -10, -10, -11, -10, -10, -11, -13, -13, 0,  -1},
        {-1,  -2,  -3,  -4,  -5,  -6,  -7,  -8,  -9,  -10, -11, -12, -13, -14,
         -15, -16, -17, -18, -19, -20, -21, -22, -23, -24, -25, -26, 0}};

/*! @brief Table used to correct wall positions in cells for the movement of one
 *  of their neighbours through a given wall. This table was constructed using
 *  the python script construct_3D_tables.py */
static const int replace_pos_table[27][27] = {
        {-26, -26, -25, -26, -26, -25, -23, -23, -22, -26, -26, -25, -26, -25,
         -23, -23, -22, -17, -17, -16, -17, -17, -16, -14, -14, 0,   -26},
        {-26, -25, -24, -26, -25, -24, -23, -22, -21, -26, -25, -24, -26, -24,
         -23, -22, -21, -17, -16, -15, -17, -16, -15, -14, 0,   -13, -25},
        {-25, -24, -24, -25, -24, -24, -22, -21, -21, -25, -24, -24, -25, -24,
         -22, -21, -21, -16, -15, -15, -16, -15, -15, 0,   -13, -13, -24},
        {-26, -26, -25, -23, -23, -22, -20, -20, -19, -26, -26, -25, -23, -22,
         -20, -20, -19, -17, -17, -16, -14, -14, 0,   -12, -12, -11, -23},
        {-26, -25, -24, -23, -22, -21, -20, -19, -18, -26, -25, -24, -23, -21,
         -20, -19, -18, -17, -16, -15, -14, 0,   -13, -12, -11, -10, -22},
        {-25, -24, -24, -22, -21, -21, -19, -18, -18, -25, -24, -24, -22, -21,
         -19, -18, -18, -16, -15, -15, 0,   -13, -13, -11, -10, -10, -21},
        {-23, -23, -22, -20, -20, -19, -20, -20, -19, -23, -23, -22, -20, -19,
         -20, -20, -19, -14, -14, 0,   -12, -12, -11, -12, -12, -11, -20},
        {-23, -22, -21, -20, -19, -18, -20, -19, -18, -23, -22, -21, -20, -18,
         -20, -19, -18, -14, 0,   -13, -12, -11, -10, -12, -11, -10, -19},
        {-22, -21, -21, -19, -18, -18, -19, -18, -18, -22, -21, -21, -19, -18,
         -19, -18, -18, 0,   -13, -13, -11, -10, -10, -11, -10, -10, -18},
        {-26, -26, -25, -26, -26, -25, -23, -23, -22, -17, -17, -16, -17, -16,
         -14, -14, 0,   -9,  -9,  -8,  -9,  -9,  -8,  -6,  -6,  -5,  -17},
        {-26, -25, -24, -26, -25, -24, -23, -22, -21, -17, -16, -15, -17, -15,
         -14, 0,   -13, -9,  -8,  -7,  -9,  -8,  -7,  -6,  -5,  -4,  -16},
        {-25, -24, -24, -25, -24, -24, -22, -21, -21, -16, -15, -15, -16, -15,
         0,   -13, -13, -8,  -7,  -7,  -8,  -7,  -7,  -5,  -4,  -4,  -15},
        {-26, -26, -25, -23, -23, -22, -20, -20, -19, -17, -17, -16, -14, 0,
         -12, -12, -11, -9,  -9,  -8,  -6,  -6,  -5,  -3,  -3,  -2,  -14},
        {-25, -24, -24, -22, -21, -21, -19, -18, -18, -16, -15, -15, 0,  -13,
         -11, -10, -10, -8,  -7,  -7,  -5,  -4,  -4,  -2,  -1,  -1,  -13},
        {-23, -23, -22, -20, -20, -19, -20, -20, -19, -14, -14, 0,  -12, -11,
         -12, -12, -11, -6,  -6,  -5,  -3,  -3,  -2,  -3,  -3,  -2, -12},
        {-23, -22, -21, -20, -19, -18, -20, -19, -18, -14, 0,  -13, -12, -10,
         -12, -11, -10, -6,  -5,  -4,  -3,  -2,  -1,  -3,  -2, -1,  -11},
        {-22, -21, -21, -19, -18, -18, -19, -18, -18, 0,  -13, -13, -11, -10,
         -11, -10, -10, -5,  -4,  -4,  -2,  -1,  -1,  -2, -1,  -1,  -10},
        {-17, -17, -16, -17, -17, -16, -14, -14, 0,  -9, -9, -8, -9, -8,
         -6,  -6,  -5,  -9,  -9,  -8,  -9,  -9,  -8, -6, -6, -5, -9},
        {-17, -16, -15, -17, -16, -15, -14, 0,  -13, -9, -8, -7, -9, -7,
         -6,  -5,  -4,  -9,  -8,  -7,  -9,  -8, -7,  -6, -5, -4, -8},
        {-16, -15, -15, -16, -15, -15, 0,  -13, -13, -8, -7, -7, -8, -7,
         -5,  -4,  -4,  -8,  -7,  -7,  -8, -7,  -7,  -5, -4, -4, -7},
        {-17, -17, -16, -14, -14, 0,  -12, -12, -11, -9, -9, -8, -6, -5,
         -3,  -3,  -2,  -9,  -9,  -8, -6,  -6,  -5,  -3, -3, -2, -6},
        {-17, -16, -15, -14, 0,  -13, -12, -11, -10, -9, -8, -7, -6, -4,
         -3,  -2,  -1,  -9,  -8, -7,  -6,  -5,  -4,  -3, -2, -1, -5},
        {-16, -15, -15, 0,  -13, -13, -11, -10, -10, -8, -7, -7, -5, -4,
         -2,  -1,  -1,  -8, -7,  -7,  -5,  -4,  -4,  -2, -1, -1, -4},
        {-14, -14, 0,  -12, -12, -11, -12, -12, -11, -6, -6, -5, -3, -2,
         -3,  -3,  -2, -6,  -6,  -5,  -3,  -3,  -2,  -3, -3, -2, -3},
        {-14, 0,  -13, -12, -11, -10, -12, -11, -10, -6, -5, -4, -3, -1,
         -3,  -2, -1,  -6,  -5,  -4,  -3,  -2,  -1,  -3, -2, -1, -2},
        {0,  -13, -13, -11, -10, -10, -11, -10, -10, -5, -4, -4, -2, -1,
         -2, -1,  -1,  -5,  -4,  -4,  -2,  -1,  -1,  -2, -1, -1, -1},
        {-26, -25, -24, -23, -22, -21, -20, -19, -18, -17, -16, -15, -14, -13,
         -12, -11, -10, -9,  -8,  -7,  -6,  -5,  -4,  -3,  -2,  -1,  0}};

/*! @brief Table used to convert two given wall positions to the wallpos of one
 *  w.r.t. the other. This table was constructed using the python script
 *  construct_3D_tables.py */
static const int wallpos_table[27][27] = {
        {0,  -13, -13, -11, -10, -10, -11, -10, -10, -5, -4, -4, -2, -1,
         -2, -1,  -1,  -5,  -4,  -4,  -2,  -1,  -1,  -2, -1, -1, -1},
        {-14, 0,  -13, -12, -11, -10, -12, -11, -10, -6, -5, -4, -3, -1,
         -3,  -2, -1,  -6,  -5,  -4,  -3,  -2,  -1,  -3, -2, -1, -2},
        {-14, -14, 0,  -12, -12, -11, -12, -12, -11, -6, -6, -5, -3, -2,
         -3,  -3,  -2, -6,  -6,  -5,  -3,  -3,  -2,  -3, -3, -2, -3},
        {-16, -15, -15, 0,  -13, -13, -11, -10, -10, -8, -7, -7, -5, -4,
         -2,  -1,  -1,  -8, -7,  -7,  -5,  -4,  -4,  -2, -1, -1, -4},
        {-17, -16, -15, -14, 0,  -13, -12, -11, -10, -9, -8, -7, -6, -4,
         -3,  -2,  -1,  -9,  -8, -7,  -6,  -5,  -4,  -3, -2, -1, -5},
        {-17, -17, -16, -14, -14, 0,  -12, -12, -11, -9, -9, -8, -6, -5,
         -3,  -3,  -2,  -9,  -9,  -8, -6,  -6,  -5,  -3, -3, -2, -6},
        {-16, -15, -15, -16, -15, -15, 0,  -13, -13, -8, -7, -7, -8, -7,
         -5,  -4,  -4,  -8,  -7,  -7,  -8, -7,  -7,  -5, -4, -4, -7},
        {-17, -16, -15, -17, -16, -15, -14, 0,  -13, -9, -8, -7, -9, -7,
         -6,  -5,  -4,  -9,  -8,  -7,  -9,  -8, -7,  -6, -5, -4, -8},
        {-17, -17, -16, -17, -17, -16, -14, -14, 0,  -9, -9, -8, -9, -8,
         -6,  -6,  -5,  -9,  -9,  -8,  -9,  -9,  -8, -6, -6, -5, -9},
        {-22, -21, -21, -19, -18, -18, -19, -18, -18, 0,  -13, -13, -11, -10,
         -11, -10, -10, -5,  -4,  -4,  -2,  -1,  -1,  -2, -1,  -1,  -10},
        {-23, -22, -21, -20, -19, -18, -20, -19, -18, -14, 0,  -13, -12, -10,
         -12, -11, -10, -6,  -5,  -4,  -3,  -2,  -1,  -3,  -2, -1,  -11},
        {-23, -23, -22, -20, -20, -19, -20, -20, -19, -14, -14, 0,  -12, -11,
         -12, -12, -11, -6,  -6,  -5,  -3,  -3,  -2,  -3,  -3,  -2, -12},
        {-25, -24, -24, -22, -21, -21, -19, -18, -18, -16, -15, -15, 0,  -13,
         -11, -10, -10, -8,  -7,  -7,  -5,  -4,  -4,  -2,  -1,  -1,  -13},
        {-26, -26, -25, -23, -23, -22, -20, -20, -19, -17, -17, -16, -14, 0,
         -12, -12, -11, -9,  -9,  -8,  -6,  -6,  -5,  -3,  -3,  -2,  -14},
        {-25, -24, -24, -25, -24, -24, -22, -21, -21, -16, -15, -15, -16, -15,
         0,   -13, -13, -8,  -7,  -7,  -8,  -7,  -7,  -5,  -4,  -4,  -15},
        {-26, -25, -24, -26, -25, -24, -23, -22, -21, -17, -16, -15, -17, -15,
         -14, 0,   -13, -9,  -8,  -7,  -9,  -8,  -7,  -6,  -5,  -4,  -16},
        {-26, -26, -25, -26, -26, -25, -23, -23, -22, -17, -17, -16, -17, -16,
         -14, -14, 0,   -9,  -9,  -8,  -9,  -9,  -8,  -6,  -6,  -5,  -17},
        {-22, -21, -21, -19, -18, -18, -19, -18, -18, -22, -21, -21, -19, -18,
         -19, -18, -18, 0,   -13, -13, -11, -10, -10, -11, -10, -10, -18},
        {-23, -22, -21, -20, -19, -18, -20, -19, -18, -23, -22, -21, -20, -18,
         -20, -19, -18, -14, 0,   -13, -12, -11, -10, -12, -11, -10, -19},
        {-23, -23, -22, -20, -20, -19, -20, -20, -19, -23, -23, -22, -20, -19,
         -20, -20, -19, -14, -14, 0,   -12, -12, -11, -12, -12, -11, -20},
        {-25, -24, -24, -22, -21, -21, -19, -18, -18, -25, -24, -24, -22, -21,
         -19, -18, -18, -16, -15, -15, 0,   -13, -13, -11, -10, -10, -21},
        {-26, -25, -24, -23, -22, -21, -20, -19, -18, -26, -25, -24, -23, -21,
         -20, -19, -18, -17, -16, -15, -14, 0,   -13, -12, -11, -10, -22},
        {-26, -26, -25, -23, -23, -22, -20, -20, -19, -26, -26, -25, -23, -22,
         -20, -20, -19, -17, -17, -16, -14, -14, 0,   -12, -12, -11, -23},
        {-25, -24, -24, -25, -24, -24, -22, -21, -21, -25, -24, -24, -25, -24,
         -22, -21, -21, -16, -15, -15, -16, -15, -15, 0,   -13, -13, -24},
        {-26, -25, -24, -26, -25, -24, -23, -22, -21, -26, -25, -24, -26, -24,
         -23, -22, -21, -17, -16, -15, -17, -16, -15, -14, 0,   -13, -25},
        {-26, -26, -25, -26, -26, -25, -23, -23, -22, -26, -26, -25, -26, -25,
         -23, -23, -22, -17, -17, -16, -17, -17, -16, -14, -14, 0,   -26},
        {-26, -25, -24, -23, -22, -21, -20, -19, -18, -17, -16, -15, -14, -13,
         -12, -11, -10, -9,  -8,  -7,  -6,  -5,  -4,  -3,  -2,  -1,  0}};
#endif

/*
 * Documentation: because an image says more than a thousand words:
 *
 * -5  -1   -6
 *    _____
 * -4 |   | -2
 *    |___|
 * -8  -3   -7
 *
 * For 3D, see the file Walls3D.pdf
 */

/**
 * @brief Get the new position by converting the given position using the given
 * wall index
 *
 * If an unknown wall index is given, the code will abort.
 *
 * @param pos Position to convert
 * @param newpos Array to store the new position in
 * @param wall Wall index
 * @param box Cuboid specifying the dimensions of the simulation box
 */
void AdaptiveMeshUtils::get_wall_position(double* pos, double* newpos, int wall,
                                          Cuboid& box) {
#if ndim_ == 3
    if(wall == -5) {
        newpos[0] = pos[0];
        newpos[1] = pos[1];
        newpos[2] = 2. * (box.get_anchor()[2] + box.get_sides()[2]) - pos[2];
        return;
    }
    if(wall == -14) {
        newpos[0] = 2. * (box.get_anchor()[0] + box.get_sides()[0]) - pos[0];
        newpos[1] = pos[1];
        newpos[2] = pos[2];
        return;
    }
    if(wall == -22) {
        newpos[0] = pos[0];
        newpos[1] = pos[1];
        newpos[2] = 2. * box.get_anchor()[2] - pos[2];
        return;
    }
    if(wall == -13) {
        newpos[0] = 2. * box.get_anchor()[0] - pos[0];
        newpos[1] = pos[1];
        newpos[2] = pos[2];
        return;
    }
    if(wall == -11) {
        newpos[0] = pos[0];
        newpos[1] = 2. * (box.get_anchor()[1] + box.get_sides()[1]) - pos[1];
        newpos[2] = pos[2];
        return;
    }
    if(wall == -16) {
        newpos[0] = pos[0];
        newpos[1] = 2. * box.get_anchor()[1] - pos[1];
        newpos[2] = pos[2];
        return;
    }
    if(wall == -4) {
        newpos[0] = 2. * box.get_anchor()[0] - pos[0];
        newpos[1] = pos[1];
        newpos[2] = 2. * (box.get_anchor()[2] + box.get_sides()[2]) - pos[2];
        return;
    }
    if(wall == -6) {
        newpos[0] = 2. * (box.get_anchor()[0] + box.get_sides()[0]) - pos[0];
        newpos[1] = pos[1];
        newpos[2] = 2. * (box.get_anchor()[2] + box.get_sides()[2]) - pos[2];
        return;
    }
    if(wall == -23) {
        newpos[0] = 2. * (box.get_anchor()[0] + box.get_sides()[0]) - pos[0];
        newpos[1] = pos[1];
        newpos[2] = 2. * box.get_anchor()[2] - pos[2];
        return;
    }
    if(wall == -21) {
        newpos[0] = 2. * box.get_anchor()[0] - pos[0];
        newpos[1] = pos[1];
        newpos[2] = 2. * box.get_anchor()[2] - pos[2];
        return;
    }
    if(wall == -10) {
        newpos[0] = 2. * box.get_anchor()[0] - pos[0];
        newpos[1] = 2. * (box.get_anchor()[1] + box.get_sides()[1]) - pos[1];
        newpos[2] = pos[2];
        return;
    }
    if(wall == -12) {
        newpos[0] = 2. * (box.get_anchor()[0] + box.get_sides()[0]) - pos[0];
        newpos[1] = 2. * (box.get_anchor()[1] + box.get_sides()[1]) - pos[1];
        newpos[2] = pos[2];
        return;
    }
    if(wall == -17) {
        newpos[0] = 2. * (box.get_anchor()[0] + box.get_sides()[0]) - pos[0];
        newpos[1] = 2. * box.get_anchor()[1] - pos[1];
        newpos[2] = pos[2];
        return;
    }
    if(wall == -15) {
        newpos[0] = 2. * box.get_anchor()[0] - pos[0];
        newpos[1] = 2. * box.get_anchor()[1] - pos[1];
        newpos[2] = pos[2];
        return;
    }
    if(wall == -2) {
        newpos[0] = pos[0];
        newpos[1] = 2. * (box.get_anchor()[1] + box.get_sides()[1]) - pos[1];
        newpos[2] = 2. * (box.get_anchor()[2] + box.get_sides()[2]) - pos[2];
        return;
    }
    if(wall == -8) {
        newpos[0] = pos[0];
        newpos[1] = 2. * box.get_anchor()[1] - pos[1];
        newpos[2] = 2. * (box.get_anchor()[2] + box.get_sides()[2]) - pos[2];
        return;
    }
    if(wall == -19) {
        newpos[0] = pos[0];
        newpos[1] = 2. * (box.get_anchor()[1] + box.get_sides()[1]) - pos[1];
        newpos[2] = 2. * box.get_anchor()[2] - pos[2];
        return;
    }
    if(wall == -25) {
        newpos[0] = pos[0];
        newpos[1] = 2. * box.get_anchor()[1] - pos[1];
        newpos[2] = 2. * box.get_anchor()[2] - pos[2];
        return;
    }
    if(wall == -1) {
        newpos[0] = 2. * box.get_anchor()[0] - pos[0];
        newpos[1] = 2. * (box.get_anchor()[1] + box.get_sides()[1]) - pos[1];
        newpos[2] = 2. * (box.get_anchor()[2] + box.get_sides()[2]) - pos[2];
        return;
    }
    if(wall == -3) {
        newpos[0] = 2. * (box.get_anchor()[0] + box.get_sides()[0]) - pos[0];
        newpos[1] = 2. * (box.get_anchor()[1] + box.get_sides()[1]) - pos[1];
        newpos[2] = 2. * (box.get_anchor()[2] + box.get_sides()[2]) - pos[2];
        return;
    }
    if(wall == -9) {
        newpos[0] = 2. * (box.get_anchor()[0] + box.get_sides()[0]) - pos[0];
        newpos[1] = 2. * box.get_anchor()[1] - pos[1];
        newpos[2] = 2. * (box.get_anchor()[2] + box.get_sides()[2]) - pos[2];
        return;
    }
    if(wall == -7) {
        newpos[0] = 2. * box.get_anchor()[0] - pos[0];
        newpos[1] = 2. * box.get_anchor()[1] - pos[1];
        newpos[2] = 2. * (box.get_anchor()[2] + box.get_sides()[2]) - pos[2];
        return;
    }
    if(wall == -18) {
        newpos[0] = 2. * box.get_anchor()[0] - pos[0];
        newpos[1] = 2. * (box.get_anchor()[1] + box.get_sides()[1]) - pos[1];
        newpos[2] = 2. * box.get_anchor()[2] - pos[2];
        return;
    }
    if(wall == -20) {
        newpos[0] = 2. * (box.get_anchor()[0] + box.get_sides()[0]) - pos[0];
        newpos[1] = 2. * (box.get_anchor()[1] + box.get_sides()[1]) - pos[1];
        newpos[2] = 2. * box.get_anchor()[2] - pos[2];
        return;
    }
    if(wall == -26) {
        newpos[0] = 2. * (box.get_anchor()[0] + box.get_sides()[0]) - pos[0];
        newpos[1] = 2. * box.get_anchor()[1] - pos[1];
        newpos[2] = 2. * box.get_anchor()[2] - pos[2];
        return;
    }
    if(wall == -24) {
        newpos[0] = 2. * box.get_anchor()[0] - pos[0];
        newpos[1] = 2. * box.get_anchor()[1] - pos[1];
        newpos[2] = 2. * box.get_anchor()[2] - pos[2];
        return;
    }
    cerr << "Impossible wall case: " << wall << endl;
    my_exit();
#else
    if(wall == -1) {
        newpos[0] = pos[0];
        newpos[1] = 2. * (box.get_anchor()[1] + box.get_sides()[1]) - pos[1];
        return;
    }
    if(wall == -2) {
        newpos[0] = 2. * (box.get_anchor()[0] + box.get_sides()[0]) - pos[0];
        newpos[1] = pos[1];
        return;
    }
    if(wall == -3) {
        newpos[0] = pos[0];
        newpos[1] = 2. * box.get_anchor()[1] - pos[1];
        return;
    }
    if(wall == -4) {
        newpos[0] = 2. * box.get_anchor()[0] - pos[0];
        newpos[1] = pos[1];
        return;
    }
    if(wall == -5) {
        newpos[0] = 2. * box.get_anchor()[0] - pos[0];
        newpos[1] = 2. * (box.get_anchor()[1] + box.get_sides()[1]) - pos[1];
        return;
    }
    if(wall == -6) {
        newpos[0] = 2. * (box.get_anchor()[0] + box.get_sides()[0]) - pos[0];
        newpos[1] = 2. * (box.get_anchor()[1] + box.get_sides()[1]) - pos[1];
        return;
    }
    if(wall == -7) {
        newpos[0] = 2. * (box.get_anchor()[0] + box.get_sides()[0]) - pos[0];
        newpos[1] = 2. * box.get_anchor()[1] - pos[1];
        return;
    }
    if(wall == -8) {
        newpos[0] = 2. * box.get_anchor()[0] - pos[0];
        newpos[1] = 2. * box.get_anchor()[1] - pos[1];
        return;
    }
    // never get here
    cerr << "Error: requested impossible wall!" << endl;
    my_exit();
#endif
}

/**
 * @brief Get the wall index of the wall that would convert the given first
 * position into the given second position
 *
 * @param pos Original position
 * @param opos Converted position
 * @return Wall index that converts the original position into the converted
 * position
 */
int AdaptiveMeshUtils::get_wall(double* pos, double* opos) {
#if ndim_ == 3
    double epsilon = 1.e-10;
    double xdiff = pos[0] - opos[0];
    double ydiff = pos[1] - opos[1];
    double zdiff = pos[2] - opos[2];
    if(ydiff > epsilon) {
        if(xdiff > epsilon) {
            if(zdiff > epsilon) {
                return -3;
            }
            if(zdiff < -epsilon) {
                return -20;
            }
            return -12;
        }
        if(xdiff < -epsilon) {
            if(zdiff > epsilon) {
                return -1;
            }
            if(zdiff < -epsilon) {
                return -18;
            }
            return -10;
        }
        if(zdiff > epsilon) {
            return -2;
        }
        if(zdiff < -epsilon) {
            return -19;
        }
        return -11;
    }
    if(ydiff < -epsilon) {
        if(xdiff > epsilon) {
            if(zdiff > epsilon) {
                return -9;
            }
            if(zdiff < -epsilon) {
                return -26;
            }
            return -17;
        }
        if(xdiff < -epsilon) {
            if(zdiff > epsilon) {
                return -7;
            }
            if(zdiff < -epsilon) {
                return -24;
            }
            return -15;
        }
        if(zdiff > epsilon) {
            return -8;
        }
        if(zdiff < -epsilon) {
            return -25;
        }
        return -16;
    }
    if(xdiff > epsilon) {
        if(zdiff > epsilon) {
            return -6;
        }
        if(zdiff < -epsilon) {
            return -23;
        }
        return -14;
    }
    if(xdiff < -epsilon) {
        if(zdiff > epsilon) {
            return -4;
        }
        if(zdiff < -epsilon) {
            return -21;
        }
        return -13;
    }
    if(zdiff > epsilon) {
        return -5;
    }
    if(zdiff < epsilon) {
        return -22;
    }
    // never get here
    cerr << "Error in wall computation!" << endl;
    my_exit();
    return 0;
#else
    double xdiff, ydiff;
    xdiff = pos[0] - opos[0];
    ydiff = pos[1] - opos[1];
    double epsilon = 1.e-13;
    if(ydiff > epsilon) {
        if(xdiff > epsilon) {
            // upper right corner
            return -6;
        }
        if(xdiff < -epsilon) {
            // upper left corner
            return -5;
        }
        // upper wall
        return -1;
    }
    if(ydiff < -epsilon) {
        if(xdiff > epsilon) {
            // lower right corner
            return -7;
        }
        if(xdiff < -epsilon) {
            // lower left corner
            return -8;
        }
        // lower wall
        return -3;
    }
    if(xdiff > epsilon) {
        // right wall
        return -2;
    }
    if(xdiff < -epsilon) {
        // left wall
        return -4;
    }
    // never get here
    cerr << "Error in wall computation!" << endl;
    my_exit();
    return 0;
#endif
}

/**
 * @brief Convert the given position for a periodic movement through the wall
 * with given index
 *
 * @param pos Position to convert (is changed in place)
 * @param wall Wall index
 * @param box Cuboid specifying the dimensions of the simulation box
 */
void AdaptiveMeshUtils::get_periodic_position(double* pos, int wall,
                                              Cuboid& box) {
#if ndim_ == 3
    if(wall == -5) {
        pos[0] = pos[0];
        pos[1] = pos[1];
        pos[2] = pos[2] + box.get_sides()[2];
        return;
    }
    if(wall == -14) {
        pos[0] = pos[0] + box.get_sides()[0];
        pos[1] = pos[1];
        pos[2] = pos[2];
        return;
    }
    if(wall == -22) {
        pos[0] = pos[0];
        pos[1] = pos[1];
        pos[2] = pos[2] - box.get_sides()[2];
        return;
    }
    if(wall == -13) {
        pos[0] = pos[0] - box.get_sides()[0];
        pos[1] = pos[1];
        pos[2] = pos[2];
        return;
    }
    if(wall == -11) {
        pos[0] = pos[0];
        pos[1] = pos[1] + box.get_sides()[1];
        pos[2] = pos[2];
        return;
    }
    if(wall == -16) {
        pos[0] = pos[0];
        pos[1] = pos[1] - box.get_sides()[1];
        pos[2] = pos[2];
        return;
    }
    if(wall == -4) {
        pos[0] = pos[0] - box.get_sides()[0];
        pos[1] = pos[1];
        pos[2] = pos[2] + box.get_sides()[2];
        return;
    }
    if(wall == -6) {
        pos[0] = pos[0] + box.get_sides()[0];
        pos[1] = pos[1];
        pos[2] = pos[2] + box.get_sides()[2];
        return;
    }
    if(wall == -23) {
        pos[0] = pos[0] + box.get_sides()[0];
        pos[1] = pos[1];
        pos[2] = pos[2] - box.get_sides()[2];
        return;
    }
    if(wall == -21) {
        pos[0] = pos[0] - box.get_sides()[0];
        pos[1] = pos[1];
        pos[2] = pos[2] - box.get_sides()[2];
        return;
    }
    if(wall == -10) {
        pos[0] = pos[0] - box.get_sides()[0];
        pos[1] = pos[1] + box.get_sides()[1];
        pos[2] = pos[2];
        return;
    }
    if(wall == -12) {
        pos[0] = pos[0] + box.get_sides()[0];
        pos[1] = pos[1] + box.get_sides()[1];
        pos[2] = pos[2];
        return;
    }
    if(wall == -17) {
        pos[0] = pos[0] + box.get_sides()[0];
        pos[1] = pos[1] - box.get_sides()[1];
        pos[2] = pos[2];
        return;
    }
    if(wall == -15) {
        pos[0] = pos[0] - box.get_sides()[0];
        pos[1] = pos[1] - box.get_sides()[1];
        pos[2] = pos[2];
        return;
    }
    if(wall == -2) {
        pos[0] = pos[0];
        pos[1] = pos[1] + box.get_sides()[1];
        pos[2] = pos[2] + box.get_sides()[2];
        return;
    }
    if(wall == -8) {
        pos[0] = pos[0];
        pos[1] = pos[1] - box.get_sides()[1];
        pos[2] = pos[2] + box.get_sides()[2];
        return;
    }
    if(wall == -19) {
        pos[0] = pos[0];
        pos[1] = pos[1] + box.get_sides()[1];
        pos[2] = pos[2] - box.get_sides()[2];
        return;
    }
    if(wall == -25) {
        pos[0] = pos[0];
        pos[1] = pos[1] - box.get_sides()[1];
        pos[2] = pos[2] - box.get_sides()[2];
        return;
    }
    if(wall == -1) {
        pos[0] = pos[0] - box.get_sides()[0];
        pos[1] = pos[1] + box.get_sides()[1];
        pos[2] = pos[2] + box.get_sides()[2];
        return;
    }
    if(wall == -3) {
        pos[0] = pos[0] + box.get_sides()[0];
        pos[1] = pos[1] + box.get_sides()[1];
        pos[2] = pos[2] + box.get_sides()[2];
        return;
    }
    if(wall == -9) {
        pos[0] = pos[0] + box.get_sides()[0];
        pos[1] = pos[1] - box.get_sides()[1];
        pos[2] = pos[2] + box.get_sides()[2];
        return;
    }
    if(wall == -7) {
        pos[0] = pos[0] - box.get_sides()[0];
        pos[1] = pos[1] - box.get_sides()[1];
        pos[2] = pos[2] + box.get_sides()[2];
        return;
    }
    if(wall == -18) {
        pos[0] = pos[0] - box.get_sides()[0];
        pos[1] = pos[1] + box.get_sides()[1];
        pos[2] = pos[2] - box.get_sides()[2];
        return;
    }
    if(wall == -20) {
        pos[0] = pos[0] + box.get_sides()[0];
        pos[1] = pos[1] + box.get_sides()[1];
        pos[2] = pos[2] - box.get_sides()[2];
        return;
    }
    if(wall == -26) {
        pos[0] = pos[0] + box.get_sides()[0];
        pos[1] = pos[1] - box.get_sides()[1];
        pos[2] = pos[2] - box.get_sides()[2];
        return;
    }
    if(wall == -24) {
        pos[0] = pos[0] - box.get_sides()[0];
        pos[1] = pos[1] - box.get_sides()[1];
        pos[2] = pos[2] - box.get_sides()[2];
        return;
    }
    cerr << "Impossible wall case: " << wall << endl;
    my_exit();
#else
    if(wall == -1) {
        pos[0] = pos[0];
        pos[1] = pos[1] + box.get_sides()[1];
        return;
    }
    if(wall == -2) {
        pos[0] = pos[0] + box.get_sides()[0];
        pos[1] = pos[1];
        return;
    }
    if(wall == -3) {
        pos[0] = pos[0];
        pos[1] = pos[1] - box.get_sides()[1];
        return;
    }
    if(wall == -4) {
        pos[0] = pos[0] - box.get_sides()[0];
        pos[1] = pos[1];
        return;
    }
    if(wall == -5) {
        pos[0] = pos[0] - box.get_sides()[0];
        pos[1] = pos[1] + box.get_sides()[1];
        return;
    }
    if(wall == -6) {
        pos[0] = pos[0] + box.get_sides()[0];
        pos[1] = pos[1] + box.get_sides()[1];
        return;
    }
    if(wall == -7) {
        pos[0] = pos[0] + box.get_sides()[0];
        pos[1] = pos[1] - box.get_sides()[1];
        return;
    }
    if(wall == -8) {
        pos[0] = pos[0] - box.get_sides()[0];
        pos[1] = pos[1] - box.get_sides()[1];
        return;
    }
    // never get here
    cerr << "Error: requested impossible wall (" << wall << ")!" << endl;
    my_exit();
#endif
}

/**
 * @brief Calculate the wallpos for left given the wallpos of left and right.
 *
 * The wallpos tells us in which region we can locate the periodic copy of a
 * given ngb. If it is 0, we deal with the actual location of the particle.
 * All other values indicate cells that reside on another location and that
 * have to be periodically translated to yield the copy we need.
 *
 * The input values are given with respect to a given cell and hold for two
 * neighbours of this cell that will become mutual neighbours. We here calculate
 * the wallpos of the right neighbour with respect to the left neighbour, to
 * put this value in the internal list of the left neighbour.
 *
 * If both left and right are in the same region (have the same wallpos), then
 * right has wallpos 0 with respect to left. If left has wallpos 0, then right
 * will have the same wallpos w.r.t. left as w.r.t. the original cell. In all
 * other cases, we have to calculate a new wallpos for right w.r.t. left.
 *
 * @param left Wallpos of the receiving neighbour w.r.t. a cell
 * @param right Wallpos of the other neighbour w.r.t. the same cell
 * @return Wallpos of right w.r.t. left.
 */
int AdaptiveMeshUtils::get_wallpos(int left, int right) {
#if ndim_ == 3
    return wallpos_table[left + 26][right + 26];
#else
    if(left == right) {
        return 0;
    }
    if(!left) {
        return right;
    }
    if(!right) {
        if(left == -1) {
            return -3;
        }
        if(left == -2) {
            return -4;
        }
        if(left == -3) {
            return -1;
        }
        if(left == -4) {
            return -2;
        }
        if(left == -5) {
            return -7;
        }
        if(left == -6) {
            return -8;
        }
        if(left == -7) {
            return -5;
        }
        if(left == -8) {
            return -6;
        }
    }
    if(left == -1) {
        if(right == -5) {
            return -4;
        }
        if(right == -6) {
            return -2;
        }
        if(right == -2) {
            return -7;
        }
        if(right == -4) {
            return -8;
        }
    }
    if(left == -2) {
        if(right == -6) {
            return -1;
        }
        if(right == -7) {
            return -3;
        }
        if(right == -1) {
            return -5;
        }
        if(right == -3) {
            return -8;
        }
    }
    if(left == -3) {
        if(right == -7) {
            return -2;
        }
        if(right == -8) {
            return -4;
        }
        if(right == -4) {
            return -5;
        }
        if(right == -2) {
            return -6;
        }
    }
    if(left == -4) {
        if(right == -8) {
            return -3;
        }
        if(right == -5) {
            return -1;
        }
        if(right == -3) {
            return -7;
        }
        if(right == -1) {
            return -6;
        }
    }
    if(left == -5) {
        if(right == -4) {
            return -3;
        }
        if(right == -1) {
            return -2;
        }
    }
    if(left == -6) {
        if(right == -1) {
            return -4;
        }
        if(right == -2) {
            return -3;
        }
    }
    if(left == -7) {
        if(right == -2) {
            return -1;
        }
        if(right == -3) {
            return -4;
        }
    }
    if(left == -8) {
        if(right == -3) {
            return -2;
        }
        if(right == -4) {
            return -1;
        }
    }
    cerr << "get_wallpos: Case not covered (left: " << left
         << ", right: " << right << ")" << endl;
    my_exit();
    return 0;
#endif
}

/**
 * @brief Keep the given position inside the box and return the index of the
 * wall it passes through (if applicable)
 *
 * @param pos Position inside or outside of the simulation box
 * @param box Cuboid specifying the dimensions of the simulation box
 * @return Wall index corresponding to a possible wall position passes through
 * to stay inside the box
 */
int AdaptiveMeshUtils::trim_wallpos(double* pos, Cuboid& box) {
#if ndim_ == 3
    if(pos[0] < 0.) {
        pos[0] += box.get_sides()[0];
        if(pos[1] < 0.) {
            pos[1] += box.get_sides()[1];
            if(pos[2] < 0.) {
                pos[2] += box.get_sides()[2];
                return -24;
            }
            if(pos[2] > 1.) {
                pos[2] -= box.get_sides()[2];
                return -7;
            }
            return -15;
        }
        if(pos[1] > 1.) {
            pos[1] -= box.get_sides()[1];
            if(pos[2] < 0.) {
                pos[2] += box.get_sides()[2];
                return -18;
            }
            if(pos[2] > 1.) {
                pos[2] -= box.get_sides()[2];
                return -1;
            }
            return -10;
        }
        if(pos[2] < 0.) {
            pos[2] += box.get_sides()[2];
            return -21;
        }
        if(pos[2] > 1.) {
            pos[2] -= box.get_sides()[2];
            return -4;
        }
        return -13;
    }
    if(pos[0] > 1.) {
        pos[0] -= box.get_sides()[0];
        if(pos[1] < 0.) {
            pos[1] += box.get_sides()[1];
            if(pos[2] < 0.) {
                pos[2] += box.get_sides()[2];
                return -26;
            }
            if(pos[2] > 1.) {
                pos[2] -= box.get_sides()[2];
                return -9;
            }
            return -17;
        }
        if(pos[1] > 1.) {
            pos[1] -= box.get_sides()[1];
            if(pos[2] < 0.) {
                pos[2] += box.get_sides()[2];
                return -20;
            }
            if(pos[2] > 1.) {
                pos[2] -= box.get_sides()[2];
                return -3;
            }
            return -12;
        }
        if(pos[2] < 0.) {
            pos[2] += box.get_sides()[2];
            return -23;
        }
        if(pos[2] > 1.) {
            pos[2] -= box.get_sides()[2];
            return -6;
        }
        return -14;
    }
    if(pos[1] < 0.) {
        pos[1] += box.get_sides()[1];
        if(pos[2] < 0.) {
            pos[2] += box.get_sides()[2];
            return -25;
        }
        if(pos[2] > 1.) {
            pos[2] -= box.get_sides()[2];
            return -8;
        }
        return -16;
    }
    if(pos[1] > 1.) {
        pos[1] -= box.get_sides()[1];
        if(pos[2] < 0.) {
            pos[2] += box.get_sides()[2];
            return -19;
        }
        if(pos[2] > 1.) {
            pos[2] -= box.get_sides()[2];
            return -2;
        }
        return -11;
    }
    if(pos[2] < 0.) {
        pos[2] += box.get_sides()[2];
        return -22;
    }
    if(pos[2] > 1.) {
        pos[2] -= box.get_sides()[2];
        return -5;
    }
    return 0;
#else
    if(pos[0] < box.get_anchor()[0]) {
        pos[0] += box.get_sides()[0];
        if(pos[1] < box.get_anchor()[1]) {
            pos[1] += box.get_sides()[1];
            return -8;
        }
        if(pos[1] > box.get_anchor()[1] + box.get_sides()[1]) {
            pos[1] -= box.get_sides()[1];
            return -5;
        }
        return -4;
    }
    if(pos[0] > box.get_anchor()[0] + box.get_sides()[0]) {
        pos[0] -= box.get_sides()[0];
        if(pos[1] < box.get_anchor()[1]) {
            pos[1] += box.get_sides()[1];
            return -7;
        }
        if(pos[1] > box.get_anchor()[1] + box.get_sides()[1]) {
            pos[1] -= box.get_sides()[1];
            return -6;
        }
        return -2;
    }
    if(pos[1] < box.get_anchor()[1]) {
        pos[1] += box.get_sides()[1];
        return -3;
    }
    if(pos[1] > box.get_anchor()[1] + box.get_sides()[1]) {
        pos[1] -= box.get_sides()[1];
        return -1;
    }
    return 0;
#endif
}

/**
 * @brief Change the given wall index to a new wall index to correct for the
 * movement of the cell that stores the index through the given wall
 *
 * Coverage might be incomplete for 2D, in which the code will abort.
 *
 * @param oldvalue Old wall index
 * @param wallpos Wall index of the wall the cell passes through
 * @return New wall index
 */
int AdaptiveMeshUtils::swap_wallpos(int oldvalue, int wallpos) {
#if ndim_ == 3
    return swap_pos_table[oldvalue + 26][wallpos + 26];
#else
    if(!oldvalue) {
        if(wallpos == -1) {
            return -3;
        }
        if(wallpos == -2) {
            return -4;
        }
        if(wallpos == -3) {
            return -1;
        }
        if(wallpos == -4) {
            return -2;
        }
        if(wallpos == -5) {
            return -7;
        }
        if(wallpos == -6) {
            return -8;
        }
        if(wallpos == -7) {
            return -5;
        }
        if(wallpos == -8) {
            return -6;
        }
    }
    if(wallpos == oldvalue) {
        return 0;
    }
    if(wallpos == -1) {
        if(oldvalue == -5) {
            return -4;
        }
        if(oldvalue == -6) {
            return -2;
        }
        if(oldvalue == -4) {
            return -8;
        }
        if(oldvalue == -2) {
            return -7;
        }
    }
    if(wallpos == -2) {
        if(oldvalue == -6) {
            return -1;
        }
        if(oldvalue == -7) {
            return -3;
        }
        if(oldvalue == -3) {
            return -8;
        }
        if(oldvalue == -1) {
            return -5;
        }
    }
    if(wallpos == -3) {
        if(oldvalue == -7) {
            return -2;
        }
        if(oldvalue == -8) {
            return -4;
        }
        if(oldvalue == -4) {
            return -5;
        }
        if(oldvalue == -2) {
            return -6;
        }
    }
    if(wallpos == -4) {
        if(oldvalue == -8) {
            return -3;
        }
        if(oldvalue == -5) {
            return -1;
        }
        if(oldvalue == -1) {
            return -6;
        }
        if(oldvalue == -3) {
            return -7;
        }
    }
    if(wallpos == -5) {
        if(oldvalue == -4) {
            return -3;
        }
        if(oldvalue == -1) {
            return -2;
        }
    }
    if(wallpos == -6) {
        if(oldvalue == -1) {
            return -4;
        }
        if(oldvalue == -2) {
            return -3;
        }
    }
    if(wallpos == -7) {
        if(oldvalue == -2) {
            return -1;
        }
        if(oldvalue == -3) {
            return -4;
        }
    }
    if(wallpos == -8) {
        if(oldvalue == -3) {
            return -2;
        }
        if(oldvalue == -4) {
            return -1;
        }
    }
    cerr << "swap_wallpos: case not covered (" << wallpos << "," << oldvalue
         << ")" << endl;
    my_exit();
    return 0;
#endif
}

/**
 * @brief Replace the wall index of a neighbour due to a movement of that
 * neighbour through the wall with given index
 *
 * Coverage might be incomplete in 2D, in which case the code will abort.
 *
 * @param oldvalue Old wall index of the neighbour
 * @param wallpos Wall index of the wall the neighbour passes through
 * @return New wall index of the neighbour
 */
int AdaptiveMeshUtils::replace_wallpos(int oldvalue, int wallpos) {
#if ndim_ == 3
    return replace_pos_table[oldvalue + 26][wallpos + 26];
#else
    if(!oldvalue) {
        return wallpos;
    }
    if(oldvalue == -1) {
        if(wallpos == -3) {
            return 0;
        }
        if(wallpos == -7) {
            return -2;
        }
        if(wallpos == -8) {
            return -4;
        }
        if(wallpos == -2) {
            return -6;
        }
        if(wallpos == -4) {
            return -5;
        }
    }
    if(oldvalue == -2) {
        if(wallpos == -4) {
            return 0;
        }
        if(wallpos == -5) {
            return -1;
        }
        if(wallpos == -8) {
            return -3;
        }
        if(wallpos == -1) {
            return -6;
        }
        if(wallpos == -3) {
            return -7;
        }
    }
    if(oldvalue == -3) {
        if(wallpos == -1) {
            return 0;
        }
        if(wallpos == -5) {
            return -4;
        }
        if(wallpos == -6) {
            return -2;
        }
        if(wallpos == -4) {
            return -8;
        }
        if(wallpos == -2) {
            return -7;
        }
    }
    if(oldvalue == -4) {
        if(wallpos == -2) {
            return 0;
        }
        if(wallpos == -6) {
            return -1;
        }
        if(wallpos == -7) {
            return -3;
        }
        if(wallpos == -3) {
            return -8;
        }
        if(wallpos == -1) {
            return -5;
        }
    }
    if(oldvalue == -5) {
        if(wallpos == -7) {
            return 0;
        }
        if(wallpos == -2) {
            return -1;
        }
        if(wallpos == -3) {
            return -4;
        }
    }
    if(oldvalue == -6) {
        if(wallpos == -8) {
            return 0;
        }
        if(wallpos == -3) {
            return -2;
        }
        if(wallpos == -4) {
            return -1;
        }
    }
    if(oldvalue == -7) {
        if(wallpos == -5) {
            return 0;
        }
        if(wallpos == -4) {
            return -3;
        }
        if(wallpos == -1) {
            return -2;
        }
    }
    if(oldvalue == -8) {
        if(wallpos == -6) {
            return 0;
        }
        if(wallpos == -1) {
            return -4;
        }
        if(wallpos == -2) {
            return -3;
        }
    }
    cerr << "Case not covered:" << endl;
    cerr << oldvalue << "\t" << wallpos << endl;
    my_exit();
    return 0;
#endif
}

/**
 * @brief Get the Hilbert key for the given position inside the given box
 *
 * To mimick the particle sorting as closely as possible, we use the center of
 * the box as a reference position and then add half the maximum sidelength of
 * all sides to that to make sure all positions are inside a cubic box with the
 * same center as the original box.
 *
 * @param pos Position
 * @param box Box
 * @return The Hilbert key for the position in the box
 */
unsigned long AdaptiveMeshUtils::get_key(double* pos, Cuboid& box) {
    unsigned long bits[ndim_] = {0};
    double maxside = std::max(box.get_sides()[0], box.get_sides()[1]);
#if ndim_ == 3
    maxside = std::max(maxside, box.get_sides()[2]);
#endif
    unsigned long bitwidth = 1;
    bitwidth <<= (60 / ndim_);
    for(unsigned int i = ndim_; i--;) {
        bits[i] = ((pos[i] - (box.get_anchor()[i] + 0.5 * box.get_sides()[i]) +
                    0.5 * maxside) /
                   maxside) *
                  bitwidth;
    }
    return HB::get_key(bits, 60);
}
