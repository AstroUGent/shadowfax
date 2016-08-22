/*******************************************************************************
 * This file is part of Shadowfax
 * Copyright (C) 2015 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file ExArith.h
 *
 * @brief Arbitrary exact arithmetics functions
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef HEAD_EXARITH
#define HEAD_EXARITH

#include <cfloat>  // for DBL_MANT_DIG

class Vec;

#define MANTISSAMASK ((long)1 << DBL_MANT_DIG)

/**
 * @brief Union used to convert a double precision floating point value to a
 * 64-bit binary representation
 */
union binaryDouble {
    /*! @brief Double precision floating point value */
    double dval;
    /*! @brief 64-bit binary representation */
    unsigned long bval;
};

/**
  * @brief Exact geometric tests
  *
  * A collection of exact geometric tests needed during the Delaunay
  * tesselation.
  */
namespace ExactArithmetic {
double orient3d_quick(Vec& pa, Vec& pb, Vec& pc, Vec& pd);
double orient3d_exact(Vec& a, Vec& b, Vec& c, Vec& d);
double orient3d_adaptive(Vec& a, Vec& b, Vec& c, Vec& d);

double insphere_quick(Vec& a, Vec& b, Vec& c, Vec& d, Vec& e);
double insphere_exact(Vec& a, Vec& b, Vec& c, Vec& d, Vec& e);
double insphere_adaptive(Vec& a, Vec& b, Vec& c, Vec& d, Vec& e);

double orient2d_quick(Vec& a, Vec& b, Vec& c);
double orient2d_exact(Vec& a, Vec& b, Vec& c);
double orient2d_adaptive(Vec& a, Vec& b, Vec& c);

double incircle_quick(Vec& a, Vec& b, Vec& c, Vec& d);
double incircle_exact(Vec& a, Vec& b, Vec& c, Vec& d);
double incircle_adaptive(Vec& a, Vec& b, Vec& c, Vec& d);

void print_value(const char* str, double val);
double get_value(unsigned long bval);

/**
 * @brief Wrapper around the orient3d test
 *
 * @param a Coordinates of the first point
 * @param b Coordinates of the second point
 * @param c Coordinates of the third point
 * @param d Coordinates of the fourth point
 * @return Value indicating the orientation of the four points
 */
inline double orient3d(Vec& a, Vec& b, Vec& c, Vec& d) {
    return orient3d_adaptive(a, b, c, d);
}

/**
 * @brief Wrapper around the insphere test
 *
 * @param a Coordinates of the first point
 * @param b Coordinates of the second point
 * @param c Coordinates of the third point
 * @param d Coordinates of the fourth point
 * @param e Coordinates of the fifth point
 * @return Value indicating if the fifth point is inside the sphere through
 * the four other ones
 */
inline double insphere(Vec& a, Vec& b, Vec& c, Vec& d, Vec& e) {
    return insphere_adaptive(a, b, c, d, e);
}

/**
 * @brief Wrapper around the 2D orientation test
 *
 * @param a Coordinates of the first point
 * @param b Coordinates of the second point
 * @param c Coordinates of the third point
 * @return Value of the orientation of the three points
 */
inline double orient2d(Vec& a, Vec& b, Vec& c) {
    return orient2d_adaptive(a, b, c);
}

/**
 * @brief Wrapper around the incircle test
 *
 * @param a Coordinates of the first point
 * @param b Coordinates of the second point
 * @param c Coordinates of the third point
 * @param d Coordinates of the fourth point
 * @return Value indicating if the fourth point lies inside the circle
 * through the three other points
 */
inline double incircle(Vec& a, Vec& b, Vec& c, Vec& d) {
    return incircle_adaptive(a, b, c, d);
}

/**
 * @brief Get the 53-bit integer mantissa of a given double precision
 * floating point value
 *
 * @param val Double precision floating point value
 * @return 53-bit long integer mantissa
 */
inline long get_mantissa(double val) {
    // the first expression should in principle be the fastest way to
    // multiply val with 2^DBL_MANT_DIG, however, on my computer, the second
    // version is significantly faster...
    //        return ldexp(val, DBL_MANT_DIG);
    return val * MANTISSAMASK;
}
}

#endif
