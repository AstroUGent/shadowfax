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
 * @file ExArith.hpp
 *
 * @brief Wrapper around Jonathan R. Shewchuk's arbitrary precision predicates
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef HEAD_EXARITH
#define HEAD_EXARITH

#include <cfloat>  // for DBL_MANT_DIG

/**
  * @brief Exact geometric tests
  *
  * A collection of exact geometric tests needed during the Delaunay
  * tesselation. These tests were written by Jonathan Shewchuk. Full credits can
  * be found in the file predicates.cpp
  */
namespace predicates {
/**
 * @brief Initialize some internal variables used by the exact geometric
 * tests
 */
void exactinit();

/**
  * @brief Check whether the points with coordinates pa, pb and pc are
  * counterclockwise or clockwise oriented
  *
  * @param pa,pb,pc Coordinate doublets specifying three distinct points in
  * space
  * @return A positive value if pa, pb and pc are in counterclockwise order,
  * a negative value if they are in clockwise order. A value 0 means the
  * three points are colinear
  */
double orient2d_old(double* pa, double* pb, double* pc);

/**
  * @brief Check whether the point with coordinates pd is above or below the
  * plane through the points with coordinates pa, pb and pc
  *
  * Above and below are defined by the orientation of the points pa, pb and
  * pc. Above is the direction from which pa, pb and pc are seen in
  * counterclockwise order.
  *
  * This function should not be used if pa, pb and pc are colinear, since
  * then there is no unique plane through these three points.
  *
  * @param pa,pb,pc,pd Coordinate triplets specifying four distinct points
  * in space
  * @return A positive value if pd is below the plane (as defined above),
  * negative if it is above the plane. Returns 0 if the four points are
  * coplanar
  */
double orient3d_old(double* pa, double* pb, double* pc, double* pd);

/**
  * @brief Test if the point with coordinates pd is inside the circle
  * through the points with coordinates pa, pb and pc
  *
  * The points pa, pb and pc should be ordered counterclockwise, as defined
  * by a positive value of predicates::orient2d. If not, the result of
  * this function will have the opposite sign.
  *
  * @param pa,pb,pc,pd Coordinate doublets specifying four distinct points
  * in space
  * @return A positive value if pd is inside the circle, a negative value if
  * it is outside. Returns 0 if the four points are cocircular
  */
double incircle_old(double* pa, double* pb, double* pc, double* pd);

/**
  * @brief Test if the point with coordinates pe is inside the sphere
  * through the points with coordinates pa, pb, pc and pd
  *
  * The points pa, pb, pc and pd should be ordered positively, as defined by
  * a positive value of predicates::orient3d. If not, the result of
  * this function will have the opposite sign.
  *
  * @param pa,pb,pc,pd,pe Coordinate triplets specifying five distinct
  * points in space
  * @return A positive value if pe is inside the sphere, a negative value if
  * it is outside. Returns 0 if the five points are cospherical
  */
double insphere_old(double* pa, double* pb, double* pc, double* pd, double* pe);
}

/**
  * @brief Dummy class wrapper around exact geometric tests
  *
  * This class is only used to correctly initialize the exact geometric tests
  * provided by predicates
  */
class PredicatesInitializer {
  public:
    /**
     * @brief Constructor
     *
     * Calls predicates::exactinit()
     */
    PredicatesInitializer() {
        predicates::exactinit();
    }

    ~PredicatesInitializer() {}

    /*! @brief Machine precision, used to discriminate between exact and inexact
     *  tests where needed */
    static double machine_epsilon;
};

#endif
