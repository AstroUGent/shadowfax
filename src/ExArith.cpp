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
 * @file ExArith.cpp
 *
 * @brief Custom implementations of the arbitrary precision arithmetics
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "ExArith.hpp"
#include "Vec.hpp"                           // for Vec, operator-
#include "utilities/HelperFunctions.hpp"     // for rand_double, sign
#include "utilities/Timer.hpp"               // for Timer
#include <algorithm>                         // for copy
#include <boost/multiprecision/cpp_int.hpp>  // for cpp_int_backend, etc
#include <boost/multiprecision/cpp_int/cpp_int_config.hpp>
#include <boost/multiprecision/detail/no_et_ops.hpp>  // for operator*, etc
#include <boost/multiprecision/detail/number_base.hpp>
#include <boost/multiprecision/detail/number_compare.hpp>
#include <boost/multiprecision/number.hpp>  // for number
#include <cmath>                            // for fabs
#include <iostream>  // for operator<<, cout, ostream, etc
#include <stdlib.h>  // for srand
// if you want to use this, you have to link to the GMP library, which is not
// done in the basic CMake setup
//#include <boost/multiprecision/gmp.hpp>
using namespace std;

// it turns out the predefined int256_t is faster than a type with the correct
// size
/**
 * @brief Big integer used for the exact orient3D test
 *
 * Detailed breakdown is done in the function. Upper bound is 165, 256 is fast
 * enough.
 */
typedef boost::multiprecision::int256_t int_orient3d;
// tests indicate that using boost::multiprecision::int256_t is more than 3
// times faster than using an unlimited mpz_int
// typedef boost::multiprecision::mpz_int int_orient3d;

/**
 * @brief Big integer used for the exact insphere test
 *
 * Breakdown in bits is done in the function itself. 277 is a safe upper bound.
 * 512 is too slow, so we define a custom type.
 */
typedef boost::multiprecision::number<boost::multiprecision::cpp_int_backend<
        278, 278, boost::multiprecision::signed_magnitude,
        boost::multiprecision::unchecked, void> >
        int_insphere;

/**
 * @brief Big integer used for exact orient2D test
 *
 * The orientation test is the difference of two multiplications, the factors of
 * which are differences themselves. For 53-bit terms, each factor is at most 54
 * bits. Each term in the overall difference is then at most 108 bits, so that
 * the result can be at most 109 bits.
 *
 * Detailed breakdown is done in the function. Upper bound is 109, 128 is fast
 * enough.
 */
typedef boost::multiprecision::int128_t int_orient2d;

/**
 * @brief Big integer used for exact incircle test
 *
 * Detailed breakdown of the upper bounds on lengths is done in the function.
 * Upper bound is 220, 256 is fast enough.
 */
typedef boost::multiprecision::int256_t int_incircle;

#if ndim_ == 3
/**
 * @brief Quick 3D orientation test
 *
 * Returns a positive value if the fourth point is below the plane through the
 * three other points, a negative value if it is above and 0 if the four points
 * are coplanar.
 *
 * Above is defined as the direction from which the three points in the plane
 * are seen in counterclockwise order, e.g. if the points are A (0,0,0), B
 * (0,0,1), C (0,1,0) and D (1,0,0), then this function returns 1.
 *
 * @warning Might fail due to roundoff error. Use orient3d_adaptive for a test
 * that always gives the correct sign.
 *
 * @param a Coordinates of the first point
 * @param b Coordinates of the second point
 * @param c Coordinates of the third point
 * @param d Coordinates of the fourth point
 * @return A positive, negative or zero value, depending on the outcome of the
 * test
 */
double ExactArithmetic::orient3d_quick(Vec& a, Vec& b, Vec& c, Vec& d) {
    Vec ad = a - d;
    Vec bd = b - d;
    Vec cd = c - d;

    return ad.z() * (bd.x() * cd.y() - cd.x() * bd.y()) +
           bd.z() * (cd.x() * ad.y() - ad.x() * cd.y()) +
           cd.z() * (ad.x() * bd.y() - bd.x() * ad.y());
}

/**
 * @brief Exact 3D orientation test
 *
 * Returns 1 if the fourth point is below the plane through the three other
 * points, a -1 if it is above and 0 if the four points are coplanar.
 *
 * Above is defined as the direction from which the three points in the plane
 * are seen in counterclockwise order, e.g. if the points are A (0,0,0), B
 * (0,0,1), C (0,1,0) and D (1,0,0), then this function returns 1.
 *
 * This function translates the given coordinates to 53-bit integers and then
 * uses 165-bit arithmetic to exactly determine the sign of the determinant.
 * For this to work, all coordinates that are passed on to this function have
 * to lie in the interval [1,2].
 *
 * @param a Coordinates of the first point, have to be in the interval [1,2]
 * @param b Coordinates of the second point, have to be in the interval [1,2]
 * @param c Coordinates of the third point, have to be in the interval [1,2]
 * @param d Coordinates of the fourth point, have to be in the interval [1,2]
 * @return 1, -1 or 0, depending on the outcome of the test
 */
double ExactArithmetic::orient3d_exact(Vec& a, Vec& b, Vec& c, Vec& d) {
    // 53-bit numbers
    int_orient3d axp = get_mantissa(a.x());
    int_orient3d ayp = get_mantissa(a.y());
    int_orient3d azp = get_mantissa(a.z());

    int_orient3d bxp = get_mantissa(b.x());
    int_orient3d byp = get_mantissa(b.y());
    int_orient3d bzp = get_mantissa(b.z());

    int_orient3d cxp = get_mantissa(c.x());
    int_orient3d cyp = get_mantissa(c.y());
    int_orient3d czp = get_mantissa(c.z());

    int_orient3d dxp = get_mantissa(d.x());
    int_orient3d dyp = get_mantissa(d.y());
    int_orient3d dzp = get_mantissa(d.z());

    // differences are at most 54 bits
    int_orient3d adx = axp - dxp;
    int_orient3d ady = ayp - dyp;
    int_orient3d adz = azp - dzp;

    int_orient3d bdx = bxp - dxp;
    int_orient3d bdy = byp - dyp;
    int_orient3d bdz = bzp - dzp;

    int_orient3d cdx = cxp - dxp;
    int_orient3d cdy = cyp - dyp;
    int_orient3d cdz = czp - dzp;

    // multiplications are at most 108 bits
    int_orient3d bdxcdy = bdx * cdy;
    int_orient3d cdxbdy = cdx * bdy;

    int_orient3d cdxady = cdx * ady;
    int_orient3d adxcdy = adx * cdy;

    int_orient3d adxbdy = adx * bdy;
    int_orient3d bdxady = bdx * ady;

    // the factors between brackets are at most 109 bits
    // the terms are at most 163 bits
    // the total sum is at most 165 bits
    int_orient3d result = adz * (bdxcdy - cdxbdy) + bdz * (cdxady - adxcdy) +
                          cdz * (adxbdy - bdxady);

    return HelperFunctions::sign(result);
}

/**
 * @brief Adaptive 3D orientation test
 *
 * Returns a positive value if the fourth point is below the plane through the
 * three other points, a negative value if it is above and 0 if the four points
 * are coplanar.
 *
 * Above is defined as the direction from which the three points in the plane
 * are seen in counterclockwise order, e.g. if the points are A (0,0,0), B
 * (0,0,1), C (0,1,0) and D (1,0,0), then this function returns 1.
 *
 * This function calculates a maximal error bound on the result due to numerical
 * roundoff error. If the result is smaller than this error bound, the function
 * calls orient3d_exact to determine the sign of the determinant using exact
 * integer arithmetics. For this to work, the coordinates passed on to this
 * function have to lie in the interval [1,2].
 *
 * @param a Coordinates of the first point, have to be in the interval [1,2]
 * @param b Coordinates of the second point, have to be in the interval [1,2]
 * @param c Coordinates of the third point, have to be in the interval [1,2]
 * @param d Coordinates of the fourth point, have to be in the interval [1,2]
 * @return A positive, negative or zero value, depending on the outcome of the
 * test
 */
double ExactArithmetic::orient3d_adaptive(Vec& a, Vec& b, Vec& c, Vec& d) {
    Vec ad = a - d;
    Vec bd = b - d;
    Vec cd = c - d;

    double bdxcdy = bd.x() * cd.y();
    double cdxbdy = cd.x() * bd.y();

    double cdxady = cd.x() * ad.y();
    double adxcdy = ad.x() * cd.y();

    double adxbdy = ad.x() * bd.y();
    double bdxady = bd.x() * ad.y();

    double errbound = (fabs(bdxcdy) + fabs(cdxbdy)) * fabs(ad.z()) +
                      (fabs(cdxady) + fabs(adxcdy)) * fabs(bd.z()) +
                      (fabs(adxbdy) + fabs(bdxady)) * fabs(cd.z());
    // not really the right factor (which is 7.77156e-16 on my local machine),
    // but this will do
    errbound *= 1.e-10;

    double result = ad.z() * (bdxcdy - cdxbdy) + bd.z() * (cdxady - adxcdy) +
                    cd.z() * (adxbdy - bdxady);
    if(result < -errbound || result > errbound) {
        return result;
    }

    return orient3d_exact(a, b, c, d);
}

/**
 * @brief Quick in sphere test
 *
 * Returns a positive value if the fifth point is inside the sphere through the
 * first four points, a negative value if it is outside, and 0 if it is on the
 * sphere.
 *
 * This function only works if the first four points are positively oriented, as
 * defined by a positive return value of orient3d_quick, orient3d_exact or
 * orient3d_adaptive. If the points are negatively oriented, the sign of this
 * test flips and a positive value denotes a point outside the sphere.
 * If the four points are coplanar, then this test will result in undefined
 * behaviour.
 *
 * @warning Might fail due to roundoff error. Use insphere_adaptive for a test
 * that always gives the correct sign.
 *
 * @param a Coordinates of the first point
 * @param b Coordinates of the second point
 * @param c Coordinates of the third point
 * @param d Coordinates of the fourth point
 * @param e Coordinates of the fifth point
 * @return A positive, negative or zero value, depending on the outcome of the
 * test
 */
double ExactArithmetic::insphere_quick(Vec& a, Vec& b, Vec& c, Vec& d, Vec& e) {
    Vec ae = a - e;
    Vec be = b - e;
    Vec ce = c - e;
    Vec de = d - e;

    double ab = ae.x() * be.y() - be.x() * ae.y();
    double bc = be.x() * ce.y() - ce.x() * be.y();
    double cd = ce.x() * de.y() - de.x() * ce.y();
    double da = de.x() * ae.y() - ae.x() * de.y();
    double ac = ae.x() * ce.y() - ce.x() * ae.y();
    double bd = be.x() * de.y() - de.x() * be.y();

    double abc = ae.z() * bc - be.z() * ac + ce.z() * ab;
    double bcd = be.z() * cd - ce.z() * bd + de.z() * bc;
    double cda = ce.z() * da + de.z() * ac + ae.z() * cd;
    double dab = de.z() * ab + ae.z() * bd + be.z() * da;

    return (de.norm2() * abc - ce.norm2() * dab) +
           (be.norm2() * cda - ae.norm2() * bcd);
}

/**
 * @brief Exact in sphere test
 *
 * Returns 1 if the fifth point is inside the sphere through the first four
 * points, a -1 if it is outside, and 0 if it is on the sphere.
 *
 * This function only works if the first four points are positively oriented, as
 * defined by a positive return value of orient3d_quick, orient3d_exact or
 * orient3d_adaptive. If the points are negatively oriented, the sign of this
 * test flips and a positive value denotes a point outside the sphere.
 * If the four points are coplanar, then this test will result in undefined
 * behaviour.
 *
 * This function translates the given coordinates to 53-bit integers and then
 * uses 277-bit arithmetic to exactly determine the sign of the determinant.
 * For this to work, all coordinates that are passed on to this function have
 * to lie in the interval [1,2].
 *
 * @param a Coordinates of the first point, have to be in the interval [1,2]
 * @param b Coordinates of the second point, have to be in the interval [1,2]
 * @param c Coordinates of the third point, have to be in the interval [1,2]
 * @param d Coordinates of the fourth point, have to be in the interval [1,2]
 * @param e Coordinates of the fifth point, have to be in the interval [1,2]
 * @return 1, -1 or 0 value, depending on the outcome of the test
 */
double ExactArithmetic::insphere_exact(Vec& a, Vec& b, Vec& c, Vec& d, Vec& e) {
    // 53-bit numbers
    int_insphere axp = get_mantissa(a.x());
    int_insphere ayp = get_mantissa(a.y());
    int_insphere azp = get_mantissa(a.z());

    int_insphere bxp = get_mantissa(b.x());
    int_insphere byp = get_mantissa(b.y());
    int_insphere bzp = get_mantissa(b.z());

    int_insphere cxp = get_mantissa(c.x());
    int_insphere cyp = get_mantissa(c.y());
    int_insphere czp = get_mantissa(c.z());

    int_insphere dxp = get_mantissa(d.x());
    int_insphere dyp = get_mantissa(d.y());
    int_insphere dzp = get_mantissa(d.z());

    int_insphere exp = get_mantissa(e.x());
    int_insphere eyp = get_mantissa(e.y());
    int_insphere ezp = get_mantissa(e.z());

    // at most 54 significant bits
    int_insphere aex = axp - exp;
    int_insphere aey = ayp - eyp;
    int_insphere aez = azp - ezp;
    int_insphere bex = bxp - exp;
    int_insphere bey = byp - eyp;
    int_insphere bez = bzp - ezp;
    int_insphere cex = cxp - exp;
    int_insphere cey = cyp - eyp;
    int_insphere cez = czp - ezp;
    int_insphere dex = dxp - exp;
    int_insphere dey = dyp - eyp;
    int_insphere dez = dzp - ezp;

    // every term has at most 108 significant bits
    // the difference has at most 109 significant bits
    int_insphere ab = aex * bey - bex * aey;
    int_insphere bc = bex * cey - cex * bey;
    int_insphere cd = cex * dey - dex * cey;
    int_insphere da = dex * aey - aex * dey;
    int_insphere ac = aex * cey - cex * aey;
    int_insphere bd = bex * dey - dex * bey;

    // every term has at most 163 significant bits
    // the total sum has at most 165 significant bits
    int_insphere abc = aez * bc - bez * ac + cez * ab;
    int_insphere bcd = bez * cd - cez * bd + dez * bc;
    int_insphere cda = cez * da + dez * ac + aez * cd;
    int_insphere dab = dez * ab + aez * bd + bez * da;

    // every term has at most 108 significant bits
    // the total sum has at most 110 significant bits
    int_insphere aenrm2 = aex * aex + aey * aey + aez * aez;
    int_insphere benrm2 = bex * bex + bey * bey + bez * bez;
    int_insphere cenrm2 = cex * cex + cey * cey + cez * cez;
    int_insphere denrm2 = dex * dex + dey * dey + dez * dez;

    // every term has at most 275 significant bits
    // the total sum has at most 277 significant bits
    int_insphere result =
            (denrm2 * abc - cenrm2 * dab) + (benrm2 * cda - aenrm2 * bcd);

    return HelperFunctions::sign(result);
}

/**
 * @brief Adaptive in sphere test
 *
 * Returns a positive value if the fifth point is inside the sphere through the
 * first four points, a negative value if it is outside, and 0 if it is on the
 * sphere.
 *
 * This function only works if the first four points are positively oriented, as
 * defined by a positive return value of orient3d_quick, orient3d_exact or
 * orient3d_adaptive. If the points are negatively oriented, the sign of this
 * test flips and a positive value denotes a point outside the sphere.
 * If the four points are coplanar, then this test will result in undefined
 * behaviour.
 *
 * This function calculates a maximal error bound on the result due to numerical
 * roundoff error. If the result is smaller than this error bound, the function
 * calls insphere_exact to determine the sign of the determinant using exact
 * integer arithmetics. For this to work, the coordinates passed on to this
 * function have to lie in the interval [1,2].
 *
 * @param a Coordinates of the first point, have to be in the interval [1,2]
 * @param b Coordinates of the second point, have to be in the interval [1,2]
 * @param c Coordinates of the third point, have to be in the interval [1,2]
 * @param d Coordinates of the fourth point, have to be in the interval [1,2]
 * @param e Coordinates of the fifth point, have to be in the interval [1,2]
 * @return A positive, negative or zero value, depending on the outcome of the
 * test
 */
double ExactArithmetic::insphere_adaptive(Vec& a, Vec& b, Vec& c, Vec& d,
                                          Vec& e) {
    Vec ae = a - e;
    Vec be = b - e;
    Vec ce = c - e;
    Vec de = d - e;

    double aexbey = ae.x() * be.y();
    double bexaey = be.x() * ae.y();
    double ab = aexbey - bexaey;
    double bexcey = be.x() * ce.y();
    double cexbey = ce.x() * be.y();
    double bc = bexcey - cexbey;
    double cexdey = ce.x() * de.y();
    double dexcey = de.x() * ce.y();
    double cd = cexdey - dexcey;
    double dexaey = de.x() * ae.y();
    double aexdey = ae.x() * de.y();
    double da = dexaey - aexdey;
    double aexcey = ae.x() * ce.y();
    double cexaey = ce.x() * ae.y();
    double ac = aexcey - cexaey;
    double bexdey = be.x() * de.y();
    double dexbey = de.x() * be.y();
    double bd = bexdey - dexbey;

    double abc = ae.z() * bc - be.z() * ac + ce.z() * ab;
    double bcd = be.z() * cd - ce.z() * bd + de.z() * bc;
    double cda = ce.z() * da + de.z() * ac + ae.z() * cd;
    double dab = de.z() * ab + ae.z() * bd + be.z() * da;

    double aenrm2 = ae.norm2();
    double benrm2 = be.norm2();
    double cenrm2 = ce.norm2();
    double denrm2 = de.norm2();

    double aezplus = fabs(ae.z());
    double bezplus = fabs(be.z());
    double cezplus = fabs(ce.z());
    double dezplus = fabs(de.z());
    double aexbeyplus = fabs(aexbey);
    double bexaeyplus = fabs(bexaey);
    double bexceyplus = fabs(bexcey);
    double cexbeyplus = fabs(cexbey);
    double cexdeyplus = fabs(cexdey);
    double dexceyplus = fabs(dexcey);
    double dexaeyplus = fabs(dexaey);
    double aexdeyplus = fabs(aexdey);
    double aexceyplus = fabs(aexcey);
    double cexaeyplus = fabs(cexaey);
    double bexdeyplus = fabs(bexdey);
    double dexbeyplus = fabs(dexbey);
    double errbound = ((cexdeyplus + dexceyplus) * bezplus +
                       (dexbeyplus + bexdeyplus) * cezplus +
                       (bexceyplus + cexbeyplus) * dezplus) *
                              aenrm2 +
                      ((dexaeyplus + aexdeyplus) * cezplus +
                       (aexceyplus + cexaeyplus) * dezplus +
                       (cexdeyplus + dexceyplus) * aezplus) *
                              benrm2 +
                      ((aexbeyplus + bexaeyplus) * dezplus +
                       (bexdeyplus + dexbeyplus) * aezplus +
                       (dexaeyplus + aexdeyplus) * bezplus) *
                              cenrm2 +
                      ((bexceyplus + cexbeyplus) * aezplus +
                       (cexaeyplus + aexceyplus) * bezplus +
                       (aexbeyplus + bexaeyplus) * cezplus) *
                              denrm2;
    // not really the right factor (which is 1.77636e-15 on my local machine),
    // but this will do
    errbound *= 1.e-10;

    double result =
            (denrm2 * abc - cenrm2 * dab) + (benrm2 * cda - aenrm2 * bcd);

    if(result < -errbound || result > errbound) {
        return result;
    }

    return insphere_exact(a, b, c, d, e);
}
#endif

#if ndim_ == 2
/**
 * @brief Quick 2D orientation test
 *
 * Returns a positive value if the three given points are counterclockwise
 * ordered, a negative value if they are in clockwise order, and zero if the
 * three points are colinear.
 *
 * This is assuming a positively oriented 2D cartesian reference frame.
 *
 * @warning Might fail due to roundoff error. Use orient2d_adaptive for a test
 * that always gives the correct sign.
 *
 * @param a Coordinates of the first point
 * @param b Coordinates of the second point
 * @param c Coordinates of the third point
 * @return A positive, negative or zero value, depending on the outcome of the
 * test
 */
double ExactArithmetic::orient2d_quick(Vec& a, Vec& b, Vec& c) {
    return (a.x() - c.x()) * (b.y() - c.y()) -
           (a.y() - c.y()) * (b.x() - c.x());
}

/**
 * @brief Exact 2D orientation test
 *
 * Returns 1 if the three given points are counterclockwise ordered, -1 if they
 * are in clockwise order, and 0 if the three points are colinear.
 *
 * This is assuming a positively oriented 2D cartesian reference frame.
 *
 * This function translates the given coordinates to 53-bit integers and then
 * uses 109-bit arithmetic to exactly determine the sign of the determinant.
 * For this to work, all coordinates that are passed on to this function have
 * to lie in the interval [1,2].
 *
 * @param a Coordinates of the first point, have to be in the interval [1,2]
 * @param b Coordinates of the second point, have to be in the interval [1,2]
 * @param c Coordinates of the third point, have to be in the interval [1,2]
 * @return 1, -1 or 0, depending on the outcome of the test
 */
double ExactArithmetic::orient2d_exact(Vec& a, Vec& b, Vec& c) {
    // 53-bit numbers
    int_orient2d axp = get_mantissa(a.x());
    int_orient2d ayp = get_mantissa(a.y());

    int_orient2d bxp = get_mantissa(b.x());
    int_orient2d byp = get_mantissa(b.y());

    int_orient2d cxp = get_mantissa(c.x());
    int_orient2d cyp = get_mantissa(c.y());

    // every factor between brackets is at most 54 bits
    // every term is at most 108 bits
    // the difference is at most 109 bits
    int_orient2d result = (axp - cxp) * (byp - cyp) - (ayp - cyp) * (bxp - cxp);

    return HelperFunctions::sign(result);
}

/**
 * @brief Adaptive 2D orientation test
 *
 * Returns a positive value if the three given points are counterclockwise
 * ordered, a negative value if they are in clockwise order, and zero if the
 * three points are colinear.
 *
 * This is assuming a positively oriented 2D cartesian reference frame.
 *
 * This function calculates a maximal error bound on the result due to numerical
 * roundoff error. If the result is smaller than this error bound, the function
 * calls orient2d_exact to determine the sign of the determinant using exact
 * integer arithmetics. For this to work, the coordinates passed on to this
 * function have to lie in the interval [1,2].
 *
 * @param a Coordinates of the first point, have to be in the interval [1,2]
 * @param b Coordinates of the second point, have to be in the interval [1,2]
 * @param c Coordinates of the third point, have to be in the interval [1,2]
 * @return A positive, negative or zero value, depending on the outcome of the
 * test
 */
double ExactArithmetic::orient2d_adaptive(Vec& a, Vec& b, Vec& c) {
    double detleft = (a.x() - c.x()) * (b.y() - c.y());
    double detright = (a.y() - c.y()) * (b.x() - c.x());

    double errbound = 0.;
    double result = detleft - detright;

    if(detleft > 0.0) {
        if(detright <= 0.0) {
            return result;
        } else {
            errbound = detleft + detright;
        }
    } else if(detleft < 0.0) {
        if(detright >= 0.0) {
            return result;
        } else {
            errbound = -detleft - detright;
        }
    } else {
        return result;
    }

    // actual value is smaller (3.33067e-16 on my machine), but this will do
    errbound *= 1.e-10;

    if((result >= errbound) || (-result >= errbound)) {
        return result;
    }

    return orient2d_exact(a, b, c);
}

/**
 * @brief Quick in circle test
 *
 * Returns a positive value if the fourth point is inside the circle through the
 * first three points, a negative value if it is outside, and zero if it is on
 * the circle.
 *
 * This function only works if the first three points are positively oriented,
 * as defined by a positive return value of orient2d_quick, orient2d_exact or
 * orient2d_adaptive. If the points are negatively oriented, the sign of this
 * test flips and a positive value denotes a point outside the circle.
 * If the three points are colinear, then this test will result in undefined
 * behaviour.
 *
 * @warning Might fail due to roundoff error. Use incircle_adaptive for a test
 * that always gives the correct sign.
 *
 * @param a Coordinates of the first point
 * @param b Coordinates of the second point
 * @param c Coordinates of the third point
 * @param d Coordinates of the fourth point
 * @return A positive, negative or zero value, depending on the outcome of the
 * test
 */
double ExactArithmetic::incircle_quick(Vec& a, Vec& b, Vec& c, Vec& d) {
    Vec ad = a - d;
    Vec bd = b - d;
    Vec cd = c - d;

    double bdxcdy = bd.x() * cd.y();
    double cdxbdy = cd.x() * bd.y();
    double adnrm2 = ad.norm2();

    double cdxady = cd.x() * ad.y();
    double adxcdy = ad.x() * cd.y();
    double bdnrm2 = bd.norm2();

    double adxbdy = ad.x() * bd.y();
    double bdxady = bd.x() * ad.y();
    double cdnrm2 = cd.norm2();

    double result = adnrm2 * (bdxcdy - cdxbdy) + bdnrm2 * (cdxady - adxcdy) +
                    cdnrm2 * (adxbdy - bdxady);

    return result;
}

/**
 * @brief Exact in circle test
 *
 * Returns 1 if the fourth point is inside the circle through the first three
 * points, -1 if it is outside, and 0 if it is on the circle.
 *
 * This function only works if the first three points are positively oriented,
 * as defined by a positive return value of orient2d_quick, orient2d_exact or
 * orient2d_adaptive. If the points are negatively oriented, the sign of this
 * test flips and a positive value denotes a point outside the circle.
 * If the three points are colinear, then this test will result in undefined
 * behaviour.
 *
 * This function translates the given coordinates to 53-bit integers and then
 * uses 220-bit arithmetic to exactly determine the sign of the determinant.
 * For this to work, all coordinates that are passed on to this function have
 * to lie in the interval [1,2].
 *
 * @param a Coordinates of the first point, have to be in the interval [1,2]
 * @param b Coordinates of the second point, have to be in the interval [1,2]
 * @param c Coordinates of the third point, have to be in the interval [1,2]
 * @param d Coordinates of the fourth point, have to be in the interval [1,2]
 * @return 1, -1 or 0, depending on the outcome of the test
 */
double ExactArithmetic::incircle_exact(Vec& a, Vec& b, Vec& c, Vec& d) {
    // 53-bit numbers
    int_incircle axp = get_mantissa(a.x());
    int_incircle ayp = get_mantissa(a.y());

    int_incircle bxp = get_mantissa(b.x());
    int_incircle byp = get_mantissa(b.y());

    int_incircle cxp = get_mantissa(c.x());
    int_incircle cyp = get_mantissa(c.y());

    int_incircle dxp = get_mantissa(d.x());
    int_incircle dyp = get_mantissa(d.y());

    // every difference is at most 54 bits
    int_incircle adx = axp - dxp;
    int_incircle ady = ayp - dyp;
    int_incircle bdx = bxp - dxp;
    int_incircle bdy = byp - dyp;
    int_incircle cdx = cxp - dxp;
    int_incircle cdy = cyp - dyp;

    // multiplication is at most 108 bits
    int_incircle bdxcdy = bdx * cdy;
    int_incircle cdxbdy = cdx * bdy;
    // every term is at most 108 bits
    // the sum is at most 109 bits
    int_incircle adnrm2 = adx * adx + ady * ady;

    int_incircle cdxady = cdx * ady;
    int_incircle adxcdy = adx * cdy;
    int_incircle bdnrm2 = bdx * bdx + bdy * bdy;

    int_incircle adxbdy = adx * bdy;
    int_incircle bdxady = bdx * ady;
    int_incircle cdnrm2 = cdx * cdx + cdy * cdy;

    // the factors between brackets are at most 109 bits
    // the terms are at most 218 bits
    // the total sum is at most 220 bits
    int_incircle result = adnrm2 * (bdxcdy - cdxbdy) +
                          bdnrm2 * (cdxady - adxcdy) +
                          cdnrm2 * (adxbdy - bdxady);

    return HelperFunctions::sign(result);
}

/**
 * @brief Adaptive in circle test
 *
 * Returns a positive value if the fourth point is inside the circle through the
 * first three points, a negative value if it is outside, and zero if it is on
 * the circle.
 *
 * This function only works if the first three points are positively oriented,
 * as defined by a positive return value of orient2d_quick, orient2d_exact or
 * orient2d_adaptive. If the points are negatively oriented, the sign of this
 * test flips and a positive value denotes a point outside the circle.
 * If the three points are colinear, then this test will result in undefined
 * behaviour.
 *
 * This function calculates a maximal error bound on the result due to numerical
 * roundoff error. If the result is smaller than this error bound, the function
 * calls incircle_exact to determine the sign of the determinant using exact
 * integer arithmetics. For this to work, the coordinates passed on to this
 * function have to lie in the interval [1,2].
 *
 * @param a Coordinates of the first point, have to be in the interval [1,2]
 * @param b Coordinates of the second point, have to be in the interval [1,2]
 * @param c Coordinates of the third point, have to be in the interval [1,2]
 * @param d Coordinates of the fourth point, have to be in the interval [1,2]
 * @return A positive, negative or zero value, depending on the outcome of the
 * test
 */
double ExactArithmetic::incircle_adaptive(Vec& a, Vec& b, Vec& c, Vec& d) {
    Vec ad = a - d;
    Vec bd = b - d;
    Vec cd = c - d;

    double bdxcdy = bd.x() * cd.y();
    double cdxbdy = cd.x() * bd.y();
    double adnrm2 = ad.norm2();

    double cdxady = cd.x() * ad.y();
    double adxcdy = ad.x() * cd.y();
    double bdnrm2 = bd.norm2();

    double adxbdy = ad.x() * bd.y();
    double bdxady = bd.x() * ad.y();
    double cdnrm2 = cd.norm2();

    double result = adnrm2 * (bdxcdy - cdxbdy) + bdnrm2 * (cdxady - adxcdy) +
                    cdnrm2 * (adxbdy - bdxady);

    double errbound = (fabs(bdxcdy) + fabs(cdxbdy)) * adnrm2 +
                      (fabs(cdxady) + fabs(adxcdy)) * bdnrm2 +
                      (fabs(adxbdy) + fabs(bdxady)) * cdnrm2;

    // actual value is smaller (1.11022e-15 on my local machine), but this will
    // do
    errbound *= 1.e-10;
    if(result < -errbound || result > errbound) {
        return result;
    }

    return incircle_exact(a, b, c, d);
}
#endif

/**
 * @brief Print a binary representation of the given value
 *
 * @param str Label for the printed value
 * @param val Double precision value to print
 */
void ExactArithmetic::print_value(const char* str, double val) {
    binaryDouble cval;
    cval.dval = val;
    cout << str << ": " << std::hex << cval.bval << endl;
    cout << str << ": " << cval.dval << endl;
}

/**
 * @brief Get the double precision floating point corresponding to the given
 * 64-bit binary representation
 *
 * @param bval 64-bit binary representation of a floating point value
 * @return Double precision floating point value that is being represented
 */
double ExactArithmetic::get_value(unsigned long bval) {
    binaryDouble cval;
    cval.bval = bval;
    return cval.dval;
}
