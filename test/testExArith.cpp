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
 * @file testExArith.cpp
 *
 * @brief Unit test for the exact arithmetics algorithms
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "../src/ExArith.hpp"                    // for get_value, etc
#include "../src/utilities/HelperFunctions.hpp"  // for rand_double, sign
#include "../src/utilities/Timer.hpp"            // for Timer
#include "ShadowfaxVersion.hpp"
#include "Vec.hpp"       // for Vec
#include "myAssert.hpp"  // for my_assert
#include <cstdlib>       // for srand
#include <iostream>      // for basic_ostream, ostream, etc

#ifdef RELEASE_BUILD
#define TESTEXARITH_NUMSPEED 100000
#else
#define TESTEXARITH_NUMSPEED 1
#endif

/**
 * @brief Exact arithmetics unit test
 *
 * Tests the quick, exact and adaptive orientation and in circle/sphere
 * functions and times the performance of these functions. Outputs a detailed
 * report to stdout.
 *
 * @param argc Number of command line arguments (ignored)
 * @param argv Command line arguments (ignored)
 * @return 0 on succes. Aborts otherwise
 */
int main(int argc, char** argv) {
#if ndim_ == 3
    cout << "###### Orient3d test ######" << endl;
    // Test 1: should return 1
    {
        cout << "## Standard positive orientation ##" << endl;
        Vec a(1., 1., 1.);
        Vec b(1., 1., 2.);
        Vec c(1., 2., 1.);
        Vec d(2., 1., 1.);
        double result;

        cout << "orient3d_quick: ";
        result = ExactArithmetic::orient3d_quick(a, b, c, d);
        cout << result << endl;
        my_assert(result > 0., "Quick orient3d failed!");

        cout << "orient3d_exact: ";
        result = ExactArithmetic::orient3d_exact(a, b, c, d);
        cout << result << endl;
        my_assert(result > 0., "Exact orient3d failed!");

        cout << "orient3d_adaptive: ";
        result = ExactArithmetic::orient3d_adaptive(a, b, c, d);
        cout << result << endl;
        my_assert(result > 0., "Adaptive orient3d failed!");
    }
    // Test 2: should return -1
    {
        cout << "## Standard negative orientation ##" << endl;
        Vec a(1., 1., 1.);
        Vec b(1., 1., 2.);
        Vec c(2., 1., 1.);
        Vec d(1., 2., 1.);
        double result;

        cout << "orient3d_quick: ";
        result = ExactArithmetic::orient3d_quick(a, b, c, d);
        cout << result << endl;
        my_assert(result < 0., "Quick orient3d failed!");

        cout << "orient3d_exact: ";
        result = ExactArithmetic::orient3d_exact(a, b, c, d);
        cout << result << endl;
        my_assert(result < 0., "Exact orient3d failed!");

        cout << "orient3d_adaptive: ";
        result = ExactArithmetic::orient3d_adaptive(a, b, c, d);
        cout << result << endl;
        my_assert(result < 0., "Adaptive orient3d failed!");
    }
    // Test 3: borderline case
    {
        cout << "## Borderline case ##" << endl;
        // since the floating point values are stored in a binary
        // representation, we cannot simply initialize the vectors for a given
        // case from decimal floating points
        // the hex values below correspond to a case that is known to fail and
        // produce the wrong sign
        // I checked that this case indeed fails with the Intel compilers on the
        // UGent HPC when using Schewchuk's predicates
        Vec a(ExactArithmetic::get_value(0x3ff4427c75c7d120),
              ExactArithmetic::get_value(0x3ff31624fd618368),
              ExactArithmetic::get_value(0x3ff7f481c6bbdc8e));
        Vec b(ExactArithmetic::get_value(0x3ff436fe3c83adae),
              ExactArithmetic::get_value(0x3ff30aa6c41d5ff6),
              ExactArithmetic::get_value(0x3ff7f481c6bbdc8e));
        Vec c(ExactArithmetic::get_value(0x3ff44dfaaf0bf493),
              ExactArithmetic::get_value(0x3ff321a336a5a6da),
              ExactArithmetic::get_value(0x3ff7f481c6bbdc8e));
        Vec d(ExactArithmetic::get_value(0x3ff0000a7c4d0614),
              ExactArithmetic::get_value(0x3ff0000a7c4d0614),
              ExactArithmetic::get_value(0x3ffffff583b2f9ed));
        double result;

        cout << "orient3d_quick: ";
        result = ExactArithmetic::orient3d_quick(a, b, c, d);
        cout << result << endl;
        // no assertion on this result, since it could be wrong

        cout << "orient3d_exact: ";
        result = ExactArithmetic::orient3d_exact(a, b, c, d);
        cout << result << endl;
        my_assert(result < 0., "Exact orient3d failed!");

        cout << "orient3d_adaptive: ";
        result = ExactArithmetic::orient3d_adaptive(a, b, c, d);
        cout << result << endl;
        my_assert(result < 0., "Adaptive orient3d failed!");
    }
    // Test 4: speed test
    {
        cout << "## Speed test ##" << endl;
        // quick
        Timer timerA;
        double sum_quick = 0.;
        srand(42);
        for(unsigned int i = 0; i < TESTEXARITH_NUMSPEED; i++) {
            Vec a(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec b(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec c(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec d(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            sum_quick += HelperFunctions::sign(
                    ExactArithmetic::orient3d_quick(a, b, c, d));
        }
        cout << "orient3d_quick took " << timerA.stop() << "s" << endl;
        cout << "Checksum: " << sum_quick << endl;
        // exact
        Timer timerB;
        double sum_exact = 0.;
        srand(42);
        for(unsigned int i = 0; i < TESTEXARITH_NUMSPEED; i++) {
            Vec a(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec b(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec c(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec d(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            sum_exact += HelperFunctions::sign(
                    ExactArithmetic::orient3d_exact(a, b, c, d));
        }
        cout << "orient3d_exact took " << timerB.stop() << "s" << endl;
        cout << "Checksum: " << sum_exact << endl;
        // adaptive
        Timer timerC;
        double sum_adaptive = 0.;
        srand(42);
        for(unsigned int i = 0; i < TESTEXARITH_NUMSPEED; i++) {
            Vec a(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec b(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec c(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec d(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            sum_adaptive += HelperFunctions::sign(
                    ExactArithmetic::orient3d_adaptive(a, b, c, d));
        }
        cout << "orient3d_adaptive took " << timerC.stop() << "s" << endl;
        cout << "Checksum: " << sum_adaptive << endl;

        my_assert(sum_exact == sum_adaptive,
                  "Adaptive results do not match with exact results!");
    }

    cout << "###### InSphere test ######" << endl;
    // Test 1: should return 0.75
    {
        cout << "## Standard point in sphere ##" << endl;
        Vec a(1., 1., 1.);
        Vec b(1., 1., 2.);
        Vec c(1., 2., 1.);
        Vec d(2., 1., 1.);
        Vec e(1.5, 1.5, 1.5);
        double result;

        cout << "insphere_quick: ";
        result = ExactArithmetic::insphere_quick(a, b, c, d, e);
        cout << result << endl;
        my_assert(result > 0., "Quick insphere failed!");

        cout << "insphere_exact: ";
        result = ExactArithmetic::insphere_exact(a, b, c, d, e);
        cout << result << endl;
        my_assert(result > 0., "Exact insphere failed!");

        cout << "insphere_adaptive: ";
        result = ExactArithmetic::insphere_adaptive(a, b, c, d, e);
        cout << result << endl;
        my_assert(result > 0., "Adaptive insphere failed!");
    }
    // Test 2: should return -0.0625
    {
        cout << "## Standard point outside sphere ##" << endl;
        Vec a(1., 1., 1.);
        Vec b(1., 1., 1.5);
        Vec c(1., 1.5, 1.);
        Vec d(1.5, 1., 1.);
        Vec e(2., 1., 1.);
        double result;

        cout << "insphere_quick: ";
        result = ExactArithmetic::insphere_quick(a, b, c, d, e);
        cout << result << endl;
        my_assert(result < 0., "Quick insphere failed!");

        cout << "insphere_exact: ";
        result = ExactArithmetic::insphere_exact(a, b, c, d, e);
        cout << result << endl;
        my_assert(result < 0., "Exact insphere failed!");

        cout << "insphere_adaptive: ";
        result = ExactArithmetic::insphere_adaptive(a, b, c, d, e);
        cout << result << endl;
        my_assert(result < 0., "Adaptive insphere failed!");
    }
    // Test 3: borderline case
    {
        cout << "## Borderline case ##" << endl;
        Vec a(ExactArithmetic::get_value(0x3ff5175d9c6aad1d),
              ExactArithmetic::get_value(0x3ff378315baa5e27),
              ExactArithmetic::get_value(0x3ff884188b736f90));
        Vec b(ExactArithmetic::get_value(0x3ff51891c1b81df8),
              ExactArithmetic::get_value(0x3ff30d426295616f),
              ExactArithmetic::get_value(0x3ff88648dfb41982));
        Vec c(ExactArithmetic::get_value(0x3ff0000a7c4d0614),
              ExactArithmetic::get_value(0x3ff0000a7c4d0614),
              ExactArithmetic::get_value(0x3ffffff583b2f9ed));
        Vec d(ExactArithmetic::get_value(0x3ff516d995fe14ae),
              ExactArithmetic::get_value(0x3ff378315baa5e27),
              ExactArithmetic::get_value(0x3ff884188b736f90));
        Vec e(ExactArithmetic::get_value(0x3ff515a570b0a3d4),
              ExactArithmetic::get_value(0x3ff30d426295616f),
              ExactArithmetic::get_value(0x3ff88648dfb41982));
        double result;

        cout << "insphere_quick: ";
        result = ExactArithmetic::insphere_quick(a, b, c, d, e);
        cout << result << endl;
        // no assertion, since the quick result will be wrong

        cout << "insphere_exact: ";
        result = ExactArithmetic::insphere_exact(a, b, c, d, e);
        cout << result << endl;
        my_assert(result < 0., "Exact insphere failed!");

        cout << "insphere_adaptive: ";
        result = ExactArithmetic::insphere_adaptive(a, b, c, d, e);
        cout << result << endl;
        my_assert(result < 0., "Adaptive insphere failed!");
    }
    // Test 4: speed test
    {
        cout << "## Speed test ##" << endl;
        // quick
        Timer timerA;
        double sum_quick = 0.;
        srand(42);
        for(unsigned int i = 0; i < TESTEXARITH_NUMSPEED; i++) {
            Vec a(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec b(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec c(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec d(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec e(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            sum_quick += HelperFunctions::sign(
                    ExactArithmetic::insphere_quick(a, b, c, d, e));
        }
        cout << "insphere_quick took " << timerA.stop() << "s" << endl;
        cout << "Checksum: " << sum_quick << endl;
        // exact
        Timer timerB;
        double sum_exact = 0.;
        srand(42);
        for(unsigned int i = 0; i < TESTEXARITH_NUMSPEED; i++) {
            Vec a(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec b(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec c(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec d(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec e(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            sum_exact += HelperFunctions::sign(
                    ExactArithmetic::insphere_exact(a, b, c, d, e));
        }
        cout << "insphere_exact took " << timerB.stop() << "s" << endl;
        cout << "Checksum: " << sum_exact << endl;
        // adaptive
        Timer timerC;
        double sum_adaptive = 0.;
        srand(42);
        for(unsigned int i = 0; i < TESTEXARITH_NUMSPEED; i++) {
            Vec a(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec b(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec c(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec d(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec e(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            sum_adaptive += HelperFunctions::sign(
                    ExactArithmetic::insphere_adaptive(a, b, c, d, e));
        }
        cout << "insphere_adaptive took " << timerC.stop() << "s" << endl;
        cout << "Checksum: " << sum_adaptive << endl;

        my_assert(sum_exact == sum_adaptive,
                  "Adaptive results do not match with exact results!");
    }
#else
    cout << "###### Orient2d test ######" << endl;
    // Test 1: should return 1
    {
        cout << "## Standard positive orientation ##" << endl;
        Vec a(1., 1.);
        Vec b(2., 1.);
        Vec c(1., 2.);
        double result;

        cout << "orient2d_quick: ";
        result = ExactArithmetic::orient2d_quick(a, b, c);
        cout << result << endl;
        my_assert(result > 0., "Quick orient2d failed!");

        cout << "orient2d_exact: ";
        result = ExactArithmetic::orient2d_exact(a, b, c);
        cout << result << endl;
        my_assert(result > 0., "Exact orient2d failed!");

        cout << "orient2d_adaptive: ";
        result = ExactArithmetic::orient2d_adaptive(a, b, c);
        cout << result << endl;
        my_assert(result > 0., "Adaptive orient2d failed!");
    }
    // Test 2: should return -1
    {
        cout << "## Standard negative orientation ##" << endl;
        Vec a(1., 1.);
        Vec b(1., 2.);
        Vec c(2., 1.);
        double result;

        cout << "orient2d_quick: ";
        result = ExactArithmetic::orient2d_quick(a, b, c);
        cout << result << endl;
        my_assert(result < 0., "Quick orient2d failed!");

        cout << "orient2d_exact: ";
        result = ExactArithmetic::orient2d_exact(a, b, c);
        cout << result << endl;
        my_assert(result < 0., "Exact orient2d failed!");

        cout << "orient2d_adaptive: ";
        result = ExactArithmetic::orient2d_adaptive(a, b, c);
        cout << result << endl;
        my_assert(result < 0., "Adaptive orient2d failed!");
    }
    // Test 3: borderline case
    {
        cout << "## Borderline case ##" << endl;
        Vec a(ExactArithmetic::get_value(0x3ff7c3626efe9827),
              ExactArithmetic::get_value(0x3ff3d7dcd7a686dc));
        Vec b(ExactArithmetic::get_value(0x3ff7c098fb743b1b),
              ExactArithmetic::get_value(0x3ff3d6781de15856));
        Vec c(ExactArithmetic::get_value(0x3ff85184738f21ba),
              ExactArithmetic::get_value(0x3ff41eedd9eecba6));
        double result;

        cout << "orient2d_quick: ";
        result = ExactArithmetic::orient2d_quick(a, b, c);
        cout << result << endl;
        // no assertion, because we expect this result to be wrong

        cout << "orient2d_exact: ";
        result = ExactArithmetic::orient2d_exact(a, b, c);
        cout << result << endl;
        my_assert(result < 0., "Exact orient2d failed!");

        cout << "orient2d_adaptive: ";
        result = ExactArithmetic::orient2d_adaptive(a, b, c);
        cout << result << endl;
        my_assert(result < 0., "Adaptive orient2d failed!");
    }
    // Test 4: speed test
    {
        cout << "## Speed test ##" << endl;
        // quick
        Timer timerA;
        double sum_quick = 0.;
        srand(42);
        for(unsigned int i = 0; i < TESTEXARITH_NUMSPEED; i++) {
            Vec a(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec b(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec c(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            sum_quick += HelperFunctions::sign(
                    ExactArithmetic::orient2d_quick(a, b, c));
        }
        cout << "orient2d_quick took " << timerA.stop() << "s" << endl;
        cout << "Checksum: " << sum_quick << endl;
        // exact
        Timer timerB;
        double sum_exact = 0.;
        srand(42);
        for(unsigned int i = 0; i < TESTEXARITH_NUMSPEED; i++) {
            Vec a(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec b(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec c(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            sum_exact += HelperFunctions::sign(
                    ExactArithmetic::orient2d_exact(a, b, c));
        }
        cout << "orient2d_exact took " << timerB.stop() << "s" << endl;
        cout << "Checksum: " << sum_exact << endl;
        // adaptive
        Timer timerC;
        double sum_adaptive = 0.;
        srand(42);
        for(unsigned int i = 0; i < TESTEXARITH_NUMSPEED; i++) {
            Vec a(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec b(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec c(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            sum_adaptive += HelperFunctions::sign(
                    ExactArithmetic::orient2d_adaptive(a, b, c));
        }
        cout << "orient2d_adaptive took " << timerC.stop() << "s" << endl;
        cout << "Checksum: " << sum_adaptive << endl;

        my_assert(sum_exact == sum_adaptive,
                  "Adaptive results do not match with exact results!");
    }

    cout << "###### InCircle test ######" << endl;
    // Test 1: should return 0.5
    {
        cout << "## Standard point in circle ##" << endl;
        Vec a(1., 1.);
        Vec b(2., 1.);
        Vec c(1., 2.);
        Vec d(1.5, 1.5);
        double result;

        cout << "incircle_quick: ";
        result = ExactArithmetic::incircle_quick(a, b, c, d);
        cout << result << endl;
        my_assert(result > 0., "Quick incircle failed!");

        cout << "incircle_exact: ";
        result = ExactArithmetic::incircle_exact(a, b, c, d);
        cout << result << endl;
        my_assert(result > 0., "Exact incircle failed!");

        cout << "incircle_adaptive: ";
        result = ExactArithmetic::incircle_adaptive(a, b, c, d);
        cout << result << endl;
        my_assert(result > 0., "Adaptive incircle failed!");
    }
    // Test 2: should return -0.125
    {
        cout << "## Standard point outside circle ##" << endl;
        Vec a(1., 1.);
        Vec b(1.5, 1.);
        Vec c(1., 1.5);
        Vec d(2., 1.);
        double result;

        cout << "incircle_quick: ";
        result = ExactArithmetic::incircle_quick(a, b, c, d);
        cout << result << endl;
        my_assert(result < 0., "Quick incircle failed!");

        cout << "incircle_exact: ";
        result = ExactArithmetic::incircle_exact(a, b, c, d);
        cout << result << endl;
        my_assert(result < 0., "Exact incircle failed!");

        cout << "incircle_adaptive: ";
        result = ExactArithmetic::incircle_adaptive(a, b, c, d);
        cout << result << endl;
        my_assert(result < 0., "Adaptive incircle failed!");
    }
    // Test 3: borderline case
    {
        cout << "## Borderline case ##" << endl;
        Vec a(ExactArithmetic::get_value(0x3ff894674887daee),
              ExactArithmetic::get_value(0x3ff3fd7c6f726f0b));
        Vec b(ExactArithmetic::get_value(0x3ff8666b561cdb9a),
              ExactArithmetic::get_value(0x3ff42b7861dd6e60));
        Vec c(ExactArithmetic::get_value(0x3ff890391b384f5c),
              ExactArithmetic::get_value(0x3ff3fd7c6f726f0b));
        Vec d(ExactArithmetic::get_value(0x3ff8666b561cdb9a),
              ExactArithmetic::get_value(0x3ff4274a348de2cc));
        double result;

        cout << "incircle_quick: ";
        result = ExactArithmetic::incircle_quick(a, b, c, d);
        cout << result << endl;
        // no assertion, because we expect this result to be wrong

        cout << "incircle_exact: ";
        result = ExactArithmetic::incircle_exact(a, b, c, d);
        cout << result << endl;
        my_assert(result < 0., "Exact incircle failed!");

        cout << "incircle_adaptive: ";
        result = ExactArithmetic::incircle_adaptive(a, b, c, d);
        cout << result << endl;
        my_assert(result < 0., "Adaptive incircle failed!");
    }
    // Test 4: speed test
    {
        cout << "## Speed test ##" << endl;
        // quick
        Timer timerA;
        double sum_quick = 0.;
        srand(42);
        for(unsigned int i = 0; i < TESTEXARITH_NUMSPEED; i++) {
            Vec a(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec b(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec c(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec d(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            sum_quick += HelperFunctions::sign(
                    ExactArithmetic::incircle_quick(a, b, c, d));
        }
        cout << "incircle_quick took " << timerA.stop() << "s" << endl;
        cout << "Checksum: " << sum_quick << endl;
        // exact
        Timer timerB;
        double sum_exact = 0.;
        srand(42);
        for(unsigned int i = 0; i < TESTEXARITH_NUMSPEED; i++) {
            Vec a(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec b(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec c(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec d(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            sum_exact += HelperFunctions::sign(
                    ExactArithmetic::incircle_exact(a, b, c, d));
        }
        cout << "incircle_exact took " << timerB.stop() << "s" << endl;
        cout << "Checksum: " << sum_exact << endl;
        // adaptive
        Timer timerC;
        double sum_adaptive = 0.;
        srand(42);
        for(unsigned int i = 0; i < TESTEXARITH_NUMSPEED; i++) {
            Vec a(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec b(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec c(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            Vec d(1. + HelperFunctions::rand_double(),
                  1. + HelperFunctions::rand_double());
            sum_adaptive += HelperFunctions::sign(
                    ExactArithmetic::incircle_adaptive(a, b, c, d));
        }
        cout << "incircle_adaptive took " << timerC.stop() << "s" << endl;
        cout << "Checksum: " << sum_adaptive << endl;

        my_assert(sum_exact == sum_adaptive,
                  "Adaptive results do not match with exact results!");
    }
#endif
}
