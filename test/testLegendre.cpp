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
 * @file testLegendre.cpp
 *
 * @brief Unit test for the Legendre polynomials
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "../src/utilities/HelperFunctions.hpp"
#include "../src/utilities/Timer.hpp"
#include "Legendre.hpp"
#include "ShadowfaxVersion.hpp"
#include "myAssert.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>

#ifdef RELEASE_BUILD
#define TESTLEGENDRE_NSPEED 10000
#else
#define TESTLEGENDRE_NSPEED 1
#endif

/**
 * @brief 1D function used to test the Gauss-Legendre quadrature
 *
 * @param x Real value
 * @return x squared
 */
double test_function1d(double x) {
    return x * x;
}

/**
 * @brief 2D function used to test the Gauss-Legendre quadrature
 *
 * @param x Real value
 * @param y Real value
 * @return x squared
 */
double test_function2d(double x, double y) {
    return x * x + y * y;
}

/**
 * @brief 3D function used to test the Gauss-Legendre quadrature
 *
 * @param x Real value
 * @param y Real value
 * @param z Real value
 * @return x squared
 */
double test_function3d(double x, double y, double z) {
    return x * x + y * y + z * z;
}

/**
 * @brief Test the Legendre polynomials
 *
 * @param argc Number of command line arguments (ignored)
 * @param argv Command line arguments (ignored)
 * @return 0 on succes. Aborts otherwise
 */
int main(int argc, char** argv) {

    // test zeros and weights
    // we use an order 3 polynomial, since this is the lowest order polynomial
    // with non-trivial weights
    Legendre legendre3(3);

    std::vector<double> zeros = legendre3.get_zeros();
    assert_values_equal(zeros[0], sqrt(3. / 5.), "Wrong root!");
    assert_values_equal(zeros[1], 0., "Wrong root!");
    assert_values_equal(zeros[2], -sqrt(3. / 5.), "Wrong root!");

    std::vector<double> weights = legendre3.get_weights();
    assert_values_equal(weights[0], 5. / 9., "Wrong weight!");
    assert_values_equal(weights[1], 8. / 9., "Wrong weight!");
    assert_values_equal(weights[2], 5. / 9., "Wrong weight!");

    double quad1d = legendre3.quadrature1d(test_function1d);
    assert_values_equal(quad1d, 2. / 3., "Quadrature error!");

    double quad2d = legendre3.quadrature2d(test_function2d);
    assert_values_equal(quad2d, 8. / 3., "Quadrature error!");

    double quad3d = legendre3.quadrature3d(test_function3d);
    assert_values_equal(quad3d, 8., "Quadrature error!");

    // test fast representation
    double prec = legendre3.polynomial1d_recursive(0.4);
    double pdir = legendre3.polynomial1d_direct(0.4);
    assert_values_equal(prec, pdir, "Polynomial representations not equal!");

    // time representations
    std::cout << "Speed tests:" << std::endl;
    for(unsigned char order = 1; order < 10; order++) {
        std::cout << "Order " << (unsigned int)order << std::endl;
        Legendre legendre(order);
        double rsum = 0.;
        double dsum = 0.;
        // recursive definition
        {
            Timer timer;
            srand(42);
            for(unsigned int i = 0; i < TESTLEGENDRE_NSPEED; i++) {
                rsum += legendre.polynomial1d_recursive(
                        2. * HelperFunctions::rand_double() - 1.);
            }
            std::cout << "Recursive: " << timer.stop() << "s" << std::endl;
        }
        // direct definition
        {
            Timer timer;
            srand(42);
            for(unsigned int i = 0; i < TESTLEGENDRE_NSPEED; i++) {
                dsum += legendre.polynomial1d_direct(
                        2. * HelperFunctions::rand_double() - 1.);
            }
            std::cout << "Direct: " << timer.stop() << "s" << std::endl;
        }
        assert_values_equal(rsum, dsum, "Representations not equal!");
    }

    return 0;
}
