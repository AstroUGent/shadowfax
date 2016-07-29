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
 * @file testGSL.cpp
 *
 * @brief Unit test for the GSL module
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "GSL.hpp"
#include <cmath>
#include <iostream>

/**
 * @brief Test function to integrate
 *
 * @param x Evaluation point
 * @param args Extra arguments (not used)
 * @return Value of the function
 */
double f(double x, void* args) {
    return sqrt(x);
}

/**
 * @brief Test function to find the zeropoint
 *
 * @param x Evaluation point
 * @param args Extra arguments (not used)
 * @return Value of the function
 */
double f_brent(double x, void* args) {
    return 3. * x - 2.;
}

/**
 * @brief Test function for interpolation
 *
 * @param x Evaluation point
 * @param args Extra arguments (not used)
 * @return Value of the function
 */
double f_interpol(double x, void* args) {
    return 3. * x * x + 4. * x - 2.;
}

/**
 * @brief Unit test for the GSL module
 *
 * @param argc Number of command line arguments
 * @param argv Command line arguments
 * @return Exit code
 */
int main(int argc, char** argv) {
    double result = GSL::qag(f, 0., 1., 1.e-8, NULL);

    std::cout << "Result: " << result << " (" << (2. / 3.) << ")" << std::endl;

    result = GSL::brent(f_brent, 0., 1., 1.e-8, NULL);

    std::cout << "Result: " << result << " (" << (2. / 3.) << ")" << std::endl;

    double x[5] = {0., 0.25, 0.5, 0.75, 1.};
    double y[5];
    for(unsigned int i = 0; i < 5; i++) {
        y[i] = f_interpol(x[i], NULL);
    }

    GSL::GSLInterpolator* lin = GSL::GSLInterpolator::create(
            GSL::TypeGSLLinearInterpolator, x, y, 5);
    GSL::GSLInterpolator* spline = GSL::GSLInterpolator::create(
            GSL::TypeGSLCubicSplineInterpolator, x, y, 5);
    for(unsigned int i = 0; i < 10; i++) {
        double xtest = i * 0.1;
        std::cout << lin->eval(xtest) << "\t" << spline->eval(xtest) << " ("
                  << f_interpol(xtest, NULL) << ")" << std::endl;
    }
    delete lin;
    delete spline;

    return 0;
}
