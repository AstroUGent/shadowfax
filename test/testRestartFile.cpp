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
 * @file testRestartFile.cpp
 *
 * @brief Unit test for the RestartFile class
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "CoolingLocation.hpp"
#include "CoolingTable.hpp"
#include "GSL.hpp"
#include "RestartFile.hpp"
#include "io/UnitSetGenerator.hpp"
#include "myAssert.hpp"
#include "utilities/HelperFunctions.hpp"
using namespace std;

#define TESTRESTARTFILE_MODE_WRITE true
#define TESTRESTARTFILE_MODE_READ false

/**
 * @brief RestartFile test
 *
 * Tests the functionality of the RestartFile class
 *
 * @param argc Number of command line arguments
 * @param argv Command line arguments
 * @return Exit code
 */
int main(int argc, char** argv) {

    vector<pair<double, double>> vecpair1;
    double x[100], y[100];
    for(unsigned int i = 0; i < 100; ++i) {
        vecpair1.push_back(make_pair(HelperFunctions::rand_double(),
                                     HelperFunctions::rand_double()));
        x[i] = i * 0.01;
        y[i] = 3. * x[i] * x[i] + 2. * x[i] + 2.;
    }

    UnitSet* units = UnitSetGenerator::generate("SI");
    CoolingTable cooling(COOLING_LOCATION, units);
    GSL::GSLInterpolator* interpolator_spline = GSL::GSLInterpolator::create(
            GSL::TypeGSLCubicSplineInterpolator, x, y, 100);
    GSL::GSLInterpolator* interpolator_linear = GSL::GSLInterpolator::create(
            GSL::TypeGSLLinearInterpolator, x, y, 100);

    // we put this in a separate block to make sure rfile is closed before
    // attempting to read it
    {
        RestartFile rfile(TESTRESTARTFILE_MODE_WRITE);

        cout << "Writing vecpair1 to restartfile" << endl;
        rfile.write(vecpair1);

        cout << "Writing cooling to restartfile" << endl;
        cooling.dump(rfile);

        cout << "Writing cubic spline interpolator to restartfile" << endl;
        interpolator_spline->dump(rfile);
        cout << "Writing linear interpolator to restartfile" << endl;
        interpolator_linear->dump(rfile);
    }

    RestartFile rfile2(TESTRESTARTFILE_MODE_READ);

    cout << "Reading vecpair2 from restartfile" << endl;
    vector<pair<double, double>> vecpair2;
    rfile2.read(vecpair2);

    cout << "Checking vecpair..." << endl;
    for(unsigned int i = 0; i < 100; ++i) {
        my_assert(vecpair1[i].first == vecpair2[i].first,
                  "Elements not equal!");
        my_assert(vecpair1[i].second == vecpair2[i].second,
                  "Elements not equal!");
    }
    cout << "Done." << endl;

    cout << "Reading cooling2 from restartfile" << endl;
    CoolingTable cooling2(rfile2);

    // compare cooling and cooling2
    // we randomly draw some parameter values and compare the results
    cout << "Checking cooling..." << endl;
    for(unsigned int i = 0; i < 100; i++) {
        double Fe = HelperFunctions::rand_value(-4., 0.5);
        double Mg = HelperFunctions::rand_value(-0.55, 0.47);
        double z = HelperFunctions::rand_value(0., 12.);
        double n = HelperFunctions::rand_value(2.e-33, 2.e-22);
        double T = HelperFunctions::rand_value(10., 1.e9);
        double c1 = cooling.get_value(Fe, Mg, z, n, T);
        double c2 = cooling2.get_value(Fe, Mg, z, n, T);
        my_assert(c1 == c2, "Cooling tables are not the same!");
    }
    cout << "Done." << endl;

    cout << "Reading cubic spline interpolator from restartfile" << endl;
    GSL::GSLInterpolator* interpolator_spline2 =
            GSL::GSLInterpolator::create(rfile2);
    cout << "Checking cubic spline interpolator..." << endl;
    if(!interpolator_spline->equal(interpolator_spline2)) {
        cerr << "Cubic spline interpolator is not the same after restart!"
             << endl;
        abort();
    }
    cout << "Done." << endl;

    cout << "Reading linear interpolator from restartfile" << endl;
    GSL::GSLInterpolator* interpolator_linear2 =
            GSL::GSLInterpolator::create(rfile2);
    cout << "Checking linear interpolator..." << endl;
    if(!interpolator_linear->equal(interpolator_linear2)) {
        cerr << "Linear interpolator is not the same after restart!" << endl;
        abort();
    }
    cout << "Done." << endl;

    delete units;
    delete interpolator_linear;
    delete interpolator_linear2;
    delete interpolator_spline;
    delete interpolator_spline2;

    return 0;
}
