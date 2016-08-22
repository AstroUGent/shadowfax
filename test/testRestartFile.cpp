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
    for(unsigned int i = 0; i < 100; ++i) {
        vecpair1.push_back(make_pair(HelperFunctions::rand_double(),
                                     HelperFunctions::rand_double()));
    }

    UnitSet* units = UnitSetGenerator::generate("SI");
    CoolingTable cooling(COOLING_LOCATION, units);

    // we put this in a separate block to make sure rfile is closed before
    // attempting to read it
    {
        RestartFile rfile(TESTRESTARTFILE_MODE_WRITE);

        cout << "Writing vecpair1 to restartfile" << endl;
        rfile.write(vecpair1);

        cooling.dump(rfile);
    }

    RestartFile rfile2(TESTRESTARTFILE_MODE_READ);

    vector<pair<double, double>> vecpair2;
    rfile2.read(vecpair2);

    for(unsigned int i = 0; i < 100; ++i) {
        my_assert(vecpair1[i].first == vecpair2[i].first,
                  "Elements not equal!");
        my_assert(vecpair1[i].second == vecpair2[i].second,
                  "Elements not equal!");
    }

    CoolingTable cooling2(rfile2);

    // compare cooling and cooling2
    // we randomly draw some parameter values and compare the results
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

    delete units;

    return 0;
}
