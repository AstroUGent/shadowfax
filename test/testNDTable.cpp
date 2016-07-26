/*******************************************************************************
 * This file is part of Shadowfax
 * Copyright (C) 2016 Yorick Van Den Bossche (yorick.vandenbossche@ugent.be)
 *                    Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file testNDTable.cpp
 *
 * @brief Unit test for the 5D regular table
 *
 * @author Yorick Van Den Bossche (yorick.vandenbossche@ugent.be)
 */
#include "NDTable.hpp"
#include "io/UnitSet.hpp"
#include "myAssert.hpp"
#include <dirent.h>
#include <math.h>
#include <vector>

/**
 * @brief Test the 5D table
 *
 * @param table FiveDTable to test
 */
void test_NDTable(FiveDTable& table) {
    for(double i = -5.; i < 10.; i += 4.7) {
        for(double j = -5.; j < 10.; j += 4.7) {
            for(double k = 0.; k < 3.; k += 1.4) {
                for(double l = 0.; l < 3.; l += 1.4) {
                    for(double m = 0.; m < 3.; m += 1.4) {
                        double b = table.get_value(
                                vector<double>({i, j, k, l, m}));
                        double a = pow(10., i + j + k + l + m);
                        cout << b << "\t" << a << endl;
                        assert_values_equal(b, a, "Interpolation error");
                    }
                }
            }
        }
    }
}

/**
 * @brief Main test routine
 *
 * @param argc Number of command line arguments
 * @param argv Command line arguments
 * @return Exit code
 */
int main(int argc, char** argv) {
    DIR* dir;
    struct dirent* ent;
    vector<string> result;
    cout << endl << "Looking for cooling tables..." << endl;
    int n = 0;
    if((dir = opendir("coolingtables/")) != NULL) {
        /* print all the files and directories within directory */
        while((ent = readdir(dir)) != NULL) {
            string str(ent->d_name);
            if(str.find(".rates") != string::npos) {
                result.push_back(string("coolingtables/") + str);
                n++;
            }
        }
        closedir(dir);
        cout << n << " tables found." << endl;
    } else {
        cout << "None found." << endl;
    }
    cout << endl;

    Unit _unit_mass("mass", "g", 0.001);
    Unit _unit_length("length", "cm", 0.01);
    Unit _unit_time("time", "s", 1);
    UnitSet units(_unit_length, _unit_mass, _unit_time);
    FiveDTable table(result, &units);
    cout << "Starting test" << endl;
    test_NDTable(table);
    cout << "Successfully finished test" << endl;

    return 0;
}
