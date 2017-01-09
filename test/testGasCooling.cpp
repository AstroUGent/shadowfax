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
 * @file testGasCooling.cpp
 *
 * @brief Unit test for the GasCooling
 *
 * @author Yorick Van Den Bossche (yorick.vandenbossche@ugent.be)
 */
#include "GasCooling.hpp"
#include "StateVector.hpp"
#include "io/UnitSet.hpp"
#include "myAssert.hpp"
#include "utilities/GasParticle.hpp"
#include <boost/algorithm/string.hpp>
#include <dirent.h>
#include <math.h>
#include <vector>

#define TESTGASCOOLING_NSTEPS 100

/**
 * @brief Test gas cooling routines
 *
 * @param gascooling GasCooling object to use
 */
void test_GasCooling(GasCooling& gascooling) {
    GasParticle* gasparticle = new GasParticle();
    // density, vx, vx, vz, pressure
    double rho = 1.e-19;
    double p = 1.e-11;
    double vol = pow(1.e21, 3.);
    double gamma = 5. / 3.;
    StateVector W(rho, 0., 0., 0., p);
    StateVector Q;
    Q.set_m(W.rho() * vol);
    Q.set_e((0.5 * W.rho() *
                     (W.vx() * W.vx() + W.vy() * W.vy() + W.vz() * W.vz()) +
             W.p() / (gamma - 1.)) *
            vol);
    double E0 = Q.e();
    gasparticle->set_W(W);
    gasparticle->set_Q(Q);
    double dt = 0.1;
    gasparticle->set_real_timestep(dt);
    StateVector state;
    double dE_tot = 0.;
    double dE;
    for(unsigned int i = 0; i < TESTGASCOOLING_NSTEPS; ++i) {
        dE = gascooling.calc_cooling(gasparticle);
        dE_tot += dE;

        state.set_e(dE);
        gasparticle->increase_dQ(state);
        gasparticle->update_Q();
        Q = gasparticle->get_Qvec();
        W.set_rho(Q.m() / vol);
        W.set_p((gamma - 1.) * (Q.e()) / vol);
        gasparticle->set_W(W);
    }
    double E = gasparticle->get_Qvec().e();
    assert_values_approx(E, E0 * exp(-0.01 * TESTGASCOOLING_NSTEPS * dt), 0.005,
                         "Error interpolating");

    delete gasparticle;
}

/**
 * @brief Main test program
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

    Unit _unit_mass("mass", "kg", 1.);
    Unit _unit_length("length", "m", 1.);
    Unit _unit_time("time", "s", 1.);
    Unit unit_temperature("temperature", "K", 1.);
    UnitSet units(_unit_length, _unit_mass, _unit_time, unit_temperature);
    GasCooling gascooling("coolingtables/", &units, NULL, NULL);
    cout << "Starting test" << endl;
    test_GasCooling(gascooling);
    cout << "Test successfully finished" << endl;

    return 0;
}
