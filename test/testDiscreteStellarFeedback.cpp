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
 * @file testDiscreteStellarFeedback.cpp
 *
 * @brief Unit test for discrete stellar feedback
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "ClosestNgbSearch.hpp"
#include "DiscreteStellarFeedback.hpp"
#include "DiscreteStellarFeedbackData.hpp"
#include "MPIMethods.hpp"
#include "io/UnitConverter.hpp"
#include "io/UnitSet.hpp"
#include "io/UnitSetGenerator.hpp"
#include "myAssert.hpp"
#include "utilities/ParticleVector.hpp"
#include "utilities/StarParticle.hpp"
#include "utilities/Tree.hpp"
#include <iostream>
using namespace std;

/**
 * @brief Test the feedback given by a single star with given mass and
 * metallicity  over the given time interval
 *
 * @param feedback DiscreteStellarFeedback instance used to give feedback
 * @param FeH Metallicity of the star
 * @param mass Mass of the star
 * @param PopII_SNII_expect Expected number of PopII SNII
 * @param PopII_SNIa_expect Expected number of PopII SNIa
 * @param PopIII_SN_expect Expected number of PopIII SN
 * @param tstart Start time of the feedback
 * @param tend End time of the feedback
 * @param dm_expect Expected change in mass of the surrounding gas
 * @param dE_expect Expected change in energy of the surrounding gas
 */
void test_feedback(DiscreteStellarFeedback& feedback, double FeH, double mass,
                   unsigned int PopII_SNII_expect,
                   unsigned int PopII_SNIa_expect,
                   unsigned int PopIII_SN_expect, double tstart, double tend,
                   double dm_expect, double dE_expect) {

    Vec position(0.5, 0.5, 0.5);
    StarParticle* star = new StarParticle(position);
    star->set_initial_mass(mass);
    star->set_birthtime(0.);

    // normal PopII star
    star->set_FeH(FeH);
    DiscreteStellarFeedbackData* data =
            (DiscreteStellarFeedbackData*)feedback.initialize_data(star);

    // check numbers:
    my_assert(data->get_PopII_SNII_number() == PopII_SNII_expect,
              "Wrong number of PopII SNII explosions!");
    my_assert(data->get_PopII_SNIa_number() == PopII_SNIa_expect,
              "Wrong number of PopII SNIa explosions!");
    my_assert(data->get_PopIII_SN_number() == PopIII_SN_expect,
              "Wrong number of PopIII SN explosions!");

    star->set_feedback_data(data);

    Vec center(0.5, 0.5, 0.5);
    Vec sides(1., 1., 1.);
    RectangularBox treebox(center, sides);
    StellarFeedbackTreeFilter filter(&feedback, tstart, tend);
    ParticleVector particles(false, treebox, true);
    particles.add_star_particle(star);

    // add a single GasParticle to receive all feedback
    Vec gasposition(0.6, 0.6, 0.6);
    particles.add_gas_particle(new GasParticle(gasposition));
    particles.sort();
    particles.get_tree().walk_tree<ClosestNgbSearch>(particles, false, false,
                                                     true, filter);

    my_assert(star->get_closest_gasparticle() == particles.gasback(),
              "Wrong closest GasParticle!");

    feedback.do_feedback(star, tstart, tend);

    StateVector dQ = particles.gas(0)->get_dQvec();
    assert_values_equal(dQ.m(), dm_expect, "Unexpected mass change!");
    assert_values_equal(dQ.e(), dE_expect, "Unexpected energy change!");
}

/**
 * @brief Unit test for discrete stellar feedback
 *
 * @param argc Number of command line arguments (ignored)
 * @param argv Command line arguments (ignored)
 * @return 0 on succes. Aborts otherwise
 */
int main(int argc, char** argv) {
    // galactic units: kpc, Msol, Gyr
    UnitSet* units = UnitSetGenerator::generate("galactic");
    DiscreteStellarFeedback feedback(*units);

    Unit energy_unit("length*length*mass/time/time", "erg", 1.e-7);
    UnitConverter energy_converter(energy_unit, units->get_energy_unit());

    // PopII star
    double mass_expected = 0.;
    double energy_expected = 0.;
    // SNII
    mass_expected += 0.191445322565 * 2500.;
    energy_expected += 28.7646 * energy_converter.convert(1.e51) * 0.7;
    // SNIa
    mass_expected += 0.00655147325196 * 2500.;
    energy_expected += 4.31469 * energy_converter.convert(1.e51) * 0.7;
    // SW
    energy_expected += 28.7646 * energy_converter.convert(1.e50) * 0.7;
    test_feedback(feedback, -2., 2500., 28, 4, 0, 0., 13.8, -mass_expected,
                  -energy_expected);

    // PopIII star
    mass_expected = 0.;
    energy_expected = 0.;
    // SN
    mass_expected += 0.45 * 2500.;
    energy_expected += 5.31551e+09;
    // SW
    energy_expected += 45.4619 * energy_converter.convert(1.e51) * 0.7;
    test_feedback(feedback, -6., 2500., 0, 0, 45, 0., 13.8, -mass_expected,
                  -energy_expected);

    delete units;

    return 0;
}
