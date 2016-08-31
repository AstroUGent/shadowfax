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
#include "myAssert.hpp"
#include "utilities/ParticleVector.hpp"
#include "utilities/StarParticle.hpp"
#include "utilities/Tree.hpp"
#include <iostream>
using namespace std;

/**
 * @brief Unit test for discrete stellar feedback
 *
 * @param argc Number of command line arguments (ignored)
 * @param argv Command line arguments (ignored)
 * @return 0 on succes. Aborts otherwise
 */
int main(int argc, char** argv) {
    MyMPI_Init(&argc, &argv);

    DiscreteStellarFeedback feedback;

    Vec position(0.5, 0.5, 0.5);
    StarParticle* star = new StarParticle(position);
    star->set_initial_mass(2500.);  // in solar masses

    // normal PopII star
    star->set_FeH(-2.);
    DiscreteStellarFeedbackData* data =
            (DiscreteStellarFeedbackData*)feedback.initialize_data(star);

    // check numbers:
    my_assert(data->get_PopII_SNII_number() == 28,
              "Wrong number of PopII SNII explosions!");
    my_assert(data->get_PopII_SNIa_number() == 4,
              "Wrong number of PopII SNIa explosions!");
    my_assert(data->get_PopIII_SN_number() == 0,
              "Wrong number of PopIII SN explosions!");

    // PopIII star
    star->set_FeH(-6.);
    data = (DiscreteStellarFeedbackData*)feedback.initialize_data(star);
    star->set_feedback_data(data);

    // check numbers:
    my_assert(data->get_PopII_SNII_number() == 0,
              "Wrong number of PopII SNII explosions!");
    my_assert(data->get_PopII_SNIa_number() == 0,
              "Wrong number of PopII SNIa explosions!");
    my_assert(data->get_PopIII_SN_number() == 45,
              "Wrong number of PopIII SN explosions!");

    Vec center(0.5, 0.5, 0.5);
    Vec sides(1., 1., 1.);
    RectangularBox treebox(center, sides);
    StellarFeedbackTreeFilter filter(&feedback, 0., 13.8);
    ParticleVector particles(false, treebox, true);
    particles.add_star_particle(star);

    // add a single GasParticle to receive all feedback
    Vec gasposition(0.6, 0.6, 0.6);
    particles.add_gas_particle(new GasParticle(gasposition));
    particles.finalize();
    particles.sort();
    particles.get_tree().walk_tree<ClosestNgbSearch>(particles, false, false,
                                                     true, filter);

    my_assert(star->get_closest_gasparticle() == particles.gasback(),
              "Wrong closest GasParticle!");

    feedback.do_feedback(star, 0., 13.8);

    return MyMPI_Finalize();
}
