/*******************************************************************************
 * This file is part of Shadowfax
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *                    Yorick Van Den Bossche (yorick.vandenbossche@ugent.be)
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
 * @file ContinuousStellarFeedback.cpp
 *
 * @brief StellarFeedback implementation that spreads out the stellar feedback
 * uniformly over the lifetime of the various feedback mechanisms:
 * implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "ContinuousStellarFeedback.hpp"
#include "NgbSearch.hpp"                 // for ClosestNgbSearch
#include "ParameterFile.hpp"             // for ParameterFile
#include "StateVector.hpp"               // for StateVector
#include "utilities/GasParticle.hpp"     // for GasParticle
#include "utilities/ParticleVector.hpp"  // for ParticleVector
#include "utilities/StarParticle.hpp"    // for StarParticle
#include "utilities/Tree.hpp"            // for Tree
#include <algorithm>                     // for max, min
#include <cstddef>                       // for NULL
#include <sstream>                       // for basic_stringbuf<>::int_type, etc
using namespace std;

/**
 * @brief Do the feedback for the given StarParticle during the given timestep
 *
 * @param star StarParticle that does feedback
 * @param particles Reference to the ParticleVector holding all particles
 * @param dt Timestep over which the feedback is given
 */
void ContinuousStellarFeedback::do_feedback(StarParticle* star,
                                            ParticleVector& particles,
                                            double dt) {
    double age = star->get_age();
    // stellar wind period
    double dtsw = std::max(
            0., (std::min(age + dt, _sw_end) - std::max(age, _sw_start)));
    // Supernova Ia period
    double dtsia = std::max(
            0., (std::min(age + dt, _sia_end) - std::max(age, _sia_start)));
    // Supernova II period
    double dtsii = std::max(
            0., (std::min(age + dt, _sii_end) - std::max(age, _sii_start)));
    // fetch gasparticle if anything needs doing, otherwise return
    // should turn this condition around...
    if(dtsw > 0. || dtsia > 0. || dtsii > 0.) {
    } else {
        return;
    }
    ClosestNgbSearch ngbsearch(star->get_position(), _sradius);
    particles.get_tree().walk_tree(ngbsearch);
    GasParticle* part = ngbsearch.get_closest();
    while(part == NULL) {
        ngbsearch.increase_radius();
        particles.get_tree().walk_tree(ngbsearch);
        part = ngbsearch.get_closest();
    }
    StateVector dQ;
    double mass = star->get_mass();
    // stellar wind
    dQ.set_e(dQ.e() -
             dtsw * mass * _massfac_snii * _sw_E / (_sw_end - _sw_start));

    // supernova Ia
    dQ.set_e(dQ.e() -
             dtsia * mass * _massfac_snii * _nsnia_per_snii * _sia_E /
                     (_sia_end - _sia_start));
    dQ.set_m(dQ.m() -
             dtsia * mass * _massfac_snii * _nsnia_per_snii * _sia_m /
                     (_sia_end - _sia_start));
    dQ.set_Fe(dQ.Fe() -
              dtsia * mass * _massfac_snii * _nsnia_per_snii * _sia_Fe /
                      (_sia_end - _sia_start));
    dQ.set_Mg(dQ.Mg() -
              dtsia * mass * _massfac_snii * _nsnia_per_snii * _sia_Mg /
                      (_sia_end - _sia_start));

    // supernova II
    dQ.set_e(dQ.e() -
             dtsii * mass * _massfac_snii * _sii_E / (_sii_end - _sii_start));
    dQ.set_m(dQ.m() -
             dtsii * mass * _massfac_snii * _sii_m / (_sii_end - _sii_start));
    dQ.set_Fe(dQ.Fe() -
              dtsii * mass * _massfac_snii * _sii_Fe / (_sii_end - _sii_start));
    dQ.set_Mg(dQ.Mg() -
              dtsii * mass * _massfac_snii * _sii_Mg / (_sii_end - _sii_start));

    star->set_mass(star->get_mass() + dQ.m());
    part->increase_dQ(dQ);
}

/**
 * @brief Constructor
 *
 * @warning Unit conversion is not yet implemented!!
 *
 * @param parameters The parameterfile used
 */
ContinuousStellarFeedback::ContinuousStellarFeedback(
        ParameterFile* parameters) {
    _sw_start = parameters->get_parameter<double>(
            "StellarFeedback.SWStart", STELLARFEEDBACK_DEFAULT_SWSTART);
    _sw_end = parameters->get_parameter<double>("StellarFeedback.SWEnd",
                                                STELLARFEEDBACK_DEFAULT_SWEND);
    _sii_start = parameters->get_parameter<double>(
            "StellarFeedback.SNIIStart", STELLARFEEDBACK_DEFAULT_SNIISTART);
    _sii_end = parameters->get_parameter<double>(
            "StellarFeedback.SNIIEnd", STELLARFEEDBACK_DEFAULT_SNIIEND);
    _sia_start = parameters->get_parameter<double>(
            "StellarFeedback.SNIaStart", STELLARFEEDBACK_DEFAULT_SNIASTART);
    _sia_end = parameters->get_parameter<double>(
            "StellarFeedback.SNIaEnd", STELLARFEEDBACK_DEFAULT_SNIAEND);

    _sii_Fe = parameters->get_parameter<double>("StellarFeedback.SNII_Fe",
                                                STELLARFEEDBACK_DEFAULT_SNIIFE);
    _sii_Mg = parameters->get_parameter<double>("StellarFeedback.SNII_Mg",
                                                STELLARFEEDBACK_DEFAULT_SNIIMG);
    _sia_Fe = parameters->get_parameter<double>("StellarFeedback.SNIa_Fe",
                                                STELLARFEEDBACK_DEFAULT_SNIAFE);
    _sia_Mg = parameters->get_parameter<double>("StellarFeedback.SNIa_Mg",
                                                STELLARFEEDBACK_DEFAULT_SNIAMG);
    _sii_E = parameters->get_parameter<double>("StellarFeedback.SNII_E",
                                               STELLARFEEDBACK_DEFAULT_SNIIE);
    _sw_E = parameters->get_parameter<double>("StellarFeedback.SW_E",
                                              STELLARFEEDBACK_DEFAULT_SWE);
    _sia_E = parameters->get_parameter<double>("StellarFeedback.SNIa_E",
                                               STELLARFEEDBACK_DEFAULT_SNIAE);
    _sii_m = parameters->get_parameter<double>("StellarFeedback.SNII_m",
                                               STELLARFEEDBACK_DEFAULT_SNIIM);
    _sia_m = parameters->get_parameter<double>("StellarFeedback.SNIa_m",
                                               STELLARFEEDBACK_DEFAULT_SNIAM);
    _nsnia_per_snii = parameters->get_parameter<double>(
            "StellarFeedback.nSNIa_per_SNII",
            STELLARFEEDBACK_DEFAULT_NSNIAPERSNII);
    _massfac_snii = parameters->get_parameter<double>(
            "StellarFeedback.MassFac_SNII",
            STELLARFEEDBACK_DEFAULT_MASSFACSNII);

    _sradius = 1.;
}

/**
 * @brief Set the radius for the search algorithm
 *
 * @param rad Value for the radius
 */
void ContinuousStellarFeedback::set_radius(double rad) {
    _sradius = rad;
}

/**
 * @brief Dump the object to the given RestartFile
 *
 * @param rfile RestartFile to write to
 */
void ContinuousStellarFeedback::dump(RestartFile& rfile) {
    rfile.write(_sw_start);
    rfile.write(_sw_end);
    rfile.write(_sii_start);
    rfile.write(_sii_end);
    rfile.write(_sia_start);
    rfile.write(_sia_end);

    rfile.write(_sii_Fe);
    rfile.write(_sii_Mg);
    rfile.write(_sia_Fe);
    rfile.write(_sia_Mg);
    rfile.write(_sii_E);
    rfile.write(_sw_E);
    rfile.write(_sia_E);
    rfile.write(_sia_m);
    rfile.write(_sii_m);
    rfile.write(_nsnia_per_snii);
    rfile.write(_massfac_snii);
    rfile.write(_sradius);
}

/**
 * @brief Restart constructor
 *
 * @param rfile RestartFile to read from
 */
ContinuousStellarFeedback::ContinuousStellarFeedback(RestartFile& rfile) {
    rfile.read(_sw_start);
    rfile.read(_sw_end);
    rfile.read(_sii_start);
    rfile.read(_sii_end);
    rfile.read(_sia_start);
    rfile.read(_sia_end);

    rfile.read(_sii_Fe);
    rfile.read(_sii_Mg);
    rfile.read(_sia_Fe);
    rfile.read(_sia_Mg);
    rfile.read(_sii_E);
    rfile.read(_sw_E);
    rfile.read(_sia_E);
    rfile.read(_sia_m);
    rfile.read(_sii_m);
    rfile.read(_nsnia_per_snii);
    rfile.read(_massfac_snii);
    rfile.read(_sradius);
}
