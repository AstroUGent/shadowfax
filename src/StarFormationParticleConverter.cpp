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
 * @file StarFormationParticleConverter.cpp
 *
 * @brief ParticleConverter implementation that converts a GasParticle to a
 * StarParticle if star formation criteria are met: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "StarFormationParticleConverter.hpp"
#include "ParameterFile.hpp"            // for ParameterFile
#include "Physics.hpp"                  // for Physics
#include "StateVector.hpp"              // for StateVector
#include "Vec.hpp"                      // for Vec
#include "io/Unit.hpp"                  // for Unit, operator/, operator*
#include "io/UnitConverter.hpp"         // for UnitConverter
#include "io/UnitSet.hpp"               // for UnitSet, etc
#include "utilities/GasParticle.hpp"    // for GasParticle
#include "utilities/Particle.hpp"       // for Particle
#include "utilities/ParticleTypes.hpp"  // for ParticleType, etc
#include "utilities/StarParticle.hpp"   // for StarParticle
#include <sstream>                      // for basic_stringbuf<>::int_type, etc

/**
 * @brief Check if the given ParticleType can be converted
 *
 * Only GasParticles can be converted by this ParticleConverter
 *
 * @param type ParticleType
 * @return True if the ParticleType is PARTTYPE_GAS
 */
bool StarFormationParticleConverter::do_conversion(ParticleType type) {
    return type == PARTTYPE_GAS;
}

/**
 * @brief Convert the given GasParticle to a StarParticle if the star formation
 * criteria are met
 *
 * @param particle GasParticle
 * @return A new StarParticle or the original GasParticle if no star is formed
 */
Particle* StarFormationParticleConverter::convert(Particle* particle) {
    StateVector state = ((GasParticle*)particle)->get_Wvec();
    double pressure = state.p();
    double density = state.rho();
    // no need to calculate divergence if density is too low
    if(density < _densitylimit) {
        return particle;
    }
    double T = pressure / density / _k_per_mH;
    // likewise, no need to go on if T is too high
    if(T >= _templimit) {
        return particle;
    }
    // calculate divergence
    StateVector state2[ndim_];
    ((GasParticle*)particle)->get_gradients(state2);
    double divergence = state2[0].vx() + state2[1].vy();
#if ndim_ == 3
    divergence += state2[2].vz();
#endif
    if(divergence >= 0) {
        return particle;
    } else {
        // star formation happens
        Vec pos = particle->get_position();
        Vec vel = particle->get_velocity();
        Vec a = particle->get_gravitational_acceleration();
        double old_a = particle->get_old_acceleration();
        unsigned long id = particle->id();
        unsigned long dtstart = particle->get_starttime();
        unsigned long dtend = particle->get_endtime();
        unsigned int compcost = particle->get_comp_cost();
        double epot = particle->get_gravitational_potential();
        double hsoft = particle->get_hsoft();
        double mass = particle->get_mass();
        StarParticle* part = new StarParticle(pos);
        part->set_velocity(vel);
        part->set_gravitational_acceleration(a);
        part->set_old_acceleration(old_a);
        part->set_id(id);
        part->set_starttime(dtstart);
        part->set_endtime(dtend);
        part->add_comp_cost(compcost);
        part->set_gravitational_potential(epot);
        part->set_hsoft(hsoft);
        part->set_mass(mass);
        // flag this StarParticle as new; it's birthtime will be set by the
        // main simulation loop
        part->set_birthtime(-1.);
        delete particle;
        return part;
    }
}

/**
 * @brief Constructor
 *
 * @param parameters ParameterFile holding the parameters of the simulation
 * @param simulation_units Internal units
 * @param physics Physical constants
 */
StarFormationParticleConverter::StarFormationParticleConverter(
        ParameterFile* parameters, UnitSet* simulation_units,
        Physics* physics) {
    UnitSet SI;
    _densitylimit = parameters->get_parameter<double>(
            "StarFormation.MinSFDensity",
            STARFORMATIONPARTICLECONVERTER_DEFAULT_MINSFDENSITY);
    _templimit = parameters->get_parameter<double>(
            "StarFormation.MaxSFTemp",
            STARFORMATIONPARTICLECONVERTER_DEFAULT_MAXSFTEMP);

    // boltzman is in energy per Kelvin but UnitSet doesn't seem to have
    // temperature
    Unit unit_kperm = SI.get_length_unit() * SI.get_length_unit() /
                      SI.get_time_unit() / SI.get_time_unit();
    Unit sim_kperm = simulation_units->get_length_unit() *
                     simulation_units->get_length_unit() /
                     simulation_units->get_time_unit() /
                     simulation_units->get_time_unit();
    UnitConverter _kperm_conv(unit_kperm, sim_kperm);
    // k/amu = 1.3806488e-23 / 1.67372e-27
    _k_per_mH =
            _kperm_conv.convert(8248.98310351 / physics->get_mean_mol_weight());
}

/**
 * @brief Dump the object to the given RestartFile
 *
 * @param rfile RestartFile to write to
 */
void StarFormationParticleConverter::dump(RestartFile& rfile) {
    rfile.write(_densitylimit);
    rfile.write(_templimit);
    rfile.write(_k_per_mH);
}

/**
 * @brief Restart constructor
 *
 * @param rfile RestartFile to read from
 */
StarFormationParticleConverter::StarFormationParticleConverter(
        RestartFile& rfile) {
    rfile.read(_densitylimit);
    rfile.read(_templimit);
    rfile.read(_k_per_mH);
}
