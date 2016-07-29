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
 * @file StarFormationParticleConverter.hpp
 *
 * @brief ParticleConverter implementation that converts a GasParticle to a
 * StarParticle if star formation criteria are met: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef STARFORMATIONPARTICLECONVERTER_HPP
#define STARFORMATIONPARTICLECONVERTER_HPP

#include "ParameterFile.hpp"
#include "StateVector.hpp"
#include "io/Unit.hpp"
#include "io/UnitConverter.hpp"
#include "io/UnitSet.hpp"
#include "io/UnitSetGenerator.hpp"
#include "utilities/ParticleConverter.hpp"

class Physics;
class RestartFile;

#define STARFORMATIONPARTICLECONVERTER_DEFAULT_MINSFDENSITY 2.e-22
#define STARFORMATIONPARTICLECONVERTER_DEFAULT_MAXSFTEMP 15000.

/**
 * @brief ParticleConverter implementation that converts a GasParticle to a
 * StarParticle if star formation criteria are met
 */
class StarFormationParticleConverter : public ParticleConverter {
  private:
    /*! @brief Density criterion for star formation */
    double _densitylimit;
    /*! @brief Temperature criterion */
    double _templimit;
    /*! @brief Factor used to convert between pressure and temperature */
    double _k_per_mH;

  public:
    StarFormationParticleConverter(ParameterFile* parameters,
                                   UnitSet* simulation_units, Physics* physics);
    virtual ~StarFormationParticleConverter() {}

    virtual bool do_conversion(ParticleType type);
    virtual Particle* convert(Particle* particle);

    void dump(RestartFile& rfile);
    StarFormationParticleConverter(RestartFile& rfile);
};

#endif  // STARFORMATIONPARTICLECONVERTER_HPP
