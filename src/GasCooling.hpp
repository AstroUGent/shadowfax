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
 * @file GasCooling.hpp
 *
 * @brief Gas cooling: header
 *
 * @author Yorick Van Den Bossche (yorick.vandenbossche@ugent.be)
 */
#ifndef GASCOOLING_HPP
#define GASCOOLING_HPP
#include "CoolingTable.hpp"
#include "LogFiles.hpp"
#include "ParameterFile.hpp"
#include "io/Unit.hpp"
#include "io/UnitSet.hpp"
#include "utilities/GasParticle.hpp"
#include <vector>
using namespace std;

class Physics;
class RestartFile;

#define GASCOOLING_DEFAULT_REQPREC 0.001
#define GASCOOLING_DEFAULT_MAXIT 128
#define GASCOOLING_DEFAULT_TIMESPLIT 10.

/**
  * @brief A class that calculates the radiative cooling experienced by a
  * gas particle.
  *
  * The current implementation uses a linear cooling, first order accurate,
  * for testing purposes.
  */
class GasCooling {
  private:
    /*! @brief Table of cooling rates */
    CoolingTable* _cooling_table;
    /*! @brief Boltzman constant/mass of hydrogen in simulation units*/
    double _k_per_mH;
    /*! @brief Required precision of the integration */
    double _req_prec;
    /*! @brief Maximum number of iterations */
    int _max_it;
    /*! @brief No idea */
    int _max_fit;
    /*! @brief No idea */
    double _timesplit;
    /*! @brief Hydrogen fraction */
    double _Hfrac;
    /*! @brief Results? */
    vector<double> _results;
    /*! @brief Calc A ? */
    vector<double> _calca;
    /*! @brief Calc B ? */
    vector<double> _calcb;

    double cooling(double Fe, double Mg, double density, double temp);
    double calc_cooling_mmm(GasParticle* particle);
    double calc_cooling_bs(GasParticle* particle);
    vector<double>& calcvec(int i);

  public:
    GasCooling(string directory, UnitSet* simulation_units,
               ParameterFile* parameters, Physics* physics);
    ~GasCooling();

    double calc_cooling(GasParticle* particle);

    void dump(RestartFile& rfile);
    GasCooling(RestartFile& rfile);
};

#endif
