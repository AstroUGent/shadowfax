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
 * @file GasCooling.cpp
 *
 * @brief Gas cooling: implementation
 *
 * @author Yorick Van Den Bossche (yorick.vandenbossche@ugent.be)
 */
#include "GasCooling.hpp"
#include "CoolingTable.hpp"           // for CoolingTable
#include "ParameterFile.hpp"          // for ParameterFile
#include "Physics.hpp"                // for Physics
#include "StateVector.hpp"            // for StateVector
#include "io/Unit.hpp"                // for Unit, operator/, operator*
#include "io/UnitConverter.hpp"       // for UnitConverter
#include "io/UnitSet.hpp"             // for UnitSet, etc
#include "utilities/GasParticle.hpp"  // for GasParticle
#include <algorithm>                  // for min, max
#include <cmath>                      // for abs, log10
#include <iostream>                   // for operator<<, basic_ostream, etc
#include <sstream>                    // for basic_stringbuf<>::int_type, etc
#include <string>                     // for string
using namespace std;

/**
 * @brief Cooling Table constructor
 *
 * Reads the files and makes the table.
 *
 * @param directory Folder where the cooling tables are stored on the system
 * @param simulation_units The UnitSet used in the simulation, for correct
 * energy units
 * @param parameters ParameterFile containing the parameters of the simulation
 * @param physics Physical constants
 */
GasCooling::GasCooling(std::string directory, UnitSet* simulation_units,
                       ParameterFile* parameters, Physics* physics) {
    _cooling_table = new CoolingTable(directory, simulation_units);
    UnitSet SI(UNITS_DEFAULT);
    // boltzman is in energy per Kelvin but UnitSet doesn't seem to have
    // temperature
    Unit unit_kperm = SI.get_length_unit() * SI.get_length_unit() /
                      SI.get_time_unit() / SI.get_time_unit();
    Unit sim_kperm = simulation_units->get_length_unit() *
                     simulation_units->get_length_unit() /
                     simulation_units->get_time_unit() /
                     simulation_units->get_time_unit();
    UnitConverter kperm_conv(unit_kperm, sim_kperm);
    // k/amu = 1.3806488e-23 / 1.67372e-27
    double mean_mol_weight;
    if(parameters) {
        // need to define this variable in Physics, since it is used in both
        // GasCooling and StarFormationParticleConverter
        mean_mol_weight = physics->get_mean_mol_weight();
        _k_per_mH = kperm_conv.convert(8248.98310351 / mean_mol_weight);
        _req_prec = parameters->get_parameter<double>(
                "Cooling.ReqPrec", GASCOOLING_DEFAULT_REQPREC);
        _max_it = parameters->get_parameter<int>("Cooling.MaxIt",
                                                 GASCOOLING_DEFAULT_MAXIT);
        // used for fitting only last n points, currently unused
        _max_fit = _max_it;
        _timesplit = parameters->get_parameter<double>(
                "Cooling.TimeSplitFactor", GASCOOLING_DEFAULT_TIMESPLIT);
    } else {
        // resort to default parameters
        _k_per_mH = kperm_conv.convert(8248.98310351 / 1.219512195);
        _req_prec = GASCOOLING_DEFAULT_REQPREC;
        _max_it = GASCOOLING_DEFAULT_MAXIT;
        _max_fit = GASCOOLING_DEFAULT_MAXIT;
        _timesplit = GASCOOLING_DEFAULT_TIMESPLIT;
        mean_mol_weight = 1.219512195;
    }
    _Hfrac = (4. - mean_mol_weight) / (3. * mean_mol_weight);
    _results = vector<double>(_max_it, 0.);
    _calca = vector<double>(_max_it + 1, 0.);
    _calcb = vector<double>(_max_it + 1, 0.);
}

/**
 * @brief Destructor.
 */
GasCooling::~GasCooling() {
    delete _cooling_table;
}

/**
 * @brief Get cooling rate
 *
 * Gets the cooling rate from the table
 *
 * @warning The redshift is currently hardcoded to a value of 11!
 *
 * @param Fe Iron abundance
 * @param Mg Magnesium abundance
 * @param density gas density in SI units
 * @param temp temperature in Kelvin
 *
 * @return Cooling rate for the given values of Fe, Mg, n, z and T
 */
double GasCooling::cooling(double Fe, double Mg, double density, double temp) {
    return -_cooling_table->get_value(Fe, Mg, 11., density, temp);
}

/**
 * @brief Calculate cooling for a particle for the current timestep
 *
 * Calculates the cooling amount for a gas particle
 *
 * @param particle the gas particle to be cooled
 *
 * @return The energy loss due to cooling for the given particle and timestep
 */
double GasCooling::calc_cooling(GasParticle* particle) {
    return calc_cooling_bs(particle);
}

/**
 * @brief No idea
 *
 * technically the interpolation uses an iterations x iterations array,
 * but with the right approach we only need to use the last two rows
 * so we can use 2 vectors and use this function to pretend it's a 2D array
 *
 * @param i No idea
 * @return No idea
 */
inline vector<double>& GasCooling::calcvec(int i) {
    switch(i % 2) {
        case 0:
            return _calca;
        case 1:
            return _calcb;
        default:
            return _calca;
    }
}

/**
 * @brief Return (reduced) timestep (?)
 *
 * @param j No idea
 * @return No idea
 */
inline double x_(int j) {
    return 1. / j;
}

/**
 * @brief Calculate cooling for a particle for the current timestep
 *
 * Calculates the cooling amount for a gas particle using the Bulirsch-Stoer
 * method.
 * Should be accurate for fairly large timesteps.
 *
 * @param particle the gas particle to be cooled
 *
 * @return The energy loss due to cooling for the given particle and timestep
 */
double GasCooling::calc_cooling_bs(GasParticle* particle) {
    // use WVec for P and rho, so as to calculate T
    StateVector state = particle->get_Wvec();
    double pressure = state.p();
    double density = state.rho();
    if(!density) {
        // vacuum, no cooling
        return 0.;
    }
    // use Qvec for mass (for total energy change at the end)
    state = particle->get_Qvec();
    double m = state.m();
    double ecomp = state.e();
    // for checking overcooling
    // bad idea to calculate the specific energy like this; better use the
    // corresponding RiemannSolver method...
    // also, the energy also depends on the density...
    double energy = 1.5 * pressure;
    double mH = m * _Hfrac;
    // 2.937794191 = log10(1.154e-3)
    double Fe = log10(state.Fe() / mH / 56.) + 2.937794191;
    Fe = std::max(-99., Fe);
    // log10(M_mgsol/M_fesol) = -0.281304, source: Asplund et al, 2005
    // (=log10((10^7.53)×24.3050÷((10^7.45)×55.845)))
    double Mg;
    if(Fe <= -99.) {
        Mg = 0.;
    } else {
        Mg = log10(state.Mg() / state.Fe() * 56. / 24.) + 0.281304;
    }

    // calculate starting temperature
    // should be a RiemannSolver method as well
    double T = pressure / density / _k_per_mH;
    // get total integration time
    double dt_glob = particle->get_real_timestep();
    // factor often used in calculation, for converting thermal energy to T
    double per_rho_k = 2. / 3. / _k_per_mH / density;

    /*Debugging output (could use ProgramLog for this...)
    cout<<"--------------------------------------"<<endl;
    cout<<"T="<<T<<endl;
    cout<<"dttot="<<dt_glob<<endl;*/

    // initialize variables
    double t_left = dt_glob;
    double dt, T_prev, TI, T0, T1, T2, DT_initial;
    bool accurate;
    int it = 0;
    TI = T;

    /* General idea: take a part of the timestep based on the first-order
     * temperature change DT_initial, such that the first-order estimate of
     * the temperature change is no larger than T_initial / timesplit_factor.
     * Do a BS integation over this part.
     * Repeat with the resulting energy until the timestep is
     * completely integrated.
     * When there is strong cooling often the first tenth of the timestep needs
     * hundreds of iterations while the remaining 90% can be done in just 2,
     * so this seems to work quite well. */
    while(t_left > 0.) {
        /*calculate part timestep*/
        // Initial DT/dt
        DT_initial = this->cooling(Fe, Mg, density, TI) * per_rho_k;
        dt = -TI / (_timesplit * DT_initial);
        if(dt > t_left) {
            dt = t_left;
        }
        /*Debugging output
        cout<<"doing "<<dt<<" of "<<dt_glob<<"with "<<t_left<<" left which is
        1/"<<pow(2,n)<<endl;
        cout<<"Initial: "<<this->cooling(density, T)<<" so
        "<<DT_initial<<"K/s"<<endl;
        */

        // The actual BS integration
        accurate = false;
        T_prev = 0.;
        double T_new = 0.;
        for(int j = 1; (j <= _max_it && !accurate); ++j) {
            // get timestep for this try
            double h = dt / (2. * j);
            // integrate the temperature by means of 2*j midpoint steps
            T0 = TI;
            T1 = T0 + h * DT_initial;
            for(int i = 0; i < 2 * j - 1; ++i) {
                it++;
                T2 = T0 +
                     2. * h * this->cooling(Fe, Mg, density, T1) * per_rho_k;
                T0 = T1;
                T1 = T2;
            }
            _results[j - 1] =
                    0.5 * (T2 + T0 +
                           h * this->cooling(Fe, Mg, density, T2) * per_rho_k);

            // make a new estimate for the total deltaT by interpolating
            // data for different h to h=0
            T_prev = T_new;
            T_new = 0.;
            // we use rational interpolation because it works well for even
            // lots of datapoints (polynomial breaks past 8)
            calcvec(j)[1] = _results[j - 1];
            for(int k = 2; k <= min(j, _max_fit); ++k) {
                calcvec(j)[k] =
                        calcvec(j)[k - 1] +
                        (calcvec(j)[k - 1] - calcvec(j - 1)[k - 1]) /
                                (x_(j - k + 1) / x_(j) *
                                         (1. -
                                          (calcvec(j)[k - 1] -
                                           calcvec(j - 1)[k - 1]) /
                                                  (calcvec(j)[k - 1] -
                                                   calcvec(j - 1)[k - 2])) -
                                 1.);
            }
            T_new = calcvec(j)[min(j, _max_fit)];

            /*Debugging output
            cout<<"attempt "<<j<<", "<<results[j-1]<<" K"<<endl;
            cout<<"which gives "<<T_new<<" K interpolated, "<<abs(( T_new-T_prev
            ))<<endl;
            */
            if(j > 1 && T_new > 0. &&
               abs((T_new - T_prev)) < abs(_req_prec * T_new)) {
                // accurate enough?
                accurate = true;
            }
        }
        if(!accurate) {
            cout << "Cooling did not converge" << endl;
        }
        // set new starting T, susbstract done time
        TI = T_new;
        t_left -= dt;

        // cout<<"T_i= "<<TI<<endl;
    }
    // cout<<"Cooling finished in "<<it<<" iterations"<<endl;

    double dE = 1.5 * (TI - T) * _k_per_mH * m;
    // prevent energy from going below zero - this should never happen
    if(dE >= ecomp) {
        // say something about overcooling in the ProgramLog
        return 0.99 * energy * m / density;
    } else {
        // dQ uses the opposite sign, so cooling needs to be positive
        return -dE;
    }
}

/**
 * @brief Calculate cooling for a particle for the current timestep
 *
 * Calculates the cooling amount for a gas particle using the modified midpoint
 * method
 *
 * @param particle the gas particle to be cooled
 *
 * @return The energy loss due to cooling for the given particle and timestep
 */
double GasCooling::calc_cooling_mmm(GasParticle* particle) {
    // must use WVec and calculate
    StateVector state = particle->get_Wvec();
    double pressure = state.p();
    double density = state.rho();
    state = particle->get_Qvec();
    double Fe = -99.;
    double Mg = 0.;
    double m = state.m();
    double energy = 2. / 3. * pressure;
    double T = pressure / density / _k_per_mH;
    double dt = particle->get_real_timestep();
    // propagate the energy through 1 timestep, currently uses the second order
    // midpoint method
    double a = this->cooling(Fe, Mg, density, T);
    double k1 = dt * a / density * 0.75 / _k_per_mH;
    if(k1 < -T) {
        // say something about overcooling in the ProgramLog
        return 0.99 * energy * m / density;
    }
    double b = this->cooling(Fe, Mg, density, T + k1);
    double dE = dt * b * m / density;
    // prevent energy from going below zero
    if(dE >= energy * m / density) {
        // say something about overcooling in the ProgramLog
        return 0.99 * energy * m / density;
    } else {
        // dQ uses the opposite sign, so cooling needs to be positive
        return -dE;
    }
}

/**
 * @brief Dump the object to the given RestartFile
 *
 * @param rfile RestartFile to write to
 */
void GasCooling::dump(RestartFile& rfile) {
    _cooling_table->dump(rfile);
    rfile.write(_k_per_mH);
    rfile.write(_req_prec);
    rfile.write(_max_it);
    rfile.write(_max_fit);
    rfile.write(_timesplit);
    rfile.write(_Hfrac);
    rfile.write(_results);
    rfile.write(_calca);
    rfile.write(_calcb);
}

/**
 * @brief Restart constructor
 *
 * @param rfile RestartFile to read from
 */
GasCooling::GasCooling(RestartFile& rfile) {
    _cooling_table = new CoolingTable(rfile);
    rfile.read(_k_per_mH);
    rfile.read(_req_prec);
    rfile.read(_max_it);
    rfile.read(_max_fit);
    rfile.read(_timesplit);
    rfile.read(_Hfrac);
    rfile.read(_results);
    rfile.read(_calca);
    rfile.read(_calcb);
}
