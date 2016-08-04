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
 * @file Cosmology.cpp
 *
 * @brief Cosmological constants and parameters: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "Cosmology.hpp"
#include "GSL.hpp"  // for qag
#include "ParameterFile.hpp"
#include "RestartFile.hpp"  // for RestartFile
#include <cmath>            // for log, sqrt, exp

/**
 * @brief Integrand for the correction factor used in the force equation
 *
 * @param a Value of the scale factor
 * @return Correction factor used in the force equation
 */
double Cosmology::integrand_acceleration_factor(double a) {
    double a2 = a * a;
    double a3 = a2 * a;
    double H = _H0 * sqrt(_omega0 / a3 + _omegaComplement / a2 + _omegaLambda);
    return 1. / (a2 * H);
}

/**
 * @brief Integrand for the correction factor used in the velocity equation
 *
 * @param a Value of the scale factor
 * @return Correction factor used in the velocity equation
 */
double Cosmology::integrand_velocity_factor(double a) {
    double a2 = a * a;
    double a3 = a2 * a;
    double H = _H0 * sqrt(_omega0 / a3 + _omegaComplement / a2 + _omegaLambda);
    return 1. / (a3 * H);
}

/**
 * @brief Precompute look-up tables for the cosmological correction factors
 */
void Cosmology::tabulate_factors() {
    double logBegin = log(_a0);
    _da = log(_a0 + _da) - logBegin;
    double da = _da / COSMOLOGY_TABLE_SIZE;

    double result;
    for(unsigned int i = 0; i < COSMOLOGY_TABLE_SIZE; i++) {
        double a1 = exp(logBegin + da * (i + 1));
        result = GSL::qag(static_integrand_acceleration_factor, _a0, a1, 1.e-8,
                          this);
        _acceleration_table[i] = result;
    }

    for(unsigned int i = 0; i < COSMOLOGY_TABLE_SIZE; i++) {
        double a1 = exp(logBegin + da * (i + 1));
        result = GSL::qag(static_integrand_velocity_factor, _a0, a1, 1.e-8,
                          this);
        _velocity_table[i] = result;
    }

    _a0 = logBegin;
}

/**
 * @brief Constructor
 *
 * @param H0 Current value of the Hubble parameter
 * @param omega0 Current value of \f$\Omega{}_0\f$
 * @param omegaLambda Current value of \f$\Omega{}_\Lambda{}\f$
 * @param a0 Lower bound of the scale factor range
 * @param a1 Upper bound of the scale factor range
 */
Cosmology::Cosmology(double H0, double omega0, double omegaLambda, double a0,
                     double a1) {
    _H0 = H0;
    _omega0 = omega0;
    _omegaLambda = omegaLambda;
    _omegaComplement = 1. - _omega0 - _omegaLambda;
    _a0 = a0;
    _da = a1 - a0;
    tabulate_factors();
}

/**
 * @brief ParameterFile constructor
 *
 * @param parameters ParameterFile containing the desired parameter values
 * @param hubble_factor Factor to convert the Hubble constant to internal units
 */
Cosmology::Cosmology(ParameterFile* parameters, double hubble_factor)
        : Cosmology(hubble_factor *
                            parameters->get_parameter<double>(
                                    "Cosmology.Hubble0",
                                    COSMOLOGY_DEFAULT_HUBBLE0),
                    parameters->get_parameter<double>("Cosmology.OmegaMatter",
                                                      COSMOLOGY_DEFAULT_OMEGA0),
                    parameters->get_parameter<double>(
                            "Cosmology.OmegaLambda",
                            COSMOLOGY_DEFAULT_OMEGALAMBDA),
                    1. / (1. +
                          parameters->get_parameter<double>(
                                  "Cosmology.InitialRedshift",
                                  COSMOLOGY_DEFAULT_INITIALREDSHIFT))) {}

/**
 * @brief Hubble parameter
 *
 * @param a Value of the scale factor
 * @return Hubble parameter at that moment in the expansion history of the
 * Universe
 */
double Cosmology::H(double a) {
    double a2 = a * a;
    double a3 = a2 * a;
    return _H0 * sqrt(_omega0 / a3 + _omegaComplement / a2 + _omegaLambda);
}

/**
 * @brief Get the initial scale factor
 *
 * @return Initial scale factor
 */
double Cosmology::get_a0() {
    // for efficiency reasons, we store log(a0) internally, need to convert
    return exp(_a0);
}

/**
 * @brief Convert the given physical position to a comoving position
 *
 * @param physical_position Physical position
 * @param a Scale factor
 * @return Comoving position
 */
Vec Cosmology::comoving_position(const Vec& physical_position, double a) {
    return physical_position / a;
}

/**
 * @brief Convert the given comoving position to a physical position
 *
 * @param comoving_position Comoving position
 * @param a Scale factor
 * @return Physical position
 */
Vec Cosmology::physical_position(const Vec& comoving_position, double a) {
    return a * comoving_position;
}

/**
 * @brief Convert the given velocity to comoving velocity
 *
 * @param physical_position Physical position
 * @param physical_velocity Physical velocity
 * @param a Scale factor
 * @return Comoving velocity
 */
Vec Cosmology::comoving_velocity(const Vec& physical_position,
                                 const Vec& physical_velocity, double a) {
    return a * (physical_velocity - H(a) * physical_position);
}

/**
 * @brief Convert the given comoving velocity to physical velocity
 *
 * @param physical_position Physical position
 * @param comoving_velocity Comoving velocity
 * @param a Scale factor
 * @return Physical velocity
 */
Vec Cosmology::physical_velocity(const Vec& physical_position,
                                 const Vec& comoving_velocity, double a) {
    return comoving_velocity / a + H(a) * physical_position;
}

/**
 * @brief Convert the given comoving velocity to the Gadget peculiar velocity
 *
 * The Gadget peculiar velocity is the velocity as generated for Gadget by
 * commonly used IC generators (e.g. MUSIC), and is defined as
 * \f[
 *   \vec{v} = \sqrt{a} \frac{{\rm d}}{{\rm d}t} \vec{x},
 * \f]
 * with \f$\vec{x}\f$ the comoving position. This definition differs from the
 * usual definition of the peculiar velocity by a factor \f$\sqrt{a}\f$.
 *
 * @param comoving_velocity Comoving velocity
 * @param a Scale factor
 * @return Gadget peculiar velocity
 */
Vec Cosmology::peculiar_velocity(const Vec& comoving_velocity, double a) {
    return comoving_velocity / a / sqrt(a);
}

/**
 * @brief Convert the given Gadget peculiar velocity to the comoving velocity
 *
 * The Gadget peculiar velocity is the velocity as generated for Gadget by
 * commonly used IC generators (e.g. MUSIC), and is defined as
 * \f[
 *   \vec{v} = \sqrt{a} \frac{{\rm d}}{{\rm d}t} \vec{x},
 * \f]
 * with \f$\vec{x}\f$ the comoving position. This definition differs from the
 * usual definition of the peculiar velocity by a factor \f$\sqrt{a}\f$.
 *
 * @param peculiar_velocity Gadget peculiar velocity
 * @param a Scale factor
 * @return Comoving velocity
 */
Vec Cosmology::comoving_velocity(const Vec& peculiar_velocity, double a) {
    return a * sqrt(a) * peculiar_velocity;
}

/**
 * @brief Get the correction factor for comoving coordinates in the
 * gravitational force equation
 *
 * @param a0 Scale factor at the beginning of the kick step
 * @param a1 Scale factor at the end of the kick step
 * @return Correction factor that is multiplied with the acceleration before the
 * kick
 */
double Cosmology::get_acceleration_factor(double a0, double a1) {
    // this code comes straight from Gadget2
    double df1, df2, u1, u2;
    int i1, i2;

    u1 = (log(a0) - _a0) / _da * COSMOLOGY_TABLE_SIZE;
    i1 = (int)u1;
    if(i1 >= COSMOLOGY_TABLE_SIZE) {
        i1 = COSMOLOGY_TABLE_SIZE - 1;
    }

    if(i1 <= 1) {
        df1 = u1 * _acceleration_table[0];
    } else {
        df1 = _acceleration_table[i1 - 1] +
              (_acceleration_table[i1] - _acceleration_table[i1 - 1]) *
                      (u1 - i1);
    }

    u2 = (log(a1) - _a0) / _da * COSMOLOGY_TABLE_SIZE;
    i2 = (int)u2;
    if(i2 >= COSMOLOGY_TABLE_SIZE) {
        i2 = COSMOLOGY_TABLE_SIZE - 1;
    }

    if(i2 <= 1) {
        df2 = u2 * _acceleration_table[0];
    } else {
        df2 = _acceleration_table[i2 - 1] +
              (_acceleration_table[i2] - _acceleration_table[i2 - 1]) *
                      (u2 - i2);
    }

    return df2 - df1;
}

/**
 * @brief Get the correction factor for comoving coordinates in the velocity
 * equation
 *
 * @param a0 Scale factor at the beginning of the drift step
 * @param a1 Scale factor at the end of the drift step
 * @return Correction factor that is multiplied with the velocity before the
 * drift
 */
double Cosmology::get_velocity_factor(double a0, double a1) {
    // this code comes straight from Gadget2
    double df1, df2, u1, u2;
    int i1, i2;

    u1 = (log(a0) - _a0) / _da * COSMOLOGY_TABLE_SIZE;
    i1 = (int)u1;
    if(i1 >= COSMOLOGY_TABLE_SIZE) {
        i1 = COSMOLOGY_TABLE_SIZE - 1;
    }

    if(i1 <= 1) {
        df1 = u1 * _velocity_table[0];
    } else {
        df1 = _velocity_table[i1 - 1] +
              (_velocity_table[i1] - _velocity_table[i1 - 1]) * (u1 - i1);
    }

    u2 = (log(a1) - _a0) / _da * COSMOLOGY_TABLE_SIZE;
    i2 = (int)u2;
    if(i2 >= COSMOLOGY_TABLE_SIZE) {
        i2 = COSMOLOGY_TABLE_SIZE - 1;
    }

    if(i2 <= 1) {
        df2 = u2 * _velocity_table[0];
    } else {
        df2 = _velocity_table[i2 - 1] +
              (_velocity_table[i2] - _velocity_table[i2 - 1]) * (u2 - i2);
    }

    return df2 - df1;
}

/**
 * @brief Get the scale factor interval corresponding to the given time interval
 *
 * @param dt Physical time interval
 * @param a Current value of the scale factor
 * @return Scale factor interval
 */
double Cosmology::get_da(double dt, double a) {
    return a * H(a) * dt;
}

/**
 * @brief Dump the Cosmology to the given RestartFile
 *
 * @param rfile RestartFile to write to
 */
void Cosmology::dump(RestartFile& rfile) {
    rfile.write(_H0);
    rfile.write(_omega0);
    rfile.write(_omegaLambda);
    rfile.write(_omegaComplement);
    rfile.write(_acceleration_table, COSMOLOGY_TABLE_SIZE);
    rfile.write(_velocity_table, COSMOLOGY_TABLE_SIZE);
    rfile.write(_a0);
    rfile.write(_da);
}

/**
 * @brief Restart constructor
 *
 * @param rfile RestartFile to read from
 */
Cosmology::Cosmology(RestartFile& rfile) {
    rfile.read(_H0);
    rfile.read(_omega0);
    rfile.read(_omegaLambda);
    rfile.read(_omegaComplement);
    rfile.read(_acceleration_table, COSMOLOGY_TABLE_SIZE);
    rfile.read(_velocity_table, COSMOLOGY_TABLE_SIZE);
    rfile.read(_a0);
    rfile.read(_da);
}
