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
 * @file Cosmology.hpp
 *
 * @brief Cosmological constants and parameters: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef COSMOLOGY_HPP
#define COSMOLOGY_HPP

// size of the look-up tables used to calculate the correction factors
#define COSMOLOGY_TABLE_SIZE 1000

#include "Vec.hpp"  // for Vec

class ParameterFile;
class RestartFile;

#define COSMOLOGY_DEFAULT_HUBBLE0 0.7
#define COSMOLOGY_DEFAULT_OMEGA0 0.28
#define COSMOLOGY_DEFAULT_OMEGALAMBDA 0.72
#define COSMOLOGY_DEFAULT_INITIALREDSHIFT 99.

/**
 * @brief Class containing the cosmological constants and parameters used by the
 * code
 */
class Cosmology {
  private:
    /*! @brief Current value of the Hubble parameter */
    double _H0;

    /*! @brief Current value of \f$\Omega{}_0\f$ */
    double _omega0;

    /*! @brief Current value of \f$\Omega{}_\Lambda{}\f$ */
    double _omegaLambda;

    /*! @brief \f$1 - \Omega{}_0 - \Omega{}_\Lambda{}\f$ */
    double _omegaComplement;

    /*! @brief Look-up table for acceleration factors */
    double _acceleration_table[COSMOLOGY_TABLE_SIZE];

    /*! @brief Look-up table for velocity factors */
    double _velocity_table[COSMOLOGY_TABLE_SIZE];

    /*! @brief Scale factor for the first acceleration factor in the look-up
     * table */
    double _a0;

    /*! @brief Total scale factor range tabulated in the look-up table */
    double _da;

    double integrand_acceleration_factor(double a);

    double integrand_velocity_factor(double a);

    /**
     * @brief Static exposure of the integrand_acceleration_factor function,
     * needed for usage with the GSL library
     *
     * @param a Value of the scale factor
     * @param param Pointer to a Cosmology instance
     * @return The return value of the integrand_acceleration_factor for the
     * given Cosmology instance
     */
    static double static_integrand_acceleration_factor(double a, void* param) {
        Cosmology* cosmo = (Cosmology*)param;
        return cosmo->integrand_acceleration_factor(a);
    }

    /**
     * @brief Static exposure of the integrand_velocity_factor function,
     * needed for usage with the GSL library
     *
     * @param a Value of the scale factor
     * @param param Pointer to a Cosmology instance
     * @return The return value of the integrand_velocity_factor for the
     * given Cosmology instance
     */
    static double static_integrand_velocity_factor(double a, void* param) {
        Cosmology* cosmo = (Cosmology*)param;
        return cosmo->integrand_velocity_factor(a);
    }

    void tabulate_factors();

  public:
    Cosmology(double H0, double omega0 = COSMOLOGY_DEFAULT_OMEGA0,
              double omegaLambda = COSMOLOGY_DEFAULT_OMEGALAMBDA,
              double a0 = 1.e-9, double a1 = 1.);
    Cosmology(ParameterFile* parameters, double hubble_factor);

    double H(double a);

    double get_a0();

    Vec comoving_position(const Vec& physical_position, double a);
    Vec physical_position(const Vec& comoving_position, double a);

    Vec comoving_velocity(const Vec& physical_position,
                          const Vec& physical_velocity, double a);
    Vec physical_velocity(const Vec& physical_position,
                          const Vec& comoving_velocity, double a);

    Vec peculiar_velocity(const Vec& comoving_velocity, double a);
    Vec comoving_velocity(const Vec& peculiar_velocity, double a);

    double get_acceleration_factor(double a0, double a1);
    double get_velocity_factor(double a0, double a1);

    double get_da(double dt, double a);

    void dump(RestartFile& rfile);
    Cosmology(RestartFile& rfile);
};

#endif  // COSMOLOGY_HPP
