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
 * @file Legendre.hpp
 *
 * @brief Legendre polynomials: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef LEGENDRE_HPP
#define LEGENDRE_HPP

#include <vector>

/**
 * @brief Legendre polynomials and Gauss-Legendre quadrature for a given order
 */
class Legendre {
  private:
    /*! @brief Order of the polynomials */
    unsigned char _order;

    /*! @brief Zeros of the polynomial */
    std::vector<double> _zeros;

    /*! @brief Weights used for Gauss-Legendre quadrature */
    std::vector<double> _weights;

    /*! @brief Highest order coefficient of the polynomial */
    double _coefficient;

    /*! @brief Scale factors for 1D, 2D and 3D */
    double _scale_factor[3];

    static double polynomial1d(double x, int n);
    static double polynomial1d_prime(double x, int n);

    void calculate_zeros();
    void calculate_weights();

  public:
    Legendre(unsigned char order);

    unsigned char get_order();

    static double polynomial1d_scaled(double x, int n);
    static double polynomial1d_scaled_prime(double x, int n);

    double polynomial1d_recursive(double x);
    double polynomial1d_prime(double x);

    double polynomial1d_direct(double x);

    double polynomial1d_scaled(double x);
    double polynomial2d_scaled(double x, double y);
    double polynomial3d_scaled(double x, double y, double z);

    double quadrature1d(double (*f)(double));
    double quadrature2d(double (*f)(double, double));
    double quadrature3d(double (*f)(double, double, double));

    std::vector<double> get_zeros();
    std::vector<double> get_weights();
};

#endif  // LEGENDRE_HPP
