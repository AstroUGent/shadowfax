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
 * @file Legendre.cpp
 *
 * @brief Legendre polynomials: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "Legendre.hpp"
#include <cmath>

/**
 * @brief Static Legendre polynomial of a given order
 *
 * Recursively defined using
 * \f{align*}{
 *   P_0(x) &= 1 \\
 *   P_1(x) &= x \\
 *   P_n(x) &= \frac{2n-1}{n}xP_{n-1}(x) - \frac{n-1}{n}P_{n-2}(x)
 * \f}
 *
 * @param x Real variable
 * @param n Order of the polynomial
 * @return Value of the Legendre polynomial of order n at coordinate x
 */
double Legendre::polynomial1d(double x, int n) {
    if(!n) {
        return 1.;
    }
    if(n == 1) {
        return x;
    }
    double Pnm1 = polynomial1d(x, n - 1);
    double Pnm2 = polynomial1d(x, n - 2);
    return (2. * n - 1.) / n * x * Pnm1 - (n - 1.) / n * Pnm2;
}

/**
 * @brief Static first derivative of the Legendre polynomial of given order
 *
 * Defined using
 * \f[
 *   \frac{\text{d}}{\text{d}x}P_n(x) = \frac{nx}{x^2-1}P_n(x) -
 *                                      \frac{n}{x^2-1}P_{n-1}(x)
 * \f]
 *
 * @param x Real variable
 * @param n Order of the polynomial
 * @return Value of the first derivative of the Legendre polynomial of order n
 * at coordinate x
 */
double Legendre::polynomial1d_prime(double x, int n) {
    if(!n) {
        return 0.;
    }
    double Pn = polynomial1d(x, n);
    double Pnm1 = polynomial1d(x, n - 1);
    return n * x / (x * x - 1.) * Pn - n * Pnm1 / (x * x - 1.);
}

/**
 * @brief Pre-calculate the roots of the Legendre polynomial
 *
 * Roots are determined using a Newton-Raphson method.
 */
void Legendre::calculate_zeros() {
    _zeros.resize(_order);
    double o2 = _order * _order;
    double o3 = o2 * _order;
    for(unsigned char q = 1; q < _order + 1; q++) {
        double guess = (1. - 0.125 / o2 + 0.125 / o3) *
                       cos((4. * q - 1.) * M_PI / (4. * _order + 2.));
        double fguess = polynomial1d(guess, _order);
        double fguess_prime = polynomial1d_prime(guess, _order);
        double guess2 = guess - fguess / fguess_prime;
        while(fabs(guess - guess2) > 1.e-7 * fabs(guess + guess2)) {
            guess = guess2;
            fguess = polynomial1d(guess, _order);
            fguess_prime = polynomial1d_prime(guess, _order);
            guess2 = guess - fguess / fguess_prime;
        }
        _zeros[q - 1] = guess2;
    }
}

/**
 * @brief 1D scaled Legendre polynomial
 *
 * @param x Real value
 * @param n Order
 * @return Value of the scaled 1D Legendre polynomial of order n at coordinate x
 */
double Legendre::polynomial1d_scaled(double x, int n) {
    return sqrt(2. * n + 1.) * polynomial1d(x, n);
}

/**
 * @brief First derivative of the 1D scaled Legendre polynomial
 *
 * @param x Real value
 * @param n Order
 * @return Value of the first derivative of the scaled 1D Legendre polynomial of
 * order n at coordinate x
 */
double Legendre::polynomial1d_scaled_prime(double x, int n) {
    return sqrt(2. * n + 1.) * polynomial1d_prime(x, n);
}

/**
 * @brief Pre-calculate the weights used for Gauss-Legendre quadrature
 */
void Legendre::calculate_weights() {
    _weights.resize(_order);
    for(unsigned char i = 0; i < _order; i++) {
        double lprime = polynomial1d_prime(_zeros[i], _order);
        _weights[i] = 2. / (1. - _zeros[i] * _zeros[i]) / lprime / lprime;
    }
}

/**
 * @brief Constructor
 *
 * Pre-calculates roots and weights for Gauss-Legendre quadrature.
 *
 * @param order Order of the polynomial
 */
Legendre::Legendre(unsigned char order) {
    _order = order;
    calculate_zeros();
    calculate_weights();
    _coefficient = 1.;
    _coefficient = polynomial1d_recursive(2.) / polynomial1d_direct(2.);

    _scale_factor[1] = 2. * _order + 1.;
    _scale_factor[0] = sqrt(_scale_factor[1]);
    _scale_factor[2] = _scale_factor[0] * _scale_factor[1];
}

/**
 * @brief Get the order of the Legendre polynomials
 *
 * @return Order
 */
unsigned char Legendre::get_order() {
    return _order;
}

/**
 * @brief Member wrapper around static Legendre function
 *
 * @param x Real value
 * @return Value of the Legendre polynomial at coordinate x
 */
double Legendre::polynomial1d_recursive(double x) {
    return polynomial1d(x, _order);
}

/**
 * @brief Member wrapper around static Legendre derivative function
 *
 * @param x Real value
 * @return Value of the first derivative of the Legendre polynomial at
 * coordinate x
 */
double Legendre::polynomial1d_prime(double x) {
    return polynomial1d_prime(x, _order);
}

/**
 * @brief Calculate the polynomial directly by using the coefficient of the
 * highest order term and its roots
 *
 * @param x Real value
 * @return Value of the Legende polynomial at coordinate x
 */
double Legendre::polynomial1d_direct(double x) {
    double ans = _coefficient;
    for(unsigned char i = 0; i < _order; i++) {
        ans *= (x - _zeros[i]);
    }
    return ans;
}

/**
 * @brief 1D scaled Legendre polynomial
 *
 * @param x Real value
 * @return Value of the scaled 1D Legendre polynomial at coordinate x
 */
double Legendre::polynomial1d_scaled(double x) {
    return _scale_factor[0] * polynomial1d_direct(x);
}

/**
 * @brief 2D scaled Legendre polynomial
 *
 * @param x Real value
 * @param y Real value
 * @return Value of the scaled 2D Legendre polynomial at coordinates x and y
 */
double Legendre::polynomial2d_scaled(double x, double y) {
    return _scale_factor[1] * polynomial1d_direct(x) * polynomial1d_direct(y);
}

/**
 * @brief 3D scaled Legendre polynomial
 *
 * @param x Real value
 * @param y Real value
 * @param z Real value
 * @return Value of the scaled 3D Legendre polynomial at coordinates x, y and z
 */
double Legendre::polynomial3d_scaled(double x, double y, double z) {
    return _scale_factor[2] * polynomial1d_direct(x) * polynomial1d_direct(y) *
           polynomial1d_direct(z);
}

/**
 * @brief Calculate the 1D quadrature for the given function
 *
 * @param f 1D function to integrate
 * @return Integral of the function over the interval \f$[-1,1]\f$
 */
double Legendre::quadrature1d(double (*f)(double)) {
    double ans = 0.;
    for(unsigned char i = 0; i < _order; i++) {
        ans += _weights[i] * f(_zeros[i]);
    }
    return ans;
}

/**
 * @brief Calculate the 2D quadrature for the given function
 *
 * @param f 2D function to integrate
 * @return Integral of the function over the interval \f$[-1,1]\times{}[-1,
 * 1]\f$
 */
double Legendre::quadrature2d(double (*f)(double, double)) {
    double ans = 0.;
    for(unsigned char i = 0; i < _order; i++) {
        for(unsigned char j = 0; j < _order; j++) {
            ans += _weights[i] * _weights[j] * f(_zeros[i], _zeros[j]);
        }
    }
    return ans;
}

/**
 * @brief Calculate the 3D quadrature for the given function
 *
 * @param f 3D function to integrate
 * @return Integral of the function over the interval \f$[-1,1]\times{}[-1,
 * 1]\times{}[-1, 1]\f$
 */
double Legendre::quadrature3d(double (*f)(double, double, double)) {
    double ans = 0.;
    for(unsigned char i = 0; i < _order; i++) {
        for(unsigned char j = 0; j < _order; j++) {
            for(unsigned char k = 0; k < _order; k++) {
                ans += _weights[i] * _weights[j] * _weights[k] *
                       f(_zeros[i], _zeros[j], _zeros[k]);
            }
        }
    }
    return ans;
}

/**
 * @brief Access the roots of the polynomial
 *
 * For debugging purposes.
 *
 * @return The roots of the Legendre polynomial
 */
std::vector<double> Legendre::get_zeros() {
    return _zeros;
}

/**
 * @brief Access the weights of the polynomial
 *
 * For debugging purposes.
 *
 * @return The weights of the Legendre polynomial
 */
std::vector<double> Legendre::get_weights() {
    return _weights;
}
