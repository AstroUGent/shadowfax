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
 * @file testRungeKutta.cpp
 *
 * @brief Unit test for the Runge-Kutta integator
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "RungeKutta.hpp"
#include "RungeKuttaFlux.hpp"
#include <cmath>
#include <iostream>

/**
 * @brief Basic implementation of a RungeKuttaFlux for a two body system
 */
template <class SystemState>
class TwoBodyRungeKuttaFlux : public RungeKuttaFlux<SystemState> {
  public:
    /**
     * @brief Get the two body flux
     *
     * The flux in this case is given by the velocities and accelerations of the
     * two bodies.
     *
     * @param t Current time of the system (not used)
     * @param x Current state of the system
     * @return Flux of the system
     */
    virtual SystemState get_flux(double t, SystemState& x) {
        double force[2];
        double r[4];
        r[0] = x[0] - x[2];
        r[1] = x[1] - x[3];
        r[2] = r[0] * r[0] + r[1] * r[1];
        r[3] = sqrt(r[2]);
        force[0] = r[0] / r[2] / r[3];
        force[1] = r[1] / r[2] / r[3];
        SystemState dx(x);
        dx[0] = x[4];
        dx[1] = x[5];
        dx[2] = x[6];
        dx[3] = x[7];
        dx[4] = -force[0] * x.m(1);
        dx[5] = -force[1] * x.m(1);
        dx[6] = force[0] * x.m(0);
        dx[7] = force[1] * x.m(0);
        return dx;
    }
};

/**
 * @brief SystemState for the two body system
 */
class TwoBodySystem {
  private:
    /*! @brief Positions and velocities of the two bodies */
    double _vars[8];
    /*! @brief Masses of the two bodies */
    double _m[2];

  public:
    /**
     * @brief Constructor
     *
     * @param m1 Mass of body 1
     * @param m2 Mass of body 2
     * @param x1 2D position of body 1
     * @param x2 2D position of body 2
     * @param v1 2D velocity of body 1
     * @param v2 2D velocity of body 2
     */
    TwoBodySystem(double m1, double m2, double* x1, double* x2, double* v1,
                  double* v2) {
        _m[0] = m1;
        _m[1] = m2;
        _vars[0] = x1[0];
        _vars[1] = x1[1];
        _vars[2] = x2[0];
        _vars[3] = x2[1];
        _vars[4] = v1[0];
        _vars[5] = v1[1];
        _vars[6] = v2[0];
        _vars[7] = v2[1];
    }

    /**
     * @brief Dereference operator for positions and velocities
     *
     * @param i Index of the component to access
     * @return Reference to the requested component
     */
    double& operator[](unsigned int i) {
        return _vars[i];
    }

    /**
     * @brief Get the mass of one of the bodies
     *
     * @param i Index of one of the bodies
     * @return Mass of the body
     */
    double m(unsigned int i) {
        return _m[i];
    }

    /**
     * @brief Multiplication operator
     *
     * Only affects the positions and velocities
     *
     * @param s Scalar to multiply with
     * @return Reference to the resulting TwoBodySystem
     */
    TwoBodySystem& operator*=(double s) {
        _vars[0] *= s;
        _vars[1] *= s;
        _vars[2] *= s;
        _vars[3] *= s;
        _vars[4] *= s;
        _vars[5] *= s;
        _vars[6] *= s;
        _vars[7] *= s;
        return *this;
    }

    /**
     * @brief Addition operator
     *
     * Only affects the positions and velocities
     *
     * @param s Scalar to add to all components
     * @return Reference to the resulting TwoBodySystem
     */
    TwoBodySystem& operator+=(const TwoBodySystem& s) {
        _vars[0] += s._vars[0];
        _vars[1] += s._vars[1];
        _vars[2] += s._vars[2];
        _vars[3] += s._vars[3];
        _vars[4] += s._vars[4];
        _vars[5] += s._vars[5];
        _vars[6] += s._vars[6];
        _vars[7] += s._vars[7];
        return *this;
    }

    /**
     * @brief Negation operator
     *
     * Only affects the positions and velocities
     *
     * @return New TwoBodySystem with positions and velocities that are opposite
     * to those of the current instance
     */
    TwoBodySystem operator-() {
        TwoBodySystem s(*this);
        s._vars[0] = -_vars[0];
        s._vars[1] = -_vars[1];
        s._vars[2] = -_vars[2];
        s._vars[3] = -_vars[3];
        s._vars[4] = -_vars[4];
        s._vars[5] = -_vars[5];
        s._vars[6] = -_vars[6];
        s._vars[7] = -_vars[7];
        return s;
    }
};

/**
 * @brief Free multiplication operator for a scalar and a TwoBodySystem
 *
 * @param s Scalar
 * @param b TwoBodySystem
 * @return Result of the multiplication
 */
TwoBodySystem operator*(double s, TwoBodySystem b) {
    return b *= s;
}

/**
 * @brief Free addition operator for two TwoBodySystem instances
 *
 * @param a TwoBodySystem
 * @param b TwoBodySystem
 * @return Result of the addition
 */
TwoBodySystem operator+(TwoBodySystem a, TwoBodySystem b) {
    return a += b;
}

/**
 * @brief Runge-Kutta integrator test
 *
 * Applies the Runge-Kutta integrator to a simple two-body problem.
 *
 * @param argc Number of command line arguments
 * @param argv Command line arguments
 * @return Exit code
 */
int main(int argc, char** argv) {
    TwoBodyRungeKuttaFlux<TwoBodySystem> flux;

    double x1[2] = {0., 0.};
    double x2[2] = {10., 0.};
    double v1[2] = {0., 0.};
    double v2[2] = {0., 10.};
    TwoBodySystem system(1000., 1., x1, x2, v1, v2);

    RungeKutta<TwoBodySystem> rk(4, flux, system);

    std::cout << "x1\ty1\tx2\ty2" << std::endl;
    for(unsigned int i = 0; i < 100; i++) {
        system = rk.integrate(0.1);
        std::cout << system[0] << "\t" << system[1] << "\t" << system[2] << "\t"
                  << system[3] << std::endl;
    }

    return 0;
}
