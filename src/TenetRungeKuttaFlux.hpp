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
 * @file TenetRungeKuttaFlux.hpp
 *
 * @brief Runge-Kutta flux used for Tenet
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef TENETRUNGEKUTTAFLUX_HPP
#define TENETRUNGEKUTTAFLUX_HPP

#include "Legendre.hpp"
#include "RungeKuttaFlux.hpp"
#include "StateVector.hpp"
#include "TenetBasisFunctions.hpp"
#include "TenetCellWeights.hpp"
#include "TenetGrid.hpp"
#include "Vec.hpp"
#include "riemann/RiemannSolver.hpp"
#include <vector>

/**
 * @brief Runge-Kutta flux for Tenet
 */
template <class SystemState>
class TenetRungeKuttaFlux : public RungeKuttaFlux<SystemState> {
  private:
    /*! @brief Reference to the TenetGrid on which we integrate */
    TenetGrid& _grid;

    /*! @brief Reference to the Legendre polynomials used as basis for the
     * integration */
    Legendre& _legendre;

    /*! @brief TenetBasisFunctions used as a basis for the integration */
    TenetBasisFunctions& _basis;

    /*! @brief Pointer to a RiemannSolver used for flux calculations */
    RiemannSolver* _solver;

    /*! @brief Roots of the Legendre polynomials, used for quadrature */
    std::vector<double> _zeros;

    /*! @brief Weights used for Gauss-Legendre quadrature */
    std::vector<double> _weights;

  public:
    /**
     * @brief Constructor
     *
     * @param grid Reference to the TenetGrid
     * @param legendre Reference to the Legendre polynomials
     * @param basis Reference to the TenetBasisFunctions
     * @param solver Pointer to a RiemannSolver
     */
    TenetRungeKuttaFlux(TenetGrid& grid, Legendre& legendre,
                        TenetBasisFunctions& basis, RiemannSolver* solver)
            : _grid(grid), _legendre(legendre), _basis(basis) {
        _solver = solver;
        _zeros = legendre.get_zeros();
        _weights = legendre.get_weights();
    }

    /**
     * @brief Get the Runge-Kutta flux corresponding to the Tenet weight flux
     *
     * @param t Current system time (not used)
     * @param x Current state of the system
     * @return Flux
     */
    virtual SystemState get_flux(double t, SystemState& x) {
        SystemState dweights(x.size());
        double gamma = _solver->get_gamma();
        double cellwidth = _grid.get_cellwidth();
        double cellvolume = _grid.get_cellvolume();
        Vec vface;
        for(auto it = _grid.begin(); it != _grid.end(); ++it) {
            unsigned int i = it.index();
            dweights[i] = _basis.get_empty_weights();
            for(auto bit = _basis.begin(); bit != _basis.end(); ++bit) {
                for(auto cit = _basis.cell_quadrature_begin();
                    cit != _basis.cell_quadrature_end(); ++cit) {
                    Vec ksi = cit.get_ksi();
                    StateVector U = it->get_variables(ksi, x[i], _basis);

#if ndim_ == 3
                    StateVector fksi[3];
                    double phij[3];
                    phij[0] = bit.prime_x(ksi);
                    phij[1] = bit.prime_y(ksi);
                    phij[2] = bit.prime_z(ksi);

                    double p =
                            (gamma - 1.) *
                            (U.e() -
                             0.5 * (U[1] * U[1] + U[2] * U[2] + U[3] * U[3]) /
                                     U[0]);
                    // rho*v_x = U[1]
                    fksi[0][0] = U[1];
                    // rho*v_x^2 + p = U[1]^2/U[0] + p
                    fksi[0][1] = U[1] * U[1] / U[0] + p;
                    fksi[0][2] = U[1] * U[2] / U[0];
                    fksi[0][3] = U[1] * U[3] / U[0];
                    fksi[0][4] = (U[4] + p) * U[1] / U[0];

                    fksi[1][0] = U[2];
                    fksi[1][1] = U[1] * U[2] / U[0];
                    fksi[1][2] = U[2] * U[2] / U[0] + p;
                    fksi[1][3] = U[3] * U[3] / U[0];
                    fksi[1][4] = (U[4] + p) * U[2] / U[0];

                    fksi[2][0] = U[3];
                    fksi[2][1] = U[1] * U[3] / U[0];
                    fksi[2][2] = U[2] * U[3] / U[0];
                    fksi[2][3] = U[3] * U[3] / U[0] + p;
                    fksi[2][4] = (U[4] + p) * U[3] / U[0];

                    StateVector phij_fksi = fksi[0] * phij[0] +
                                            fksi[1] * phij[1] +
                                            fksi[2] * phij[2];
#else
                    StateVector fksi[2];
                    double phij[2];
                    phij[0] = bit.prime_x(ksi);
                    phij[1] = bit.prime_y(ksi);

                    double p =
                            (gamma - 1.) *
                            (U.e() - 0.5 * (U[1] * U[1] + U[2] * U[2]) / U[0]);
                    // rho*v_x = U[1]
                    fksi[0][0] = U[1];
                    // rho*v_x^2 + p = U[1]^2/U[0] + p
                    fksi[0][1] = U[1] * U[1] / U[0] + p;
                    fksi[0][2] = U[1] * U[2] / U[0];
                    fksi[0][3] = (U[3] + p) * U[1] / U[0];

                    fksi[1][0] = U[2];
                    fksi[1][1] = U[1] * U[2] / U[0];
                    fksi[1][2] = U[2] * U[2] / U[0] + p;
                    fksi[1][3] = (U[3] + p) * U[2] / U[0];

                    StateVector phij_fksi =
                            fksi[0] * phij[0] + fksi[1] * phij[1];
#endif

                    dweights[i][bit.index()] += phij_fksi * cit.get_weight();
                }
            }
        }

        for(auto it = _grid.faces_begin(); it != _grid.faces_end(); ++it) {
            unsigned int iL = it.get_index_L();
            unsigned int iR = it.get_index_R();
            Vec midpoint = it.get_midpoint();
            Vec normal = it.get_normal();
            for(auto fqit = _basis.face_quadrature_begin(midpoint, normal);
                fqit != _basis.face_quadrature_end(midpoint, normal); ++fqit) {
                Vec ksiL = fqit.get_ksiL();
                Vec ksiR = fqit.get_ksiR();
                StateVector UL = _grid[iL]->get_variables(ksiL, x[iL], _basis);
                StateVector UR = _grid[iR]->get_variables(ksiR, x[iR], _basis);
                StateVector WL, WR;

#if ndim_ == 3
                WL[0] = UL[0];
                WL[1] = UL[1] / UL[0];
                WL[2] = UL[2] / UL[0];
                WL[3] = UL[3] / UL[0];
                WL[4] = (gamma - 1.) *
                        (UL[4] -
                         0.5 * (UL[1] * UL[1] + UL[2] * UL[2] + UL[3] * UL[3]) /
                                 UL[0]);

                WR[0] = UR[0];
                WR[1] = UR[1] / UR[0];
                WR[2] = UR[2] / UR[0];
                WR[3] = UR[3] / UR[0];
                WR[4] = (gamma - 1.) *
                        (UR[4] -
                         0.5 * (UR[1] * UR[1] + UR[2] * UR[2] + UR[3] * UR[3]) /
                                 UR[0]);

                if(WL[0] < 0. || WL[4] < 0. || WR[0] < 0. || WR[4] < 0.) {
                    cerr << "Negative positive quantity!" << endl;
                    my_exit();
                }
#else
                WL[0] = UL[0];
                WL[1] = UL[1] / UL[0];
                WL[2] = UL[2] / UL[0];
                WL[3] = (gamma - 1.) *
                        (UL[3] - 0.5 * (UL[1] * UL[1] + UL[2] * UL[2]) / UL[0]);

                WR[0] = UR[0];
                WR[1] = UR[1] / UR[0];
                WR[2] = UR[2] / UR[0];
                WR[3] = (gamma - 1.) *
                        (UR[3] - 0.5 * (UR[1] * UR[1] + UR[2] * UR[2]) / UR[0]);

                if(WL[0] < 0. || WL[3] < 0. || WR[0] < 0. || WR[3] < 0.) {
                    cerr << "Negative positive quantity!" << endl;
                    my_exit();
                }
#endif

                StateVector flux =
                        _solver->solve_for_flux(WL, WR, normal, vface);
                for(auto bit = _basis.begin(); bit != _basis.end(); ++bit) {
                    dweights[iL][bit.index()] -=
                            flux * bit(ksiL) * fqit.get_weight();
                    dweights[iR][bit.index()] +=
                            flux * bit(ksiR) * fqit.get_weight();
                }
            }
        }

        for(auto it = _grid.begin(); it != _grid.end(); ++it) {
            unsigned int i = it.index();
            for(auto bit = _basis.begin(); bit != _basis.end(); ++bit) {
#if ndim_ == 3
                dweights[i][bit.index()] *=
                        0.25 * cellwidth * cellwidth / cellvolume;
#else
                dweights[i][bit.index()] *= 0.5 * cellwidth / cellvolume;
#endif
            }
        }

        return dweights;
    }
};

#endif  // TENETRUNGEKUTTAFLUX_HPP
