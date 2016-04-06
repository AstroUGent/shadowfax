/*******************************************************************************
 * This file is part of Shadowfax
 * Copyright (C) 2015 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file ICRegion.cpp
 *
 * @brief Region in a BlockICGenerator setup: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "ICRegion.hpp"
#include "SymbolicFunction.hpp"
#include <cmath>
#include <iostream>
using namespace std;

/**
 * @brief Constructor
 *
 * @param origin Origin of the region
 * @param sides Sidelengths of the region
 * @param exponent Exponent of the region
 * @param hydrofunctions Hydrodynamical SymbolicFunctions in the region
 * @param dmfunction Dark matter SymbolicFunctions in the region
 */
ICRegion::ICRegion(Vec origin, Vec sides, double exponent,
                   vector<string> hydrofunctions, vector<string> dmfunction)
        : _origin(origin), _sides(sides) {
    _exponent = exponent;
    if(hydrofunctions.size()) {
        _hydrofunctions.resize(ndim_ + 2);
        for(unsigned int i = 0; i < ndim_ + 2; i++) {
            _hydrofunctions[i] = new SymbolicFunction(hydrofunctions[i]);
        }
    }
    if(dmfunction.size()) {
        _dmfunction.resize(1);
        _dmfunction[0] = new SymbolicFunction(dmfunction[0]);
    }
}

/**
 * @brief Destructor
 *
 * Clean up the symbolic functions.
 */
ICRegion::~ICRegion() {
    if(_hydrofunctions.size()) {
        for(unsigned int i = 0; i < ndim_ + 2; i++) {
            delete _hydrofunctions[i];
        }
    }
    if(_dmfunction.size()) {
        delete _dmfunction[0];
    }
}

/**
 * @brief Test whether the point with given coordinates is inside the region
 *
 * If we call \f$\vec{o} = (o_x, o_y, o_z)\f$ the origin of the region,
 * \f$\vec{l} = (l_x, l_y, l_z)\f$ its sidelength and \f$e\f$ its exponents,
 * then a point with coordinates \f$\vec{x} = (x_x, x_y, x_z)\f$ lies inside the
 * region if
 * \f[\sqrt[1/e]{\left(2\frac{x_x-o_x}{l_x}\right)^e +
 *               \left(2\frac{x_y-o_y}{l_y}\right)^e +
 *               \left(2\frac{x_z-o_z}{l_z}\right)^e} \leq 1.
 * \f]
 *
 * @param position Position of the point
 * @return True if the point lies inside the region, false otherwise
 */
bool ICRegion::inside(Vec position) {
    double r = 0.;
    for(unsigned int k = ndim_; k--;) {
        double x = 2. * fabs(position[k] - _origin[k]) / _sides[k];
        if(_exponent < 10.) {
            r += pow(x, _exponent);
        } else {
            r = std::max(r, x);
        }
    }
    if(_exponent < 10.) {
        r = pow(r, 1. / _exponent);
    }
    return (r <= 1.);
}

/**
 * @brief Get the primitive variables at the given position in the region
 *
 * @warning This function does not check if the given position lies inside the
 * region!
 *
 * @param position Position in the region
 * @return Primitive variables Statevector
 */
StateVector ICRegion::get_hydro(Vec position) {
    StateVector hydro;
    if(_hydrofunctions.size()) {
        double r = 0.;
        for(unsigned int i = 0; i < ndim_; i++) {
            double x = position[i] - _origin[i];
            r += x * x;
        }
        r = sqrt(r);
        for(unsigned int i = 0; i < ndim_ + 2; i++) {
#if ndim_ == 3
            hydro[i] = (*_hydrofunctions[i])(r, position.x(), position.y(),
                                             position.z());
#else
            hydro[i] = (*_hydrofunctions[i])(r, position.x(), position.y());
#endif
        }
    }
    return hydro;
}

/**
 * @brief Integrate the given SymbolicFunction over this region
 *
 * The region is divided into at most 1,000,000 cells (10,000 in 2D) and
 * evaluated at the center of every cell (if the center is inside the region).
 * The integral is then given by the volume weigthed sum of all contributions.
 *
 * @warning This approach is not entirely correct...
 * @deprecated This function is not used
 *
 * @param function SymbolicFunction to evaluate
 * @param origin Origin of the region
 * @param sides Sidelengths of the region
 * @param exponent Exponent of the region
 * @return Integral of the SymbolicFunction over the region
 */
double ICRegion::integrate(SymbolicFunction* function, Vec origin, Vec sides,
                           double exponent) {
#if ndim_ == 3
    double integral = 0.;
    for(unsigned int i = 0; i < 100; i++) {
        for(unsigned int j = 0; j < 100; j++) {
            for(unsigned int z = 0; z < 100; z++) {
                double pos[3] = {(i + 0.5) * 0.01 * sides[0] - origin[0],
                                 (j + 0.5) * 0.01 * sides[1] - origin[1],
                                 (z + 0.5) * 0.01 * sides[2] - origin[2]};
                double re = 0.;
                for(unsigned int k = 0; k < ndim_; k++) {
                    double x = 2. * fabs(pos[k]) / sides[k];
                    if(exponent < 10.) {
                        re += pow(x, exponent);
                    } else {
                        re = std::max(re, x);
                    }
                }
                if(exponent < 10.) {
                    re = pow(re, 1. / exponent);
                }
                if(re <= 1.) {
                    re = sqrt(pos[0] * pos[0] + pos[1] * pos[1] +
                              pos[2] * pos[2]);
                    integral +=
                            1.e-6 * ((*function)(re, pos[0], pos[1], pos[2]));
                    _volume += 1.e-6;
                }
            }
        }
    }
    return integral;
#else
    double integral = 0.;
    for(unsigned int i = 0; i < 100; i++) {
        for(unsigned int j = 0; j < 100; j++) {
            double pos[2] = {(i + 0.5) * 0.01 * sides[0] - origin[0],
                             (j + 0.5) * 0.01 * sides[1] - origin[1]};
            double re = 0.;
            for(unsigned int k = 0; k < ndim_; k++) {
                double x = 2. * fabs(pos[k]) / sides[k];
                if(exponent < 10.) {
                    re += pow(x, exponent);
                } else {
                    re = std::max(re, x);
                }
            }
            if(exponent < 10.) {
                re = pow(re, 1. / exponent);
            }
            if(re <= 1.) {
                re = sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
                integral += 1.e-4 * ((*function)(re, pos[0], pos[1]));
                _volume += 1.e-4;
            }
        }
    }
    return integral;
#endif
}

/**
 * @brief Calculate the maximal hydrodynamical density in the region
 *
 * We divide the region in at most 1,000,000 cells (10,000 in 2D) and calculate
 * the maximal value of the density of all centers of cells that lie inside the
 * region, but not inside the given cut out region.
 *
 * @param cut_out_region ICRegion that (partially) overlaps with this region and
 * that should not be considered when calculating the maximum
 * @return The maximal hydrodynamical density inside the region
 */
double ICRegion::get_max_value_hydro(ICRegion* cut_out_region) {
    if(!_hydrofunctions.size()) {
        return 0;
    }
#if ndim_ == 3
    _max_value_hydro = 0.;
    for(unsigned int i = 0; i < 100; i++) {
        for(unsigned int j = 0; j < 100; j++) {
            for(unsigned int z = 0; z < 100; z++) {
                Vec pos((i + 0.5) * 0.01 * _sides[0] - _origin[0],
                        (j + 0.5) * 0.01 * _sides[1] - _origin[1],
                        (z + 0.5) * 0.01 * _sides[2] - _origin[2]);
                if(!cut_out_region || !cut_out_region->inside(pos)) {
                    double re = 0.;
                    for(unsigned int k = 0; k < ndim_; k++) {
                        double x = 2. * fabs(pos[k]) / _sides[k];
                        if(_exponent < 10.) {
                            re += pow(x, _exponent);
                        } else {
                            re = std::max(re, x);
                        }
                    }
                    if(_exponent < 10.) {
                        re = pow(re, 1. / _exponent);
                    }
                    if(re <= 1.) {
                        re = sqrt(pos[0] * pos[0] + pos[1] * pos[1] +
                                  pos[2] * pos[2]);
                        _max_value_hydro =
                                std::max(_max_value_hydro,
                                         (*_hydrofunctions[0])(re, pos[0],
                                                               pos[1], pos[2]));
                    }
                }
            }
        }
    }
    return _max_value_hydro;
#else
    _max_value_hydro = 0.;
    for(unsigned int i = 0; i < 100; i++) {
        for(unsigned int j = 0; j < 100; j++) {
            Vec pos((i + 0.5) * 0.01 * _sides[0] - _origin[0],
                    (j + 0.5) * 0.01 * _sides[1] - _origin[1]);
            if(!cut_out_region || !cut_out_region->inside(pos)) {
                double re = 0.;
                for(unsigned int k = 0; k < ndim_; k++) {
                    double x = 2. * fabs(pos[k]) / _sides[k];
                    if(_exponent < 10.) {
                        re += pow(x, _exponent);
                    } else {
                        re = std::max(re, x);
                    }
                }
                if(_exponent < 10.) {
                    re = pow(re, 1. / _exponent);
                }
                if(re <= 1.) {
                    re = sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
                    _max_value_hydro =
                            std::max(_max_value_hydro,
                                     (*_hydrofunctions[0])(re, pos[0], pos[1]));
                }
            }
        }
    }
    return _max_value_hydro;
#endif
}

/**
 * @brief Calculate the maximal dark matter density in the region
 *
 * We divide the region in at most 1,000,000 cells (10,000 in 2D) and calculate
 * the maximal value of the density of all centers of cells that lie inside the
 * region, but not inside the given cut out region.
 *
 * @param cut_out_region ICRegion that (partially) overlaps with this region and
 * that should not be considered when calculating the maximum
 * @return The maximal dark matter density inside the region
 */
double ICRegion::get_max_value_dm(ICRegion* cut_out_region) {
    if(!_dmfunction.size()) {
        return 0;
    }
#if ndim_ == 3
    _max_value_dm = 0.;
    for(unsigned int i = 0; i < 100; i++) {
        for(unsigned int j = 0; j < 100; j++) {
            for(unsigned int z = 0; z < 100; z++) {
                Vec pos((i + 0.5) * 0.01 * _sides[0] - _origin[0],
                        (j + 0.5) * 0.01 * _sides[1] - _origin[1],
                        (z + 0.5) * 0.01 * _sides[2] - _origin[2]);
                if(!cut_out_region || !cut_out_region->inside(pos)) {
                    double re = 0.;
                    for(unsigned int k = 0; k < ndim_; k++) {
                        double x = 2. * fabs(pos[k]) / _sides[k];
                        if(_exponent < 10.) {
                            re += pow(x, _exponent);
                        } else {
                            re = std::max(re, x);
                        }
                    }
                    if(_exponent < 10.) {
                        re = pow(re, 1. / _exponent);
                    }
                    if(re <= 1.) {
                        re = sqrt(pos[0] * pos[0] + pos[1] * pos[1] +
                                  pos[2] * pos[2]);
                        _max_value_dm = std::max(
                                _max_value_dm,
                                (*_dmfunction[0])(re, pos[0], pos[1], pos[2]));
                    }
                }
            }
        }
    }
    return _max_value_dm;
#else
    _max_value_dm = 0.;
    for(unsigned int i = 0; i < 100; i++) {
        for(unsigned int j = 0; j < 100; j++) {
            Vec pos((i + 0.5) * 0.01 * _sides[0] - _origin[0],
                    (j + 0.5) * 0.01 * _sides[1] - _origin[1]);
            if(!cut_out_region || !cut_out_region->inside(pos)) {
                double re = 0.;
                for(unsigned int k = 0; k < ndim_; k++) {
                    double x = 2. * fabs(pos[k]) / _sides[k];
                    if(_exponent < 10.) {
                        re += pow(x, _exponent);
                    } else {
                        re = std::max(re, x);
                    }
                }
                if(_exponent < 10.) {
                    re = pow(re, 1. / _exponent);
                }
                if(re <= 1.) {
                    re = sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
                    _max_value_dm =
                            std::max(_max_value_dm,
                                     (*_dmfunction[0])(re, pos[0], pos[1]));
                }
            }
        }
    }
    return _max_value_dm;
#endif
}

/**
 * @brief Check whether the given position inside this region is acceptable in a
 * Monte Carlo rejection sampling of the density field for the hydro
 *
 * @param position Position inside the region
 * @return True if the position is accepted, false if it is rejected
 */
bool ICRegion::accept_hydro(Vec position) {
    Vec p = position - _origin;
#if ndim_ == 3
    return ((double)rand()) / ((double)RAND_MAX) <
           ((*_hydrofunctions[0])(p.norm(), p.x(), p.y(), p.z())) /
                   _max_value_hydro;
#else
    return ((double)rand()) / ((double)RAND_MAX) <
           ((*_hydrofunctions[0])(p.norm(), p.x(), p.y())) / _max_value_hydro;
#endif
}

/**
 * @brief Check whether the given position inside this region is acceptable in a
 * Monte Carlo rejection sampling of the density field for the dark matter
 *
 * @param position Position inside the region
 * @return True if the position is accepted, false if it is rejected
 */
bool ICRegion::accept_dm(Vec position) {
    Vec p = position - _origin;
#if ndim_ == 3
    return ((double)rand()) / ((double)RAND_MAX) <
           ((*_dmfunction[0])(p.norm(), p.x(), p.y(), p.z())) / _max_value_dm;
#else
    return ((double)rand()) / ((double)RAND_MAX) <
           ((*_dmfunction[0])(p.norm(), p.x(), p.y())) / _max_value_dm;
#endif
}
