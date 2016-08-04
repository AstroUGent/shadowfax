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
 * @file EwaldTable.cpp
 *
 * @brief Ewald tables used for periodic force corrections: header
 *
 * This code is largely based on similar code in Gadget2 (Springel 2005).
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "EwaldTable.hpp"
#include "EwaldTableLocation.hpp"
#include <fstream>
#include <iostream>  // for operator<<, basic_ostream, etc
#include <string>    // for string
using namespace std;

/**
  * @brief Constructor
  *
  * Initialize parameters and allocate memory to store the internal table.
  * Read in or calculate the table (if it does not yet exist) and initialize the
  * conversion factor.
  *
  * @param alpha The Ewald alpha parameter, which specifies the cutoff between
  * long range and short range parts of the force
  * @param size The integer size of one row in the multidimensional table. The
  * actual size is size+1, since both lower and upper limit are inclusive.
  */
EwaldTable::EwaldTable(double alpha, unsigned int size) {
    _size = size;
    _alpha = alpha;
#if ndim_ == 3
    _f = new Vec**[size + 1];
    for(unsigned int i = 0; i < size + 1; i++) {
        _f[i] = new Vec*[size + 1];
        for(unsigned int j = 0; j < size + 1; j++) {
            _f[i][j] = new Vec[size + 1];
        }
    }
#else
    _f = new Vec*[size + 1];
    for(unsigned int i = 0; i < size + 1; i++) {
        _f[i] = new Vec[size + 1];
    }
#endif
    if(!read_table()) {
        // Show some information to the user. This way, the user does not think
        // the program crashed or hangs.
        // We're not the Belgian railway company, you know...
        cout << "Pre-calculated Ewald tables not found. Calculating them."
             << endl;
        cout << "This can take some time (~1 min)" << endl;
        construct();
    }
    _ewald_intp_fac = 2. * _size;
}

/**
  * @brief Garbage collection. Free memory associated with the internal table.
  */
EwaldTable::~EwaldTable() {
#if ndim_ == 3
    for(unsigned int i = 0; i < _size + 1; i++) {
        for(unsigned int j = 0; j < _size + 1; j++) {
            delete[] _f[i][j];
        }
        delete[] _f[i];
    }
    delete[] _f;
#else
    for(unsigned int i = 0; i < _size + 1; i++) {
        delete[] _f[i];
    }
    delete[] _f;
#endif
}

/**
  * @brief Reset the boxsize from generic side 1 to the given length
  *
  * Internally, this means we divide the conversion factor by the new length and
  * we divide all force corrections by the square of the new length.
  *
  * @param L The side of the periodic box.
  */
void EwaldTable::set_boxsize(double L) {
    _ewald_intp_fac /= L;
#if ndim_ == 3
    for(unsigned int i = _size + 1; i--;) {
        for(unsigned int j = _size + 1; j--;) {
            for(unsigned int k = _size + 1; k--;) {
                for(unsigned int n = 3; n--;) {
                    _f[i][j][k][n] /= (L * L);
                }
            }
        }
    }
#else
    for(unsigned int i = _size + 1; i--;) {
        for(unsigned int j = _size + 1; j--;) {
            for(unsigned int n = 2; n--;) {
                _f[i][j][n] /= (L * L);
            }
        }
    }
#endif
}

/**
  * @brief Read a dump of previously calculated Ewald corrections if it exists
  *
  * @return True if the dump exists and the Ewald table was succesfully filled,
  * false otherwise
  */
bool EwaldTable::read_table() {
#if ndim_ == 3
    string path = string(EWALDTABLE_LOCATION) + string("ewaldtable3d.dat");
    ifstream file(path, ios::in | ios::binary);
    if(file) {
        for(unsigned int i = _size + 1; i--;) {
            for(unsigned int j = _size + 1; j--;) {
                for(unsigned int k = _size + 1; k--;) {
                    file.read(reinterpret_cast<char*>(&_f[i][j][k][0]),
                              sizeof(_f[i][j][k][0]));
                    file.read(reinterpret_cast<char*>(&_f[i][j][k][1]),
                              sizeof(_f[i][j][k][1]));
                    file.read(reinterpret_cast<char*>(&_f[i][j][k][2]),
                              sizeof(_f[i][j][k][2]));
                }
            }
        }
        return true;
    } else {
        return false;
    }
#else
    string path = string(EWALDTABLE_LOCATION) + string("ewaldtable2d.dat");
    ifstream file(path, ios::in | ios::binary);
    if(file) {
        for(unsigned int i = _size + 1; i--;) {
            for(unsigned int j = _size + 1; j--;) {
                file.read(reinterpret_cast<char*>(&_f[i][j][0]),
                          sizeof(_f[i][j][0]));
                file.read(reinterpret_cast<char*>(&_f[i][j][1]),
                          sizeof(_f[i][j][1]));
            }
        }
        return true;
    } else {
        return false;
    }
#endif
}

/**
  * @brief Calculate Ewald force corrections and store them in the internal
  * table
  *
  * @warning This scales as the third (or second) power of the table size and
  * hence can take quite some time.
  *
  * The resulting force corrections are dumped in a file so they may be used
  * again in a future run.
  */
void EwaldTable::construct() {
#if ndim_ == 3
    Vec n;
    Vec h;
    for(unsigned int i = _size + 1; i--;) {
        for(unsigned int j = _size + 1; j--;) {
            for(unsigned int k = _size + 1; k--;) {
                if(!i && !j && !k) {
                    // all forces are initialized to zero vectors, no need to do
                    // it here
                    continue;
                }
                Vec pos(0.5 * i / (double)_size, 0.5 * j / (double)_size,
                        0.5 * k / (double)_size);
                _f[i][j][k] = pos / (pos.norm() * pos.norm2());
                for(n[0] = -4.; n[0] <= 4.; n[0]++) {
                    for(n[1] = -4.; n[1] <= 4.; n[1]++) {
                        for(n[2] = -4.; n[2] <= 4.; n[2]++) {
                            Vec dx = pos - n;
                            double r = dx.norm();
                            double r2 = dx.norm2();
                            _f[i][j][k] -= dx / (r * r2) *
                                           (erfc(_alpha * r) +
                                            2. * _alpha * r / sqrt(M_PI) *
                                                    exp(-_alpha * _alpha * r2));
                        }
                    }
                }
                for(h[0] = -4.; h[0] <= 4.; h[0]++) {
                    for(h[1] = -4.; h[1] <= 4.; h[1]++) {
                        for(h[2] = -4.; h[2] <= 4.; h[2]++) {
                            double h2 = h.norm2();
                            if(h2) {
                                _f[i][j][k] -=
                                        2. * h / h2 * exp(-M_PI * M_PI * h2 /
                                                          (_alpha * _alpha)) *
                                        sin(2. * M_PI * Vec::dot(h, pos));
                            }
                        }
                    }
                }
            }
        }
    }
    string path = string(EWALDTABLE_LOCATION) + string("ewaldtable3d.dat");
    ofstream file(path, ios::out | ios::binary);
    for(unsigned int i = _size + 1; i--;) {
        for(unsigned int j = _size + 1; j--;) {
            for(unsigned int k = _size + 1; k--;) {
                file.write(reinterpret_cast<char*>(&_f[i][j][k][0]),
                           sizeof(_f[i][j][k][0]));
                file.write(reinterpret_cast<char*>(&_f[i][j][k][1]),
                           sizeof(_f[i][j][k][1]));
                file.write(reinterpret_cast<char*>(&_f[i][j][k][2]),
                           sizeof(_f[i][j][k][2]));
            }
        }
    }
#else
    Vec n;
    Vec h;
    for(unsigned int i = _size + 1; i--;) {
        for(unsigned int j = _size + 1; j--;) {
            if(!i && !j) {
                // all forces are initialized to zero vectors, no need to do it
                // here
                continue;
            }
            Vec pos(0.5 * i / (double)_size, 0.5 * j / (double)_size);
            _f[i][j] = pos / (pos.norm() * pos.norm2());
            for(n[0] = -4.; n[0] <= 4.; n[0]++) {
                for(n[1] = -4.; n[1] <= 4.; n[1]++) {
                    Vec dx = pos - n;
                    double r = dx.norm();
                    double r2 = dx.norm2();
                    _f[i][j] -= dx / (r * r2) *
                                (erfc(_alpha * r) +
                                 2. * _alpha * r / sqrt(M_PI) *
                                         exp(-_alpha * _alpha * r2));
                }
            }
            for(h[0] = -4.; h[0] <= 4.; h[0]++) {
                for(h[1] = -4.; h[1] <= 4.; h[1]++) {
                    double h2 = h.norm2();
                    if(h2) {
                        _f[i][j] -= 2. * h / h2 *
                                    exp(-M_PI * M_PI * h2 / (_alpha * _alpha)) *
                                    sin(2. * M_PI * Vec::dot(h, pos));
                    }
                }
            }
        }
    }
    string path = string(EWALDTABLE_LOCATION) + string("ewaldtable2d.dat");
    ofstream file(path, ios::out | ios::binary);
    for(unsigned int i = _size + 1; i--;) {
        for(unsigned int j = _size + 1; j--;) {
            file.write(reinterpret_cast<char*>(&_f[i][j][0]),
                       sizeof(_f[i][j][0]));
            file.write(reinterpret_cast<char*>(&_f[i][j][1]),
                       sizeof(_f[i][j][1]));
        }
    }
#endif
}

/**
  * @brief Get the Ewald correction at a given position
  *
  * The correction is calculated using linear interpolation on the internal
  * table.
  *
  * @param position A Vec specifying a position inside the periodic box
  * @return A Vec containing the force correction due to periodic copies of the
  * particles in the box
  */
Vec EwaldTable::get_correction(Vec& position) {
    Vec correction;
    Vec sign;
    for(unsigned int i = ndim_; i--;) {
        if(position[i] < 0) {
            position[i] = -position[i];
            sign[i] = 1.;
        } else {
            sign[i] = -1.;
        }
    }
    Vec u = position * _ewald_intp_fac;
    int integers[ndim_] = {0};
    for(unsigned int i = ndim_; i--;) {
        integers[i] = (int)u[i];
        if(integers[i] >= static_cast<int>(_size)) {
            integers[i] = _size - 1;
        }
        u[i] -= integers[i];
    }
#if ndim_ == 3
    double f[8] = {0.};
    f[0] = (1. - u[0]) * (1. - u[1]) * (1. - u[2]);
    f[1] = (1. - u[0]) * (1. - u[1]) * u[2];
    f[2] = (1. - u[0]) * u[1] * (1. - u[2]);
    f[3] = (1. - u[0]) * u[1] * u[2];
    f[4] = u[0] * (1. - u[1]) * (1. - u[2]);
    f[5] = u[0] * (1. - u[1]) * u[2];
    f[6] = u[0] * u[1] * (1. - u[2]);
    f[7] = u[0] * u[1] * u[2];
    correction = _f[integers[0]][integers[1]][integers[2]] * f[0] +
                 _f[integers[0]][integers[1]][integers[2] + 1] * f[1] +
                 _f[integers[0]][integers[1] + 1][integers[2]] * f[2] +
                 _f[integers[0]][integers[1] + 1][integers[2] + 1] * f[3] +
                 _f[integers[0] + 1][integers[1]][integers[2]] * f[4] +
                 _f[integers[0] + 1][integers[1]][integers[2] + 1] * f[5] +
                 _f[integers[0] + 1][integers[1] + 1][integers[2]] * f[6] +
                 _f[integers[0] + 1][integers[1] + 1][integers[2] + 1] * f[7];
#else
    double f[4] = {0.};
    f[0] = (1. - u[0]) * (1. - u[1]);
    f[1] = (1. - u[0]) * u[1];
    f[2] = u[0] * (1. - u[1]);
    f[3] = u[0] * u[1];
    correction = _f[integers[0]][integers[1]] * f[0] +
                 _f[integers[0]][integers[1] + 1] * f[1] +
                 _f[integers[0] + 1][integers[1]] * f[2] +
                 _f[integers[0] + 1][integers[1] + 1] * f[3];
#endif
    for(unsigned int i = ndim_; i--;) {
        correction[i] *= sign[i];
    }
    return correction;
}
