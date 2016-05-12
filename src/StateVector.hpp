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
 * @file StateVector.hpp
 *
 * @brief Hydrodynamical state vector
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef STATEVECTOR_HPP
#define STATEVECTOR_HPP

#define GAMMA 1.66667

#include "Vec.hpp"
// for std::max and std::min
#include <algorithm>

#if ndim_ == 3
/**
 * @brief 2D or 3D hydrodynamical state vector
 *
 * The StateVector can either contain primitive variables (density, velocity,
 * pressure) or conserved quantities (mass, momentum, total energy). It also
 * contains an extra field for Passively Advected Quantities (PAQs).
 * We store them in a single container because this is how the hydrodynamical
 * equations are expressed in literature, e.g. Toro (2009).
 */
class StateVector {
  protected:
    /**
     * \brief Contents of the state vector
     *
     * We use a union with an array to be able to index the elements.
     *
     * We use separate structs for the primitive and conserved variables to be
     * able to use logical names for the quantities everywhere.
     *
     * Notice that we only have 5 variables (due to the union): it is not
     * possible to store e.g. a density and a mass at the same time, because
     * they use the same memory!
     */
    union {
        /*! \brief Auxiliary array to help indexing the contents of this
         *  vector */
        double _c[5];
        /*! \brief Primitive variables */
        struct {
            /*! \brief Density */
            double _rho;
            /*! \brief x-component of the velocity */
            double _vx;
            /*! \brief y-component of the velocity */
            double _vy;
            /*! \brief z-component of the velocity */
            double _vz;
            /*! \brief Pressure */
            double _p;
        };
        /*! \brief Conserved variables */
        struct {
            /*! \brief Mass */
            double _m;
            /*! \brief x-component of the momentum */
            double _px;
            /*! \brief y-component of the momentum */
            double _py;
            /*! \brief z-component of the momentum */
            double _pz;
            /*! \brief Energy */
            double _e;
        };
    };

    /*! \brief Passively Advected Quantity that is advected with the flow */
    double _paq;

  public:
    /**
     * @brief Empty constructor
     */
    inline StateVector() : _rho(0), _vx(0), _vy(0), _vz(0), _p(0), _paq(0) {}

    /**
     * @brief Single value constructor
     *
     * Initialize a vector where all elements have the same given value.
     *
     * @param singleVal Value for all 5 components of the vector
     */
    inline StateVector(double singleVal)
            : _rho(singleVal), _vx(singleVal), _vy(singleVal), _vz(singleVal),
              _p(singleVal), _paq(singleVal) {}

    /**
     * @brief Constructor
     *
     * @param rho Density/mass
     * @param vx x-component of the velocity/momentum
     * @param vy y-component of the velocity/momentum
     * @param vz z-component of the velocity/momentum
     * @param p Pressure/energy
     */
    inline StateVector(double rho, double vx, double vy, double vz, double p)
            : _rho(rho), _vx(vx), _vy(vy), _vz(vz), _p(p), _paq(0) {}

    /**
     * @brief Set the elements of the state vector to the given values
     *
     * @param rho Density/mass
     * @param vx x-component of the velocity/momentum
     * @param vy y-component of the velocity/momentum
     * @param vz z-component of the velocity/momentum
     * @param p Pressure/energy
     */
    inline void set(double rho, double vx, double vy, double vz, double p) {
        _rho = rho;
        _vx = vx;
        _vy = vy;
        _vz = vz;
        _p = p;
    }

    /**
     * @brief Set the density
     *
     * @param rho Density
     */
    inline void set_rho(double rho) {
        _rho = rho;
    }

    /**
     * @brief Set the x-component of the velocity
     *
     * @param u x-component of the velocity
     */
    inline void set_vx(double u) {
        _vx = u;
    }

    /**
     * @brief Set the y-component of the velocity
     *
     * @param u y-component of the velocity
     */
    inline void set_vy(double u) {
        _vy = u;
    }

    /**
     * @brief Set the z-component of the velocity
     *
     * @param u z-component of the velocity
     */
    inline void set_vz(double u) {
        _vz = u;
    }

    /**
     * @brief Set the pressure
     *
     * @param p Pressure
     */
    inline void set_p(double p) {
        _p = p;
    }

    /**
     * @brief Set the mass
     *
     * @param m Mass
     */
    inline void set_m(double m) {
        _m = m;
    }

    /**
     * @brief Set the x-component of the momentum
     *
     * @param px x-component of the momentum
     */
    inline void set_px(double px) {
        _px = px;
    }

    /**
     * @brief Set the y-component of the momentum
     *
     * @param py y-component of the momentum
     */
    inline void set_py(double py) {
        _py = py;
    }

    /**
     * @brief Set the z-component of the momentum
     *
     * @param pz z-component of the momentum
     */
    inline void set_pz(double pz) {
        _pz = pz;
    }

    /**
     * @brief Set the energy
     *
     * @param e Energy
     */
    inline void set_e(double e) {
        _e = e;
    }

    /**
     * @brief Set the Passively Advected Quantity (PAQ)
     *
     * @param paq PAQ
     */
    inline void set_paq(double paq) {
        _paq = paq;
    }

    /**
     * @brief Get the density
     *
     * @return Density
     */
    inline double rho() const {
        return _rho;
    }

    /**
     * @brief Get the x-component of the velocity
     *
     * @return x-component of the velocity
     */
    inline double vx() const {
        return _vx;
    }

    /**
     * @brief Get the y-component of the velocity
     *
     * @return y-component of the velocity
     */
    inline double vy() const {
        return _vy;
    }

    /**
     * @brief Get the z-component of the velocity
     *
     * @return z-component of the velocity
     */
    inline double vz() const {
        return _vz;
    }

    /**
     * @brief Get the pressure
     *
     * @return Pressure
     */
    inline double p() const {
        return _p;
    }

    /**
     * @brief Get the mass
     *
     * @return Mass
     */
    inline double m() const {
        return _m;
    }

    /**
     * @brief Get the x-component of the momentum
     *
     * @return x-component of the momentum
     */
    inline double px() const {
        return _px;
    }

    /**
     * @brief Get the y-component of the momentum
     *
     * @return y-component of the momentum
     */
    inline double py() const {
        return _py;
    }

    /**
     * @brief Get the z-component of the momentum
     *
     * @return z-component of the momentum
     */
    inline double pz() const {
        return _pz;
    }

    /**
     * @brief Get the energy
     *
     * @return Energy
     */
    inline double e() const {
        return _e;
    }

    /**
     * @brief Get the Passively Advected Quantity (PAQ)
     *
     * @return PAQ
     */
    inline double paq() const {
        return _paq;
    }

    /**
     * @brief Add another state vector to this one
     *
     * @param v StateVector to add
     * @return Reference to this state vector
     */
    inline StateVector& operator+=(StateVector v) {
        _m += v._m;
        _px += v._px;
        _py += v._py;
        _pz += v._pz;
        _e += v._e;
        _paq += v._paq;
        return *this;
    }

    /**
     * @brief Subtract another state vector from this one
     *
     * @param v StateVector to subtract
     * @return Reference to this state vector
     */
    inline StateVector& operator-=(StateVector v) {
        _m -= v._m;
        _px -= v._px;
        _py -= v._py;
        _pz -= v._pz;
        _e -= v._e;
        _paq -= v._paq;
        return *this;
    }

    /**
     * @brief Add a vector to the velocity components of this state vector
     *
     * @param v Vec to add
     * @return Reference to this state vector
     */
    inline StateVector& operator+=(Vec v) {
        _vx += v.x();
        _vy += v.y();
        _vz += v.z();
        return *this;
    }

    /**
     * @brief Subtract a vector from the velocity components of this state
     * vector
     *
     * @param v Vec to subtract
     * @return Reference to this state vector
     */
    inline StateVector& operator-=(Vec v) {
        _vx -= v.x();
        _vy -= v.y();
        _vz -= v.z();
        return *this;
    }

    /**
     * @brief Multiply all elements of the state vector with a scalar
     *
     * @param s Scalar to multiply with
     * @return Reference to this state vector
     */
    inline StateVector& operator*=(double s) {
        _m *= s;
        _px *= s;
        _py *= s;
        _pz *= s;
        _e *= s;
        _paq *= s;
        return *this;
    }

    /**
     * @brief Divide all elements of this state vector by a scalar
     *
     * @param s Scalar to divide by
     * @return Reference to this state vector
     */
    inline StateVector& operator/=(double s) {
        _m /= s;
        _px /= s;
        _py /= s;
        _pz /= s;
        _e /= s;
        _paq /= s;
        return *this;
    }

    /**
     * @brief Multiply the elements of this state vector one-by-one with the
     * elements of another state vector
     *
     * @param v StateVector to multiply with
     * @return Reference to this state vector
     */
    inline StateVector& operator*=(StateVector v) {
        _m *= v._m;
        _px *= v._px;
        _py *= v._py;
        _pz *= v._pz;
        _e *= v._e;
        _paq *= v._paq;
        return *this;
    }

    /**
     * @brief Divide the elements of this state vector one-by-one by the
     * elements of another state vector
     *
     * If an element of the second state vector is zero, we set the
     * corresponding component of the state vector to one.
     *
     * @param v StateVector to divide by
     * @return Reference to this state vector
     */
    inline StateVector& operator/=(StateVector v) {
        if(v._m) {
            _m /= v._m;
        } else {
            _m = 1;
        }
        if(v._px) {
            _px /= v._px;
        } else {
            _px = 1;
        }
        if(v._py) {
            _py /= v._py;
        } else {
            _py = 1;
        }
        if(v._pz) {
            _pz /= v._pz;
        } else {
            _pz = 1;
        }
        if(v._e) {
            _e /= v._e;
        } else {
            _e = 1;
        }
        if(v._paq) {
            _paq /= v._paq;
        } else {
            _paq = 1;
        }
        return *this;
    }

    /**
     * @brief Set the components of this state vector to the maximum of the
     * components of this state vector and another state vector
     *
     * @param v Second StateVector
     * @return Reference to this state vector
     */
    inline StateVector& max(StateVector v) {
        _m = std::max(_m, v._m);
        _px = std::max(_px, v._px);
        _py = std::max(_py, v._py);
        _pz = std::max(_pz, v._pz);
        _e = std::max(_e, v._e);
        _paq = std::max(_paq, v._paq);
        return *this;
    }

    /**
     * @brief Set the components of this state vector to the minimum of the
     * components of this state vector and another state vector
     *
     * @param v Second StateVector
     * @return Reference to this state vector
     */
    inline StateVector& min(StateVector v) {
        _m = std::min(_m, v._m);
        _px = std::min(_px, v._px);
        _py = std::min(_py, v._py);
        _pz = std::min(_pz, v._pz);
        _e = std::min(_e, v._e);
        _paq = std::min(_paq, v._paq);
        return *this;
    }

    /**
     * @brief Clear the state vector by setting all of its elements to zero
     */
    inline void reset() {
        _rho = 0;
        _vx = 0;
        _vy = 0;
        _vz = 0;
        _p = 0;
        _paq = 0;
    }

    /**
     * @brief Access individual elements by index
     *
     * @param i Index of the element in the internal array
     * @return Element at the given index
     */
    inline double operator[](int i) const {
        return _c[i];
    }

    /**
     * @brief Access individual elements by index
     *
     * @param i Index of the element in the internal array
     * @return Element at the given index
     */
    inline double& operator[](int i) {
        return _c[i];
    }
};
#else
/**
 * @brief 2D or 3D hydrodynamical state vector
 *
 * The StateVector can either contain primitive variables (density, velocity,
 * pressure) or conserved quantities (mass, momentum, total energy). It also
 * contains an extra field for Passively Advected Quantities (PAQs).
 * We store them in a single container because this is how the hydrodynamical
 * equations are expressed in literature, e.g. Toro (2009).
 */
class StateVector {
  protected:
    /**
     * \brief Contents of the state vector
     *
     * We use a union with an array to be able to index the elements.
     *
     * We use separate structs for the primitive and conserved variables to be
     * able to use logical names for the quantities everywhere.
     *
     * Notice that we only have 4 variables (due to the union): it is not
     * possible to store e.g. a density and a mass at the same time, because
     * they use the same memory!
     */
    union {
        /*! \brief Auxiliary array to help indexing the contents of this
         *  vector */
        double _c[4];
        /*! \brief Primitive variables */
        struct {
            /*! \brief Density */
            double _rho;
            /*! \brief x-component of the velocity */
            double _vx;
            /*! \brief y-component of the velocity */
            double _vy;
            /*! \brief Pressure */
            double _p;
        };
        /*! \brief Conserved variables */
        struct {
            /*! \brief Mass */
            double _m;
            /*! \brief x-component of the momentum */
            double _px;
            /*! \brief y-component of the momentum */
            double _py;
            /*! \brief Energy */
            double _e;
        };
    };

    /*! \brief Passively Advected Quantity that is advected with the flow */
    double _paq;

  public:
    /**
     * @brief Empty constructor
     */
    inline StateVector() : _rho(0), _vx(0), _vy(0), _p(0), _paq(0) {}

    /**
     * @brief Single value constructor
     *
     * Initialize a vector where all elements have the same given value.
     *
     * @param singleVal Value for all 4 components of the vector
     */
    inline StateVector(double singleVal)
            : _rho(singleVal), _vx(singleVal), _vy(singleVal), _p(singleVal),
              _paq(singleVal) {}

    /**
     * @brief Constructor
     *
     * @param rho Density/mass
     * @param vx x-component of the velocity/momentum
     * @param vy y-component of the velocity/momentum
     * @param p Pressure/energy
     */
    inline StateVector(double rho, double vx, double vy, double p)
            : _rho(rho), _vx(vx), _vy(vy), _p(p), _paq(0) {}

    /**
     * @brief Set the elements of the state vector to the given values
     *
     * @param rho Density/mass
     * @param vx x-component of the velocity/momentum
     * @param vy y-component of the velocity/momentum
     * @param p Pressure/energy
     */
    inline void set(double rho, double vx, double vy, double p) {
        _rho = rho;
        _vx = vx;
        _vy = vy;
        _p = p;
    }

    /**
     * @brief Set the density
     *
     * @param rho Density
     */
    inline void set_rho(double rho) {
        _rho = rho;
    }

    /**
     * @brief Set the x-component of the velocity
     *
     * @param u x-component of the velocity
     */
    inline void set_vx(double u) {
        _vx = u;
    }

    /**
     * @brief Set the y-component of the velocity
     *
     * @param u y-component of the velocity
     */
    inline void set_vy(double u) {
        _vy = u;
    }

    /**
     * @brief Set the pressure
     *
     * @param p Pressure
     */
    inline void set_p(double p) {
        _p = p;
    }

    /**
     * @brief Set the mass
     *
     * @param m Mass
     */
    inline void set_m(double m) {
        _m = m;
    }

    /**
     * @brief Set the x-component of the momentum
     *
     * @param px x-component of the momentum
     */
    inline void set_px(double px) {
        _px = px;
    }

    /**
     * @brief Set the y-component of the momentum
     *
     * @param py y-component of the momentum
     */
    inline void set_py(double py) {
        _py = py;
    }

    /**
     * @brief Set the energy
     *
     * @param e Energy
     */
    inline void set_e(double e) {
        _e = e;
    }

    /**
     * @brief Set the Passively Advected Quantity (PAQ)
     *
     * @param paq PAQ
     */
    inline void set_paq(double paq) {
        _paq = paq;
    }

    /**
     * @brief Get the density
     *
     * @return Density
     */
    inline double rho() const {
        return _rho;
    }

    /**
     * @brief Get the x-component of the velocity
     *
     * @return x-component of the velocity
     */
    inline double vx() const {
        return _vx;
    }

    /**
     * @brief Get the y-component of the velocity
     *
     * @return y-component of the velocity
     */
    inline double vy() const {
        return _vy;
    }

    /**
     * @brief Get the pressure
     *
     * @return Pressure
     */
    inline double p() const {
        return _p;
    }

    /**
     * @brief Get the mass
     *
     * @return Mass
     */
    inline double m() const {
        return _m;
    }

    /**
     * @brief Get the x-component of the momentum
     *
     * @return x-component of the momentum
     */
    inline double px() const {
        return _px;
    }

    /**
     * @brief Get the y-component of the momentum
     *
     * @return y-component of the momentum
     */
    inline double py() const {
        return _py;
    }

    /**
     * @brief Get the energy
     *
     * @return Energy
     */
    inline double e() const {
        return _e;
    }

    /**
     * @brief Get the Passively Advected Quantity (PAQ)
     *
     * @return PAQ
     */
    inline double paq() const {
        return _paq;
    }

    /**
     * @brief Add another state vector to this one
     *
     * @param v StateVector to add
     * @return Reference to this state vector
     */
    inline StateVector& operator+=(StateVector v) {
        _m += v._m;
        _px += v._px;
        _py += v._py;
        _e += v._e;
        _paq += v._paq;
        return *this;
    }

    /**
     * @brief Subtract another state vector from this one
     *
     * @param v StateVector to subtract
     * @return Reference to this state vector
     */
    inline StateVector& operator-=(StateVector v) {
        _m -= v._m;
        _px -= v._px;
        _py -= v._py;
        _e -= v._e;
        _paq -= v._paq;
        return *this;
    }

    /**
     * @brief Add a vector to the velocity components of this state vector
     *
     * @param v Vec to add
     * @return Reference to this state vector
     */
    inline StateVector& operator+=(Vec v) {
        _vx += v.x();
        _vy += v.y();
        return *this;
    }

    /**
     * @brief Subtract a vector from the velocity components of this state
     * vector
     *
     * @param v Vec to subtract
     * @return Reference to this state vector
     */
    inline StateVector& operator-=(Vec v) {
        _vx -= v.x();
        _vy -= v.y();
        return *this;
    }

    /**
     * @brief Multiply all elements of the state vector with a scalar
     *
     * @param s Scalar to multiply with
     * @return Reference to this state vector
     */
    inline StateVector& operator*=(double s) {
        _m *= s;
        _px *= s;
        _py *= s;
        _e *= s;
        _paq *= s;
        return *this;
    }

    /**
     * @brief Divide all elements of this state vector by a scalar
     *
     * @param s Scalar to divide by
     * @return Reference to this state vector
     */
    inline StateVector& operator/=(double s) {
        _m /= s;
        _px /= s;
        _py /= s;
        _e /= s;
        _paq /= s;
        return *this;
    }

    /**
     * @brief Multiply the elements of this state vector one-by-one with the
     * elements of another state vector
     *
     * @param v StateVector to multiply with
     * @return Reference to this state vector
     */
    inline StateVector& operator*=(StateVector v) {
        _m *= v._m;
        _px *= v._px;
        _py *= v._py;
        _e *= v._e;
        _paq *= v._paq;
        return *this;
    }

    /**
     * @brief Divide the elements of this state vector one-by-one by the
     * elements of another state vector
     *
     * If an element of the second state vector is zero, we set the
     * corresponding component of the state vector to one.
     *
     * @param v StateVector to divide by
     * @return Reference to this state vector
     */
    inline StateVector& operator/=(StateVector v) {
        if(v._m) {
            _m /= v._m;
        } else {
            _m = 1;
        }
        if(v._px) {
            _px /= v._px;
        } else {
            _px = 1;
        }
        if(v._py) {
            _py /= v._py;
        } else {
            _py = 1;
        }
        if(v._e) {
            _e /= v._e;
        } else {
            _e = 1;
        }
        if(v._paq) {
            _paq /= v._paq;
        } else {
            _paq = 1;
        }
        return *this;
    }

    /**
     * @brief Set the components of this state vector to the maximum of the
     * components of this state vector and another state vector
     *
     * @param v Second StateVector
     * @return Reference to this state vector
     */
    inline StateVector& max(StateVector v) {
        _m = std::max(_m, v._m);
        _px = std::max(_px, v._px);
        _py = std::max(_py, v._py);
        _e = std::max(_e, v._e);
        _paq = std::max(_paq, v._paq);
        return *this;
    }

    /**
     * @brief Set the components of this state vector to the minimum of the
     * components of this state vector and another state vector
     *
     * @param v Second StateVector
     * @return Reference to this state vector
     */
    inline StateVector& min(StateVector v) {
        _m = std::min(_m, v._m);
        _px = std::min(_px, v._px);
        _py = std::min(_py, v._py);
        _e = std::min(_e, v._e);
        _paq = std::min(_paq, v._paq);
        return *this;
    }

    /**
     * @brief Clear the state vector by setting all of its elements to zero
     */
    inline void reset() {
        _rho = 0;
        _vx = 0;
        _vy = 0;
        _p = 0;
        _paq = 0;
    }

    /**
     * @brief Access individual elements by index
     *
     * @param i Index of the element in the internal array
     * @return Element at the given index
     */
    inline double operator[](int i) const {
        return _c[i];
    }

    /**
     * @brief Access individual elements by index
     *
     * @param i Index of the element in the internal array
     * @return Element at the given index
     */
    inline double& operator[](int i) {
        return _c[i];
    }
};
#endif

/**
 * @brief This operator adds two state vectors
 *
 * @param a First StateVector
 * @param b Second StateVector
 * @return Resulting StateVector
 */
inline StateVector operator+(StateVector a, StateVector b) {
    return a += b;
}

/**
 * @brief This operator subtracts a state vector from another state vector
 *
 * @param a First StateVector
 * @param b Second StateVector that is subtracted from the first
 * @return Resulting StateVector
 */
inline StateVector operator-(StateVector a, StateVector b) {
    return a -= b;
}

/**
 * @brief This operator multiplies two state vectors element-wise
 *
 * @param a First StateVector
 * @param b Second StateVector
 * @return Resulting StateVector
 */
inline StateVector operator*(StateVector a, StateVector b) {
    return a *= b;
}

/**
 * @brief This operator multiplies a state vector with a scalar
 *
 * @param a StateVector
 * @param s Scalar to multiply with
 * @return Resulting StateVector
 */
inline StateVector operator*(StateVector a, double s) {
    return a *= s;
}

/**
 * @brief This operator multiplies a scalar with a state vector
 *
 * @param s Scalar
 * @param b StateVector to multiply with
 * @return Resulting StateVector
 */
inline StateVector operator*(double s, StateVector b) {
    return b *= s;
}

/**
 * @brief This operator divides a state vector by a scalar
 *
 * @param a StateVector
 * @param s Scalar to divide by
 * @return Resulting StateVector
 */
inline StateVector operator/(StateVector a, double s) {
    return a /= s;
}

/**
 * @brief Get a state vector with components that are the maximum of the
 * components of the given state vectors
 *
 * @param a First StateVector
 * @param b Second StateVector
 * @return Resulting StateVector
 */
inline StateVector max(StateVector a, StateVector b) {
    return a.max(b);
}

/**
 * @brief Get a state vector with components that are the minimum of the
 * components of the given state vectors
 *
 * @param a First StateVector
 * @param b Second StateVector
 * @return Resulting StateVector
 */
inline StateVector min(StateVector a, StateVector b) {
    return a.min(b);
}

/**
 * @brief Operator that divides a state vector element-wise by another state
 * vector
 *
 * @param a First StateVector
 * @param b Second StateVector to divide by
 * @return Resulting StateVector
 */
inline StateVector operator/(StateVector a, StateVector b) {
    return a /= b;
}

#endif  // STATEVECTOR_HPP
