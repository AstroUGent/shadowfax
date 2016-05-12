/*******************************************************************************
 * This file is part of Shadowfax
 * Copyright (C) 2015 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *               2013 Peter Camps (peter.camps@ugent.be)
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
 * @file Vec.hpp
 *
 * @brief A 2D or 3D vector
 *
 * Code mostly copied from the Vec class in the radiative transfer code SKIRT,
 * that was kindly provided to us by Peter Camps. We added some small
 * modifications to be able to index the elements of the Vec.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 * @author Peter Camps (peter.camps@ugent.be)
 */
#ifndef VEC_HPP
#define VEC_HPP

#include <cmath>

using namespace std;

//////////////////////////////////////////////////////////////////////

#if ndim_ == 3
/**
  * \brief A 2D or 3D vector
  *
  * Vec is a low-level class for working with three-dimensional vectors: each
  * instance represents a vector with three cartesian components called \em x,
  * \em y, and \em z. The class defines operators for adding two vectors and for
  * multiplying a vector by a scalar. It offers functions to retrieve the vector
  * components, to get the norm or the squared norm of a vector, and for
  * calculating the dot product and cross product of two vectors. The Vec class
  * is fully implemented inline (in this header file). Most compilers optimize
  * away all overhead so that using this class is just as efficient as directly
  * writing the code in terms of the vector components.
  *
  * We acknowledge Peter Camps (Ghent University) for making this code available
  * to us.
  */
class Vec {
  protected:
    /*! \brief These data members represent the cartesian vector components.
     *  The anonymous union is used to permit indexing the components */
    union {
        /*! \brief Array storing the actual vector contents */
        double _c[3];
        struct {
            /*! \brief x-component of the vector */
            double _x;
            /*! \brief y-component of the vector */
            double _y;
            /*! \brief z-component of the vector */
            double _z;
        };
    };

  public:
    /**
     * @brief Default constructor
     *
     * All vector components are initialized to zero.
     */
    inline Vec() : _x(0), _y(0), _z(0) {}

    /**
     * @brief Constructor
     *
     * Initializes the vector components to the values provided as arguments.
     *
     * @param x x-component of the vector
     * @param y y-component of the vector
     * @param z z-component of the vector
     */
    inline Vec(double x, double y, double z) : _x(x), _y(y), _z(z) {}

    /**
     * @brief This function sets the vector components to the values provided
     * as arguments
     *
     * @param x x-component of the vector
     * @param y y-component of the vector
     * @param z z-component of the vector
     */
    inline void set(double x, double y, double z) {
        _x = x;
        _y = y;
        _z = z;
    }

    /**
     * @brief This function returns the \em x component of the vector
     *
     * @return x-component of the vector
     */
    inline double x() const {
        return _x;
    }

    /**
     * @brief This function returns the \em y component of the vector
     *
     * @return y-component of the vector
     */
    inline double y() const {
        return _y;
    }

    /**
     * @brief This function returns the \em z component of the vector
     *
     * @return z-component of the vector
     */
    inline double z() const {
        return _z;
    }

    /**
     * @brief This function returns the norm of the vector
     *
     * @return Norm of the vector
     */
    inline double norm() const {
        return sqrt(_x * _x + _y * _y + _z * _z);
    }

    /**
     * @brief This function returns the squared norm of the vector
     *
     * @return Norm of the vector squared
     */
    inline double norm2() const {
        return _x * _x + _y * _y + _z * _z;
    }

    /**
     * @brief This static function returns the dot product (inner product) of
     * two vectors
     *
     * @param a First vector
     * @param b Second vector
     * @return Dot product of the two given vectors
     */
    inline static double dot(Vec a, Vec b) {
        return a._x * b._x + a._y * b._y + a._z * b._z;
    }

    /**
     * @brief This static function returns the vector product (outer product) of
     * two vectors
     *
     * @param a First vector
     * @param b Second vector
     * @return Vec that is the vector product of the two given vectors
     */
    inline static Vec cross(Vec a, Vec b) {
        return Vec(a._y * b._z - a._z * b._y, a._z * b._x - a._x * b._z,
                   a._x * b._y - a._y * b._x);
    }

    /**
     * @brief This operator adds another vector to this one
     *
     * @param v Vector to add to this one
     * @return Reference to this vector
     */
    inline Vec& operator+=(Vec v) {
        _x += v._x;
        _y += v._y;
        _z += v._z;
        return *this;
    }

    /**
     * @brief This operator adds a scalar to this vector
     *
     * @param s Scalar to add to every component of this vector
     * @return Reference to this vector
     */
    inline Vec& operator+=(double s) {
        _x += s;
        _y += s;
        _z += s;
        return *this;
    }

    /**
     * @brief This operator subtracts another vector from this one
     *
     * @param v Vector to subtract from this one
     * @return Reference to this vector
     */
    inline Vec& operator-=(Vec v) {
        _x -= v._x;
        _y -= v._y;
        _z -= v._z;
        return *this;
    }

    /**
     * @brief This operator piecewise multiplies this vector with another vector
     *
     * @param v Vector to multiply piecewise with this one
     * @return Reference to this vector
     */
    inline Vec& operator*=(Vec v) {
        _x *= v._x;
        _y *= v._y;
        _z *= v._z;
        return *this;
    }

    /**
     * @brief This operator subtracts a scalar from this vector
     *
     * @param s Scalar to subtract from every component of this vector
     * @return Reference to this vector
     */
    inline Vec& operator-=(double s) {
        _x -= s;
        _y -= s;
        _z -= s;
        return *this;
    }

    /**
     * @brief This operator multiplies this vector by a scalar
     *
     * @param s Scalar to multiply with every component of this vector
     * @return Reference to this vector
     */
    inline Vec& operator*=(double s) {
        _x *= s;
        _y *= s;
        _z *= s;
        return *this;
    }

    /**
     * @brief This operator divides this vector by a scalar
     *
     * @param s Scalar to divide every component of this vector by
     * @return Reference to this vector
     */
    inline Vec& operator/=(double s) {
        _x /= s;
        _y /= s;
        _z /= s;
        return *this;
    }

    /**
     * @brief This operator returns the ith element of the vector, where i = 0
     * stands for _x and so on
     *
     * @param i Index of requested component
     * @return ith element of the vector
     */
    inline double operator[](int i) const {
        return _c[i];
    }

    /**
     * @brief This operator return a reference to the ith element, where i = 0
     * stands for _x and so on
     *
     * This reference can be used to change the contents of the ith element.
     *
     * @param i Index of the requested component
     * @return Reference to the ith element of the vector
     */
    inline double& operator[](int i) {
        return _c[i];
    }
};
#else
/**
  * \brief A 2D or 3D vector
  *
  * Vec is a low-level class for working with three-dimensional vectors: each
  * instance represents a vector with three cartesian components called \em x,
  * \em y, and \em z. The class defines operators for adding two vectors and for
  * multiplying a vector by a scalar. It offers functions to retrieve the vector
  * components, to get the norm or the squared norm of a vector, and for
  * calculating the dot product and cross product of two vectors. The Vec class
  * is fully implemented inline (in this header file). Most compilers optimize
  * away all overhead so that using this class is just as efficient as directly
  * writing the code in terms of the vector components.
  *
  * We acknowledge Peter Camps (Ghent University) for making this code available
  * to us.
  */
class Vec {
  protected:
    /*! \brief These data members represent the cartesian vector components.
     *  The anonymous union is used to permit indexing the components */
    union {
        /*! \brief Array storing the actual vector contents */
        double _c[2];
        struct {
            /*! \brief x-component of the vector */
            double _x;
            /*! \brief y-component of the vector */
            double _y;
        };
    };

  public:
    /**
     * @brief Default constructor
     *
     * All vector components are initialized to zero.
     */
    inline Vec() : _x(0), _y(0) {}

    /**
     * @brief Constructor
     *
     * Initializes the vector components to the values provided as arguments.
     *
     * @param x x-component of the vector
     * @param y y-component of the vector
     */
    inline Vec(double x, double y) : _x(x), _y(y) {}

    /**
     * @brief This function sets the vector components to the values provided
     * as arguments
     *
     * @param x x-component of the vector
     * @param y y-component of the vector
     */
    inline void set(double x, double y) {
        _x = x;
        _y = y;
    }

    /**
     * @brief This function returns the \em x component of the vector
     *
     * @return x-component of the vector
     */
    inline double x() const {
        return _x;
    }

    /**
     * @brief This function returns the \em y component of the vector
     *
     * @return y-component of the vector
     */
    inline double y() const {
        return _y;
    }

    /**
     * @brief This function returns the norm of the vector
     *
     * @return Norm of the vector
     */
    inline double norm() const {
        return sqrt(_x * _x + _y * _y);
    }

    /**
     * @brief This function returns the squared norm of the vector
     *
     * @return Norm of the vector squared
     */
    inline double norm2() const {
        return _x * _x + _y * _y;
    }

    /**
     * @brief This static function returns the dot product (inner product) of
     * two vectors
     *
     * @param a First vector
     * @param b Second vector
     * @return Dot product of the two given vectors
     */
    inline static double dot(Vec a, Vec b) {
        return a._x * b._x + a._y * b._y;
    }

    /**
     * @brief This operator adds another vector to this one
     *
     * @param v Vector to add to this one
     * @return Reference to this vector
     */
    inline Vec& operator+=(Vec v) {
        _x += v._x;
        _y += v._y;
        return *this;
    }

    /**
     * @brief This operator adds a scalar to this vector
     *
     * @param s Scalar to add to every component of this vector
     * @return Reference to this vector
     */
    inline Vec& operator+=(double s) {
        _x += s;
        _y += s;
        return *this;
    }

    /**
     * @brief This operator subtracts another vector from this one
     *
     * @param v Vector to subtract from this one
     * @return Reference to this vector
     */
    inline Vec& operator-=(Vec v) {
        _x -= v._x;
        _y -= v._y;
        return *this;
    }

    /**
     * @brief This operator piecewise multiplies this vector with another vector
     *
     * @param v Vector to multiply piecewise with this one
     * @return Reference to this vector
     */
    inline Vec& operator*=(Vec v) {
        _x *= v._x;
        _y *= v._y;
        return *this;
    }

    /**
     * @brief This operator subtracts a scalar from this vector
     *
     * @param s Scalar to subtract from every component of this vector
     * @return Reference to this vector
     */
    inline Vec& operator-=(double s) {
        _x -= s;
        _y -= s;
        return *this;
    }

    /**
     * @brief This operator multiplies this vector by a scalar
     *
     * @param s Scalar to multiply with every component of this vector
     * @return Reference to this vector
     */
    inline Vec& operator*=(double s) {
        _x *= s;
        _y *= s;
        return *this;
    }

    /**
     * @brief This operator divides this vector by a scalar
     *
     * @param s Scalar to divide every component of this vector by
     * @return Reference to this vector
     */
    inline Vec& operator/=(double s) {
        _x /= s;
        _y /= s;
        return *this;
    }

    /**
     * @brief This operator returns the ith element of the vector, where i = 0
     * stands for _x and so on
     *
     * @param i Index of requested component
     * @return ith element of the vector
     */
    inline double operator[](int i) const {
        return _c[i];
    }

    /**
     * @brief This operator return a reference to the ith element, where i = 0
     * stands for _x and so on
     *
     * This reference can be used to change the contents of the ith element.
     *
     * @param i Index of the requested component
     * @return Reference to the ith element of the vector
     */
    inline double& operator[](int i) {
        return _c[i];
    }
};
#endif
//////////////////////////////////////////////////////////////////////

/**
 * @brief This free operator adds two vectors
 *
 * @param a First vector
 * @param b Second vector
 * @return Vec that is the sum of the two given vectors
 */
inline Vec operator+(Vec a, Vec b) {
    return a += b;
}

/**
 * @brief This free operator subtracts two vectors
 *
 * @param a First vector
 * @param b Second vector
 * @return Vec that is the first vector minus the second one
 */
inline Vec operator-(Vec a, Vec b) {
    return a -= b;
}

/**
 * @brief This free operator piecewise multiplies two vectors
 *
 * @param a First vector
 * @param b Second vector
 * @return Vec that is the piecewise multiplication of the two given vectors
 */
inline Vec operator*(Vec a, Vec b) {
    return a *= b;
}

/**
 * @brief This free operator multiplies a vector by a scalar
 *
 * @param a Vector
 * @param s Scalar
 * @return Vec that is the result of the multiplication
 */
inline Vec operator*(Vec a, double s) {
    return a *= s;
}

/**
 * @brief This free operator multiplies a scalar by a vector
 *
 * @param s Scalar
 * @param b Vector
 * @return Vec that is the result of the multiplication
 */
inline Vec operator*(double s, Vec b) {
    return b *= s;
}

/**
 * @brief This free operator adds a scalar to a vector
 *
 * @param a Vector
 * @param s Scalar
 * @return Vec that is the result of the addition
 */
inline Vec operator+(Vec a, double s) {
    return a += s;
}

/**
 * @brief This free operator adds a vector to a scalar
 *
 * @param s Scalar
 * @param a Vector
 * @return Vec that is the result of the addition
 */
inline Vec operator+(double s, Vec a) {
    return a += s;
}

/**
 * @brief This free operator subtracts a scalar from a vector
 *
 * @param a Vector
 * @param s Scalar
 * @return Vec that is the given vector minus the given scalar
 */
inline Vec operator-(Vec a, double s) {
    return a -= s;
}

/**
 * @brief This free operator divides a vector by a scalar
 *
 * @param a Vector
 * @param s Scalar
 * @return Vec that is the result of the division
 */
inline Vec operator/(Vec a, double s) {
    return a /= s;
}

//////////////////////////////////////////////////////////////////////

#endif  // VEC_HPP
