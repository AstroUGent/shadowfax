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
 * @file GSL.hpp
 *
 * @brief C++ version of basic GSL functions: header
 *
 * We use this approach, since (a) GSL is not supported by the default CMake
 * version we use, and (b) GSL is an ancient C library.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef GSL_HPP
#define GSL_HPP

#include <ostream>

class RestartFile;

namespace GSL {
// quadrature
double rescale_error(double err, const double result_abs,
                     const double result_asc);
double qk(double (*f)(double, void*), double x0, double x1, void* params,
          double& abserr, double& resabs, double& resasc);
double qag(double (*f)(double, void*), double x0, double x1, double tolerance,
           void* params);

// root finding
double brent(double (*f)(double, void*), double xlow, double xupp,
             double tolerance, void* params);

// interpolation

unsigned int bsearch(const double* xvals, double x, unsigned int ilow,
                     unsigned int iupp);

void solve_tridiag(const double* diag, const double* offdiag, const double* b,
                   double* x, unsigned int N);

/**
 * @brief Supported GSL interpolator types
 */
enum GSLInterpolatorType {
    /*! @brief Linear interpolation */
    TypeGSLLinearInterpolator = 0,
    /*! @brief Cubic spline interpolation */
    TypeGSLCubicSplineInterpolator
};

/**
 * @brief General interface for GSL interpolators
 */
class GSLInterpolator {
  protected:
    /*! @brief Copy of the data array underlying the interpolation: x-values */
    double* _x;
    /*! @brief Copy of the data array underlying the interpolation: y-values */
    double* _y;

    /*! @brief Size of the underlying data array */
    unsigned int _size;

  public:
    GSLInterpolator(double* x, double* y, unsigned int size);
    virtual ~GSLInterpolator();

    static GSLInterpolator* create(GSLInterpolatorType type, double* x,
                                   double* y, unsigned int size);
    static GSLInterpolator* create(RestartFile& rfile);

    virtual void dump(RestartFile& rfile);
    GSLInterpolator(RestartFile& rfile);

    virtual void dump(std::ostream& stream);
    virtual bool equal(GSLInterpolator* interpolator);

    /**
     * @brief Evaluate the interpolating function at the given position
     *
     * @param x Position to evaluate at
     * @return Value of the interpolating function
     */
    virtual double eval(double x) = 0;
};

/**
 * @brief Linear GSLInterpolator
 */
class GSLLinearInterpolator : public GSLInterpolator {
  private:
    /*! @brief Index of the last bin that was encountered */
    unsigned int _cache;

  public:
    GSLLinearInterpolator(double* x, double* y, unsigned int size);
    virtual ~GSLLinearInterpolator();

    virtual double eval(double x);

    virtual void dump(RestartFile& rfile);
    GSLLinearInterpolator(RestartFile& rfile);

    virtual void dump(std::ostream& stream);
    virtual bool equal(GSLInterpolator* interpolator);
};

/**
 * @brief Cubic spline GSLInterpolator
 */
class GSLCubicSplineInterpolator : public GSLInterpolator {
  private:
    /*! @brief Index of the last bin that was encountered */
    unsigned int _cache;

    /*! @brief Extra coefficients for spline */
    double* _c;
    /*! @brief Auxiliary variable for coefficient calculation */
    double* _g;
    /*! @brief Diagonal elements for coefficient calculation */
    double* _diag;
    /*! @brief Off-diagonal elements for coefficient calculation */
    double* _offdiag;

  public:
    GSLCubicSplineInterpolator(double* x, double* y, unsigned int size);
    virtual ~GSLCubicSplineInterpolator();

    virtual double eval(double x);

    virtual void dump(RestartFile& rfile);
    GSLCubicSplineInterpolator(RestartFile& rfile);

    virtual void dump(std::ostream& stream);
    virtual bool equal(GSLInterpolator* interpolator);
};
}

#endif  // GSL_HPP
