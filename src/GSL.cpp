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
 * @file GSL.cpp
 *
 * @brief C++ version of basic GSL functions: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "GSL.hpp"
#include "Error.hpp"  // for my_exit
#include <algorithm>  // for min
#include <cmath>      // for fabs, pow
#include <iostream>   // for operator<<, basic_ostream, etc
#include <stddef.h>   // for NULL

/*
 * The code in this file is mostly copied from the GSL library, which was made
 * public under the GNU General Public License. We have used code from
 *  - integration/gsl_integration.h
 *  - integration/qag.c
 *  - integration/qak41.c
 *  - integration/qak.c
 *  - integration/workspace.c
 *  - integration/util.c
 *  - integration/set_initial.c
 *  - integration/initialise.c
 *  - integration/qpsrt.c
 *  - integration/err.c
 *
 *  - roots/brent.c
 *  - roots/convergence.c
 *  - roots/fsolver.c
 *
 *  - interpolation/linear.c
 *  - interpolation/cspline.c
 *  - interpolation/accel.c
 *  - interpolation/gsl_interp.h
 *  - linalg/tridiag.c
 *
 * Some values from gsl_machine.h were also used.
 *
 * The original license info for this code is below.
 * The original code can be found on
 * http://git.savannah.gnu.org/cgit/gsl.git/tree/
 */

/*
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2002, 2004 2007 Reid Priedhorsky,
 * Brian Gough, Gerard Jungman, David Necas
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
 * USA.
 */

/* Gauss quadrature weights and kronrod quadrature abscissae and
   weights as evaluated with 80 decimal digit arithmetic by
   L. W. Fullerton, Bell Labs, Nov. 1981. */

/*! @brief Abscissae of the 41-point kronrod rule */
static const double xgk[21] = {0.998859031588277663838315576545863,
                               0.993128599185094924786122388471320,
                               0.981507877450250259193342994720217,
                               0.963971927277913791267666131197277,
                               0.940822633831754753519982722212443,
                               0.912234428251325905867752441203298,
                               0.878276811252281976077442995113078,
                               0.839116971822218823394529061701521,
                               0.795041428837551198350638833272788,
                               0.746331906460150792614305070355642,
                               0.693237656334751384805490711845932,
                               0.636053680726515025452836696226286,
                               0.575140446819710315342946036586425,
                               0.510867001950827098004364050955251,
                               0.443593175238725103199992213492640,
                               0.373706088715419560672548177024927,
                               0.301627868114913004320555356858592,
                               0.227785851141645078080496195368575,
                               0.152605465240922675505220241022678,
                               0.076526521133497333754640409398838,
                               0.000000000000000000000000000000000};

/* xgk[1], xgk[3], ... abscissae of the 20-point gauss rule.
   xgk[0], xgk[2], ... abscissae to optimally extend the 20-point gauss rule */

/*! @brief Weights of the 20-point gauss rule */
static const double wg[10] = {0.017614007139152118311861962351853,
                              0.040601429800386941331039952274932,
                              0.062672048334109063569506535187042,
                              0.083276741576704748724758143222046,
                              0.101930119817240435036750135480350,
                              0.118194531961518417312377377711382,
                              0.131688638449176626898494499748163,
                              0.142096109318382051329298325067165,
                              0.149172986472603746787828737001969,
                              0.152753387130725850698084331955098};

/*! @brief Weights of the 41-point kronrod rule */
static const double wgk[21] = {0.003073583718520531501218293246031,
                               0.008600269855642942198661787950102,
                               0.014626169256971252983787960308868,
                               0.020388373461266523598010231432755,
                               0.025882133604951158834505067096153,
                               0.031287306777032798958543119323801,
                               0.036600169758200798030557240707211,
                               0.041668873327973686263788305936895,
                               0.046434821867497674720231880926108,
                               0.050944573923728691932707670050345,
                               0.055195105348285994744832372419777,
                               0.059111400880639572374967220648594,
                               0.062653237554781168025870122174255,
                               0.065834597133618422111563556969398,
                               0.068648672928521619345623411885368,
                               0.071054423553444068305790361723210,
                               0.073030690332786667495189417658913,
                               0.074582875400499188986581418362488,
                               0.075704497684556674659542775376617,
                               0.076377867672080736705502835038061,
                               0.076600711917999656445049901530102};

/**
 * @brief GSL workspace abstraction
 */
class GSL_workspace {
  private:
    /*! @brief Maximal number of subintervals */
    unsigned int _limit;

    /*! @brief Current number of subintervals */
    unsigned int _size;

    /*! @brief Index of the interval with the largest absolute error */
    unsigned int _nrmax;

    /*! @brief Index of the current subinterval */
    unsigned int _i;

    /*! @brief Maximum subdivision level of all intervals */
    unsigned int _maximum_level;

    /*! @brief Starting points of the intervals */
    double* _alist;

    /*! @brief End points of the intervals */
    double* _blist;

    /*! @brief Results for the different intervals */
    double* _rlist;

    /*! @brief Absolute errors for the different intervals */
    double* _elist;

    /*! @brief Ordered list of intervals according to the size of the current
     * absolute error on their results */
    unsigned int* _order;

    /*! @brief Subdivision level of the intervals as a power of 2 of the initial
     * interval */
    unsigned int* _level;

  public:
    /**
     * @brief Constructor
     *
     * @param n Size of the internal arrays, maximal number of subintervals
     */
    GSL_workspace(unsigned int n) {
        _limit = n;
        _size = 0;
        _maximum_level = 0;

        _alist = new double[n];
        _blist = new double[n];
        _rlist = new double[n];
        _elist = new double[n];
        _order = new unsigned int[n];
        _level = new unsigned int[n];
    }

    /**
     * @brief Destructor
     */
    ~GSL_workspace() {
        delete[] _alist;
        delete[] _blist;
        delete[] _rlist;
        delete[] _elist;
        delete[] _order;
        delete[] _level;
    }

    /**
     * @brief Initialise the first interval
     *
     * @param a Lower bound of the interval
     * @param b Upper bound of the interval
     */
    void initialise(double a, double b) {
        _alist[0] = a;
        _blist[0] = b;
        _size = 0;
        _nrmax = 0;
        _i = 0;
        _rlist[0] = 0.;
        _elist[0] = 0.;
        _order[0] = 0;
        _level[0] = 0;
        _maximum_level = 0;
    }

    /**
     * @brief Set the result for the first interval
     *
     * @param result Result
     * @param error Absolute error on the result
     */
    void set_initial_result(double result, double error) {
        _size = 1;
        _rlist[0] = result;
        _elist[0] = error;
    }

    /**
     * @brief Retrieve the latest values for the current interval
     *
     * @param a Lower bound of the interval
     * @param b Upper bound of the interval
     * @param r Latest result of the interval
     * @param e Absolute error on the latest result for the interval
     */
    void retrieve(double& a, double& b, double& r, double& e) {
        a = _alist[_i];
        b = _blist[_i];
        r = _rlist[_i];
        e = _elist[_i];
    }

    /**
     * @brief Set the next interval to subdivide based on the largest absolute
     * error
     */
    void qpsrt() {
        const unsigned int last = _size - 1;
        const unsigned int limit = _limit;

        double* elist = _elist;
        unsigned int* order = _order;

        double errmax;
        double errmin;
        int i, k, top;

        unsigned int i_nrmax = _nrmax;
        unsigned int i_maxerr = order[i_nrmax];

        /* Check whether the list contains more than two error estimates */

        if(last < 2) {
            order[0] = 0;
            order[1] = 1;
            _i = i_maxerr;
            return;
        }

        errmax = elist[i_maxerr];

        /* This part of the routine is only executed if, due to a difficult
           integrand, subdivision increased the error estimate. In the normal
           case the insert procedure should start after the nrmax-th largest
           error estimate. */

        while(i_nrmax > 0 && errmax > elist[order[i_nrmax - 1]]) {
            order[i_nrmax] = order[i_nrmax - 1];
            i_nrmax--;
        }

        /* Compute the number of elements in the list to be maintained in
           descending order. This number depends on the number of
           subdivisions still allowed. */

        if(last < (limit / 2 + 2)) {
            top = last;
        } else {
            top = limit - last + 1;
        }

        /* Insert errmax by traversing the list top-down, starting
           comparison from the element elist(order(i_nrmax+1)). */

        i = i_nrmax + 1;

        /* The order of the tests in the following line is important to
           prevent a segmentation fault */

        while(i < top && errmax < elist[order[i]]) {
            order[i - 1] = order[i];
            i++;
        }

        order[i - 1] = i_maxerr;

        /* Insert errmin by traversing the list bottom-up */

        errmin = elist[last];

        k = top - 1;

        while(k > i - 2 && errmin >= elist[order[k]]) {
            order[k + 1] = order[k];
            k--;
        }

        order[k + 1] = last;

        /* Set i_max and e_max */

        i_maxerr = order[i_nrmax];

        _i = i_maxerr;
        _nrmax = i_nrmax;
    }

    /**
     * @brief Subdivide the current interval and set the values of the two new
     * intervals
     *
     * @param a1 Lower bound of the first interval
     * @param b1 Upper bound of the first interval
     * @param area1 Result for the first interval
     * @param error1 Absolute error for the first interval
     * @param a2 Lower bound of the second interval
     * @param b2 Upper bound of the second interval
     * @param area2 Result for the second interval
     * @param error2 Absolute error for the second interval
     */
    void update(double a1, double b1, double area1, double error1, double a2,
                double b2, double area2, double error2) {
        double* alist = _alist;
        double* blist = _blist;
        double* rlist = _rlist;
        double* elist = _elist;
        unsigned int* level = _level;

        const unsigned int i_max = _i;
        const unsigned int i_new = _size;

        const unsigned int new_level = _level[i_max] + 1;

        /* append the newly-created intervals to the list */

        if(error2 > error1) {
            alist[i_max] = a2; /* blist[maxerr] is already == b2 */
            rlist[i_max] = area2;
            elist[i_max] = error2;
            level[i_max] = new_level;

            alist[i_new] = a1;
            blist[i_new] = b1;
            rlist[i_new] = area1;
            elist[i_new] = error1;
            level[i_new] = new_level;
        } else {
            blist[i_max] = b1; /* alist[maxerr] is already == a1 */
            rlist[i_max] = area1;
            elist[i_max] = error1;
            level[i_max] = new_level;

            alist[i_new] = a2;
            blist[i_new] = b2;
            rlist[i_new] = area2;
            elist[i_new] = error2;
            level[i_new] = new_level;
        }

        _size++;

        if(new_level > _maximum_level) {
            _maximum_level = new_level;
        }

        qpsrt();
    }

    /**
     * @brief Sum the results of all subintervals to get the total result
     *
     * @return Total result
     */
    double sum_results() {
        const double* const rlist = _rlist;
        const unsigned int n = _size;

        unsigned int k;
        double result_sum = 0;

        for(k = 0; k < n; k++) {
            result_sum += rlist[k];
        }

        return result_sum;
    }
};

/**
 * @brief Rescale the given error to account for round off error
 *
 * @param err Error estimate on the result
 * @param result_abs No idea
 * @param result_asc No idea
 * @return New error estimate
 */
double GSL::rescale_error(double err, const double result_abs,
                          const double result_asc) {
    err = fabs(err);

    if(result_asc != 0 && err != 0) {
        double scale = pow((200 * err / result_asc), 1.5);

        if(scale < 1) {
            err = result_asc * scale;
        } else {
            err = result_asc;
        }
    }
    if(result_abs > 2.2250738585072014e-308 / (50 * 2.2204460492503131e-16)) {
        double min_err = 50 * 2.2204460492503131e-16 * result_abs;

        if(min_err > err) {
            err = min_err;
        }
    }

    return err;
}

/**
 * @brief Get the quadrature of the given function over the given interval using
 * a combined 20 point Gauss quadrature and a 41 point Kronrod rule
 *
 * @param f Function to integrate
 * @param x0 Lower bound of the interval
 * @param x1 Upper bound of the interval
 * @param params Extra parameters for the function
 * @param abserr Absolute error on the result, estimated from the difference
 * between the Gauss and the Kronrod result
 * @param resabs No idea
 * @param resasc No idea
 * @return Quadrature of the function over the interval
 */
double GSL::qk(double (*f)(double, void*), double x0, double x1, void* params,
               double& abserr, double& resabs, double& resasc) {
    double fv1[21];
    double fv2[21];

    const double center = 0.5 * (x0 + x1);
    const double half_length = 0.5 * (x1 - x0);
    const double abs_half_length = fabs(half_length);
    const double f_center = f(center, params);

    double result_gauss = 0.;
    double result_kronrod = f_center * wgk[20];

    double result_abs = fabs(result_kronrod);
    double result_asc = 0.;
    double mean = 0., err = 0.;

    int j;

    for(j = 0; j < 10; j++) {
        const int jtw = j * 2 + 1;
        const double abscissa = half_length * xgk[jtw];
        const double fval1 = f(center - abscissa, params);
        const double fval2 = f(center + abscissa, params);
        const double fsum = fval1 + fval2;
        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        result_gauss += wg[j] * fsum;
        result_kronrod += wgk[jtw] * fsum;
        result_abs += wgk[jtw] * (fabs(fval1) + fabs(fval2));
    }

    for(j = 0; j < 10; j++) {
        int jtwm1 = j * 2;
        const double abscissa = half_length * xgk[jtwm1];
        const double fval1 = f(center - abscissa, params);
        const double fval2 = f(center + abscissa, params);
        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        result_kronrod += wgk[jtwm1] * (fval1 + fval2);
        result_abs += wgk[jtwm1] * (fabs(fval1) + fabs(fval2));
    };

    mean = result_kronrod * 0.5;

    result_asc = wgk[20] * fabs(f_center - mean);

    for(j = 0; j < 20; j++) {
        result_asc += wgk[j] * (fabs(fv1[j] - mean) + fabs(fv2[j] - mean));
    }

    /* scale by the width of the integration region */

    err = (result_kronrod - result_gauss) * half_length;

    result_kronrod *= half_length;
    result_abs *= abs_half_length;
    result_asc *= abs_half_length;

    resabs = result_abs;
    resasc = result_asc;
    abserr = rescale_error(err, result_abs, result_asc);
    return result_kronrod;
}

/**
 * @brief Perform a numerical quadrature using gsl_integration_qag, with the
 * GSL_INTEG_GAUSS41 integration rule
 *
 * @param x0 Lower limit of the integral
 * @param x1 Upper limit of the integral
 * @param tolerance Relative tolerance for the result
 * @param params Parameters for the integrand function
 * @return Value of the integral
 */
double GSL::qag(double (*f)(double, void*), double x0, double x1,
                double tolerance, void* params) {
    GSL_workspace workspace(100000);

    double area, errsum;
    double result0, abserr0, resabs0, resasc0;
    double ltolerance;
    unsigned int iteration = 0;
    int roundoff_type1 = 0, roundoff_type2 = 0, error_type = 0;

    /* perform the first integration */
    result0 = qk(f, x0, x1, params, abserr0, resabs0, resasc0);

    workspace.initialise(x0, x1);
    workspace.set_initial_result(result0, abserr0);

    /* Test on accuracy */
    ltolerance = tolerance * fabs(result0);

    if((abserr0 <= ltolerance && abserr0 != resasc0) || abserr0 == 0.0) {
        return result0;
    }

    area = result0;
    errsum = abserr0;

    iteration = 1;

    do {
        double a1, b1, a2, b2;
        double a_i, b_i, r_i, e_i;
        double area1 = 0, area2 = 0, area12 = 0;
        double error1 = 0, error2 = 0, error12 = 0;
        double resasc1, resasc2;
        double resabs1, resabs2;

        /* Bisect the subinterval with the largest error estimate */

        workspace.retrieve(a_i, b_i, r_i, e_i);

        a1 = a_i;
        b1 = 0.5 * (a_i + b_i);
        a2 = b1;
        b2 = b_i;

        area1 = qk(f, a1, b1, params, error1, resabs1, resasc1);
        area2 = qk(f, a2, b2, params, error2, resabs2, resasc2);

        area12 = area1 + area2;
        error12 = error1 + error2;

        errsum += (error12 - e_i);
        area += area12 - r_i;

        if(resasc1 != error1 && resasc2 != error2) {
            double delta = r_i - area12;

            if(fabs(delta) <= 1.0e-5 * fabs(area12) && error12 >= 0.99 * e_i) {
                roundoff_type1++;
            }
            if(iteration >= 10 && error12 > e_i) {
                roundoff_type2++;
            }
        }

        ltolerance = tolerance * fabs(area);

        // we don't implement error checks

        workspace.update(a1, b1, area1, error1, a2, b2, area2, error2);

        iteration++;

    } while(iteration < 100000 && !error_type && errsum > ltolerance);

    return workspace.sum_results();
}

/**
 * @brief Find the root of the given function in the given interval using
 * Brent's method
 *
 * @param xlow Lower bound of the search interval
 * @param xupp Upper bound of the search interval
 * @param tolerance Desired relative accuracy for the result
 * @param params Extra parameters passed on to the function
 * @return Root of the function
 */
double GSL::brent(double (*f)(double, void*), double xlow, double xupp,
                  double tolerance, void* params) {
    // initial estimate
    double root = 0.5 * (xlow + xupp);

    // brent_init
    double a, b, c, d, e;
    double fa, fb, fc;

    double tol, m;

    a = xlow;
    b = xupp;
    fa = f(xlow, params);
    fb = f(xupp, params);

    if((fa < 0. && fb < 0.) || (fa > 0. && fb > 0.)) {
        std::cerr << "Lower and upper bound for Brent's method have same sign!"
                  << std::endl;
        my_exit();
        return std::numeric_limits<double>::quiet_NaN();
    }

    c = xupp;
    fc = fb;
    d = xupp - xlow;
    e = d;

    double tolerance2 = fabs(d);
    while(fabs(xupp - xlow) >= tolerance2) {
        bool ac_equal = false;
        // brent_iterate
        if((fb < 0. && fc < 0.) || (fb > 0. && fc > 0.)) {
            ac_equal = true;
            c = a;
            fc = fa;
            d = b - a;
            e = b - a;
        }

        if(fabs(fc) < fabs(fb)) {
            ac_equal = true;
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }

        tol = 0.5 * 2.2204460492503131e-16 * fabs(b);
        m = 0.5 * (c - b);

        if(!fb) {
            root = b;
            xlow = b;
            xupp = b;

            if((xlow > 0. && xupp > 0.) || (xlow < 0. && xupp < 0.)) {
                tolerance2 = tolerance * std::min(fabs(xupp), fabs(xlow));
            } else {
                tolerance2 = 0.;
            }
            continue;
        }

        if(fabs(m) <= tol) {
            root = b;

            if(b < c) {
                xlow = b;
                xupp = c;
            } else {
                xlow = c;
                xupp = b;
            }

            if((xlow > 0. && xupp > 0.) || (xlow < 0. && xupp < 0.)) {
                tolerance2 = tolerance * std::min(fabs(xupp), fabs(xlow));
            } else {
                tolerance2 = 0.;
            }
            continue;
        }

        if(fabs(e) < tol || fabs(fa) <= fabs(fb)) {
            // use bisection
            d = m;
            e = m;
        } else {
            // use inverse cubic interpolation
            double p, q, r;
            double s = fb / fa;

            if(ac_equal) {
                p = 2 * m * s;
                q = 1 - s;
            } else {
                q = fa / fc;
                r = fb / fc;
                p = s * (2 * m * q * (q - r) - (b - a) * (r - 1));
                q = (q - 1) * (r - 1) * (s - 1);
            }

            if(p > 0) {
                q = -q;
            } else {
                p = -p;
            }

            if(2 * p < std::min(3 * m * q - fabs(tol * q), fabs(e * q))) {
                e = d;
                d = p / q;
            } else {
                // interpolation failed, fall back to bisection
                d = m;
                e = m;
            }
        }

        a = b;
        fa = fb;

        if(fabs(d) > tol) {
            b += d;
        } else {
            if(m > 0.) {
                b += tol;
            } else {
                b -= tol;
            }
        }

        fb = f(b, params);

        // Update the best estimate of the root and bounds on each iteration

        root = b;

        if((fb < 0. && fc < 0.) || (fb > 0. && fc > 0.)) {
            c = a;
        }

        if(b < c) {
            xlow = b;
            xupp = c;
        } else {
            xlow = c;
            xupp = b;
        }

        if((xlow > 0. && xupp > 0.) || (xlow < 0. && xupp < 0.)) {
            tolerance2 = tolerance * std::min(fabs(xupp), fabs(xlow));
        } else {
            tolerance2 = 0.;
        }
    }

    return root;
}

/**
 * @brief Find the lower index of the bin containing the given value in the
 * given array, using bisection search
 *
 * @param xvals Array of X-values
 * @param x X-value to find
 * @param ilow Initial guess for the lower bin index
 * @param iupp Initial guess for the upper bin index
 * @return Lower index of the bin containing x
 */
unsigned int GSL::bsearch(const double* xvals, double x, unsigned int ilow,
                          unsigned int iupp) {
    while(iupp > ilow + 1) {
        unsigned int i = (iupp + ilow) / 2;
        if(xvals[i] > x) {
            iupp = i;
        } else {
            ilow = i;
        }
    }

    return ilow;
}

/**
 * @brief Solve the given set of linear equations
 *
 * @param diag Diagonal elements
 * @param offdiag Off-diagonal elemetns
 * @param b Coefficients
 * @param x Solution array
 * @param N Number of elements
 */
void GSL::solve_tridiag(const double* diag, const double* offdiag,
                        const double* b, double* x, unsigned int N) {
    double* gamma = new double[N];
    double* alpha = new double[N];
    double* c = new double[N];
    double* z = new double[N];

    /* Cholesky decomposition
           A = L.D.L^t
           lower_diag(L) = gamma
           diag(D) = alpha
         */
    alpha[0] = diag[0];
    gamma[0] = offdiag[0] / alpha[0];

    if(alpha[0] == 0) {
        std::cerr << "Divide by zero!" << std::endl;
        my_exit();
    }

    for(unsigned int i = 1; i < N - 1; i++) {
        alpha[i] = diag[i] - offdiag[i - 1] * gamma[i - 1];
        gamma[i] = offdiag[i] / alpha[i];
        if(alpha[i] == 0) {
            std::cerr << "Divide by zero!" << std::endl;
            my_exit();
        }
    }

    if(N > 1) {
        alpha[N - 1] = diag[N - 1] - offdiag[N - 2] * gamma[N - 2];
    }

    /* update RHS */
    z[0] = b[0];
    for(unsigned int i = 1; i < N; i++) {
        z[i] = b[i] - gamma[i - 1] * z[i - 1];
    }
    for(unsigned int i = 0; i < N; i++) {
        c[i] = z[i] / alpha[i];
    }

    /* backsubstitution */
    x[N - 1] = c[N - 1];
    if(N >= 2) {
        for(unsigned int i = N - 2, j = 0; j <= N - 2; j++, i--) {
            x[i] = c[i] - gamma[i] * x[i + 1];
        }
    }

    delete[] z;
    delete[] c;
    delete[] alpha;
    delete[] gamma;
}

/**
 * @brief Constructor
 *
 * Allocate the internal arrays and copy the original data
 *
 * @param x X-values
 * @param y Y-values
 * @param size Size of the data arrays
 */
GSL::GSLInterpolator::GSLInterpolator(double* x, double* y, unsigned int size) {
    _x = new double[size];
    for(unsigned int i = 0; i < size; i++) {
        _x[i] = x[i];
    }

    _y = new double[size];
    for(unsigned int i = 0; i < size; i++) {
        _y[i] = y[i];
    }

    _size = size;
}

/**
 * @brief Destructor
 *
 * Deallocate internal arrays
 */
GSL::GSLInterpolator::~GSLInterpolator() {
    delete[] _x;
    delete[] _y;
}

/**
 * @brief Factory function for GSLInterpolators
 *
 * @param type GSLInterpolatorType of the desired GSLInterpolator
 * @param x X-values
 * @param y Y-values
 * @param size Size of the data arrays
 * @return Pointer to a new GSLInterpolator implementation instance
 */
GSL::GSLInterpolator* GSL::GSLInterpolator::create(GSLInterpolatorType type,
                                                   double* x, double* y,
                                                   unsigned int size) {
    if(type == TypeGSLLinearInterpolator) {
        return new GSLLinearInterpolator(x, y, size);
    }

    if(type == TypeGSLCubicSplineInterpolator) {
        return new GSLCubicSplineInterpolator(x, y, size);
    }

    return NULL;
}

/**
 * @brief Constructor
 *
 * @param x X-values
 * @param y Y-values
 * @param size Size of the data arrays
 */
GSL::GSLLinearInterpolator::GSLLinearInterpolator(double* x, double* y,
                                                  unsigned int size)
        : GSLInterpolator(x, y, size) {
    _cache = 0;
}

/**
 * @brief Destructor
 */
GSL::GSLLinearInterpolator::~GSLLinearInterpolator() {}

/**
 * @brief Evaluate the interpolating function at the given position
 *
 * @param x X-value
 * @return Value of the linear interpolation
 */
double GSL::GSLLinearInterpolator::eval(double x) {
    double x_lo, x_hi;
    double y_lo, y_hi;
    double dx;

    if(x < _x[_cache]) {
        _cache = bsearch(_x, x, 0, _cache);
    } else {
        if(x >= _x[_cache + 1]) {
            _cache = bsearch(_x, x, _cache, _size - 1);
        }
    }

    // evaluate
    x_lo = _x[_cache];
    x_hi = _x[_cache + 1];
    y_lo = _y[_cache];
    y_hi = _y[_cache + 1];
    dx = x_hi - x_lo;
    if(dx > 0.0) {
        return y_lo + (x - x_lo) / dx * (y_hi - y_lo);
    } else {
        // error
        return 0.;
    }
}

/**
 * @brief Constructor
 *
 * @param x X-values
 * @param y Y-values
 * @param size Size of the data arrays
 */
GSL::GSLCubicSplineInterpolator::GSLCubicSplineInterpolator(double* x,
                                                            double* y,
                                                            unsigned int size)
        : GSLInterpolator(x, y, size) {
    _c = new double[size];
    _g = new double[size];
    _diag = new double[size];
    _offdiag = new double[size];

    _cache = 0;

    unsigned int max_index = size - 1;
    unsigned int sys_size = max_index - 1;

    _c[0] = 0.;
    _c[max_index] = 0.;

    for(unsigned int i = 0; i < sys_size; i++) {
        double h_i = _x[i + 1] - _x[i];
        double h_ip1 = _x[i + 2] - _x[i + 1];
        double ydiff_i = _y[i + 1] - _y[i];
        double ydiff_ip1 = _y[i + 2] - _y[i + 1];
        double g_i;
        if(h_i) {
            g_i = 1. / h_i;
        } else {
            g_i = 0.;
        }
        double g_ip1;
        if(h_ip1) {
            g_ip1 = 1. / h_ip1;
        } else {
            g_ip1 = 0.;
        }
        _offdiag[i] = h_ip1;
        _diag[i] = 2.0 * (h_ip1 + h_i);
        _g[i] = 3.0 * (ydiff_ip1 * g_ip1 - ydiff_i * g_i);
    }

    if(sys_size == 1) {
        _c[1] = _g[0] / _diag[0];
    } else {
        solve_tridiag(_diag, _offdiag, _g, _c, sys_size);
    }
}

/**
 * @brief Destructor
 */
GSL::GSLCubicSplineInterpolator::~GSLCubicSplineInterpolator() {
    delete[] _c;
    delete[] _g;
    delete[] _diag;
    delete[] _offdiag;
}

/**
 * @brief Evaluate the cubic spline at the given position
 *
 * @param x X-value
 * @return Value of the interpolating cubic spline
 */
double GSL::GSLCubicSplineInterpolator::eval(double x) {
    double x_lo, x_hi;
    double dx;

    if(x < _x[_cache]) {
        _cache = bsearch(_x, x, 0, _cache);
    } else {
        if(x >= _x[_cache + 1]) {
            _cache = bsearch(_x, x, _cache, _size - 1);
        }
    }

    /* evaluate */
    x_hi = _x[_cache + 1];
    x_lo = _x[_cache];
    dx = x_hi - x_lo;
    if(dx > 0.) {
        double y_lo = _y[_cache];
        double y_hi = _y[_cache + 1];
        double dy = y_hi - y_lo;
        double delx = x - x_lo;

        double c_i = _c[_cache];
        double c_ip1 = _c[_cache + 1];
        double b_i = (dy / dx) - dx * (c_ip1 + 2. * c_i) / 3.;
        double d_i = (c_ip1 - c_i) / (3. * dx);

        return y_lo + delx * (b_i + delx * (c_i + delx * d_i));
    } else {
        return 0.;
    }
}
