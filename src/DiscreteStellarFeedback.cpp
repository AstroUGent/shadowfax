/*******************************************************************************
 * This file is part of Shadowfax
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *                    Sven De Rijcke (sven.derijcke@ugent.be)
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
 * @file DiscreteStellarFeedback.cpp
 *
 * @brief StellarFeedback implementation that takes into account the discrete
 * nature of stellar feedback: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 * @author Sven De Rijcke (sven.derijcke@ugent.be)
 */
#include "DiscreteStellarFeedback.hpp"
#include "DiscreteStellarFeedbackData.hpp"
#include "RestartFile.hpp"
#include "io/UnitConverter.hpp"
#include "io/UnitSet.hpp"
#include "utilities/GasParticle.hpp"
#include "utilities/HelperFunctions.hpp"
#include "utilities/StarParticle.hpp"
#include <cmath>
using namespace std;

/**
 * @brief Low mass part of the Chabrier IMF
 *
 * @param m Mass value (in solar masses)
 * @return Value of the low-mass part of the Chabrier IMF
 */
double DiscreteStellarFeedback::PopII_IMF_sub1(double m) {
    double a = log10(m) + 1.1024;
    return exp(-0.5 * (a * a / 0.4761)) / m;
}

/**
 * @brief Chabrier IMF
 *
 * @param m Mass value (in solar masses)
 * @return Value of the Chabrier IMF
 */
double DiscreteStellarFeedback::PopII_IMF(double m) {
    if(m > _PopII_M_low && m < _PopII_M_upp) {
        if(m < 1.) {
            return _PopII_fac_IMF * PopII_IMF_sub1(m);
        } else {
            return pow(m, -2.3);
        }
    } else {
        return 0.;
    }
}

/**
 * @brief Chabrier IMF mass integrand
 *
 * @param m Mass value (in solar masses)
 * @return Integrand for the Chabrier IMF mass integral
 */
double DiscreteStellarFeedback::PopII_mIMF(double m) {
    return m * PopII_IMF(m);
}

/**
 * @brief Get the lower mass of an interval starting at the given upper mass
 * that contains exactly one star under a Chabrier IMF
 *
 * @param pmass Mass of the stellar population (in solar masses)
 * @param mupp Upper mass of the interval (in solar masses)
 * @return Lower mass (mlow) of the interval, such that the interval [mlow,mupp]
 * contains exactly one star under a Chabrier IMF (in solar masses)
 */
double DiscreteStellarFeedback::PopII_interval_mass(double pmass, double mupp) {
    return pow(1.3 * _PopII_Mint / pmass + pow(mupp, -1.3), -10. / 13.);
}

/**
 * @brief Get the lifetime of a star with the given mass
 *
 * @param m Mass of a star (in solar masses)
 * @return Lifetime of a star with the given mass (in yr)
 */
double DiscreteStellarFeedback::PopII_lifetime(double m) {
    return pow(10., _PopII_lifetime_spline->eval(m));
}

/**
 * @brief First part of the SNIa delay time distribution
 *
 * See Mannucci et al. 2006 (MNRAS, 370, pp. 773-783)
 *
 * @param t Time value (in Gyr)
 * @return First part of the SNIa delay time distribution
 */
double DiscreteStellarFeedback::PopII_SNIa_delay1(double t) {
    double a = (t - _PopII_SNIa_delay_mu) / _PopII_SNIa_delay_sigma;
    return _PopII_SNIa_delay_norm1 * (t - 0.03) * (13.6 - t) *
           exp(-0.5 * a * a);
}

/**
 * @brief Second part of the SNIa delay time distribution
 *
 * See Mannucci et al. 2006 (MNRAS, 370, pp. 773-783)
 *
 * @param t Time value (in Gyr)
 * @return Second part of the SNIa delay time distribution
 */
double DiscreteStellarFeedback::PopII_SNIa_delay2(double t) {
    double delay;
    if(t < 0.25) {
        delay = _PopII_SNIa_delay_norm2 *
                (exp((t - 0.25) / 0.1) - exp((0.03 - 0.25) / 0.1));
    } else {
        delay = _PopII_SNIa_delay_norm2 *
                (exp((0.25 - t) / 7.) - exp((0.03 - 0.25) / 0.1));
    }

    if(delay > 0.) {
        return delay;
    } else {
        return 0.;
    }
}

/**
 * @brief SNIa delay time distribution
 *
 * See Mannucci et al. 2006 (MNRAS, 370, pp. 773-783)
 *
 * @param t Time value (in Gyr)
 * @return SNIa delay time distribution
 */
double DiscreteStellarFeedback::PopII_SNIa_delay(double t) {
    return PopII_SNIa_delay1(t) + PopII_SNIa_delay2(t);
}

/**
 * @brief Cumulative SNIa delay time distribution
 *
 * The distribution is interpolated using a cubic spline.
 *
 * @param t Time value (in Gyr)
 * @return Cumulative SNIa delay time distribution
 */
double DiscreteStellarFeedback::PopII_SNIa_delay_cumul(double t) {
    return _PopII_SNIa_delay_spline->eval(log10(t));
}

/**
 * @brief Auxiliary function used to find the upper bound of a time interval
 * that contains exactly one SNIa
 *
 * @param t Time value (in Gyr)
 * @param NIa Number of SNIa
 * @param vlow Constant that depends on the lower bound of the interval
 * @return Value that should be as close to zero as possible
 */
double DiscreteStellarFeedback::PopII_SNIa_interval_func(double t, double NIa,
                                                         double vlow) {
    return 1. - NIa * PopII_SNIa_delay_cumul(t) + vlow;
}

/**
 * @brief Get the upper limit of a time interval starting at the given lower
 * limit that contains exactly one SNIa if the total number of SNIa should be
 * equal to the given number
 *
 * @param NIa Total number of SNIa
 * @param tlow Lower limit of the time interval (in Gyr)
 * @return Upper limit of the time interval (in Gyr)
 */
double DiscreteStellarFeedback::PopII_SNIa_interval_time(double NIa,
                                                         double tlow) {
    struct PopII_SNIa_interval_func_static_params params;
    params.object = this;
    params.NIa = NIa;
    params.vlow = NIa * PopII_SNIa_delay_cumul(tlow);

    return GSL::brent(&PopII_SNIa_interval_func_static, 0.02, 13.7,
                      4.4408920985006262e-16, &params);
}

/**
 * @brief Susa PopIII IMF
 *
 * @param m Mass value (in solar masses)
 * @return Value of the Susa Pop III IMF
 */
double DiscreteStellarFeedback::PopIII_IMF(double m) {
    if(m > _PopIII_M_low && m < _PopIII_M_upp) {
        double logm = log10(m);
        double IMF;
        if(logm < _PopIII_M2) {
            IMF = 0.5 * (logm - _PopIII_M1) / (_PopIII_M2 - _PopIII_M1);
        } else {
            IMF = 0.5 * (logm + _PopIII_M3 - 2. * _PopIII_M2) /
                  (_PopIII_M3 - _PopIII_M2);
        }
        IMF = _PopIII_fac * pow(IMF * (1.0 - IMF), _PopIII_pw);
        if(IMF > 0.) {
            return IMF / m;
        } else {
            return 0.;
        }
    } else {
        return 0.;
    }
}

/**
 * @brief Susa PopIII IMF mass integrand
 *
 * @param m Mass value (in solar masses)
 * @return Integrand for the Susa PopIII IMF mass integral
 */
double DiscreteStellarFeedback::PopIII_mIMF(double m) {
    return m * PopIII_IMF(m);
}

/**
 * @brief Energy of a PopIII SN with a given mass
 *
 * @param m Mass value (in solar masses)
 * @return Energy of a PopIII SN with that mass (in erg)
 */
double DiscreteStellarFeedback::PopIII_E_SN(double m) {
    return 1.e51 * _PopIII_E_SN_spline->eval(m);
}

/**
 * @brief Integrand of the PopIII IMF energy integral
 *
 * @param m Mass value (in solar masses)
 * @return Integrand value
 */
double DiscreteStellarFeedback::PopIII_EIMF(double m) {
    return PopIII_E_SN(m) * PopIII_IMF(m);
}

/**
 * @brief Cumulative Susa PopIII IMF
 *
 * Stored internally as a gsl_spline
 *
 * @param m Mass value (in solar masses)
 * @return Value of the cumulative Susa PopIII IMF
 */
double DiscreteStellarFeedback::PopIII_IMF_cumul(double m) {
    return _PopIII_IMF_spline->eval(log10(m));
}

/**
 * @brief PopIII IMF auxiliary function, used to find a mass interval containing
 * exactly one PopIII SN
 *
 * @param m Mass value (in solar masses)
 * @param pmass Initial mass of the stellar particle (in solar masses)
 * @param vlow Constant that depends on the upper bound of the interval
 * @return Value that should be as close to zero as possible
 */
double DiscreteStellarFeedback::PopIII_IMF_interval_func(double m, double pmass,
                                                         double vlow) {
    return 1. - pmass / _PopIII_Mint * PopIII_IMF_cumul(m) + vlow;
}

/**
 * @brief Get the lower mass of an interval starting at the given upper mass
 * that contains exactly one star under the Susa PopIII IMF
 *
 * @param pmass Mass of the stellar population (in solar masses)
 * @param mupp Upper mass of the interval (in solar masses)
 * @return Lower mass (mlow) of the interval, such that the interval [mlow,mupp]
 * contains exactly one star under the Susa PopIII IMF (in solar masses)
 */
double DiscreteStellarFeedback::PopIII_interval_mass(double pmass,
                                                     double mupp) {
    struct PopIII_IMF_interval_func_static_params pars;
    pars.object = this;
    pars.pmass = pmass;
    pars.vlow = pmass / _PopIII_Mint * PopIII_IMF_cumul(mupp);

    return GSL::brent(&PopIII_IMF_interval_func_static, _PopIII_M_SN_low,
                      _PopIII_M_upp, 4.4408920985006262e-16, &pars);
}

/**
 * @brief Get the lifetime of a star with the given mass
 *
 * @param m Mass of a star (in solar masses)
 * @return Lifetime of a star with the given mass (in Gyr)
 */
double DiscreteStellarFeedback::PopIII_lifetime(double m) {
    return pow(10., _PopIII_lifetime_spline->eval(m));
}

/**
 * @brief Initialize the internal variables
 */
void DiscreteStellarFeedback::initialize() {
    // Chabrier IMF: mass interval and normalization factor for low mass part
    // internally, we store all of these in solar masses
    // these values are never exposed outside this class, so we do not need to
    // perform unit conversions
    _PopII_M_low = 0.07;     // in solar masses
    _PopII_M_upp = 100.;     // in solar masses
    _PopII_M_SNII_low = 8.;  // in solar masses
    _PopII_M_SNIa_low = 3.;  // in solar masses
    _PopII_M_SNIa_upp = 8.;  // in solar masses
    _PopII_fac_IMF = 1. / PopII_IMF_sub1(1.);

    // SNIa delay time distribution
    // again, no unit conversions are needed as long as we keep these inside
    _PopII_SNIa_delay_mu = 0.05;     // in Gyr
    _PopII_SNIa_delay_sigma = 0.01;  // in Gyr
    _PopII_SNIa_delay_norm1 = 1.;
    _PopII_SNIa_delay_norm2 = 1.;

    // Cutoff metallicity for PopIII stars
    // unit conversion might apply in the future. CHECK THIS!
    _PopIII_cutoff = -5.;  // [Fe/H]

    // PopIII IMF: mass interval and factors (come from Sven's Python script)
    // no unit conversions necessary
    _PopIII_M_low = 0.7;     // in solar masses
    _PopIII_M_upp = 500.;    // in solar masses
    _PopIII_M_SN_low = 10.;  // in solar masses
    _PopIII_M1 = log10(_PopIII_M_low);
    _PopIII_M2 = 1.51130759;
    _PopIII_M3 = log10(_PopIII_M_upp);
    _PopIII_fac = 708.92544818;
    _PopIII_pw = 2.8008394;

    // Calculate number of SNII, SNIa and PopIII SN by integrating the IMFs
    _PopII_Mint = GSL::qag(&PopII_mIMF_static, _PopII_M_low, _PopII_M_upp,
                           1.e-8, this);  // in solar masses
    _PopII_NIIint = GSL::qag(&PopII_IMF_static, _PopII_M_SNII_low, _PopII_M_upp,
                             1.e-8, this);
    _PopII_NIaint = GSL::qag(&PopII_IMF_static, _PopII_M_SNIa_low,
                             _PopII_M_SNIa_upp, 1.e-8, this);

    _PopIII_Mint = GSL::qag(&PopIII_mIMF_static, _PopIII_M_low, _PopIII_M_upp,
                            1.e-8, this);  // in solar masses
    _PopIII_Nint = GSL::qag(&PopIII_IMF_static, _PopIII_M_SN_low, _PopIII_M_upp,
                            1.e-8, this);

    // Normalize the SNIa delay time function
    double dint = GSL::qag(&PopII_SNIa_delay1_static, 0.03, 13.6, 1.e-8, this);
    _PopII_SNIa_delay_norm1 = 0.4 / dint;
    dint = GSL::qag(&PopII_SNIa_delay2_static, 0.03, 13.6, 1.e-8, this);
    _PopII_SNIa_delay_norm2 = 0.6 / dint;

    // Initialize the cumulative SNIa delay time spline
    double ts[30], cumul_delay[30];
    unsigned int i = 1;
    double t;  // in Gyr
    ts[0] = -2.;
    cumul_delay[0] = 0.;
    for(t = log10(0.03); t < log10(13.6); t += 0.1) {
        ts[i] = t;
        cumul_delay[i] = GSL::qag(&PopII_SNIa_delay_static, 0.03, pow(10., t),
                                  1.e-8, this);
        i++;
    }
    ts[28] = log10(13.6);
    cumul_delay[28] = 1.;
    ts[29] = log10(13.8);
    cumul_delay[29] = 1.;
    ts[0] += cumul_delay[1];
    _PopII_SNIa_delay_spline = GSL::GSLInterpolator::create(
            GSL::TypeGSLCubicSplineInterpolator, ts, cumul_delay, 30);

    // Initialize the cumulative PopIII IMF spline
    // We need this many samples to come reliably close to Sven's intervals
    // without needing to evaluate the cumulative distribution at runtime
    double PopIII_IMF_ms[289], PopIII_IMF_IMFs[289];
    i = 1;
    double m;  // in log10(m/solar masses)
    PopIII_IMF_ms[0] = -2.;
    PopIII_IMF_IMFs[0] =
            GSL::qag(&PopIII_IMF_static, 0., _PopIII_M_upp, 1.e-8, this);
    for(m = log10(_PopIII_M_low); m < log10(_PopIII_M_upp); m += 0.01) {
        PopIII_IMF_ms[i] = m;
        PopIII_IMF_IMFs[i] = GSL::qag(&PopIII_IMF_static, pow(10., m),
                                      _PopIII_M_upp, 1.e-8, this);
        i++;
    }
    PopIII_IMF_ms[287] = log10(_PopIII_M_upp);
    PopIII_IMF_IMFs[287] = 0.;
    PopIII_IMF_ms[288] = 3.;
    PopIII_IMF_IMFs[288] = 0.;
    _PopIII_IMF_spline =
            GSL::GSLInterpolator::create(GSL::TypeGSLLinearInterpolator,
                                         PopIII_IMF_ms, PopIII_IMF_IMFs, 289);

    // Initialize PopII lifetime spline using Sven's data
    // in solar masses:
    double PopII_lifetime_ms[87] = {
            0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95, 1,    1.05, 1.1,  1.15,
            1.2,  1.25, 1.3,  1.35, 1.4,  1.45, 1.5,  1.55, 1.6,  1.65, 1.7,
            1.75, 1.8,  1.85, 1.9,  1.95, 2,    2.05, 2.1,  2.15, 2.2,  2.25,
            2.3,  2.35, 2.4,  2.6,  2.8,  3,    3.2,  3.4,  3.6,  3.8,  4,
            4.2,  4.4,  4.6,  4.8,  5,    5.2,  5.4,  5.6,  5.8,  6,    6.2,
            6.4,  7,    8,    9,    10,   11,   12,   14,   16,   18,   20,
            24,   28,   30,   40,   45,   50,   55,   60,   65,   70,   75,
            80,   90,   95,   100,  120,  150,  200,  250,  300,  350};
    // in log10(t/yr):
    double PopII_lifetime_ts[87] = {
            10.452,  10.3415, 10.2362, 10.139,  10.0468, 9.95885, 9.87605,
            9.79803, 9.72399, 9.65217, 9.58335, 9.52138, 9.46703, 9.42591,
            9.38266, 9.33495, 9.28953, 9.24572, 9.20422, 9.16552, 9.12806,
            9.0921,  9.05728, 9.02345, 8.99101, 8.95935, 8.92926, 8.90003,
            8.87177, 8.84446, 8.81782, 8.79194, 8.7669,  8.74255, 8.7187,
            8.69608, 8.60937, 8.53033, 8.45752, 8.39064, 8.32889, 8.27122,
            8.21763, 8.16709, 8.12,    8.07545, 8.03334, 7.99336, 7.95561,
            7.91966, 7.88551, 7.85298, 7.82191, 7.79234, 7.76385, 7.73664,
            7.66115, 7.55343, 7.46282, 7.38582, 7.31927, 7.26113, 7.16375,
            7.08681, 7.02325, 6.97053, 6.88722, 6.82377, 6.79741, 6.72284,
            6.68824, 6.65918, 6.63473, 6.61362, 6.59507, 6.57899, 6.56445,
            6.55241, 6.52965, 6.52028, 6.5111,  6.48241, 6.45054, 6.41528,
            6.3914,  6.37373, 6.36056};
    _PopII_lifetime_spline = GSL::GSLInterpolator::create(
            GSL::TypeGSLCubicSplineInterpolator, PopII_lifetime_ms,
            PopII_lifetime_ts, 87);

    // Initialize PopIII lifetime spline using Sven's data
    // in solar masses:
    double PopIII_lifetime_ms[24] = {0.7, 0.8, 0.9, 1,  1.1, 1.2, 1.3, 1.4,
                                     1.5, 1.6, 1.8, 2,  2.5, 3,   4,   5,
                                     10,  15,  20,  30, 50,  70,  100, 500};
    // in log10(t/Gyr)
    double PopIII_lifetime_ts[24] = {
            1.32077,   1.13354,   0.947924,  0.78533,  0.637463,  0.516384,
            0.399708,  0.289454,  0.193159,  0.100413, -0.062413, -0.208425,
            -0.493919, -0.702468, -0.979499, -1.18585, -1.75003,  -1.96392,
            -2.08837,  -2.24245,  -2.39566,  -2.47468, -2.54276,  -2.89056};
    _PopIII_lifetime_spline = GSL::GSLInterpolator::create(
            GSL::TypeGSLCubicSplineInterpolator, PopIII_lifetime_ms,
            PopIII_lifetime_ts, 24);

    // Initialize PopIII SN energy spline using Woosley energy data
    // Calculate PopIII SN energy integral
    // in solar masses:
    double ms[18] = {0.7,  10.,  35.,  40.,  50.,  60.,
                     70.,  80.,  90.,  100., 140., 140. + 1.e-10,
                     150., 170., 200., 270., 300., 500.};
    // in 10^51 erg:
    double Es[18] = {0., 1., 1., 1.,  1.,  1.,  1.,  1.,  1.,
                     1., 1., 9., 16., 28., 44., 49., 50., 50.};
    _PopIII_E_SN_spline = GSL::GSLInterpolator::create(
            GSL::TypeGSLLinearInterpolator, ms, Es, 18);
    _PopIII_Eint = GSL::qag(&PopIII_EIMF_static, _PopIII_M_SN_low,
                            _PopIII_M_upp, 1.e-8, this);  // in erg

    // Initialize feedback values
    double feedback_efficiency = 0.7;
    _PopII_SNII_energy =
            _erg_to_internal_energy->convert(1.e51 * feedback_efficiency);
    // ratios:
    _PopII_SNII_mass = 0.191445322565;
    _PopII_SNII_metals = 0.0241439721018;
    _PopII_SNII_Fe = 0.000932719658516;
    _PopII_SNII_Mg = 0.00151412640705;

    _PopII_SNIa_energy =
            _erg_to_internal_energy->convert(1.e51 * feedback_efficiency);
    // ratios:
    _PopII_SNIa_mass = 0.00655147325196;
    _PopII_SNIa_metals = 0.00655147325196;
    _PopII_SNIa_Fe = 0.00165100587997;
    _PopII_SNIa_Mg = 0.000257789470044;

    _PopII_SW_energy =
            _erg_to_internal_energy->convert(1.e50 * feedback_efficiency);
    _PopII_SW_end_time = 31.0;  // in Myr
    _PopII_SW_end_time =
            _Gyr_to_internal_time->convert(1.e-3 * _PopII_SW_end_time);
    _PopII_SW_energy /= _PopII_SW_end_time;

    _PopIII_SN_energy = _erg_to_internal_energy->convert(_PopIII_Eint *
                                                         feedback_efficiency);
    // ratios:
    _PopIII_SN_mass = 0.45;
    _PopIII_SN_metals = 0.026;
    _PopIII_SN_Fe = 0.0000932719658516;
    _PopIII_SN_Mg = 0.000151412640705;

    _PopIII_SW_energy =
            _erg_to_internal_energy->convert(1.e51 * feedback_efficiency);
    _PopIII_SW_end_time = 16.7;  // in Myr
    _PopIII_SW_end_time =
            _Gyr_to_internal_time->convert(1.e-3 * _PopIII_SW_end_time);
    _PopIII_SW_energy /= _PopIII_SW_end_time;
}

/**
 * @brief Constructor
 *
 * Calls initialize().
 */
DiscreteStellarFeedback::DiscreteStellarFeedback(UnitSet& units) {
    // unit conversion factors
    Unit time_unit("time", "Gyr", 3.154e16);
    _Gyr_to_internal_time = new UnitConverter(time_unit, units.get_time_unit());
    Unit energy_unit("length*length*mass/time/time", "erg", 1.e-7);
    _erg_to_internal_energy =
            new UnitConverter(energy_unit, units.get_energy_unit());
    initialize();
}

/**
 * @brief Destructor
 *
 * Delete unit converters and interpolators.
 */
DiscreteStellarFeedback::~DiscreteStellarFeedback() {
    delete _Gyr_to_internal_time;
    delete _erg_to_internal_energy;

    delete _PopII_SNIa_delay_spline;
    delete _PopIII_IMF_spline;
    delete _PopII_lifetime_spline;
    delete _PopIII_lifetime_spline;
    delete _PopIII_E_SN_spline;
}

/**
 * @brief Does the given StarParticle give feedback during the next time step?
 *
 * @param star StarParticle
 * @param starttime Physical start time of the next time step
 * @param endtime Physical end time of the next time step
 * @return True if the StarParticle does feedback
 */
bool DiscreteStellarFeedback::does_feedback(StarParticle* star,
                                            double starttime, double endtime) {
    DiscreteStellarFeedbackData* data =
            (DiscreteStellarFeedbackData*)star->get_feedback_data();
    if(data->get_PopII_SNII_count() < data->get_PopII_SNII_number() &&
       endtime >= data->get_PopII_SNII_next_time()) {
        return true;
    }
    if(data->get_PopII_SNIa_count() < data->get_PopII_SNIa_number() &&
       endtime >= data->get_PopII_SNIa_next_time()) {
        return true;
    }
    if(data->get_PopII_SW_fac() &&
       (starttime - star->get_birthtime()) < _PopII_SW_end_time) {
        return true;
    }
    if(data->get_PopIII_SN_count() < data->get_PopIII_SN_number() &&
       endtime >= data->get_PopIII_SN_next_time()) {
        return true;
    }
    if(data->get_PopIII_SW_fac() &&
       (starttime - star->get_birthtime()) < _PopIII_SW_end_time) {
        return true;
    }
    return false;
}

/**
 * @brief Give discrete stellar feedback
 *
 * @param star StarParticle that does feedback
 * @param starttime Physical start time of the next time step
 * @param endtime Physical end time of the next time step
 */
void DiscreteStellarFeedback::do_feedback(StarParticle* star, double starttime,
                                          double endtime) {
    double fb_energy = 0.;
    double fb_mass = 0.;
    double fb_metals = 0.;
    double fb_fe = 0.;
    double fb_mg = 0.;

    DiscreteStellarFeedbackData* data =
            (DiscreteStellarFeedbackData*)star->get_feedback_data();

    // PopII SNII feedback
    while(data->get_PopII_SNII_count() < data->get_PopII_SNII_number() &&
          endtime >= data->get_PopII_SNII_next_time()) {
        fb_energy += data->get_PopII_SNII_fac() * _PopII_SNII_energy;
        fb_mass += data->get_PopII_SNII_fac() * _PopII_SNII_mass *
                   star->get_initial_mass();
        fb_metals += data->get_PopII_SNII_fac() * _PopII_SNII_metals *
                     star->get_initial_mass();
        fb_fe += data->get_PopII_SNII_fac() * _PopII_SNII_Fe *
                 star->get_initial_mass();
        fb_mg += data->get_PopII_SNII_fac() * _PopII_SNII_Mg *
                 star->get_initial_mass();

        // set next PopII SNII feedback time
        data->increase_PopII_SNII_count();
        if(data->get_PopII_SNII_count() < data->get_PopII_SNII_number()) {
            double mupp = data->get_PopII_SNII_interval();
            double mlow = PopII_interval_mass(star->get_initial_mass(), mupp);
            double tlow = PopII_lifetime(mupp);
            double tupp = PopII_lifetime(mlow);
            // in yr:
            double trand =
                    tlow + (tupp - tlow) * HelperFunctions::rand_double();
            data->set_PopII_SNII_next_time(
                    star->get_birthtime() +
                    _Gyr_to_internal_time->convert(1.e-9 * trand));
            data->set_PopII_SNII_interval(mlow);
        }
    }

    // PopII SNIa feedback
    while(data->get_PopII_SNIa_count() < data->get_PopII_SNIa_number() &&
          endtime >= data->get_PopII_SNIa_next_time()) {
        fb_energy += data->get_PopII_SNIa_fac() * _PopII_SNIa_energy;
        fb_mass += data->get_PopII_SNIa_fac() * _PopII_SNIa_mass *
                   star->get_initial_mass();
        fb_metals += data->get_PopII_SNIa_fac() * _PopII_SNIa_metals *
                     star->get_initial_mass();
        fb_fe += data->get_PopII_SNIa_fac() * _PopII_SNIa_Fe *
                 star->get_initial_mass();
        fb_mg += data->get_PopII_SNIa_fac() * _PopII_SNIa_Mg *
                 star->get_initial_mass();

        // set next PopII SNII feedback time
        data->increase_PopII_SNIa_count();
        if(data->get_PopII_SNIa_count() < data->get_PopII_SNIa_number()) {
            double tlow = data->get_PopII_SNIa_interval();
            double tupp = PopII_SNIa_interval_time(
                    data->get_PopII_SNIa_number(), tlow);
            // in Gyr:
            double trand =
                    tlow + (tupp - tlow) * HelperFunctions::rand_double();
            data->set_PopII_SNIa_next_time(
                    star->get_birthtime() +
                    _Gyr_to_internal_time->convert(trand));
            data->set_PopII_SNIa_interval(tupp);
        }
    }

    // PopII SW feedback
    if(data->get_PopII_SW_fac() &&
       (endtime - star->get_birthtime()) <= _PopII_SW_end_time) {
        double tmin =
                std::min(endtime, star->get_birthtime() + _PopII_SW_end_time);
        fb_energy += (tmin - starttime) * data->get_PopII_SW_fac() *
                     _PopII_SW_energy;
    }

    // PopIII SN feedback
    while(data->get_PopIII_SN_count() < data->get_PopIII_SN_number() &&
          endtime >= data->get_PopIII_SN_next_time()) {
        fb_energy += data->get_PopIII_SN_fac() * _PopIII_SN_energy;
        fb_mass += data->get_PopIII_SN_fac() * _PopIII_SN_mass *
                   star->get_initial_mass();
        fb_metals += data->get_PopIII_SN_fac() * _PopIII_SN_metals *
                     star->get_initial_mass();
        fb_fe += data->get_PopIII_SN_fac() * _PopIII_SN_Fe *
                 star->get_initial_mass();
        fb_mg += data->get_PopIII_SN_fac() * _PopIII_SN_Mg *
                 star->get_initial_mass();

        data->increase_PopIII_SN_count();
        if(data->get_PopIII_SN_count() < data->get_PopIII_SN_number()) {
            double mupp = data->get_PopIII_SN_interval();
            double mlow = PopIII_interval_mass(star->get_initial_mass(), mupp);
            double tlow = PopIII_lifetime(mupp);
            double tupp = PopIII_lifetime(mlow);
            // in Gyr:
            double trand =
                    tlow + (tupp - tlow) * HelperFunctions::rand_double();
            data->set_PopIII_SN_next_time(
                    star->get_birthtime() +
                    _Gyr_to_internal_time->convert(trand));
            data->set_PopIII_SN_interval(mlow);

            data->set_PopIII_SN_fac(
                    star->get_initial_mass() * _PopIII_Nint *
                    PopIII_E_SN(0.5 * (mlow + mupp)) / _PopIII_Mint /
                    data->get_PopIII_SN_number() / _PopIII_Eint);
        }
    }

    // PopIII SW feedback
    if(data->get_PopIII_SW_fac() &&
       (endtime - star->get_birthtime()) <= _PopIII_SW_end_time) {
        double tmin =
                std::min(endtime, star->get_birthtime() + _PopIII_SW_end_time);
        fb_energy += (tmin - starttime) * data->get_PopIII_SW_fac() *
                     _PopIII_SW_energy;
    }

    GasParticle* gas = star->get_closest_gasparticle();
    StateVector dQ;
    dQ.set_m(-fb_mass);
    dQ.set_e(-fb_energy);
    dQ.set_Fe(-fb_fe);
    dQ.set_Mg(-fb_mg);

    gas->increase_dQ(dQ);
}

/**
 * @brief Initialize the StellarFeedbackData for discrete stellar feedback for
 * the given StarParticle
 *
 * @param star StarParticle for which the feedback is done
 * @return Pointer to an initialized DiscreteStellarFeedbackData instance
 */
StellarFeedbackData* DiscreteStellarFeedback::initialize_data(
        StarParticle* star) {
    DiscreteStellarFeedbackData* data = new DiscreteStellarFeedbackData();

    double FeH = star->get_FeH();
    if(FeH <= _PopIII_cutoff) {
        // only PopIII stars
        data->set_PopII_SNII_number(0);
        data->set_PopII_SNIa_number(0);
        unsigned int PopIII_SN_number =
                star->get_initial_mass() * _PopIII_Nint / _PopIII_Mint;
        data->set_PopIII_SN_number(PopIII_SN_number);
        data->set_PopII_SW_fac(0.);
        // number fraction of PopIII SW
        double PopIII_SW_fac =
                star->get_initial_mass() * _PopIII_Nint / _PopIII_Mint;
        data->set_PopIII_SW_fac(PopIII_SW_fac);
    } else {
        unsigned int PopII_SNII_number =
                star->get_initial_mass() * _PopII_NIIint / _PopII_Mint;
        data->set_PopII_SNII_number(PopII_SNII_number);
        // we cannot simply do 0.15*PopII_SNII_number, since this has less
        // precision
        unsigned int PopII_SNIa_number =
                0.15 * star->get_initial_mass() * _PopII_NIIint / _PopII_Mint;
        data->set_PopII_SNIa_number(PopII_SNIa_number);
        data->set_PopIII_SN_number(0);
        // number fraction of PopII SW
        double PopII_SW_fac =
                star->get_initial_mass() * _PopII_NIIint / _PopII_Mint;
        data->set_PopII_SW_fac(PopII_SW_fac);
        data->set_PopIII_SW_fac(0.);
    }
    data->set_PopII_SNII_count(0);
    data->set_PopII_SNIa_count(0);
    data->set_PopIII_SN_count(0);

    if(data->get_PopII_SNII_number()) {
        double mupp = _PopII_M_upp;  // in solar masses
        // in solar masses:
        double mlow = PopII_interval_mass(star->get_initial_mass(), mupp);
        double tlow = PopII_lifetime(mupp);  // in yr
        double tupp = PopII_lifetime(mlow);  // in yr
        // in yr:
        double trand = tlow + (tupp - tlow) * HelperFunctions::rand_double();
        // UNITS. OFFSET?
        data->set_PopII_SNII_next_time(
                star->get_birthtime() +
                _Gyr_to_internal_time->convert(1.e-9 * trand));
        data->set_PopII_SNII_interval(mlow);
        // we have to divide the total SNII energy over a discrete number of SN
        // the energy per SNII will hence be a small factor larger
        double PopII_SNII_fac = star->get_initial_mass() * _PopII_NIIint /
                                _PopII_Mint / data->get_PopII_SNII_number();
        data->set_PopII_SNII_fac(PopII_SNII_fac);
    }

    if(data->get_PopII_SNIa_number()) {
        double tlow = 0.03;  // in Gyr
        // in Gyr:
        double tupp =
                PopII_SNIa_interval_time(data->get_PopII_SNIa_number(), tlow);
        // in Gyr:
        double trand = tlow + (tupp - tlow) * HelperFunctions::rand_double();
        // UNITS. OFFSET?
        data->set_PopII_SNIa_next_time(star->get_birthtime() +
                                       _Gyr_to_internal_time->convert(trand));
        data->set_PopII_SNIa_interval(tupp);
        double PopII_SNIa_fac = star->get_initial_mass() * _PopII_NIaint /
                                _PopII_Mint / data->get_PopII_SNIa_number();
        data->set_PopII_SNIa_fac(PopII_SNIa_fac);
    }

    if(data->get_PopIII_SN_number()) {
        double mupp = _PopIII_M_upp;  // in solar masses
        // in solar masses:
        double mlow = PopIII_interval_mass(star->get_initial_mass(), mupp);
        double tlow = PopIII_lifetime(mupp);  // in Gyr
        double tupp = PopIII_lifetime(mlow);  // in Gyr
        // in Gyr:
        double trand = tlow + (tupp - tlow) * HelperFunctions::rand_double();
        // UNITS. OFFSET?
        data->set_PopIII_SN_next_time(star->get_birthtime() +
                                      _Gyr_to_internal_time->convert(trand));
        data->set_PopIII_SN_interval(mlow);
        double PopIII_SN_fac = star->get_initial_mass() * _PopIII_Nint *
                               PopIII_E_SN(0.5 * (mlow + mupp)) / _PopIII_Mint /
                               _PopIII_Eint / data->get_PopIII_SN_number();
        data->set_PopIII_SN_fac(PopIII_SN_fac);
    }

    return data;
}

/**
 * @brief Dump the object to the given RestartFile
 *
 * @param rfile RestartFile to write to
 */
void DiscreteStellarFeedback::dump(RestartFile& rfile) {
    rfile.write(_PopII_M_low);
    rfile.write(_PopII_M_upp);
    rfile.write(_PopII_fac_IMF);

    rfile.write(_PopII_M_SNII_low);

    rfile.write(_PopII_M_SNIa_low);
    rfile.write(_PopII_M_SNIa_upp);
    rfile.write(_PopII_SNIa_delay_mu);
    rfile.write(_PopII_SNIa_delay_sigma);
    rfile.write(_PopII_SNIa_delay_norm1);
    rfile.write(_PopII_SNIa_delay_norm2);

    rfile.write(_PopIII_cutoff);
    rfile.write(_PopIII_M_low);
    rfile.write(_PopIII_M_upp);
    rfile.write(_PopIII_M_SN_low);
    rfile.write(_PopIII_M1);
    rfile.write(_PopIII_M2);
    rfile.write(_PopIII_M3);
    rfile.write(_PopIII_fac);
    rfile.write(_PopIII_pw);

    rfile.write(_PopII_Mint);
    rfile.write(_PopII_NIIint);
    rfile.write(_PopII_NIaint);

    rfile.write(_PopIII_Mint);
    rfile.write(_PopIII_Nint);
    rfile.write(_PopIII_Eint);

    _PopII_SNIa_delay_spline->dump(rfile);

    _PopIII_IMF_spline->dump(rfile);

    _PopII_lifetime_spline->dump(rfile);
    _PopIII_lifetime_spline->dump(rfile);

    _PopIII_E_SN_spline->dump(rfile);

    rfile.write(_PopII_SNII_energy);
    rfile.write(_PopII_SNII_mass);
    rfile.write(_PopII_SNII_metals);
    rfile.write(_PopII_SNII_Fe);
    rfile.write(_PopII_SNII_Mg);

    rfile.write(_PopII_SNIa_energy);
    rfile.write(_PopII_SNIa_mass);
    rfile.write(_PopII_SNIa_metals);
    rfile.write(_PopII_SNIa_Fe);
    rfile.write(_PopII_SNIa_Mg);

    rfile.write(_PopII_SW_energy);
    rfile.write(_PopII_SW_end_time);

    rfile.write(_PopIII_SN_energy);
    rfile.write(_PopIII_SN_mass);
    rfile.write(_PopIII_SN_metals);
    rfile.write(_PopIII_SN_Fe);
    rfile.write(_PopIII_SN_Mg);

    rfile.write(_PopIII_SW_energy);
    rfile.write(_PopIII_SW_end_time);
}

/**
 * @brief Restart constructor
 *
 * @param rfile RestartFile to read from
 */
DiscreteStellarFeedback::DiscreteStellarFeedback(RestartFile& rfile) {
    rfile.read(_PopII_M_low);
    rfile.read(_PopII_M_upp);
    rfile.read(_PopII_fac_IMF);

    rfile.read(_PopII_M_SNII_low);

    rfile.read(_PopII_M_SNIa_low);
    rfile.read(_PopII_M_SNIa_upp);
    rfile.read(_PopII_SNIa_delay_mu);
    rfile.read(_PopII_SNIa_delay_sigma);
    rfile.read(_PopII_SNIa_delay_norm1);
    rfile.read(_PopII_SNIa_delay_norm2);

    rfile.read(_PopIII_cutoff);
    rfile.read(_PopIII_M_low);
    rfile.read(_PopIII_M_upp);
    rfile.read(_PopIII_M_SN_low);
    rfile.read(_PopIII_M1);
    rfile.read(_PopIII_M2);
    rfile.read(_PopIII_M3);
    rfile.read(_PopIII_fac);
    rfile.read(_PopIII_pw);

    rfile.read(_PopII_Mint);
    rfile.read(_PopII_NIIint);
    rfile.read(_PopII_NIaint);

    rfile.read(_PopIII_Mint);
    rfile.read(_PopIII_Nint);
    rfile.read(_PopIII_Eint);

    _PopII_SNIa_delay_spline = GSL::GSLInterpolator::create(rfile);

    _PopIII_IMF_spline = GSL::GSLInterpolator::create(rfile);

    _PopII_lifetime_spline = GSL::GSLInterpolator::create(rfile);
    _PopIII_lifetime_spline = GSL::GSLInterpolator::create(rfile);

    _PopIII_E_SN_spline = GSL::GSLInterpolator::create(rfile);

    rfile.read(_PopII_SNII_energy);
    rfile.read(_PopII_SNII_mass);
    rfile.read(_PopII_SNII_metals);
    rfile.read(_PopII_SNII_Fe);
    rfile.read(_PopII_SNII_Mg);

    rfile.read(_PopII_SNIa_energy);
    rfile.read(_PopII_SNIa_mass);
    rfile.read(_PopII_SNIa_metals);
    rfile.read(_PopII_SNIa_Fe);
    rfile.read(_PopII_SNIa_Mg);

    rfile.read(_PopII_SW_energy);
    rfile.read(_PopII_SW_end_time);

    rfile.read(_PopIII_SN_energy);
    rfile.read(_PopIII_SN_mass);
    rfile.read(_PopIII_SN_metals);
    rfile.read(_PopIII_SN_Fe);
    rfile.read(_PopIII_SN_Mg);

    rfile.read(_PopIII_SW_energy);
    rfile.read(_PopIII_SW_end_time);
}
