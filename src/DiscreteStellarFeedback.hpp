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
 * @file DiscreteStellarFeedback.hpp
 *
 * @brief StellarFeedback implementation that takes into account the discrete
 * nature of stellar feedback: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 * @author Sven De Rijcke (sven.derijcke@ugent.be)
 */
#ifndef DISCRETESTELLARFEEDBACK_HPP
#define DISCRETESTELLARFEEDBACK_HPP

#include "GSL.hpp"
#include "StellarFeedback.hpp"

/**
 * @brief StellarFeedback implementation that spreads out the stellar feedback
 * uniformly over the lifetime of the various feedback mechanisms
 */
class DiscreteStellarFeedback : public StellarFeedback {
  private:
    /*! @brief Chabrier IMF: lower limit mass interval */
    double _PopII_M_low;
    /*! @brief Chabrier IMF: upper limit mass interval */
    double _PopII_M_upp;
    /*! @brief Chabier IMF: normalization factor for low mass part */
    double _PopII_fac_IMF;

    /*! @brief Lower limit for SNII mass interval */
    double _PopII_M_SNII_low;

    /*! @brief Lower limit for SNIa mass interval */
    double _PopII_M_SNIa_low;
    /*! @brief Upper limit for SNIa mass interval */
    double _PopII_M_SNIa_upp;
    /*! @brief SNIa delay time distribution: \f$\mu\f$ factor */
    double _PopII_SNIa_delay_mu;
    /*! @brief SNIa delay time distribution: \f$\sigma\f$ factor */
    double _PopII_SNIa_delay_sigma;
    /*! @brief SNIa delay time distribution: normalization of first component */
    double _PopII_SNIa_delay_norm1;
    /*! @brief SNIa delay time distribution: normalization of second
     *  component */
    double _PopII_SNIa_delay_norm2;

    /*! @brief PopIII stars: cutoff metallicity */
    double _PopIII_cutoff;
    /*! @brief PopIII IMF: lower limit mass interval */
    double _PopIII_M_low;
    /*! @brief PopIII IMF: upper limit mass interval */
    double _PopIII_M_upp;
    /*! @brief PopIII IMF: lower mass limit SN feedback */
    double _PopIII_M_SN_low;
    /*! @brief PopIII IMF: M1 constant */
    double _PopIII_M1;
    /*! @brief PopIII IMF: M2 constant */
    double _PopIII_M2;
    /*! @brief PopIII IMF: M3 constant */
    double _PopIII_M3;
    /*! @brief PopIII IMF: factor */
    double _PopIII_fac;
    /*! @brief PopIII IMF: pw? */
    double _PopIII_pw;

    /*! @brief Mass integral of the Chabrier IMF */
    double _PopII_Mint;
    /*! @brief Number of stars in the SNII mass range */
    double _PopII_NIIint;
    /*! @brief Number of stars in the SNIa mass range */
    double _PopII_NIaint;

    /*! @brief Mass integral of the Susa IMF */
    double _PopIII_Mint;
    /*! @brief Number of PopIII SN */
    double _PopIII_Nint;
    /*! @brief PopIII SN energy integral */
    double _PopIII_Eint;

    /*! @brief SNIa delay time spline */
    GSL::GSLInterpolator* _PopII_SNIa_delay_spline;

    /*! @brief PopIII IMF spline */
    GSL::GSLInterpolator* _PopIII_IMF_spline;

    /*! @brief PopII lifetime spline */
    GSL::GSLInterpolator* _PopII_lifetime_spline;
    /*! @brief PopII lifetime spline */
    GSL::GSLInterpolator* _PopIII_lifetime_spline;

    /*! @brief PopIII SN energy spline */
    GSL::GSLInterpolator* _PopIII_E_SN_spline;

    /*! @brief PopII SNII energy */
    double _PopII_SNII_energy;
    /*! @brief PopII SNII mass return */
    double _PopII_SNII_mass;
    /*! @brief PopII SNII metal mass return */
    double _PopII_SNII_metals;
    /*! @brief PopII SNII iron mass return */
    double _PopII_SNII_Fe;
    /*! @brief PopII SNII magnesium mass return */
    double _PopII_SNII_Mg;

    /*! @brief PopII SNIa energy */
    double _PopII_SNIa_energy;
    /*! @brief PopII SNIa mass return */
    double _PopII_SNIa_mass;
    /*! @brief PopII SNIa metal mass return */
    double _PopII_SNIa_metals;
    /*! @brief PopII SNIa iron mass return */
    double _PopII_SNIa_Fe;
    /*! @brief PopII SNIa magnesium mass return */
    double _PopII_SNIa_Mg;

    /*! @brief PopII SW energy */
    double _PopII_SW_energy;
    /*! @brief PopII SW upper limit time interval */
    double _PopII_SW_end_time;

    /*! @brief PopIII SN energy */
    double _PopIII_SN_energy;
    /*! @brief PopIII SN mass return */
    double _PopIII_SN_mass;
    /*! @brief PopIII SN metal mass return */
    double _PopIII_SN_metals;
    /*! @brief PopIII SN iron mass return */
    double _PopIII_SN_Fe;
    /*! @brief PopIII SN magnesium mass return */
    double _PopIII_SN_Mg;

    /*! @brief PopIII SW energy */
    double _PopIII_SW_energy;
    /*! @brief PopIII SW upper limit time interval */
    double _PopIII_SW_end_time;

    double PopII_IMF_sub1(double m);
    double PopII_IMF(double m);
    double PopII_mIMF(double m);

    double PopII_SNIa_delay1(double t);
    double PopII_SNIa_delay2(double t);
    double PopII_SNIa_delay(double t);

    /**
     * @brief Static wrapper around PopII_mIMF()
     *
     * @param m Mass value
     * @param params Void pointer containing the DiscreteStellarFeedback object
     * @return Value of PopII_mIMF(m)
     */
    static double PopII_mIMF_static(double m, void* params) {
        DiscreteStellarFeedback* object = (DiscreteStellarFeedback*)params;
        return object->PopII_mIMF(m);
    }

    /**
     * @brief Static wrapper around PopII_IMF()
     *
     * @param m Mass value
     * @param params Void pointer containing the DiscreteStellarFeedback object
     * @return Value of PopII_IMF(m)
     */
    static double PopII_IMF_static(double m, void* params) {
        DiscreteStellarFeedback* object = (DiscreteStellarFeedback*)params;
        return object->PopII_IMF(m);
    }

    /**
     * @brief Static wrapper around PopII_SNIa_delay1()
     *
     * @param t Time value
     * @param params Void pointer containing the DiscreteStellarFeedback object
     * @return Value of PopII_SNIa_delay1(t)
     */
    static double PopII_SNIa_delay1_static(double t, void* params) {
        DiscreteStellarFeedback* object = (DiscreteStellarFeedback*)params;
        return object->PopII_SNIa_delay1(t);
    }

    /**
     * @brief Static wrapper around PopII_SNIa_delay2()
     *
     * @param t Time value
     * @param params Void pointer containing the DiscreteStellarFeedback object
     * @return Value of PopII_SNIa_delay2(t)
     */
    static double PopII_SNIa_delay2_static(double t, void* params) {
        DiscreteStellarFeedback* object = (DiscreteStellarFeedback*)params;
        return object->PopII_SNIa_delay2(t);
    }

    /**
     * @brief Static wrapper around PopII_SNIa_delay()
     *
     * @param t Time value
     * @param params Void pointer containing the DiscreteStellarFeedback object
     * @return Value of PopII_SNIa_delay(t)
     */
    static double PopII_SNIa_delay_static(double t, void* params) {
        DiscreteStellarFeedback* object = (DiscreteStellarFeedback*)params;
        return object->PopII_SNIa_delay(t);
    }

    double PopIII_IMF(double m);
    double PopIII_mIMF(double m);

    double PopIII_E_SN(double m);
    double PopIII_EIMF(double m);

    /**
     * @brief Static wrapper around PopIII_mIMF()
     *
     * @param m Mass value
     * @param params Void pointer containing the DiscreteStellarFeedback object
     * @return Value of PopIII_mIMF(m)
     */
    static double PopIII_mIMF_static(double m, void* params) {
        DiscreteStellarFeedback* object = (DiscreteStellarFeedback*)params;
        return object->PopIII_mIMF(m);
    }

    /**
     * @brief Static wrapper around PopIII_IMF()
     *
     * @param m Mass value
     * @param params Void pointer containing the DiscreteStellarFeedback object
     * @return Value of PopIII_IMF(m)
     */
    static double PopIII_IMF_static(double m, void* params) {
        DiscreteStellarFeedback* object = (DiscreteStellarFeedback*)params;
        return object->PopIII_IMF(m);
    }

    /**
     * @brief Static wrapper around PopIII_EIMF()
     *
     * @param m Mass value
     * @param params Void pointer containing the DiscreteStellarFeedback object
     * @return Value of PopIII_EIMF(m)
     */
    static double PopIII_EIMF_static(double m, void* params) {
        DiscreteStellarFeedback* object = (DiscreteStellarFeedback*)params;
        return object->PopIII_EIMF(m);
    }

    void initialize();

  public:
    DiscreteStellarFeedback();

    void do_feedback(StarParticle* star, ParticleVector& particles, double dt);

    void dump(RestartFile& rfile);
    DiscreteStellarFeedback(RestartFile& rfile);
};

#endif  // DISCRETESTELLARFEEDBACK_HPP
