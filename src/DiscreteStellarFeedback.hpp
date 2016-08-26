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
 * @brief StellarFeedbackData implementation for discrete stellar feedback
 */
class DiscreteStellarFeedbackData : public StellarFeedbackData {
  private:
    /*! @brief Number of PopII SNII explosions that needs to go off */
    unsigned int _PopII_SNII_number;
    /*! @brief Number of PopII SNIa explosions that needs to go off */
    unsigned int _PopII_SNIa_number;
    /*! @brief Number of PopIII SN explosions that needs to go off */
    unsigned int _PopIII_SN_number;

    /*! @brief PopII SW factor */
    double _PopII_SW_fac;
    /*! @brief PopIII SW factor */
    double _PopIII_SW_fac;

    /*! @brief Number of PopII SNII explosions that has gone off */
    unsigned int _PopII_SNII_count;
    /*! @brief Number of PopII SNIa explosions that has gone off */
    unsigned int _PopII_SNIa_count;
    /*! @brief Number of PopIII SN explosions that has gone off */
    unsigned int _PopIII_SN_count;

    /*! @brief Next time a PopII SNII needs to go off */
    double _PopII_SNII_next_time;
    /*! @brief Upper limit for the next mass interval in which a PopII SNII
     *  should go off */
    double _PopII_SNII_interval;
    /*! @brief PopII SNII factor */
    double _PopII_SNII_fac;

    /*! @brief Next time a PopII SNIa needs to go off */
    double _PopII_SNIa_next_time;
    /*! @brief Lower limit for the next time interval in which a PopII SNIa
     *  should go off */
    double _PopII_SNIa_interval;
    /*! @brief PopII SNIa factor */
    double _PopII_SNIa_fac;

    /*! @brief Next time a PopIII SN needs to go off */
    double _PopIII_SN_next_time;
    /*! @brief Upper limit for the next mass interval in which a PopIII SN
     *  should go off */
    double _PopIII_SN_interval;
    /*! @brief PopIII SN factor */
    double _PopIII_SN_fac;

  public:
    /**
     * @brief Set the total number of PopII SNII that needs to go off
     *
     * @param PopII_SNII_number Total number of PopII SNII
     */
    inline void set_PopII_SNII_number(unsigned int PopII_SNII_number) {
        _PopII_SNII_number = PopII_SNII_number;
    }

    /**
     * @brief Get the total number of PopII SNII that needs to go off
     *
     * @return Total number of PopII SNII
     */
    inline unsigned int get_PopII_SNII_number() {
        return _PopII_SNII_number;
    }

    /**
     * @brief Set the total number of PopII SNIa that needs to go off
     *
     * @param PopII_SNIa_number Total number of PopII SNIa
     */
    inline void set_PopII_SNIa_number(unsigned int PopII_SNIa_number) {
        _PopII_SNIa_number = PopII_SNIa_number;
    }

    /**
     * @brief Get the total number of PopII SNIa that needs to go off
     *
     * @return Total number of PopII SNIa
     */
    inline unsigned int get_PopII_SNIa_number() {
        return _PopII_SNIa_number;
    }

    /**
     * @brief Set the total number of PopIII SN that needs to go off
     *
     * @param PopIII_SN_number Total number of PopIII SN
     */
    inline void set_PopIII_SN_number(unsigned int PopIII_SN_number) {
        _PopIII_SN_number = PopIII_SN_number;
    }

    /**
     * @brief Get the total number of PopIII SN that needs to go off
     *
     * @return Total number of PopIII SN
     */
    inline unsigned int get_PopIII_SN_number() {
        return _PopIII_SN_number;
    }

    /**
     * @brief Set the PopII SW factor
     *
     * @param PopII_SW_fac PopII SW factor
     */
    inline void set_PopII_SW_fac(double PopII_SW_fac) {
        _PopII_SW_fac = PopII_SW_fac;
    }

    /**
     * @brief Get the PopII SW factor
     *
     * @return PopII SW factor
     */
    inline double get_PopII_SW_fac() {
        return _PopII_SW_fac;
    }

    /**
     * @brief Set the PopIII SW factor
     *
     * @param PopIII_SW_fac PopIII SW factor
     */
    inline void set_PopIII_SW_fac(double PopIII_SW_fac) {
        _PopIII_SW_fac = PopIII_SW_fac;
    }

    /**
     * @brief Get the PopIII SW factor
     *
     * @return PopIII SW factor
     */
    inline double get_PopIII_SW_fac() {
        return _PopIII_SW_fac;
    }

    /**
     * @brief Set the number of PopII SNII explosions that has gone off
     *
     * @param PopII_SNII_count Number of PopII SNII explosions that has gone off
     */
    inline void set_PopII_SNII_count(unsigned int PopII_SNII_count) {
        _PopII_SNII_count = PopII_SNII_count;
    }

    /**
     * @brief Get the number of PopII SNII explosions that has gone off
     *
     * @return Number of PopII SNII explosions that has gone off
     */
    inline unsigned int get_PopII_SNII_count() {
        return _PopII_SNII_count;
    }

    /**
     * @brief Increase the number of PopII SNII explosions that has gone off
     */
    inline void increase_PopII_SNII_count() {
        _PopII_SNII_count++;
    }

    /**
     * @brief Set the number of PopII SNIa explosions that has gone off
     *
     * @param PopII_SNIa_count Number of PopII SNIa explosions that has gone off
     */
    inline void set_PopII_SNIa_count(unsigned int PopII_SNIa_count) {
        _PopII_SNIa_count = PopII_SNIa_count;
    }

    /**
     * @brief Get the number of PopII SNIa explosions that has gone off
     *
     * @return Number of PopII SNIa explosions that has gone off
     */
    inline unsigned int get_PopII_SNIa_count() {
        return _PopII_SNIa_count;
    }

    /**
     * @brief Increase the number of PopII SNIa explosions that has gone off
     */
    inline void increase_PopII_SNIa_count() {
        _PopII_SNIa_count++;
    }

    /**
     * @brief Set the number of PopIII SN explosions that has gone off
     *
     * @param PopIII_SN_count Number of PopIII SN explosions that has gone off
     */
    inline void set_PopIII_SN_count(unsigned int PopIII_SN_count) {
        _PopIII_SN_count = PopIII_SN_count;
    }

    /**
     * @brief Get the number of PopIII SN explosions that has gone off
     *
     * @return Number of PopIII SN explosions that has gone off
     */
    inline unsigned int get_PopIII_SN_count() {
        return _PopIII_SN_count;
    }

    /**
     * @brief Increase the number of PopIII SN explosions that has gone off
     */
    inline void increase_PopIII_SN_count() {
        _PopIII_SN_count++;
    }

    /**
     * @brief Set the next time a PopII SNII should go off
     *
     * @param PopII_SNII_next_time Next time a PopII SNII should go off
     */
    inline void set_PopII_SNII_next_time(double PopII_SNII_next_time) {
        _PopII_SNII_next_time = PopII_SNII_next_time;
    }

    /**
     * @brief Get the next time a PopII SNII should go off
     *
     * @return Next time a PopII SNII should go off
     */
    inline double get_PopII_SNII_next_time() {
        return _PopII_SNII_next_time;
    }

    /**
     * @brief Set the upper limit for the next mass interval in which a PopII
     * SNII should go off
     *
     * @param PopII_SNII_interval Upper limit for the next mass interval in
     * which a PopII SNII should go off
     */
    inline void set_PopII_SNII_interval(double PopII_SNII_interval) {
        _PopII_SNII_interval = PopII_SNII_interval;
    }

    /**
     * @brief Get the upper limit for the next mass interval in which a PopII
     * SNII should go off
     *
     * @return Upper limit for the next mass interval in which a PopII SNII
     * should go off
     */
    inline double get_PopII_SNII_interval() {
        return _PopII_SNII_interval;
    }

    /**
     * @brief Set the PopII SNII factor
     *
     * @param PopII_SNII_fac PopII SNII factor
     */
    inline void set_PopII_SNII_fac(double PopII_SNII_fac) {
        _PopII_SNII_fac = PopII_SNII_fac;
    }

    /**
     * @brief Get the PopII SNII factor
     *
     * @return PopII SNII factor
     */
    inline double get_PopII_SNII_fac() {
        return _PopII_SNII_fac;
    }

    /**
     * @brief Set the next time a PopII SNIa should go off
     *
     * @param PopII_SNIa_next_time Next time a PopII SNIa should go off
     */
    inline void set_PopII_SNIa_next_time(double PopII_SNIa_next_time) {
        _PopII_SNIa_next_time = PopII_SNIa_next_time;
    }

    /**
     * @brief Get the next time a PopII SNIa should go off
     *
     * @return Next time a PopII SNIa should go off
     */
    inline double get_PopII_SNIa_next_time() {
        return _PopII_SNIa_next_time;
    }

    /**
     * @brief Set the lower limit for the next time interval in which a PopII
     * SNIa should go off
     *
     * @param PopII_SNIa_interval Upper limit for the next mass interval in
     * which a PopII SNIa should go off
     */
    inline void set_PopII_SNIa_interval(double PopII_SNIa_interval) {
        _PopII_SNIa_interval = PopII_SNIa_interval;
    }

    /**
     * @brief Get the lower limit for the next time interval in which a PopII
     * SNIa should go off
     *
     * @return Upper limit for the next mass interval in which a PopII SNIa
     * should go off
     */
    inline double get_PopII_SNIa_interval() {
        return _PopII_SNIa_interval;
    }

    /**
     * @brief Set the PopII SNIa factor
     *
     * @param PopII_SNIa_fac PopII SNIa factor
     */
    inline void set_PopII_SNIa_fac(double PopII_SNIa_fac) {
        _PopII_SNIa_fac = PopII_SNIa_fac;
    }

    /**
     * @brief Get the PopII SNIa factor
     *
     * @return PopII SNIa factor
     */
    inline double get_PopII_SNIa_fac() {
        return _PopII_SNIa_fac;
    }

    /**
     * @brief Set the next time a PopIII SN should go off
     *
     * @param PopIII_SN_next_time Next time a PopIII SN should go off
     */
    inline void set_PopIII_SN_next_time(double PopIII_SN_next_time) {
        _PopIII_SN_next_time = PopIII_SN_next_time;
    }

    /**
     * @brief Get the next time a PopIII SN should go off
     *
     * @return Next time a PopIII SN should go off
     */
    inline double get_PopIII_SN_next_time() {
        return _PopIII_SN_next_time;
    }

    /**
     * @brief Set the upper limit for the next mass interval in which a PopIII
     * SN should go off
     *
     * @param PopIII_SN_interval Upper limit for the next mass interval in
     * which a PopIII SN should go off
     */
    inline void set_PopIII_SN_interval(double PopIII_SN_interval) {
        _PopIII_SN_interval = PopIII_SN_interval;
    }

    /**
     * @brief Get the upper limit for the next mass interval in which a PopIII
     * SN should go off
     *
     * @return Upper limit for the next mass interval in which a PopIII SN
     * should go off
     */
    inline double get_PopIII_SN_interval() {
        return _PopIII_SN_interval;
    }

    /**
     * @brief Set the PopIII SN factor
     *
     * @param PopIII_SN_fac PopIII SN factor
     */
    inline void set_PopIII_SN_fac(double PopIII_SN_fac) {
        _PopIII_SN_fac = PopIII_SN_fac;
    }

    /**
     * @brief Get the PopIII SN factor
     *
     * @return PopIII SN factor
     */
    inline double get_PopIII_SN_fac() {
        return _PopIII_SN_fac;
    }
};

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

    double PopII_interval_mass(double pmass, double mupp);
    double PopII_lifetime(double m);

    double PopII_SNIa_delay1(double t);
    double PopII_SNIa_delay2(double t);
    double PopII_SNIa_delay(double t);

    double PopII_SNIa_delay_cumul(double t);
    double PopII_SNIa_interval_func(double t, double NIa, double vlow);
    double PopII_SNIa_interval_time(double NIa, double tlow);

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

    /**
     * @brief Auxiliary structure used to pass on parameters to
     * PopII_SNIa_interval_func()
     */
    struct PopII_SNIa_interval_func_static_params {
        /*! @brief Object instance that the static wrapper should call */
        DiscreteStellarFeedback* object;
        /*! @brief Number of SNIa */
        double NIa;
        /*! @brief Extra parameter that depends on the lower bound of the
         *  interval */
        double vlow;
    };

    /**
     * @brief Static wrapper around PopII_SNIa_interval_func()
     *
     * @param t Time value
     * @param params PopII_SNIa_interval_func_static_params struct
     * @return Value of PopII_SNIa_interval_func()
     */
    static double PopII_SNIa_interval_func_static(double t, void* params) {
        struct PopII_SNIa_interval_func_static_params* pars =
                (PopII_SNIa_interval_func_static_params*)params;
        DiscreteStellarFeedback* object = pars->object;
        return object->PopII_SNIa_interval_func(t, pars->NIa, pars->vlow);
    }

    double PopIII_IMF(double m);
    double PopIII_mIMF(double m);

    double PopIII_E_SN(double m);
    double PopIII_EIMF(double m);

    double PopIII_IMF_cumul(double m);
    double PopIII_IMF_interval_func(double m, double pmass, double vlow);
    double PopIII_interval_mass(double pmass, double mupp);

    double PopIII_lifetime(double m);

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

    /**
     * @brief Extra parameters for the PopIII IMF auxiliary function
     */
    struct PopIII_IMF_interval_func_static_params {
        /*! @brief Object instance that the static wrapper should use */
        DiscreteStellarFeedback* object;
        /*! @brief Total (initial) mass of the stellar particle */
        double pmass;
        /*! @brief Value of the cumulative PopIII IMF at the upper bound of the
         * mass
         *  interval */
        double vlow;
    };

    /**
     * @brief Static wrapper around PopIII_IMF_interval_func()
     *
     * @param m Mass value
     * @param params PopIII_IMF_interval_func_static_params structure holding
     * extra parameters
     * @return Value of PopIII_IMF_interval_func()
     */
    static double PopIII_IMF_interval_func_static(double m, void* params) {
        struct PopIII_IMF_interval_func_static_params* pars =
                (PopIII_IMF_interval_func_static_params*)params;
        DiscreteStellarFeedback* object = pars->object;
        return object->PopIII_IMF_interval_func(m, pars->pmass, pars->vlow);
    }

    void initialize();

  public:
    DiscreteStellarFeedback();
    ~DiscreteStellarFeedback();

    virtual void do_feedback(StarParticle* star, ParticleVector& particles,
                             double dt);
    virtual StellarFeedbackData* initialize_data(StarParticle* star);

    void dump(RestartFile& rfile);
    DiscreteStellarFeedback(RestartFile& rfile);
};

#endif  // DISCRETESTELLARFEEDBACK_HPP
