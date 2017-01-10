/*******************************************************************************
 * This file is part of Shadowfax
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *                    Yorick Van Den Bossche (yorick.vandenbossche@ugent.be)
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
 * @file StellarFeedback.hpp
 *
 * @brief ParticleConverter implementation that converts a GasParticle to a
 * StarParticle if star formation criteria are met: header
 *
 * @author Yorick Van Den Bossche (yorick.vandenbossche@ugent.be)
 */
#ifndef STELLARFEEDBACK_HPP
#define STELLARFEEDBACK_HPP

class ParameterFile;
class ParticleVector;
class RestartFile;
class StarParticle;

#define STELLARFEEDBACK_DEFAULT_SWSTART "0. yr"
#define STELLARFEEDBACK_DEFAULT_SWEND "4.3e9 yr"
#define STELLARFEEDBACK_DEFAULT_SNIISTART "3.8e6 yr"
#define STELLARFEEDBACK_DEFAULT_SNIIEND "3.1e7 yr"
#define STELLARFEEDBACK_DEFAULT_SNIASTART "1.543e9 yr"
#define STELLARFEEDBACK_DEFAULT_SNIAEND "1.87e9 yr"

#define STELLARFEEDBACK_DEFAULT_SNIIFE 0.000932719658516
#define STELLARFEEDBACK_DEFAULT_SNIIMG 0.00151412640705
#define STELLARFEEDBACK_DEFAULT_SNIAFE 0.00165100587997
#define STELLARFEEDBACK_DEFAULT_SNIAMG 0.000257789470044
#define STELLARFEEDBACK_DEFAULT_SNIIE "1.e51 erg"
#define STELLARFEEDBACK_DEFAULT_SWE "1.e50 erg"
#define STELLARFEEDBACK_DEFAULT_SNIAE "1.e51 erg"
#define STELLARFEEDBACK_DEFAULT_SNIIM 0.00655147325196
#define STELLARFEEDBACK_DEFAULT_SNIAM 0.191445322565
#define STELLARFEEDBACK_DEFAULT_NSNIAPERSNII 0.15
#define STELLARFEEDBACK_DEFAULT_MASSFACSNII "5.93062793e-33 kg^-1"

/**
 * @brief ParticleConverter implementation that converts a GasParticle to a
 * StarParticle if star formation criteria are met
 */
class StellarFeedback {
  private:
    /*! @brief Start time of stellar wind feedback */
    double _sw_start;
    /*! @brief End time of stellar wind feedback */
    double _sw_end;

    /*! @brief Start time of SNII feedback */
    double _sii_start;
    /*! @brief End time of SNII feedback */
    double _sii_end;

    /*! @brief Start time of SNIa feedback */
    double _sia_start;
    /*! @brief End time of SNIa feedback */
    double _sia_end;

    /*! @brief Fe yield of SNII feedback */
    double _sii_Fe;
    /*! @brief Mg yield of SNII feedback */
    double _sii_Mg;

    /*! @brief Fe yield of SNIa feedback */
    double _sia_Fe;
    /*! @brief Mg yield of SNIa feedback */
    double _sia_Mg;

    /*! @brief Energy emitted by SNII feedback */
    double _sii_E;
    /*! @brief Energy emitted by stellar wind feedback */
    double _sw_E;
    /*! @brief Energy emitted by SNIa feedback */
    double _sia_E;

    /*! @brief Mass emitted by SNIa feedback */
    double _sia_m;
    /*! @brief Mass emitted by SNII feedback */
    double _sii_m;

    /*! @brief Number of SNIa explosions per SNII explosion */
    double _nsnia_per_snii;

    /*! @brief Number of SNII explosions per stellar population mass unit */
    double _massfac_snii;

    /*! @brief Initial search radius for gasparticles close to the star */
    double _sradius;

  public:
    StellarFeedback(ParameterFile* parameters);

    void set_radius(double rad);
    void do_feedback(StarParticle* star, ParticleVector& particles, double dt);

    void dump(RestartFile& rfile);
    StellarFeedback(RestartFile& rfile);
};

#endif  // STELLARFEEDBACK_HPP
