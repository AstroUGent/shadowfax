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
 * @file BlockICGenerator.hpp
 *
 * @brief Classes used to generate initial conditions for simulations based on
 * the RAMSES block structure specification: header
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef BLOCKICGENERATOR_HPP
#define BLOCKICGENERATOR_HPP

#include "ICGenerator.hpp"  // for ICGenerator, etc
#include <string>           // for string
#include <vector>           // for vector

class DelCont;
class ICRegion;
class ParticleVector;
class Vec;

/**
  * \brief Generate initial conditions from a list of geometrical blocks with
  * hydrodynamical properties
  *
  * The simulation domain is composed of geometrical blocks (cubes, sphere...),
  * each having different values for the hydrodynamical quantities. This way of
  * generating initial conditions is similar to how the initial conditions are
  * set up in the AMR code RAMSES (Teyssier, 2001)
  */
class BlockICGenerator : public ICGenerator {
  private:
    /*! \brief DelCont specifying the entire simulation box */
    DelCont* _container;

    /*! \brief Number of gas particles to generate */
    unsigned int _ngas;

    /*! \brief NUmber of dark matter particles to generate */
    unsigned int _ndark;

    /*! \brief ICMode specifying if a cartesian or random grid should be set
     *  up */
    unsigned int _mode;

    /**
      * \brief Exponents of the geometrical blocks constituting the simulation
      * box
      *
      * The exponent is used to determine the metric for measuring distances in
      * the block.
      * The length of the vector \f$\vec{x} = (x_1, x_2)\f$ under the metric
      * with exponent \f$e\f$ is given by \f$(x_1^e, x_2^e)^{\frac{1}{e}}\f$. An
      * exponent 1 corresponds to diamond shape, 2 to a sphere or circle and
      * exponents larger than 10 give a cube or square. This method is copied
      * from RAMSES (Teyssier, 2001).
      */
    std::vector<ICRegion*> _regions;

    /**
      * \brief Rejection factors for the geometrical blocks constituting the
      * simulation box
      *
      * We try to achieve a constant mass in all grid cells by using more cells
      * in denser regions. We therefore calculate the maximal density in the
      * entire domain and calculate for every block the ratio of its density to
      * this maximal density. If the ratio is high, more sampled particles are
      * accepted during the random sampling.
      */
    std::vector<double> _rejection_fac;

    /*! \brief Rejection factors for the dark matter */
    std::vector<double> _rejection_fac_dm;

    /*! \brief The maximal density of the entire simulation, used for adaptive
     *  sampling */
    double _maxdens;

    /*! \brief Maximal density for the dark matter */
    double _maxdens_dm;

    /*! \brief The desired adiabatic index for the simulation */
    double _gamma;

    /*! \brief Flag indicating if the simulation box is periodic (true) or
     *  reflective (false) */
    bool _periodic;

    /*! \brief Flag indicating if the initial condition has gas */
    bool _has_gas;

    /*! \brief Flag indicating if the initial condition has dark matter */
    bool _has_dm;

    void make_cartesian_grid(ParticleVector& plist);
    void make_random_grid(ParticleVector& plist);
    void relax_grid(ParticleVector& grid);
    void apply_regions(ParticleVector& grid);
    bool regions_accepted(Vec& p);
    bool regions_accepted_dm(Vec& p);

    void add_DM(ParticleVector& plist);

  public:
    BlockICGenerator(unsigned int npart, unsigned int mode = IC_RAND,
                     unsigned int seed = 42, double gamma = 1.66667);
    virtual ~BlockICGenerator();

    void read_xml(std::string filename);

    ParticleVector generate(bool conserved_variables = false);
};

#endif  // BLOCKICGENERATOR_HPP
