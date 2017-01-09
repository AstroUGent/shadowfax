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
 * @file SpecificICGenerator.cpp
 *
 * @brief Classes used to generate initial conditions for simulations: specific
 * IC-generator implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "SpecificICGenerator.hpp"
#include "DelCont.hpp"           // for DelCont
#include "Error.hpp"             // for my_exit
#include "MPIGlobal.hpp"         // for rank, size
#include "StateVector.hpp"       // for StateVector
#include "Vec.hpp"               // for Vec
#include "VorTessManager.hpp"    // for VorTessManager
#include "io/Unit.hpp"           // for Unit
#include "io/UnitConverter.hpp"  // for UnitConverter
#include "io/UnitSet.hpp"        // for UnitSet
#include "riemann/RiemannSolverFactory.hpp"
#include "utilities/DMParticle.hpp"      // for DMParticle
#include "utilities/GasParticle.hpp"     // for GasParticle
#include "utilities/ParticleVector.hpp"  // for ParticleVector
#include <algorithm>                     // for max
#include <cmath>                         // for sqrt, sin, M_PI, cos, log, etc
#include <cstdlib>                       // for rand, srand, RAND_MAX, NULL
#include <iostream>                      // for operator<<, basic_ostream, etc
#include <string>                        // for string
using namespace std;

/**
 * @brief Constructor
 *
 * @param ngaspart Number of gas particles to generate
 * @param ndmpart Number of dark matter particles to generate
 * @param type ICSpecType specifying for which test initial conditions should be
 * generated
 * @param seed Seed for the random generator
 * @param mode ICMode specifying which type of grid to use
 * @param gamma Adiabatic index of the gas in the generated initial conditions
 */
SpecificICGenerator::SpecificICGenerator(unsigned int ngaspart,
                                         unsigned int ndmpart,
                                         unsigned int type, unsigned int seed,
                                         unsigned int mode, double gamma) {
    _ngaspart = ngaspart;
    _ndmpart = ndmpart;
    _type = type;
    _mode = mode;
    _gamma = gamma;
    _periodic = false;
    // make sure every process has a different random seed
    srand(MPIGlobal::rank + 1 + seed);

// default values
#if ndim_ == 3
    Vec origin(0.5, 0.5, 0.5);
    Vec sides(1., 1., 1.);
#else
    Vec origin(0.5, 0.5);
    Vec sides(1., 1.);
#endif
    if(type == IC_SPEC_GRESHO) {
#if ndim_ == 3
        cerr << "The Gresho vortex is a 2D test problem. Get rid of a dimension"
                " and try again!"
             << endl;
        my_exit();
#endif
        _periodic = true;
    }
    if(type == IC_SPEC_DWARF) {
#if ndim_ == 3
        origin = Vec(0., 0., 0.);
        sides = Vec(30., 30., 30.);
#else
        origin = Vec(0., 0.);
        sides = Vec(30., 30.);
#endif
        _ndmpart = _ngaspart;
    }
    if(type == IC_SPEC_KH) {
#if ndim_ == 3
        origin = Vec(2., 2., 2.);
        sides = Vec(4., 4., 4.);
#else
        origin = Vec(2., 2.);
        sides = Vec(4., 4.);
#endif
        _periodic = true;
    }
    if(type == IC_SPEC_EVRARD) {
#if ndim_ == 3
        origin = Vec(0., 0., 0.);
        sides = Vec(4., 4., 4.);
        _evrardfrac = 0.0001;
#else
        cerr << "The Evrard collapse is a 3D test problem. Buy another "
                "dimension and try again!"
             << endl;
        my_exit();
#endif
    }
    if(type == IC_SPEC_NBODY) {
#if ndim_ == 3
        origin = Vec(0., 0., 0.);
        sides = Vec(40., 40., 40.);
#else
        cerr << "The N-body test is a 3D test problem. Buy another dimension"
                "and try again!"
             << endl;
        my_exit();
#endif
        _ndmpart = _ngaspart;
        _ngaspart = 0;
    }
    _container = new RectangularBox(origin, sides);
}

/**
 * @brief Destructor. Free the simulation box
 */
SpecificICGenerator::~SpecificICGenerator() {
    delete _container;
}

/**
  * @brief Generate initial conditions
  *
  * We first set up a ParticleVector using the DelCont previously set up. We
  * then generate the grid cells and optionally regularize them. We then set up
  * the hydrodynamical properties of the cells using the specific profile
  * specified in _type.
  *
  * @param numlloyd Number of Lloyd iterations used to relax the grid
  * @param conserved_variables Generate conserved variables?
  * @return A ParticleVector suitable as initial condition for a hydrodynamical
  * simulation
  */
ParticleVector SpecificICGenerator::generate(unsigned int numlloyd,
                                             bool conserved_variables) {
    ParticleVector particles(false, *((RectangularBox*)_container));
    if(_ngaspart) {
        if(_mode == IC_CART) {
            make_cartesian_grid(particles);
        } else {
            make_random_grid(particles);
            relax_grid(particles, numlloyd);
        }
        apply_profiles(particles);
        if(conserved_variables) {
            set_conserved_variables(particles);
        }
    }
    if(_ndmpart > 0) {
        add_DM(particles);
    }
    particles.finalize();
    particles.get_header().set_periodic(_periodic);
    return particles;
}

/**
  * @brief Set up a cartesian grid of cells
  *
  * @param plist Reference to an empty ParticleVector that will be filled
  */
void SpecificICGenerator::make_cartesian_grid(ParticleVector& plist) {
    double box[ndim_ + ndim_] = {0.};
    _container->get_bounding_box(box);
#if ndim_ == 3
    unsigned int cartnum[3];
    double C = cbrt((box[4] / box[3]) * _ngaspart);
    double B = (box[4] / box[3]) * _ngaspart / (C * C);
    cartnum[0] = int(_ngaspart / B / C);
    cartnum[1] = int(B);
    cartnum[2] = int(C);
    _ngaspart = cartnum[0] * cartnum[1] * cartnum[2];

    // divide the particles over the processes
    unsigned int numpart =
            _ngaspart / MPIGlobal::size + ((_ngaspart % MPIGlobal::size) > 0);
    // make sure the process with the highest rank does not generate too many
    // particles
    unsigned int nother = MPIGlobal::rank * numpart;
    while(nother + numpart > _ngaspart) {
        numpart--;
    }

    // 3D cartesian grid
    cout << "setting up cartesian grid of " << cartnum[0] << "x" << cartnum[1]
         << "x" << cartnum[2] << " cells" << endl;
    double cartwidth[3] = {box[3] / cartnum[0], box[4] / cartnum[1],
                           box[5] / cartnum[2]};
    for(unsigned int i = cartnum[0]; i--;) {
        for(unsigned int j = cartnum[1]; j--;) {
            for(unsigned int k = cartnum[2]; k--;) {
                unsigned int id =
                        i * cartnum[0] * cartnum[0] + j * cartnum[1] + k;
                // we only generate the particles assigned to this process
                if(id < nother || id >= nother + numpart) {
                    continue;
                }
                Vec coords(cartwidth[0] * (i + 0.5) + box[0],
                           cartwidth[1] * (j + 0.5) + box[1],
                           cartwidth[2] * (k + 0.5) + box[2]);
                plist.add_gas_particle(new GasParticle(coords));
                plist.gasback()->set_id(id);
            }
        }
    }
#else
    unsigned int cartnum[2];
    double B = sqrt((box[2] / box[3]) * _ngaspart);
    cartnum[1] = int(_ngaspart / B);
    cartnum[0] = int(B);
    _ngaspart = cartnum[0] * cartnum[1];

    // divide the particles over the processes
    unsigned int numpart =
            _ngaspart / MPIGlobal::size + ((_ngaspart % MPIGlobal::size) > 0);
    // make sure the process with the highest rank does not generate too many
    // particles
    unsigned int nother = MPIGlobal::rank * numpart;
    while(nother + numpart > _ngaspart) {
        numpart--;
    }

    // 2D cartesian grid
    cout << "setting up cartesian grid of " << cartnum[0] << "x" << cartnum[1]
         << " cells" << endl;
    double cartwidth[2] = {box[2] / cartnum[0], box[3] / cartnum[1]};
    for(unsigned int i = cartnum[0]; i--;) {
        for(unsigned int j = cartnum[1]; j--;) {
            unsigned int id = i * cartnum[0] + j;
            // we only generate the particles assigned to this process
            if(id < nother || id >= nother + numpart) {
                continue;
            }
            Vec coords(cartwidth[0] * (i + 0.5) + box[0],
                       cartwidth[1] * (j + 0.5) + box[1]);
            plist.add_gas_particle(new GasParticle(coords));
            plist.gasback()->set_id(id);
        }
    }
#endif
}

/**
  * @brief Set up a uniform random grid using uniform random sampling
  *
  * @param plist Reference to an empty ParticleVector that will be filled
  */
void SpecificICGenerator::make_random_grid(ParticleVector& plist) {
    double box[ndim_ + ndim_] = {0.};
    _container->get_bounding_box(box);
    //    ofstream ofile("particles.dat");
    unsigned int numinside = 0;
    // make sure every process has a different random seed
    srand(MPIGlobal::rank + 1);
    // divide the particles over the processes
    unsigned int numpart =
            _ngaspart / MPIGlobal::size + ((_ngaspart % MPIGlobal::size) > 0);
    // make sure the process with the highest rank does not generate too many
    // particles
    unsigned int nother = MPIGlobal::rank * numpart;
    while(nother + numpart > _ngaspart) {
        numpart--;
    }
    for(unsigned int i = 0; i < numpart; i++) {
        Vec x;
        // we first sample inside the sphere with radius 1, until the number of
        // points sampled is equal to the fraction of the total mass of the
        // system inside the sphere minus the mass that is uniformly added
        // afterwards
        if(_type == IC_SPEC_EVRARD &&
           numinside < (1. - 2. * _evrardfrac / 3.) /
                               (32. * _evrardfrac / M_PI + 1. -
                                2. * _evrardfrac / 3.) *
                               _ngaspart) {
            for(unsigned int j = ndim_; j--;) {
                x[j] = (double)rand() / ((double)RAND_MAX);
            }
#if ndim_ == 3
            double r = pow(x[0], 1.5);
            double theta = 2. * M_PI * x[1];
            double cosphi = 2. * x[2] - 1.;
            x.set(r * sqrt(1. - cosphi * cosphi) * cos(theta),
                  r * sqrt(1. - cosphi * cosphi) * sin(theta), r * cosphi);
#endif
            numinside++;
        } else {
            for(unsigned int j = ndim_; j--;) {
                x[j] = box[ndim_ + j] * (double)rand() / ((double)RAND_MAX) +
                       box[j];
            }
        }
        plist.add_gas_particle(new GasParticle(x));
        plist.gasback()->set_id(i);
    }
}

/**
  * @brief Regularize a randomly sampled grid using Lloyd's algorithm
  *
  * The VorTess for the grid is calculated and the positions of the cell
  * generators are reset to the centroids of their respective VorCell. This
  * procedure is repeated until the resulting grid contains almost equal sized
  * and regularly spaced cells.
  *
  * At the moment, no specific criterion is used to determine when the grid is
  * sufficiently relaxed. Instead, we just iterate for a fixed number of times.
  *
  * @param grid A reference to a non-empty ParticleVector containing a randomly
  * sampled, noisy grid
  * @param numlloyd Number of Lloyd iterations
  */
void SpecificICGenerator::relax_grid(ParticleVector& grid,
                                     unsigned int numlloyd) {
    grid.sort();
    VorTessManager voronoi(grid, _periodic);
    double treebox[ndim_ + ndim_];
    grid.get_container().get_bounding_box(treebox);

    for(unsigned int i = numlloyd; i--;) {
        cout << "Lloyd iteration " << (numlloyd - i) << endl;
        for(unsigned int j = grid.gassize(); j--;) {
            Vec centroid = voronoi.get_centroid(j);
            grid.gas(j)->set_x(centroid[0]);
            grid.gas(j)->set_y(centroid[1]);
#if ndim_ == 3
            grid.gas(j)->set_z(centroid[2]);
#endif
            grid.gas(j)->set_vorgen(NULL);
        }

        // order is important here!
        // since the update_positions depends on the old ordering of the
        // particles, we have to call it before updating the ordering by
        // sorting
        voronoi.update_positions(_periodic);
        grid.sort();

        voronoi.reset(&grid.get_container(), _periodic);

        voronoi.update(0);

        voronoi.set_hs(0);
    }

    // save the cell volumes; we need them for some setups
    for(unsigned int i = 0; i < grid.gassize(); i++) {
        grid.gas(i)->set_mass(voronoi.get_volume(i));
    }
}

/**
 * @brief Apply the requested profiles to the hydrodynamical quantities in the
 * given ParticleVector
 *
 * This method calls the appropriate hard-coded method corresponding to the
 * ICSpecType of this generator.
 *
 * @param grid ParticleVector containing particles without hydrodynamical
 * properties
 */
void SpecificICGenerator::apply_profiles(ParticleVector& grid) {
    if(_type == IC_SPEC_GRESHO) {
        apply_profile_gresho(grid);
    }
    if(_type == IC_SPEC_SEDOV_TAYLOR) {
        apply_profile_sedov_taylor(grid);
    }
    if(_type == IC_SPEC_DWARF) {
        apply_profile_dwarf(grid);
    }
    if(_type == IC_SPEC_KH) {
        apply_profile_kh(grid);
    }
    if(_type == IC_SPEC_EVRARD) {
        apply_profile_evrard(grid);
    }
}

/**
 * @brief Set up the Gresho vortex test
 *
 * @param grid ParticleVector with particles that should be initialized
 */
void SpecificICGenerator::apply_profile_gresho(ParticleVector& grid) {
    for(unsigned int i = grid.gassize(); i--;) {
        double r2 = (grid.gas(i)->x() - 0.5) * (grid.gas(i)->x() - 0.5) +
                    (grid.gas(i)->y() - 0.5) * (grid.gas(i)->y() - 0.5);
        double r = sqrt(r2);
        StateVector W;
        W.set_rho(1.);
        if(r < 0.2) {
            W[1] = -5. * (grid.gas(i)->y() - 0.5);
            W[2] = 5. * (grid.gas(i)->x() - 0.5);
            W[3] = 5. + 12.5 * r2;
        } else {
            if(r < 0.4) {
                W[1] = -(2. - 5. * r) * (grid.gas(i)->y() - 0.5) / r;
                W[2] = (2. - 5. * r) * (grid.gas(i)->x() - 0.5) / r;
                W[3] = 9. + 12.5 * r2 - 20. * r + 4. * log(r / 0.2);
            } else {
                W[1] = 0.;
                W[2] = 0.;
                W[3] = 3. + 4. * log(2);
            }
        }
        grid.gas(i)->set_W(W);
        grid.gas(i)->set_v(W[1], W[2], 0.);
    }
}

/**
 * @brief Set up the Sedov-Taylor blastwave test
 *
 * @param grid ParticleVector with particles to initialize
 */
void SpecificICGenerator::apply_profile_sedov_taylor(ParticleVector& grid) {
    unsigned int central = 0;
    double small_r = 1.;
    for(unsigned int i = grid.gassize(); i--;) {
#if ndim_ == 3
        double r2 = (grid.gas(i)->x() - 0.5) * (grid.gas(i)->x() - 0.5) +
                    (grid.gas(i)->y() - 0.5) * (grid.gas(i)->y() - 0.5) +
                    (grid.gas(i)->z() - 0.5) * (grid.gas(i)->z() - 0.5);
#else
        double r2 = (grid.gas(i)->x() - 0.5) * (grid.gas(i)->x() - 0.5) +
                    (grid.gas(i)->y() - 0.5) * (grid.gas(i)->y() - 0.5);
#endif
        if(r2 < small_r) {
            small_r = r2;
            central = i;
        }
        StateVector W;
        W.set_rho(1.);
        W.set_p(1.e-6);
        grid.gas(i)->set_W(W);
        grid.gas(i)->set_v(0., 0., 0.);
    }

    double central_p;
    double E = 1.;
    if(_mode == IC_CART) {
        central_p = E * (_gamma - 1.) * _ngaspart;
    } else {
        central_p = E * (_gamma - 1.) / grid.gas(central)->get_mass();
    }
    StateVector Wcentral;
    Wcentral.set_rho(1.);
    Wcentral.set_p(central_p);
    grid.gas(central)->set_W(Wcentral);
}

/**
 * @brief Set up a dwarf galaxy halo
 *
 * @param grid ParticleVector with particles to initialize
 */
void SpecificICGenerator::apply_profile_dwarf(ParticleVector& grid) {
    for(unsigned int i = grid.gassize(); i--;) {
        StateVector W;
        // values estimated from Gadget simulation values
        // assuming a gravitational constant set to 1.
        // the length scale fixes the length unit to kpc
        // the density fixes the mass unit to solar mass
        // the gravitational constant then fixes the time unit
        Unit CGS_length("length", "cm", 0.01);
        Unit CGS_mass("mass", "g", 0.001);
        Unit CGS_time("time", "s", 1.);
        Unit CGS_temperature("temperature", "K", 1.);
        UnitSet CGS_units(CGS_length, CGS_mass, CGS_time, CGS_temperature);
        Unit astro_length("length", "kpc", 3.08567785e19);
        Unit astro_mass("mass", "Msol", 1.9891e30);
        Unit astro_time("time", "Gyr", 3.1556926e16);
        Unit astro_temperature("temperature", "K", 1.);
        UnitSet astro_units(astro_length, astro_mass, astro_time,
                            astro_temperature);

        UnitConverter density_converter(CGS_units.get_density_unit(),
                                        astro_units.get_density_unit());
        UnitConverter pressure_converter(CGS_units.get_pressure_unit(),
                                         astro_units.get_pressure_unit());

        W.set_rho(density_converter.convert(1e-20));
        W.set_p(pressure_converter.convert(1.e-15));

        grid.gas(i)->set_W(W);
    }
}

/**
 * @brief Set up the Kelvin-Helmholtz test
 *
 * This is the Kelvin-Helmholtz test as described by Hendrix et al. (2014)
 *
 * @param grid ParticleVector with particles to initialize
 */
void SpecificICGenerator::apply_profile_kh(ParticleVector& grid) {
    double M = 1.;
    double d = 0.06341;
    double kx = 2. * M_PI;
    double v0 = 0.5 * sqrt(5. / 3.);
    for(unsigned int i = grid.gassize(); i--;) {
        StateVector W;
        Vec p = grid.gas(i)->get_position();
        W.set_rho(1.);
        W.set_p(1.);
        double y = p.y();
        if(y > 2.) {
            y = 4. - y;
        }
        W.set_vy(0.001 * v0 * sin(kx * p.x()) *
                 exp(-(y - M) * (y - M) / 32. / d / d));
        if(y > M + d) {
            W.set_vx(-v0);
        } else {
            if(y < M - d) {
                W.set_vx(v0);
            } else {
                W.set_vx((M - y) * v0 / d);
            }
        }
        grid.gas(i)->set_W(W);
    }
}

/**
 * @brief Set up the Evrard collapse test
 *
 * Following Springel (2010)
 *
 * @param grid ParticleVector with particles to initialize
 */
void SpecificICGenerator::apply_profile_evrard(ParticleVector& grid) {
#if ndim_ == 3
    double M = 1.;
    double R = 1.;
    double fac = 0.5 * M / M_PI / R / R;
    for(unsigned int i = grid.gassize(); i--;) {
        StateVector W(_evrardfrac * fac / R, 0., 0., 0.,
                      _evrardfrac * 2. / 3. * fac / R * 0.05);
        Vec p = grid.gas(i)->get_position();
        double r = p.norm();
        if(r <= R) {
            W.set_rho(fac / r);
            W.set_vx(0.);
            W.set_vy(0.);
            W.set_vz(0.);
            W.set_p(2. / 3. * W.rho() * 0.05);
        }
        grid.gas(i)->set_W(W);
    }
#else
    cerr << "The Evrard collapse is a 3D test problem. Buy another dimension "
            "and try again!"
         << endl;
    my_exit();
#endif
}

/**
 * @brief Set the conserved variables of the given particles
 *
 * @param grid Particles to act upon
 */
void SpecificICGenerator::set_conserved_variables(ParticleVector& grid) {
    RiemannSolver* solver =
            RiemannSolverFactory::generate("Exact", _gamma, 0., 0.);
    grid.sort();
    VorTessManager voronoi(grid, _periodic);
    for(unsigned int j = grid.gassize(); j--;) {
        double volume = voronoi.get_volume(j);
        StateVector Q = solver->get_Q(volume, grid.gas(j)->get_Wvec());
        grid.gas(j)->set_Q(Q);
    }
    cout << "Set conserved variables" << endl;

    delete solver;
}

/**
 * @brief Add dark matter particles to the given ParticleVector
 *
 * There is currently only one implemented profile: a Plummer model. However,
 * if more models are written, this function should call the appropriate
 * model method.
 *
 * @param grid ParticleVector to fill
 */
void SpecificICGenerator::add_DM(ParticleVector& grid) {
    add_DM_plummer(grid);
}

/**
 * @brief Set up a Plummer dark matter profile
 *
 * @param grid ParticleVector to fill
 */
void SpecificICGenerator::add_DM_plummer(ParticleVector& grid) {
    double M = 1000.;
    double a = 1.;
    double rcut = 15.;
    double Mmax = a / rcut;
    Mmax *= Mmax;
    Mmax = sqrt(1. + Mmax);
    Mmax *= Mmax * Mmax;
    Mmax = M / Mmax;
    for(unsigned int i = 0; i < _ndmpart; i++) {
        double rad = cbrt(M / (Mmax * rand_double()));
        rad *= rad;
        rad = a / sqrt(rad - 1.);
        // phi in [0,2*pi]
        double phi = 2. * M_PI * rand_double();
        double cosp = cos(phi);
        double sinp = sin(phi);
#if ndim_ == 3
        // cost in [-1,+1]  -> theta in [0,pi]
        double cost = 2. * rand_double() - 1.;
        double sint = sqrt(max(0.0, 1. - cost * cost));
        Vec position(rad * sint * cosp, rad * sint * sinp, rad * cost);
#else
        Vec position(rad * cosp, rad * sinp);
#endif

        // the local escape velocity
        double v_e = sqrt(2. * M / sqrt(a * a + rad * rad));
        double q = rand_double();
        double gq = gplummer(q);
        double y = rand_double();
        while(y > gq) {
            q = rand_double();
            gq = gplummer(q);
            y = rand_double();
        }
        double vrad = v_e * q;
        // phi in [0,2*pi]
        phi = 2. * M_PI * rand_double();
        cosp = cos(phi);
        sinp = sin(phi);
#if ndim_ == 3
        // cost in [-1,+1]  -> theta in [0,pi]
        cost = 2. * rand_double() - 1.;
        sint = sqrt(max(0.0, 1. - cost * cost));
        Vec velocity(vrad * sint * cosp, vrad * sint * sinp, vrad * cost);
#else
        Vec velocity(vrad * cosp, vrad * sinp);
#endif

        float mass = (float)(M / ((double)_ndmpart));
        grid.add_DM_particle(new DMParticle(position));
        // accelerating from 0 to velocity is the same as setting the velocity
        // (but this function does not exist, so yeah...)
        grid.dm(i)->accelerate(velocity);
        grid.dm(i)->set_mass(mass);
        grid.dm(i)->set_id(i);
    }
}

/**
 * @brief Return a uniform random double precision floating point value in the
 * range [0,1[
 *
 * @return A double precision floating point value in the range [0,1[
 */
double SpecificICGenerator::rand_double() {
    return ((double)rand()) / ((double)RAND_MAX);
}

/**
 * @brief Plummer g-function
 *
 * @param q q-parameter
 * @return Value of the Plummer g-function
 */
double SpecificICGenerator::gplummer(double q) {
    double gfac = sqrt(9. / 7.);
    gfac *= gfac * gfac * gfac * gfac * gfac * gfac;
    gfac *= 4.5;
    double q2 = q * q;
    double g = sqrt(1. - q2);
    double g2 = g * g;
    double g3 = g2 * g;
    double g6 = g3 * g3;
    g *= g6;
    return gfac * g * q2;
}
