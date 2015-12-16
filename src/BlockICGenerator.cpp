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
 * @file BlockICGenerator.cpp
 *
 * @brief Classes used to generate initial conditions for simulations: block
 * IC-generator implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "ICGenerator.hpp"
#include "ICRegion.hpp"
#include "VorTessManager.hpp"
#include "DelCont.hpp"
#include "VorCell.hpp"
#include "Error.hpp"
#include "utilities/DMParticle.hpp"
#include "utilities/GasParticle.hpp"
#include "utilities/ParticleVector.hpp"
#include "utilities/Tree.hpp"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/optional/optional.hpp>
#include <boost/algorithm/string.hpp>
using namespace std;

/**
  * \brief Construct a BlockICGenerator with a given number of cells, ICMode and
  * adiabatic index
  *
  * If ICMode::IC_CART is selected, the sampled number of cells can be lower
  * than the given number, since the total number of cells is a second or third
  * power of an integer number. If the square or cubic root of the given number
  * is not integer, the highest integer smaller than this root is used.
  *
  * If, for example, npart = 1000 in 2D, a cartesian grid of 31x31 cells will be
  * set up. The total number of cells will be 961 and not 1000.
  *
  * @param npart Number of cells to be set up
  * @param mode ICMode specifying if cartesian (IC_CART) or regularized uniform
  * random (IC_RAND) cells should be set up
  * @param seed Random generator seed (for parallel runs, the actual seed is
  * this seed + the MPI rank + 1)
  * @param gamma Adiabatic index of the gas in the initial condition
  */
BlockICGenerator::BlockICGenerator(unsigned int npart, unsigned int mode,
                                   unsigned int seed, double gamma){
    _ngas = npart;
    _ndark = npart;
    _mode = mode;
    _maxdens = 0.;
    _gamma = gamma;
    _periodic = false;
    _has_gas = false;
    _has_dm = false;
    // make sure every process has a different random seed
    srand(MPIGlobal::rank + seed + 1);
}

/**
 * @brief Destructor. Free the IC-regions
 */
BlockICGenerator::~BlockICGenerator(){
    for(unsigned int i = 0; i < _regions.size(); i++){
        delete _regions[i];
    }
    delete _container;
}

/**
  * \brief Read in geometrical blocks and other relevant information from a
  * specified xml file
  *
  * A basic xml file could look like
\verbatim
<box>
    <height>SIMULATION_BOX_HEIGHT</height>
    <width>SIMULATION_BOX_WIDTH</width>
   (<depth>SIMULATION_BOX_DEPTH</depth>)
    <origin>
        <x>ORIGIN_X</x>
        <y>ORIGIN_Y</y>
       (<z>ORIGIN_Z</z>)
    </origin>
    <regions>
        <region>
            <height>BLOCK_HEIGHT</height>
            <width>BLOCK_WIDTH</width>
           (<depth>BLOCK_DEPTH</depth>)
            <exponent>BLOCK_EXPONENT</exponent>
            <origin>
                <x>BLOCK_ORIGIN_X</x>
                <y>BLOCK_ORIGIN_Y</y>
               (<z>BLOCK_ORIGIN_Z</z>)
            </origin>
            <hydro>
                <rho>BLOCK_DENSITY</rho>
                <vx>BLOCK_VELOCITY_X</vx>
                <vy>BLOCK_VELOCITY_Y</vy>
               (<vz>BLOCK_VELOCITY_Z</vz>)
                <p>BLOCK_PRESSURE</p>
            </hydro>
        </region>
        ...
    </regions>
</box>
\endverbatim
  * The brackets denote tags only needed in 3D.
  *
  * @param filename std::string specifying the valid name of an xml file
  */
void BlockICGenerator::read_xml(string filename){
    ifstream file(filename.c_str());
    if(!file){
        cerr << "Cannot open XML-file \"" << filename << "\"!" << endl;
        my_exit();
    }

    boost::property_tree::ptree pt;
    boost::property_tree::xml_parser::read_xml(filename, pt);
    Vec origin;
    origin[0] = pt.get<double>("box.origin.x");
    origin[1] = pt.get<double>("box.origin.y");
#if ndim_==3
    origin[2] = pt.get<double>("box.origin.z");
#endif
    Vec width;
    width[0] = pt.get<double>("box.height");
    width[1] = pt.get<double>("box.width");
#if ndim_==3
    width[2] = pt.get<double>("box.depth");
#endif
    string boundary;
    boundary = pt.get<string>("box.boundary", "reflective");
    boost::algorithm::trim(boundary);
    _periodic = (boundary == "periodic");
    if(_periodic){
        cout << "periodic boundaries" << endl;
    } else {
        cout << "reflective boundaries" << endl;
    }
    cout << "origin: " << origin[0] << "\t" << origin[1];
#if ndim_==3
    cout << "\t" << origin[2];
#endif
    cout << endl;
    cout << "height: " << width[0] << endl;
    cout << "width: " << width[1] << endl;
#if ndim_==3
    cout << "depth: " << width[2] << endl;
#endif
    _container = new RectangularBox(origin, width);
    cout << "size: " << pt.get_child("box.regions").size() << endl;
    boost::property_tree::ptree::const_iterator end =
            pt.get_child("box.regions").end();
    for(boost::property_tree::ptree::const_iterator it =
        pt.get_child("box.regions").begin(); it != end; it++){
        cout << "region:" << endl;
        Vec regorigin;
        regorigin[0] = it->second.get("origin.x", origin[0]);
        regorigin[1] = it->second.get("origin.y", origin[1]);
#if ndim_==3
        regorigin[2] = it->second.get("origin.z", origin[2]);
#endif
        cout << "origin: " << regorigin[0] << "\t" << regorigin[1];
#if ndim_==3
        cout << "\t" << regorigin[2];
#endif
        cout << endl;
        Vec regside;
        regside[0] = it->second.get("height", width[0]);
        cout << "height: " << regside[0] << endl;
        regside[1] = it->second.get("width", width[1]);
        cout << "width: " << regside[1] << endl;
#if ndim_==3
        regside[2] = it->second.get("depth", width[2]);
        cout << "depth: " << regside[2] << endl;
#endif
        double exponent = it->second.get("exponent", 10.);
        cout << "exponent: " << exponent << endl;
        vector<string> hydrofunctions;
        vector<string> dmfunction;
        boost::optional<const boost::property_tree::ptree& > hydro =
                it->second.get_child_optional("hydro");
        if(hydro){
            _has_gas = true;
            cout << "hydro functions found:" << endl;
            hydrofunctions.push_back(it->second.get<string>("hydro.rho"));
            cout << "rho: " << hydrofunctions.back() << endl;
            hydrofunctions.push_back(it->second.get<string>("hydro.vx"));
            cout << "vx: " << hydrofunctions.back() << endl;
            hydrofunctions.push_back(it->second.get<string>("hydro.vy"));
            cout << "vy: " << hydrofunctions.back() << endl;
#if ndim_==3
            hydrofunctions.push_back(it->second.get<string>("hydro.vz"));
            cout << "vz: " << hydrofunctions.back() << endl;
#endif
            hydrofunctions.push_back(it->second.get<string>("hydro.p"));
            cout << "p: " << hydrofunctions.back() << endl;
        }
        boost::optional<const boost::property_tree::ptree& > dm =
                it->second.get_child_optional("dm");
        if(dm){
            _has_dm = true;
            cout << "DM function found:" << endl;
            dmfunction.push_back(it->second.get<string>("dm.rho"));
            cout << "rho: " << dmfunction.back() << endl;
        }
        _regions.push_back(new ICRegion(regorigin, regside, exponent,
                                        hydrofunctions, dmfunction));
    }
}

/**
  * \brief Set the hydrodynamical state vectors for the cells in the given grid
  * based on the geometrical blocks
  *
  * For every block (starting from the first one added - order is important
  * here!), we check whether or not a cell is inside the block. If it is, we set
  * its StateVector to the StateVector for that block.
  *
  * It is possible to nest geometrical blocks, so that a cell is part of
  * multiple blocks. It is important to note that only the last geometrical
  * block in the list that contains this cell will set its StateVector. All
  * previous changes will be overwritten.
  *
  * @param grid A ParticleVector containing cells that lie inside the simulation
  * box
  */
void BlockICGenerator::apply_regions(ParticleVector &grid){
    for(unsigned int i = grid.gassize(); i--;){
        Vec& pos = grid.gas(i)->get_position();
        // order is important here: the last region that contains the
        // coordinates sets its actual hydro
        for(unsigned int j = 0; j < _regions.size(); j++){
            if(_regions[j]->inside(pos)){
                grid.gas(i)->set_W(_regions[j]->get_hydro(pos));
            }
        }
    }
}

/**
  * \brief Generate initial conditions
  *
  * We first set up a ParticleVector using the DelCont previously set up. We
  * then generate the grid cells and optionally regularize them. We then set up
  * the hydrodynamical properties of the cells using the geometrical blocks.
  *
  * @return A ParticleVector suitable as initial condition for a hydrodynamical
  * simulation
  */
ParticleVector BlockICGenerator::generate(){
    ParticleVector particles(false, *((RectangularBox*)_container), _periodic);
    if(_has_gas){
        if(_mode == IC_CART){
            make_cartesian_grid(particles);
        } else {
            make_random_grid(particles);
            relax_grid(particles);
        }
        apply_regions(particles);
    }
    if(_has_dm){
        add_DM(particles);
    }
    particles.finalize();
    particles.get_header().set_periodic(_periodic);
    return particles;
}

/**
  * \brief Set up a cartesian grid of cells
  *
  * @param plist Reference to an empty ParticleVector that will be filled
  */
void BlockICGenerator::make_cartesian_grid(ParticleVector& plist){
    double box[ndim_+ndim_] = {0.};
    _container->get_bounding_box(box);
#if ndim_==3
    unsigned int cartnum[3];
    double C = cbrt((box[4]/box[3])*_ngas);
    double B = (box[4]/box[3])*_ngas/(C*C);
    // we have to use round instead of int, since the latter sometimes rounds
    // down unexpectedly
    cartnum[0] = round(_ngas/B/C);
    cartnum[1] = round(B);
    cartnum[2] = round(C);
    _ngas = cartnum[0]*cartnum[1]*cartnum[2];

    // divide the particles over the processes
    unsigned int numpart = _ngas/MPIGlobal::size +
            ((_ngas%MPIGlobal::size) > 0);
    // make sure the process with the highest rank does not generate too many
    // particles
    unsigned int nother = MPIGlobal::rank*numpart;
    while(nother + numpart > _ngas){
        numpart--;
    }

    // 3D cartesian grid
    cout << "setting up cartesian grid of " << cartnum[0] << "x" << cartnum[1]
         << "x" << cartnum[2] << " cells" << endl;
    double cartwidth[3] = {box[3]/cartnum[0], box[4]/cartnum[1],
                           box[5]/cartnum[2]};
    for(unsigned int i = cartnum[0]; i--;){
        for(unsigned int j = cartnum[1]; j--;){
            for(unsigned int k = cartnum[2]; k--;){
                unsigned int id = i*cartnum[0]*cartnum[0]+j*cartnum[1]+k;
                // we only generate the particles assigned to this process
                if(id < nother || id >= nother+numpart){
                    continue;
                }
                Vec coords(cartwidth[0]*(i+0.5)+box[0],
                           cartwidth[1]*(j+0.5)+box[1],
                           cartwidth[2]*(k+0.5)+box[2]);
                plist.add_gas_particle(new GasParticle(coords));
                plist.gasback()->set_id(id);
            }
        }
    }
#else
    unsigned int cartnum[2];
    double B = sqrt((box[2]/box[3])*_ngas);
    cartnum[1] = round(_ngas/B);
    cartnum[0] = round(B);
    _ngas = cartnum[0]*cartnum[1];

    // divide the particles over the processes
    unsigned int numpart = _ngas/MPIGlobal::size +
            ((_ngas%MPIGlobal::size) > 0);
    // make sure the process with the highest rank does not generate too many
    // particles
    unsigned int nother = MPIGlobal::rank*numpart;
    while(nother + numpart > _ngas){
        numpart--;
    }

    // 2D cartesian grid
    cout << "setting up cartesian grid of " << cartnum[0] << "x" << cartnum[1]
         << " cells" << endl;
    double cartwidth[2] = {box[2]/cartnum[0], box[3]/cartnum[1]};
    for(unsigned int i = cartnum[0]; i--;){
        for(unsigned int j = cartnum[1]; j--;){
            unsigned int id = i*cartnum[0]+j;
            // we only generate the particles assigned to this process
            if(id < nother || id >= nother+numpart){
                continue;
            }
            Vec coords(cartwidth[0]*(i+0.5)+box[0],
                    cartwidth[1]*(j+0.5)+box[1]);
            plist.add_gas_particle(new GasParticle(coords));
            plist.gasback()->set_id(id);
        }
    }
#endif
}

/**
  * \brief Set up a uniform random grid using adaptive random sampling
  *
  * @param plist Reference to an empty ParticleVector that will be filled
  */
void BlockICGenerator::make_random_grid(ParticleVector& plist){
    _rejection_fac.resize(_regions.size());
    _maxdens = 0.;
    for(unsigned int i = _regions.size(); i--;){
        ICRegion* next = NULL;
        if(i < _regions.size()-2){
            next = _regions[i+1];
        }
        _rejection_fac[i] = _regions[i]->get_max_value_hydro(next);
        _maxdens = std::max(_maxdens, _rejection_fac[i]);
    }
    for(unsigned int i = 0; i < _regions.size(); i++){
        _rejection_fac[i] /= _maxdens;
    }
    double box[ndim_+ndim_] = {0.};
    _container->get_bounding_box(box);
    // divide the particles over the processes
    unsigned int numpart = _ngas/MPIGlobal::size +
            ((_ngas%MPIGlobal::size) > 0);
    // make sure the process with the highest rank does not generate too many
    // particles
    unsigned int nother = MPIGlobal::rank*numpart;
    while(nother + numpart > _ngas){
        numpart--;
    }
    for(unsigned int i = 0; i < numpart; i++){
        Vec x;
        for(unsigned int j = ndim_; j--;){
            x[j] = box[ndim_+j]*(double)rand()/((double)RAND_MAX)+box[j];
        }
        while(!regions_accepted(x)){
            for(unsigned int j = ndim_; j--;){
                x[j] = box[ndim_+j]*(double)rand()/((double)RAND_MAX)+box[j];
            }
        }
        plist.add_gas_particle(new GasParticle(x));
        plist.gasback()->set_id(i);
        plist.gasback()->set_starttime(0);
    }
}

/**
  * \brief Determine whether a given Vec is an acceptable coordinate vector for
  * a randomly sampled particle
  *
  * We take into account the geometrical blocks and their density with respect
  * to the maximal density in the simulation.
  *
  * @param p Vec specifying a valid coordinate vector inside the simulation box
  * @return true if p is accepted by the geometrical block that contains p,
  * false if it is rejected
  */
bool BlockICGenerator::regions_accepted(Vec &p){
    unsigned int lastregion = _regions.size()-1;
    while(!_regions[lastregion]->inside(p)){
        lastregion--;
    }
    if(_rejection_fac[lastregion] < 1.){
        if((double)rand()/((double)RAND_MAX) < _rejection_fac[lastregion]){
            return _regions[lastregion]->accept_hydro(p);
        } else {
            return false;
        }
    } else {
        return _regions[lastregion]->accept_hydro(p);
    }
}

/**
 * @brief Same as BlockICGenerator::regions_accepted(), but for dark matter
 * particles
 *
 * @param p Vec specifying a valid coordinate vector inside the simulation box
 * @return True if p is accepted, false if it is rejected
 */
bool BlockICGenerator::regions_accepted_dm(Vec &p){
    unsigned int lastregion = _regions.size()-1;
    while(!_regions[lastregion]->inside(p)){
        lastregion--;
    }
    if(_rejection_fac_dm[lastregion] < 1.){
        if((double)rand()/((double)RAND_MAX) < _rejection_fac_dm[lastregion]){
            return _regions[lastregion]->accept_dm(p);
        } else {
            return false;
        }
    } else {
        return _regions[lastregion]->accept_dm(p);
    }
}

/**
  * \brief Regularize a randomly sampled grid using Lloyd's algorithm
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
  */
void BlockICGenerator::relax_grid(ParticleVector &grid){
    grid.sort();
    VorTessManager voronoi(grid, _periodic);
    double treebox[ndim_+ndim_];
    grid.get_container().get_bounding_box(treebox);

    unsigned int numlloyd = 10;
    for(unsigned int i = numlloyd; i--;){
        cout << "Lloyd iteration " << (numlloyd-i) << endl;
        for(unsigned int j = grid.gassize(); j--;){
            Vec centroid = voronoi.get_centroid(j);
            grid.gas(j)->set_x(centroid[0]);
            grid.gas(j)->set_y(centroid[1]);
#if ndim_==3
            grid.gas(j)->set_z(centroid[2]);
#endif
            grid.gas(j)->set_vorgen(NULL);
        }

        voronoi.reset(&grid.get_container(), _periodic);
        // order is important here!
        // since the update_positions depends on the old ordering of the
        // particles, we have to call it before updating the ordering by
        // sorting
        voronoi.update_positions(_periodic);
        grid.sort();

        voronoi.update(0);

        voronoi.set_hs(0);
    }
}

/**
 * @brief Add dark matter particles to the given ParticleVector based on the
 * requested grid type and possibly weighted by the requested block densities
 *
 * @param plist ParticleVector to fill
 */
void BlockICGenerator::add_DM(ParticleVector &plist){
    _rejection_fac_dm.resize(_regions.size());
    _maxdens_dm = 0.;
    for(unsigned int i = _regions.size(); i--;){
        ICRegion* next = NULL;
        if(i < _regions.size()-2){
            next = _regions[i+1];
        }
        _rejection_fac_dm[i] = _regions[i]->get_max_value_dm(next);
        _maxdens_dm = std::max(_maxdens_dm, _rejection_fac_dm[i]);
    }
    for(unsigned int i = 0; i < _regions.size(); i++){
        _rejection_fac_dm[i] /= _maxdens_dm;
    }
    double box[ndim_+ndim_] = {0.};
    _container->get_bounding_box(box);
    // make sure every process has a different random seed
    srand(MPIGlobal::rank + 42);
    // divide the particles over the processes
    unsigned int numpart = _ndark/MPIGlobal::size +
            ((_ndark%MPIGlobal::size) > 0);
    // make sure the process with the highest rank does not generate too many
    // particles
    unsigned int nother = MPIGlobal::rank*numpart;
    while(nother + numpart > _ndark){
        numpart--;
    }
    for(unsigned int i = 0; i < numpart; i++){
        Vec x;
        for(unsigned int j = ndim_; j--;){
            x[j] = box[ndim_+j]*(double)rand()/((double)RAND_MAX)+box[j];
        }
        while(!regions_accepted_dm(x)){
            for(unsigned int j = ndim_; j--;){
                x[j] = box[ndim_+j]*(double)rand()/((double)RAND_MAX)+box[j];
            }
        }
        plist.add_DM_particle(new DMParticle(x));
        plist.dmback()->set_id(i);
    }
}
