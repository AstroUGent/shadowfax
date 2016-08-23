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
 * @file VorTessManager.cpp
 *
 * @brief General interface for Voronoi tesselations: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "VorTessManager.hpp"
#include "AdaptiveMeshException.hpp"
#include "DelCont.hpp"  // for DelCont
#include "FixedGrid.hpp"
#include "ProgramLog.hpp"   // for LOGS
#include "RestartFile.hpp"  // for RestartFile
#include "TimeLine.hpp"
#include "VorCell.hpp"                   // for VorCell
#include "VorFace.hpp"                   // for VorFace
#include "VorTess.hpp"                   // for VorTess
#include "utilities/GasParticle.hpp"     // for GasParticle
#include "utilities/ParticleVector.hpp"  // for ParticleVector
#include "utilities/Timer.hpp"           // for Timer
#include <iostream>                      // for operator<<, ostream, etc
#include <sstream>
#include <vector>  // for vector
#if ndim_ == 2
#include "AdaptiveFace2d.hpp"
#include "AdaptiveVorTess2d.hpp"
#endif
#if ndim_ == 3
#include "AdaptiveFace3d.hpp"
#include "AdaptiveVorTess3d.hpp"
#endif
using namespace std;

/**
 * @brief Constructor
 *
 * @param particles ParticleVector for which the tesselation is constructed
 * @param periodic Flag indicating if the simulation box is periodic (true) or
 * reflective (false)
 * @param tolerance Tolerance value used to distinguish between standard
 * floating point arithmetics and arbitrary precision arithmetics in geometry
 * tests
 * @param evolve Flag indicating if the grid should be evolved rather than
 * reconstructed for every timestep (currently ignored; grid is always
 * reconstructed)
 */
VorTessManager::VorTessManager(ParticleVector& particles, bool periodic,
                               double tolerance, bool evolve)
        : _particles(particles), _periodic(periodic), _tolerance(tolerance) {
    _gridtimer.start();
    _tesselation = new VorTess(&particles.get_container(), particles.gassize(),
                               periodic, tolerance);
    LOGS("Start adding points to VorTess");
    for(unsigned int i = 0; i < _particles.gassize(); i++) {
        _tesselation->add_point(_particles.gas(i),
                                _particles.gas(i)->get_local_id());
        _particles.gas(i)->reset_copies();
        _particles.gas(i)->reset_export();
    }
    LOGS("Done adding points.");
    _tesselation->complete(_particles.get_tree());
    _tesselation->construct();
    LOGS("VorTess ready");

//    stringstream filename;
//    filename << "mesh_" << (MPIGlobal::rank) << ".dat";
//    ofstream ofile(filename.str());
//    _tesselation->print_tesselation_gnuplot(ofile);
//    ofile.close();
//    MyMPI_Barrier();
//    MPI_Barrier(MPI_COMM_WORLD);
//    exit(1);
//    ofstream algfile("tess3d.dat");
//    _tesselation->print_tesselation_fastvor(algfile, periodic);
//    algfile.close();
//    ofstream checkfile("tess3d_gnuplot.dat");
//    _tesselation->print_tesselation_gnuplot(checkfile);
//    checkfile.close();
//    exit(1);
#if ndim_ == 3
    // make sure the h's are set, since they are used to subdivide the steps
    // in the evolution algorithm
    for(unsigned int i = 0; i < particles.gassize(); i++) {
        particles.gas(i)->set_h(_particles.gas(i)->get_cell()->get_h());
    }
    _tess2 = new AdaptiveVorTess3d(
            _tesselation, particles.get_container().get_cuboid(), periodic);

    evolve = false;
    if(evolve) {
        for(unsigned int loop = 0; loop < 100; loop++) {
            cout << "Loop " << loop << endl;
            for(unsigned int i = 0; i < particles.gassize(); i++) {
                //                Vec centroid =
                //                        _tess2->get_centroid(particles.gas(i)->get_local_id());
                Vec centroid = particles.gas(i)->get_position();
                centroid[0] +=
                        0.00725 * (((double)rand()) / ((double)RAND_MAX) - 1.);
                centroid[1] +=
                        0.00725 * (((double)rand()) / ((double)RAND_MAX) - 1.);
                centroid[2] +=
                        0.00725 * (((double)rand()) / ((double)RAND_MAX) - 1.);
                particles.gas(i)->set_x(centroid.x());
                particles.gas(i)->set_y(centroid.y());
                particles.gas(i)->set_z(centroid.z());
            }
            _tess2->update_positions(particles);
            _tess2->update(particles, 0);
            _tess2->check_mesh(cerr);
            cout << "Done." << endl;
        }
        my_exit();
    }
//    ofstream sfile("start.dat");
//    _tess2->print_tesselation_gnuplot(sfile);
//    sfile.close();
//    exit(1);
//    unsigned long numrand = 0;
//    unsigned int startloop = 0;
//    if(true){
//        _tess2 = new AdaptiveVorTess3d(_tesselation,
//    particles.get_container().get_cuboid(), periodic);
//    } else {
//        ifstream rfile("restart.dat");
//        rfile.read(reinterpret_cast<char*>(&numrand), sizeof(unsigned long));
//        double randsum = 0.;
//        for(unsigned long i = 0; i < numrand; i++){
//            randsum += ((double)rand())/((double)RAND_MAX);
//        }
//        cerr << "Randsum: " << randsum << endl;
//        rfile.read(reinterpret_cast<char*>(&startloop), sizeof(unsigned int));
//        startloop += 1;
//        ParticleVector partdump(rfile,
//    (RectangularBox&) particles.get_container(), periodic);
//        for(unsigned int i = 0; i < partdump.gassize(); i++){
//            Vec newpos = partdump.gas(i)->get_position();
//            _particles.gas(i)->set_x(newpos.x());
//            _particles.gas(i)->set_y(newpos.y());
//            _particles.gas(i)->set_z(newpos.z());
//        }
//        _tess2 = new AdaptiveVorTess3d(rfile, _particles.get_vector());
//    }
//    set_hs(0);
//    for(unsigned int loop = startloop; loop < 1000; loop++){
//        cout << "Starting loop " << loop << endl;
//        for(unsigned int i = 0; i < _particles.gassize(); i++){
////            Vec pos = _tess2->get_centroid(_particles.gas(i)->id());
//            Vec pos = _particles.gas(i)->get_position();
//            pos[0] += 0.001*(((double)rand())/((double)RAND_MAX)-1.);
//            pos[1] += 0.001*(((double)rand())/((double)RAND_MAX)-1.);
//            pos[2] += 0.001*(((double)rand())/((double)RAND_MAX)-1.);
//            numrand += 3;
//            _particles.gas(i)->set_x(pos.x());
//            _particles.gas(i)->set_y(pos.y());
//            _particles.gas(i)->set_z(pos.z());
//        }
//        _tess2->update_positions(_particles.get_vector());
//        _tess2->update(_particles.get_vector(), 0);
//        set_hs(0);
////        ofstream efile("errors_mesh.dat");
////        cout << "Checking mesh..." << endl;
////        _tess2->check_mesh(efile);
////        cout << "Done." << endl;
//        if(!(loop%50)){
//            ofstream rfile("restart.dat");
//            rfile.write(reinterpret_cast<char*>(&numrand),
//    sizeof(unsigned long));
//            rfile.write(reinterpret_cast<char*>(&loop), sizeof(unsigned int));
//            particles.dump(rfile);
//            _tess2->dump(rfile, particles.get_vector());
//        }
//        cout << "Finished loop " << loop << endl;
////        stringstream fname;
////        fname << "tess_" << loop << ".dat";
////        ofstream ofile(fname.str().c_str());
////        _tess2->print_tesselation_gnuplot(ofile);
//    }
//    exit(1);
#else
    //    stringstream orfilename;
    //    orfilename << "refstart_" << MPIGlobal::rank << ".dat";
    //    ofstream reffile(orfilename.str());
    //    _tesselation->print_tesselation_gnuplot(reffile);
    //    reffile.close();
    //    MyMPI_Barrier();
    // note: if we ever implement mesh refinement or if for some reason the
    // maximal particle ID is larger than the number of particles, this will
    // not work...
    unsigned long maxid = particles.get_header().npart();
    //    cerr << "maxid: " << maxid << endl;
    _tess2 = new AdaptiveVorTess2d(_tesselation,
                                   particles.get_container().get_cuboid(),
                                   maxid, particles, periodic);
    delete _tesselation;
    return;
    cerr << "Build tesselation" << endl;

    evolve = true;
    if(evolve) {
        delete _tesselation;
        stringstream filename;
        filename << "start_" << MPIGlobal::rank << ".dat";
        ofstream sfile(filename.str());
        _tess2->print_tesselation_gnuplot(sfile);
        sfile.close();
        MyMPI_Barrier();
        //    _tess2->check_mesh(cerr);
        ofstream timefile("times.txt");
        timefile << "0" << endl;
        for(unsigned int loop = 0; loop < 100; loop++) {
            cout << "Starting loop " << loop << endl;
            //            cerr << "Starting loop " << loop << endl;
            timefile << loop << endl;
            cout << "Moving particles..." << endl;
            for(unsigned int i = 0; i < _particles.gassize(); i++) {
                //                unsigned int j = 0;
                //                while(_particles.gas(j)->id() != i){
                //                    j++;
                //                }
                unsigned int j = i;
                Vec pos = _particles.gas(j)->get_position();
                double dpos[2];
                dpos[0] = 0.0014 * (((double)rand()) / ((double)RAND_MAX) - 1.);
                dpos[1] = 0.0014 * (((double)rand()) / ((double)RAND_MAX) - 1.);
                pos[0] += dpos[0];
                pos[1] += dpos[1];
                //                Vec pos =
                //                        _tess2->get_centroid(_particles.gas(j)->get_local_id());

                //                cout << pos.x() << "\t" <<
                // _particles.gas(j)->x() << endl;
                //                cout << pos.y() << "\t" <<
                // _particles.gas(j)->y() << endl;
                _particles.gas(j)->set_x(pos.x());
                _particles.gas(j)->set_y(pos.y());
            }
            cout << "Done." << endl;

            cout << "Updating positions..." << endl;
            _tess2->update_positions(_particles);
            cout << "Done." << endl;
            // to compensate for the particle communications invoked by sort,
            // we should
            //  (a) create new cells for incoming particles
            //  (b) replace outgoing particle cells by orphans
            // normally, new cells will be old orphans and vice versa
            // but even in that case we need new orphans...
            //
            // so, two options:
            //  - either we sort the cells as well, as far as this is possible
            //  - or we somehow gather information about the movements during
            //    the sort and use this information to communicate cells
            //    afterwards
            //
            // in principle, we only need to know
            //  - which cells leave, they are replaced by orphans and their
            //    orphan ngbs are discarded (but probably not all of them - we
            //    better process everything and then discard dangling orphans)
            //  - which cells come in, together with their ngbs (which become
            //    orphans - this probably requires a second step too)
            //            vector<unsigned int> outgoing;
            //            _particles.sort(outgoing);
            cout << "Sorting..." << endl;
            _particles.sort();
            cout << "Done." << endl;
            //            cerr << MPIGlobal::rank << ": " << outgoing.size()
            //            << " outgoing particles:" << endl;
            //            for(unsigned int i = 0; i < outgoing.size(); i++){
            //                cerr << MPIGlobal::rank << ": " << outgoing[i] <<
            //                endl;
            //            }
            cout << "Saving tesselation before update..." << endl;
            stringstream fname;
            fname << "tess_" << MPIGlobal::rank << "_" << loop << "_new.dat";
            ofstream ofile(fname.str().c_str());
            _tess2->print_tesselation_gnuplot(ofile);
            ofile.close();
            MyMPI_Barrier();
            cout << "Done." << endl;
            cout << "Updating..." << endl;
            _tess2->update(_particles, 0);
            MyMPI_Barrier();
            cout << "Done." << endl;
            cout << "Saving tesselation after update..." << endl;
            //            stringstream conname;
            //            conname << "conn_" << MPIGlobal::rank << "_" << loop
            //            << "_new.dat";
            //            ofstream confile(conname.str());
            //            dump_connectivity(confile);
            //            stringstream oldname;
            //            oldname << "conn_" << MPIGlobal::rank << "_" << loop
            //            << "_old.dat";
            //            ofstream oldfile(oldname.str());
            //            VorTessManager oldgrid(_particles, periodic, 1.e-9,
            //            false);
            //            oldgrid.dump_connectivity(oldfile);
            stringstream oldfname;
            oldfname << "tess_" << MPIGlobal::rank << "_" << loop << "_old.dat";
            ofstream oldffile(oldfname.str());
            _tess2->print_tesselation_gnuplot(oldffile);
            MyMPI_Barrier();
            cout << "Done." << endl;
            cout << "Checking mesh validity..." << endl;
            //            _tess2->check_mesh(cerr);
            cout << "Done." << endl;
        }
        timefile.close();
        my_exit();
    }
#endif
    delete _tesselation;

    _gridtimer.stop();
}

/**
 * @brief Destructor
 *
 * Clean up. Print the timers to the stdout.
 */
VorTessManager::~VorTessManager() {
    delete _tess2;
    cout << "Spent " << _gridtimer.value() << "s in grid construction" << endl;
    cout << "Spent " << _hydrotimer.value() << "s in hydro related functions"
         << endl;
}

/**
 * @brief Reset the tesselation
 *
 * For the old algorithm, this deletes the old tesselation and initializes a new
 * one. For the evolution algorithm, this method does nothing, since we keep the
 * old tesselation.
 *
 * @param cont DelCont specifying the dimensions of the simulation box
 * @param periodic Flag indicating if the simulation box is periodic (true) or
 * reflective (false)
 * @param tolerance Tolerance value used to distinguish between standard
 * floating point arithmetics and arbitrary precision arithmetics in geometry
 * tests
 */
void VorTessManager::reset(DelCont* cont, bool periodic, double tolerance) {
    //#if ndim_==3
    //    _tesselation = new VorTess(cont, _particles.gassize(), periodic,
    //    tolerance);
    //#endif
    return;
}

/**
 * @brief Update the tesselation to the given current system time
 *
 * For the old algorithm, this adds all points to the new tesselation and
 * constructs the new grid. For the evolution algorithm, this method does the
 * actual evolution of the grid, so that at the end of this method, the
 * tesselation is a valid Voronoi tesselation again.
 *
 * @param currentTime Current integer system time
 * @return Number of active cells in the tesselation
 */
unsigned int VorTessManager::update(unsigned long currentTime) {
    _gridtimer.start();
    unsigned int count = 0;
    //#if ndim_==2
    try {
        count = _tess2->update(_particles, currentTime);
    } catch(AdaptiveMeshException e) {
        cout << "Reconstructing mesh using old algorithm" << endl;
#if ndim_ == 3
        // we have to update the positions, since this was not done yet
        vector<Vec> new_positions = _tess2->get_new_positions();
        Cuboid box = _particles.get_container().get_cuboid();
        Vec anchor = box.get_anchor();
        Vec sides = box.get_sides();
        for(unsigned int i = 0; i < _particles.gassize(); i++) {
            Vec new_pos = new_positions[_particles.gas(i)->get_local_id()];
            for(unsigned int j = 0; j < ndim_; j++) {
                if(new_pos[j] < anchor[j]) {
                    new_pos[j] += sides[j];
                }
                if(new_pos[j] >= anchor[j] + sides[j]) {
                    new_pos[j] -= sides[j];
                }
            }
            _particles.gas(i)->set_x(new_pos.x());
            _particles.gas(i)->set_y(new_pos.y());
            _particles.gas(i)->set_z(new_pos.z());
        }
#endif
        delete _tess2;
        _tesselation = new VorTess(&_particles.get_container(),
                                   _particles.gassize(), _periodic, _tolerance);
        for(unsigned int i = 0; i < _particles.gassize(); i++) {
            _tesselation->add_point(_particles.gas(i),
                                    _particles.gas(i)->get_local_id());
            _particles.gas(i)->reset_copies();
            _particles.gas(i)->reset_export();
            if(_particles.gas(i)->get_endtime() == currentTime) {
                count++;
            }
        }
        _tesselation->complete(_particles.get_tree());
        _tesselation->construct();
#if ndim_ == 3
        // make sure the h's are set, since they are used to subdivide the steps
        // in the evolution algorithm
        for(unsigned int i = 0; i < _particles.gassize(); i++) {
            _particles.gas(i)->set_h(_particles.gas(i)->get_cell()->get_h());
        }

        _tess2 = new AdaptiveVorTess3d(_tesselation,
                                       _particles.get_container().get_cuboid(),
                                       _periodic);
#else
        unsigned long maxid = _particles.get_header().npart();
        _tess2 = new AdaptiveVorTess2d(_tesselation,
                                       _particles.get_container().get_cuboid(),
                                       maxid, _particles, _periodic);
#endif
        delete _tesselation;
    }

    //    if(currentTime > 36028797018963968){
    //        ofstream ofile("last_mesh.dat");
    //        _tess2->print_tesselation_gnuplot(ofile);
    //        ofile.close();
    //        _tess2->check_mesh(cerr);
    //    }
    //#else
    //    for(unsigned int i = 0; i < _particles.gassize(); i++){
    //        if(_particles.gas(i)->get_endtime() == currentTime){
    //            _tesselation->add_point(_particles.gas(i),
    //    _particles.gas(i)->id());
    //            count++;
    //        }
    //        _particles.gas(i)->reset_copies();
    //        _particles.gas(i)->reset_export();
    //    }
    //    _tesselation->complete(_particles.get_tree());
    //    _tesselation->construct();
    //    delete _tess2;
    //    _tess2 = new AdaptiveVorTess3d(_tesselation,
    //    _particles.get_container().get_cuboid(), true);
    //    delete _tesselation;
    //#endif
    _gridtimer.stop();
    return count;
}

/**
 * @brief Communicate primitive variables between MPI processes
 *
 * Currently only done for the old algorithm.
 */
void VorTessManager::update_Ws() {
    //#if ndim_==2
    //#warning Have to fix communication stuff
    //#else
    //    _tesselation->update_Ws();
    //#endif
}

/**
 * @brief Calculate characteristic lengths for the active cells of the
 * tesselation
 *
 * The characteristic length is the radius of a sphere/circle with the same
 * volume/face area as the cell
 *
 * @param currentTime Current integer system time
 */
void VorTessManager::set_hs(unsigned long currentTime) {
    //    _tesselation->set_hs();
    for(unsigned int i = 0; i < _particles.gassize(); i++) {
        if(_particles.gas(i)->get_endtime() == currentTime) {
            //#if ndim_==2
            double h = _tess2->get_h(_particles.gas(i)->get_local_id());
            //#else
            //            double h = _particles.gas(i)->get_cell()->get_h();
            //#endif
            _particles.gas(i)->set_h(h);
        }
    }
}

/**
 * @brief Communicate gradients between MPI processes
 *
 * Currently only done for the old algorithm.
 */
void VorTessManager::update_gradients() {
    //#if ndim_==2
    //#warning Have to fix communication stuff
    //#else
    //    _tesselation->update_gradients();
    //#endif
}

/**
 * @brief Communicate timesteps between MPI processes
 *
 * Currently only done for the old algorithm.
 *
 * @param currentTime Current integer system time
 */
void VorTessManager::update_dts(unsigned long currentTime) {
    //#if ndim_==2
    //#warning Have to fix communication stuff
    //#else
    //    _tesselation->update_dts(currentTime);
    //#endif
}

/**
 * @brief Communicate gravitational correction terms between MPI processes
 *
 * Currently only done for the old algorithm.
 */
void VorTessManager::update_gravitational_corrections() {
    // have to implement this
}

/**
 * @brief Do the hydrodynamical integration by calculating the flux exchanges
 * between cells through the faces of the tesselation
 *
 * @param timeline Simulation TimeLine
 * @param solver Riemann solver used to solve the Riemann problem at the faces
 * @param particles ParticleVector containing the particles of the simulation
 */
void VorTessManager::hydro(TimeLine& timeline, RiemannSolver& solver,
                           ParticleVector& particles) {
    _gridtimer.start();
#if ndim_ == 3
    vector<AdaptiveFace3d*> faces =
            _tess2->get_faces(timeline.get_integertime());
#else
    vector<AdaptiveFace2d*> faces =
            _tess2->get_faces(timeline.get_integertime());
#endif
    cout << "Faces size: " << faces.size() << endl;
    _gridtimer.stop();
    _hydrotimer.start();
    for(unsigned int i = 0; i < faces.size(); i++) {
        faces[i]->set_v();
        faces[i]->calculate_flux(timeline, solver);
    }
    for(unsigned int i = 0; i < faces.size(); i++) {
        delete faces[i];
    }
    _hydrotimer.stop();
}

/**
 * @brief Communicate fluxes between MPI processes
 *
 * Currently only done for the old algorithm.
 */
void VorTessManager::update_dQs() {
    //#if ndim_==2
    //#warning Have to fix communication stuff
    //#else
    //    _tesselation->update_dQs();
    //#endif
}

/**
 * @brief Print the tesselation to the given stream in a format that can be
 * plotted using gnuplot
 *
 * @param stream std::ostream to write to
 */
void VorTessManager::print_tesselation_gnuplot(ostream& stream) {
    _tess2->print_tesselation_gnuplot(stream);
}

/**
 * @brief Dump the tesselation to the given RestartFile
 *
 * @param rfile RestartFile to write to
 */
void VorTessManager::dump(RestartFile& rfile) {
    _gridtimer.dump(rfile);
    _hydrotimer.dump(rfile);
    //#if ndim_==2
    //    _tess2->dump(rfile, _particles.get_vector());
    //#endif
}

/**
 * @brief Restart constructor. Initialize the tesselation from the given
 * RestartFile
 *
 * Fot the old algorithm, we reconstruct the grid, since this is easier than
 * storing it and reading it in again.
 *
 * @param rfile RestartFile to read from
 * @param particles ParticleVector containing the particles of the simulation
 */
VorTessManager::VorTessManager(RestartFile& rfile, ParticleVector& particles)
        : _particles(particles), _gridtimer(rfile), _hydrotimer(rfile) {
#if ndim_ == 3
//    _tess2 = new AdaptiveVorTess3d(rfile, _particles.get_vector());
#else
//    _tess2 = new AdaptiveVorTess2d(rfile, _particles.get_vector());
#endif
}

/**
 * @brief Get the velocity of the cell corresponding to the particle with the
 * given index
 *
 * @param index Index of a GasParticle in the simulation ParticleVector
 * @return The velocity of the cell corresponding to the given particle
 */
Vec VorTessManager::get_velocity(unsigned int index) {
    //#if ndim_==2
    return _tess2->get_velocity(_particles.gas(index)->get_local_id(),
                                _particles.gas(index));
    //#else
    //    VorCell* cell = _particles.gas(index)->get_cell();
    //    return cell->get_velocity();
    //#endif
}

/**
 * @brief Estimate slope-limited gradients for the cell corresponding to the
 * particle with the given index
 *
 * @param index Index of a GasParticle in the simulation ParticleVector
 * @param delta Array to store the resulting gradients in
 */
void VorTessManager::estimate_gradients(unsigned int index,
                                        StateVector* delta) {
    _hydrotimer.start();
    //#if ndim_==2
    _tess2->estimate_gradients(_particles.gas(index)->get_local_id(), delta,
                               _particles.gas(index));
    //#else
    //    VorCell* cell = _particles.gas(index)->get_cell();
    //    cell->estimate_gradient(delta);
    //#endif
    _hydrotimer.stop();
}

/**
 * @brief Calculate the volume of the cell corresponding to the particle with
 * the given index
 *
 * @param index Index of a GasParticle in the simulation ParticleVector
 * @return Volume of the cell corresponding to the given particle
 */
double VorTessManager::get_volume(unsigned int index) {
    //#if ndim_==2
    return _tess2->get_volume(_particles.gas(index)->get_local_id());
    //#else
    //    VorCell* cell = _particles.gas(index)->get_cell();
    //    return cell->get_volume();
    //#endif
}

/**
 * @brief Get the total surface area of the Voronoi cell with the given index
 *
 * @param index Index of a GasParticle in the simulation ParticleVector
 * @return Total surface area of the Voronoi cell corresponding to the given
 * particle
 */
double VorTessManager::get_total_area(unsigned int index) {
    //#if ndim_==2
    //    return _tess2->get_total_area(_particles.gas(index)->get_local_id());
    return 0.;
    //#else
    //    VorCell* cell = _particles.gas(index)->get_cell();
    //    return cell->get_volume();
    //#endif
}

/**
 * @brief Calculate the centroid of the cell corresponding to the particle with
 * the given index
 *
 * @param index Index of a GasParticle in the simulation ParticleVector
 * @return Centroid of the cell corresponding to the given particle
 */
Vec VorTessManager::get_centroid(unsigned int index) {
    //#if ndim_==2
    return _tess2->get_centroid(_particles.gas(index)->get_local_id());
    //#else
    //    VorCell* cell = _particles.gas(index)->get_cell();
    //    return cell->get_centroid();
    //#endif
}

/**
 * @brief Update the positions of the generators of the tesselation after a
 * movement of the simulation particles
 *
 * For the old algorithm, we only make sure all generators stay inside the
 * simulation box (by trimming the coordinates for a periodic box). For the
 * evolution algorithm, we update the positions in the tesselation, without
 * correcting for deviations that might occur. The latter is done in
 * VorTessManager::update().
 *
 * @param periodic Flag indicating if the simulation box is periodic (true) or
 * reflective
 */
void VorTessManager::update_positions(bool periodic) {
    _gridtimer.start();
    //#if ndim_==2
    _tess2->update_positions(_particles);
    _gridtimer.stop();
    return;
    //#endif
    if(periodic) {
        for(unsigned int i = 0; i < _particles.gassize(); i++) {
            // make sure the particle stays inside the periodic box
            _particles.get_container().keep_inside(_particles.gas(i));
        }
    }
    _gridtimer.stop();
}

/**
 * @brief Get the correction term to the gravitational acceleration due to the
 * use of variable softening lengths for the GasParticle with the given index
 *
 * @param index Index of a GasParticle in the simulation ParticleVector
 * @return Correction to the acceleration
 */
Vec VorTessManager::get_gravitational_correction(unsigned int index) {
    VorCell* cell = _particles.gas(index)->get_cell();
    return cell->get_gravitational_correction();
}

/**
 * @brief Write connectivity information to the given stream
 *
 * This information can be used to compare tesselations constructed with
 * different algorithms.
 *
 * @param stream std::ostream to write to
 * @param currentTime Current integer system time
 */
void VorTessManager::dump_connectivity(ostream& stream,
                                       unsigned long currentTime) {
    for(unsigned int i = 0; i < 3 * _particles.gassize(); i++) {
        unsigned int j = 0;
        while(j < _particles.gassize() && _particles.gas(j)->id() != i) {
            j++;
        }
        if(j == _particles.gassize()) {
            continue;
        }
        //        unsigned int j = i;
        if(_particles.gas(j)->get_endtime() == currentTime) {
            unsigned long id = _particles.gas(j)->id();
            char c = 'c';
            stream.write(&c, 1);
            stream.write(reinterpret_cast<char*>(&id), sizeof(unsigned long));
            vector<unsigned long> ids =
                    _tess2->get_ngb_ids(_particles.gas(j)->get_local_id());
            for(unsigned int k = 0; k < ids.size(); k++) {
                stream.write(reinterpret_cast<char*>(&ids[k]),
                             sizeof(unsigned long));
            }
        }
    }
}

/**
 * @brief Print cell statistics for the Voronoi mesh
 */
void VorTessManager::print_statistics() {}
