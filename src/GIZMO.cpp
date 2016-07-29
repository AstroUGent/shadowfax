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
 * @file GIZMO.cpp
 *
 * @brief Implementation of Phil Hopkins' GIZMO: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "GIZMO.hpp"
#include "DelCont.hpp"                // for RectangularBox
#include "Error.hpp"                  // for my_exit
#include "MPIMethods.hpp"             // for MyMPI_Finalize, MyMPI_Init
#include "NgbSearch.hpp"
#include "SnapshotHandler.hpp"        // for SnapshotReader
#include "SnapshotReaderFactory.hpp"  // for SnapshotReaderFactory
#include "StateVector.hpp"            // for StateVector
#include "Vec.hpp"                    // for Vec, operator-, operator*, etc
#include "io/Header.hpp"              // for Header
#include "io/UnitSet.hpp"             // for UnitSet
#include "io/UnitSetGenerator.hpp"    // for UnitSetGenerator
#include "riemann/RiemannSolver.hpp"  // for RiemannSolver
#include "riemann/RiemannSolverFactory.hpp"  // for RiemannSolverFactory
#include "utilities/GasParticle.hpp"         // for GasParticle
#include "utilities/ParticleVector.hpp"      // for ParticleVector
#include "utilities/Tree.hpp"                // for Tree
#include <algorithm>                         // for max, min
#include <cmath>                             // for M_PI
#include <cstddef>                           // for NULL
#include <ext/alloc_traits.h>
#include <getopt.h>  // for optarg, getopt_long, etc
#include <iostream>  // for operator<<, basic_ostream, etc
#include <sstream>
#include <string>  // for string, operator<<, etc
#include <vector>  // for vector, vector<>::reference
using namespace std;

/**
 * @brief Cubic spline kernel
 *
 * Taken from Springel 2005.
 *
 * @param r Distance between two particles
 * @param h Smoothing length of one of the particles
 * @return Value of the kernel at the given distance
 */
double GIZMO::kernel(double r, double h) {
    double u = r / h;
    if(u >= 1.) {
        return 0.;
    } else {
        if(u < 0.5) {
            return 8. / M_PI / h / h / h * (1. - 6. * u * u + 6. * u * u * u);
        } else {
            double fac = 1. - u;
            return 16. / M_PI / h / h / h * fac * fac * fac;
        }
    }
}

/**
 * @brief Derivative of the cubic spline kernel
 *
 * @param r Distance between two particles
 * @param h Smoothing length of one of the particles
 * @return Value of the derivative of the kernel at the given distance
 */
double GIZMO::kernel_d(double r, double h) {
    double u = r / h;
    if(u >= 1.) {
        return 0.;
    } else {
        if(u < 0.5) {
            return 8. / M_PI / h / h / h / h * (-12. * u + 18. * u * u);
        } else {
            double fac = 1. - u;
            return -48. / M_PI / h / h / h / h * fac * fac;
        }
    }
}

/**
 * @brief Invert the given symmetric matrix
 *
 * Given are: {M[0][0], M[0][1], M[0][2], M[1][1], M[1][2], M[2][2]}.
 *
 * Outputted is: {IM[0][0], IM[0][1], IM[0][2], IM[1][1], IM[1][2], IM[2][2]}.
 *
 * @param M Matrix
 * @param IM Array to store the inverse matrix in
 */
void GIZMO::invert_matrix(double* M, double* IM) {
#ifdef LU_DECOMPOSITION
    // obtain LU decomposition
    double l[3];
    double u[6];
    // first column
    u[0] = M[0];
    l[0] = M[1] / u[0];
    l[1] = M[2] / u[0];
    // second column
    u[1] = M[1];
    u[3] = M[3] - u[1] * l[0];
    l[2] = (M[4] - u[1] * l[1]) / u[3];
    // third column
    u[2] = M[2];
    u[4] = M[4] - u[2] * l[0];
    u[5] = M[5] - u[2] * l[1] - l[2] * u[4];

    // calculate first column of inverse matrix
    // LZ = (100)
    double z[3];
    // correct version
    //    z[0] = 1.;
    //    z[1] = -l[0]*z[0];
    //    z[2] = -l[1]*z[0] - l[2]*z[1];
    // fast version
    z[0] = 1.;
    z[1] = -l[0];
    z[2] = -l[1] - l[2] * z[1];
    // UX = Z
    IM[2] = z[2] / u[5];
    IM[1] = (z[1] - u[4] * IM[2]) / u[3];
    IM[0] = (z[0] - u[1] * IM[1] - u[2] * IM[2]) / u[0];

    // calculate second column
    // correct version
    //    z[0] = 0.;
    //    z[1] = 1. - l[0]*z[0];
    //    z[2] = -l[1]*z[0] - l[2]*z[1];
    // fast version
    z[1] = 1.;
    z[2] = -l[2];
    IM[4] = z[2] / u[5];
    IM[3] = (z[1] - u[4] * IM[4]) / u[3];

    // calculate third column (single element)
    // correct version
    //    z[0] = 0.;
    //    z[1] = -l[0]*z[0];
    //    z[2] = 1. - l[1]*z[0] - l[2]*z[1];
    // fast version
    z[2] = 1.;
    IM[5] = z[2] / u[5];
#else
    double detM = M[0] * M[3] * M[5] + 2. * M[1] * M[2] * M[4] -
                  M[2] * M[2] * M[3] - M[0] * M[4] * M[4] - M[1] * M[1] * M[5];
    if(!detM) {
        cerr << "Determinant 0, inversion won't work!" << endl;
        my_exit();
    }
    IM[0] = (M[3] * M[5] - M[4] * M[4]) / detM;
    IM[1] = (M[2] * M[4] - M[1] * M[5]) / detM;
    IM[2] = (M[1] * M[4] - M[2] * M[3]) / detM;
    IM[3] = (M[0] * M[5] - M[2] * M[2]) / detM;
    IM[4] = (M[1] * M[2] - M[0] * M[4]) / detM;
    IM[5] = (M[0] * M[3] - M[1] * M[1]) / detM;
#endif
}

/**
 * @brief Check if the tree is correct by searching the particles themselves
 *
 * @param particles ParticleVector
 */
void GIZMO::test_tree(ParticleVector& particles) {
    for(unsigned int i = 0; i < particles.gassize(); i++) {
        GasParticle* p = particles.gas(i);
        vector<bool> exports;
        NgbSearch ngbsearch(p->get_position(), 1.e-10, exports, 100);
        particles.get_tree().walk_tree(ngbsearch);
        vector<GasParticle*> ngbs = ngbsearch.get_ngbs();
        if(!ngbs.size()) {
            cerr << "Can't find particle in tree!" << endl;
            my_exit();
        } else {
            if(ngbs.size() > 1) {
                cerr << "Too many ngbs!" << endl;
                my_exit();
            } else {
                if(ngbs[0] != p) {
                    cerr << "Wrong particle!" << endl;
                    my_exit();
                }
            }
        }
    }
    cout << "Passed tree test" << endl;
}

/**
 * @brief Main program
 *
 * @param argc Number of command line arguments
 * @param argv Command line arguments
 */
GIZMO::GIZMO(int argc, char** argv) {
    // test matrix inversion
    //    double M[6] = {1., 2., 3., 2., 3., 3.};
    //    double IM[6];
    //    invert_matrix(M, IM);
    //    cout << IM[0] << "\t" << IM[1] << "\t" << IM[2] << endl;
    //    cout << IM[1] << "\t" << IM[3] << "\t" << IM[4] << endl;
    //    cout << IM[2] << "\t" << IM[4] << "\t" << IM[5] << endl;
    //    return;

    // read command line arguments
    string readertype = "Gadget";
    string filename;
    string onlyname;

    double desnumngb = 5.;
    double maxnumngbdev = 1.;

    static struct option long_options[] = {
            {"type", required_argument, NULL, 't'},
            {"filename", required_argument, NULL, 'o'},
            {0, 0, 0, 0}};

    int c;
    // force rescan of the arguments
    optind = 1;
    opterr = 0;
    while((c = getopt_long(argc, argv, ":t:o:", long_options, NULL)) != -1) {
        switch(c) {
            case 't':
                readertype = optarg;
                break;
            case 'o':
                filename = optarg;
                // drop the last 5 characters from the name by setting the end
                // of
                // string
                optarg[filename.size() - 5] = '\0';
                onlyname = optarg;
                break;
            case ':':
                cerr << "Error! Missing required argument for option " << optopt
                     << "!" << endl;
                return;
        }
    }

    // command line errors
    if(!filename.size()) {
        cerr << "Error! No filename specified!" << endl;
        return;
    }

    cout << "Reading initial conditions: " << onlyname << endl;

    RectangularBox container;
    ParticleVector particles(false, container);
    UnitSet* simulation_units = UnitSetGenerator::generate("CGS");
    SnapshotReader* reader = SnapshotReaderFactory::generate(
            readertype, filename, *simulation_units);
    RiemannSolver* solver =
            RiemannSolverFactory::generate("Exact", 5. / 3., 1.e-8, -5.);
    Header header = reader->read_snapshot(particles);
    delete reader;
    bool periodic = header.periodic();
    particles.set_periodic(periodic);

    double box[ndim_ + ndim_] = {0.};
    header.box(box);
    Vec center;
    Vec sides;
    for(unsigned int i = ndim_; i--;) {
        center[i] = box[i] + 0.5 * box[ndim_ + i];
        sides[i] = box[ndim_ + i];
    }
    container = RectangularBox(center, sides);
    particles.set_container(container);

    cout << "Reading ready. Starting initialization." << endl;

    // sort particles and build tree
    particles.sort();

    test_tree(particles);

    // initialize hydro variables
    vector<StateVector> Ws(particles.gassize());
    vector<StateVector> Qs(particles.gassize());
    vector<GasParticle*> pcopy(particles.gassize());
    vector<Vec> vs(particles.gassize());
    for(unsigned int i = 0; i < particles.gassize(); i++) {
        Ws[i] = particles.gas(i)->get_Wvec();
        pcopy[i] = particles.gas(i);
        // correct for wrong conversion done during snapshot reading
        // (InternalEnergy field actually contains Pressure)
        //        Ws[i][4] = 1.5*Ws[i][4]/Ws[i][0];
        // small velocity
        //        vs[i][0] = 1.0;
    }

    // first neighbour loop: calculate volumes
    vector<double> volumes(particles.gassize(), 0.);
    vector<double> matrix(particles.gassize() * 6, 0.);
    vector<double> h(particles.gassize(), 0.038);
    vector<bool> todo(particles.gassize(), true);
    unsigned int loopcount = 0;
    unsigned int redo = 1;
    while(redo && loopcount < 100) {
        loopcount++;
        cout << "Loop " << loopcount << ", redo = " << redo << endl;

        redo = 0;
        for(unsigned int i = 0; i < particles.gassize(); i++) {
            if(todo[i]) {
                GasParticle* pi = particles.gas(i);
                vector<GasParticle*> ngbs;
                vector<bool> exports;
                particles.get_tree().get_neighbours_periodic(
                        pi->get_position(), h[i], ngbs, exports);
                double nngb = 0.;
                double dw_du = 0.;
                matrix[6 * i] = 0.;
                matrix[6 * i + 1] = 0.;
                matrix[6 * i + 2] = 0.;
                matrix[6 * i + 3] = 0.;
                matrix[6 * i + 4] = 0.;
                matrix[6 * i + 5] = 0.;
                for(unsigned int j = 0; j < ngbs.size(); j++) {
                    GasParticle* pj = ngbs[j];
                    Vec d = pi->get_position() - pj->get_position();
                    // account for periodic boundaries
                    for(unsigned idim = 0; idim < 3; idim++) {
                        if(d[idim] > 0.5) {
                            d[idim] -= 1.;
                        }
                        if(d[idim] <= -0.5) {
                            d[idim] += 1.;
                        }
                    }
                    double kernval = kernel(d.norm(), h[i]);
                    nngb += kernval;
                    dw_du += kernel_d(d.norm(), h[i]);
                    matrix[6 * i] += d[0] * d[0] * kernval;
                    matrix[6 * i + 1] += d[0] * d[1] * kernval;
                    matrix[6 * i + 2] += d[0] * d[2] * kernval;
                    matrix[6 * i + 3] += d[1] * d[1] * kernval;
                    matrix[6 * i + 4] += d[1] * d[2] * kernval;
                    matrix[6 * i + 5] += d[2] * d[2] * kernval;
                }

                volumes[i] = nngb;
                nngb *= h[i] * h[i] * h[i];
                if(nngb > (desnumngb + maxnumngbdev) ||
                   nngb < (desnumngb - maxnumngbdev)) {
                    dw_du *= h[i] * h[i] * h[i];
                    redo++;
                    double hcorr = h[i];
                    if(dw_du) {
                        hcorr = (nngb - desnumngb) / dw_du;
                        hcorr = std::min(hcorr, h[i]);
                        hcorr = std::max(hcorr, -0.5 * h[i]);
                    }
                    h[i] += hcorr;
                } else {
                    todo[i] = false;
                    volumes[i] = 1. / volumes[i];
                    if(particles.gas(i)->id() == 0) {
                        cout << "numngb: " << ngbs.size() << endl;
                    }
                }
            }
        }
    }

    // we have volumes, calculate conserved quantities and invert matrices
    for(unsigned int i = 0; i < particles.gassize(); i++) {
        double m = Ws[i].rho() * volumes[i];
        double px = m * Ws[i].vx();
        double py = m * Ws[i].vy();
        double pz = m * Ws[i].vz();
        double e =
                0.5 * m * (Ws[i].vx() * Ws[i].vx() + Ws[i].vy() * Ws[i].vy() +
                           Ws[i].vz() * Ws[i].vz()) +
                1.5 * Ws[i].p() * volumes[i];
        Qs[i].set(m, px, py, pz, e);

        // set particle velocity
        vs[i].set(Ws[i].vx(), Ws[i].vy(), Ws[i].vz());

        double inverse_matrix[6];
        invert_matrix(&matrix[6 * i], inverse_matrix);
        matrix[6 * i] = inverse_matrix[0];
        matrix[6 * i + 1] = inverse_matrix[1];
        matrix[6 * i + 2] = inverse_matrix[2];
        matrix[6 * i + 3] = inverse_matrix[3];
        matrix[6 * i + 4] = inverse_matrix[4];
        matrix[6 * i + 5] = inverse_matrix[5];
    }

    cout << "Initialization ready. Starting integration." << endl;

    double dt = 0.0015625;
    Vec mid(0.5, 0.5, 0.5);
    for(unsigned int step = 0; step < 64; step++) {
        cout << "Step " << (step + 1) << endl;
        // exchange fluxes
        vector<double> Atot(particles.gassize(), 0.);
        for(unsigned int i = 0; i < particles.gassize(); i++) {
            GasParticle* pi = particles.gas(i);
            vector<GasParticle*> ngbs;
            vector<bool> exports;
            particles.get_tree().get_neighbours_periodic(pi->get_position(),
                                                         h[i], ngbs, exports);
            for(unsigned int j = 0; j < ngbs.size(); j++) {
                GasParticle* pj = ngbs[j];
                if(pi == pj) {
                    continue;
                }
                unsigned int k = pj->get_local_id();
                Vec d = pi->get_position() - pj->get_position();
                // account for periodic boundaries
                for(unsigned idim = 0; idim < 3; idim++) {
                    if(d[idim] > 0.5) {
                        d[idim] -= 1.;
                    }
                    if(d[idim] <= -0.5) {
                        d[idim] += 1.;
                    }
                }
                double kernval = kernel(d.norm(), h[i]);
                Vec A;
                A[0] = -volumes[i] * (matrix[6 * i] * d[0] +
                                      matrix[6 * i + 1] * d[1] +
                                      matrix[6 * i + 2] * d[2]) *
                               kernval -
                       volumes[k] * (matrix[6 * k] * d[0] +
                                     matrix[6 * k + 1] * d[1] +
                                     matrix[6 * k + 2] * d[2]) *
                               kernval;
                A[1] = -volumes[i] * (matrix[6 * i + 1] * d[0] +
                                      matrix[6 * i + 3] * d[1] +
                                      matrix[6 * i + 4] * d[2]) *
                               kernval -
                       volumes[k] * (matrix[6 * k + 1] * d[0] +
                                     matrix[6 * k + 3] * d[1] +
                                     matrix[6 * k + 4] * d[2]) *
                               kernval;
                A[2] = -volumes[i] * (matrix[6 * i + 2] * d[0] +
                                      matrix[6 * i + 4] * d[1] +
                                      matrix[6 * i + 5] * d[2]) *
                               kernval -
                       volumes[k] * (matrix[6 * k + 2] * d[0] +
                                     matrix[6 * k + 4] * d[1] +
                                     matrix[6 * k + 5] * d[2]) *
                               kernval;
                double Anrm = A.norm();
                Atot[i] += Anrm;
                if(!Anrm) {
                    continue;
                }
                Vec n = A / Anrm;

                Vec x = -h[i] / (h[i] + h[k]) * d;
                Vec v = vs[i] + (vs[i] - vs[k]) * Vec::dot(d, x) / d.norm2();

                StateVector WL = Ws[i];
                StateVector WR = Ws[k];

                WL[1] -= v[0];
                WL[2] -= v[1];
                WL[3] -= v[2];
                WR[1] -= v[0];
                WR[2] -= v[1];
                WR[3] -= v[2];
                double mach;
                StateVector Wstar = solver->solve(WL, WR, n, mach);
                Wstar[1] += v[0];
                Wstar[2] += v[1];
                Wstar[3] += v[2];
                //                if(i == 3957 &&  k == 3864){
                //                    cerr << Wstar[0] << "\t" << Wstar[1] <<
                //                    "\t" << Wstar[2]
                //                         << "\t" << Wstar[3] << "\t" <<
                //                         Wstar[4] << endl;
                //                    my_exit();
                //                }
                StateVector flux[3];
                flux[0][0] = Wstar.rho() * (Wstar.vx() - v[0]);
                flux[1][0] = Wstar.rho() * (Wstar.vy() - v[1]);
                flux[2][0] = Wstar.rho() * (Wstar.vz() - v[2]);
                flux[0][1] = Wstar.rho() * Wstar.vx() * (Wstar.vx() - v[0]) +
                             Wstar.p();
                flux[1][1] = Wstar.rho() * Wstar.vx() * (Wstar.vy() - v[1]);
                flux[2][1] = Wstar.rho() * Wstar.vx() * (Wstar.vz() - v[2]);
                flux[0][2] = Wstar.rho() * Wstar.vy() * (Wstar.vx() - v[0]);
                flux[1][2] = Wstar.rho() * Wstar.vy() * (Wstar.vy() - v[1]) +
                             Wstar.p();
                flux[2][2] = Wstar.rho() * Wstar.vy() * (Wstar.vz() - v[2]);
                flux[0][3] = Wstar.rho() * Wstar.vz() * (Wstar.vx() - v[0]);
                flux[1][3] = Wstar.rho() * Wstar.vz() * (Wstar.vy() - v[1]);
                flux[2][3] = Wstar.rho() * Wstar.vz() * (Wstar.vz() - v[2]) +
                             Wstar.p();
                double ehalf = 1.5 * Wstar.p() / Wstar.rho() +
                               0.5 * (Wstar.vx() * Wstar.vx() +
                                      Wstar.vy() * Wstar.vy() +
                                      Wstar.vz() * Wstar.vz());
                flux[0][4] = Wstar.rho() * ehalf * (Wstar.vx() - v[0]) +
                             Wstar.p() * Wstar.vx();
                flux[1][4] = Wstar.rho() * ehalf * (Wstar.vy() - v[1]) +
                             Wstar.p() * Wstar.vy();
                flux[2][4] = Wstar.rho() * ehalf * (Wstar.vz() - v[2]) +
                             Wstar.p() * Wstar.vz();

                Qs[i][0] -= dt * (flux[0][0] * A[0] + flux[1][0] * A[1] +
                                  flux[2][0] * A[2]);
                Qs[i][1] -= dt * (flux[0][1] * A[0] + flux[1][1] * A[1] +
                                  flux[2][1] * A[2]);
                Qs[i][2] -= dt * (flux[0][2] * A[0] + flux[1][2] * A[1] +
                                  flux[2][2] * A[2]);
                Qs[i][3] -= dt * (flux[0][3] * A[0] + flux[1][3] * A[1] +
                                  flux[2][3] * A[2]);
                Qs[i][4] -= dt * (flux[0][4] * A[0] + flux[1][4] * A[1] +
                                  flux[2][4] * A[2]);

                // i will not be found as a neighbour of k, so we have to do the
                // flux exchange for k here as well
                if(d.norm() > h[k]) {
                    Qs[k][0] += dt * (flux[0][0] * A[0] + flux[1][0] * A[1] +
                                      flux[2][0] * A[2]);
                    Qs[k][1] += dt * (flux[0][1] * A[0] + flux[1][1] * A[1] +
                                      flux[2][1] * A[2]);
                    Qs[k][2] += dt * (flux[0][2] * A[0] + flux[1][2] * A[1] +
                                      flux[2][2] * A[2]);
                    Qs[k][3] += dt * (flux[0][3] * A[0] + flux[1][3] * A[1] +
                                      flux[2][3] * A[2]);
                    Qs[k][4] += dt * (flux[0][4] * A[0] + flux[1][4] * A[1] +
                                      flux[2][4] * A[2]);
                }
            }
        }
        // update primitive variables
        for(unsigned int i = 0; i < particles.gassize(); i++) {
            double rho = Qs[i].m() / volumes[i];
            double vx = Qs[i].px() / Qs[i].m();
            double vy = Qs[i].py() / Qs[i].m();
            double vz = Qs[i].pz() / Qs[i].m();
            double p =
                    2. / 3. *
                    (Qs[i].e() -
                     0.5 * (Qs[i].px() * Qs[i].px() + Qs[i].py() * Qs[i].py() +
                            Qs[i].pz() * Qs[i].pz()) /
                             Qs[i].m()) /
                    volumes[i];
            Ws[i].set(rho, vx, vy, vz, p);
            todo[i] = true;
        }
        //        mid += dt*vs[0];

        // move particles
        for(unsigned int i = 0; i < particles.gassize(); i++) {
            Vec pos = particles.gas(i)->get_position();
            pos += dt * vs[i];
            // keep inside
            for(unsigned int idim = 0; idim < 3; idim++) {
                if(pos[idim] > 1.) {
                    pos[idim] -= 1.;
                }
                if(pos[idim] <= 0.) {
                    pos[idim] += 1.;
                }
            }
            particles.gas(i)->set_x(pos.x());
            particles.gas(i)->set_y(pos.y());
            particles.gas(i)->set_z(pos.z());
        }
        particles.sort();

        // reorder all lists, since the order in particles changed
        // all lists: h, Ws, Qs, vs (volumes and matrix are recomputed)
        // make copies of all lists
        vector<double> hcopy = h;
        vector<StateVector> Wcopy = Ws;
        vector<StateVector> Qcopy = Qs;
        vector<Vec> vcopy = vs;
        for(unsigned int i = 0; i < pcopy.size(); i++) {
            h[pcopy[i]->get_local_id()] = hcopy[i];
            Ws[pcopy[i]->get_local_id()] = Wcopy[i];
            Qs[pcopy[i]->get_local_id()] = Qcopy[i];
            vs[pcopy[i]->get_local_id()] = vcopy[i];
        }
        for(unsigned int i = 0; i < particles.gassize(); i++) {
            pcopy[i] = particles.gas(i);
        }

        stringstream name;
        name << onlyname;
        name.fill('0');
        name.width(3);
        name << step;
        name << ".txt";
        cout << "Writing " << name.str() << "..." << endl;
        ofstream volfile(name.str());
        cout << "Mid: " << mid[0] << " " << mid[1] << " " << mid[2] << endl;
        for(unsigned int i = 0; i < particles.gassize(); i++) {
            double r = (particles.gas(i)->get_position() - mid).norm();
            volfile << r << "\t" << Ws[i].rho() << "\t" << Ws[i].p() << "\n";
        }
        volfile.flush();

        //        test_tree(particles);

        // recalculate volumes
        cout << "Recalculating volumes..." << endl;
        loopcount = 0;
        redo = 1;
        while(redo && loopcount < 100) {
            loopcount++;
            cout << "Loop " << loopcount << ", redo = " << redo << endl;

            redo = 0;
            for(unsigned int i = 0; i < particles.gassize(); i++) {
                if(todo[i]) {
                    GasParticle* pi = particles.gas(i);
                    vector<GasParticle*> ngbs;
                    vector<bool> exports;
                    particles.get_tree().get_neighbours_periodic(
                            pi->get_position(), h[i], ngbs, exports);
                    double nngb = 0.;
                    double dw_du = 0.;
                    matrix[6 * i] = 0.;
                    matrix[6 * i + 1] = 0.;
                    matrix[6 * i + 2] = 0.;
                    matrix[6 * i + 3] = 0.;
                    matrix[6 * i + 4] = 0.;
                    matrix[6 * i + 5] = 0.;
                    for(unsigned int j = 0; j < ngbs.size(); j++) {
                        GasParticle* pj = ngbs[j];
                        Vec d = pi->get_position() - pj->get_position();
                        // account for periodic boundaries
                        for(unsigned idim = 0; idim < 3; idim++) {
                            if(d[idim] > 0.5) {
                                d[idim] -= 1.;
                            }
                            if(d[idim] <= -0.5) {
                                d[idim] += 1.;
                            }
                        }
                        double kernval = kernel(d.norm(), h[i]);
                        nngb += kernval;
                        dw_du += kernel_d(d.norm(), h[i]);
                        matrix[6 * i] += d[0] * d[0] * kernval;
                        matrix[6 * i + 1] += d[0] * d[1] * kernval;
                        matrix[6 * i + 2] += d[0] * d[2] * kernval;
                        matrix[6 * i + 3] += d[1] * d[1] * kernval;
                        matrix[6 * i + 4] += d[1] * d[2] * kernval;
                        matrix[6 * i + 5] += d[2] * d[2] * kernval;
                    }

                    volumes[i] = nngb;
                    nngb *= h[i] * h[i] * h[i];
                    if(nngb > (desnumngb + maxnumngbdev) ||
                       nngb < (desnumngb - maxnumngbdev)) {
                        dw_du *= h[i] * h[i] * h[i];
                        redo++;
                        double hcorr = h[i];
                        if(dw_du) {
                            hcorr = (nngb - desnumngb) / dw_du;
                            hcorr = std::min(hcorr, h[i]);
                            hcorr = std::max(hcorr, -0.5 * h[i]);
                        }
                        h[i] += hcorr;
                    } else {
                        todo[i] = false;
                        volumes[i] = 1. / volumes[i];
                        if(particles.gas(i)->id() == 0) {
                            cout << "numngb: " << ngbs.size() << endl;
                        }
                    }
                }
            }
        }

        // recalculate primitive variables
        for(unsigned int i = 0; i < particles.gassize(); i++) {
            double rho = Qs[i].m() / volumes[i];
            double vx = Qs[i].px() / Qs[i].m();
            double vy = Qs[i].py() / Qs[i].m();
            double vz = Qs[i].pz() / Qs[i].m();
            double p =
                    2. / 3. *
                    (Qs[i].e() -
                     0.5 * (Qs[i].px() * Qs[i].px() + Qs[i].py() * Qs[i].py() +
                            Qs[i].pz() * Qs[i].pz()) /
                             Qs[i].m()) /
                    volumes[i];
            Ws[i].set(rho, vx, vy, vz, p);

            // set particle velocity
            vs[i].set(Ws[i].vx(), Ws[i].vy(), Ws[i].vz());

            // invert matrix
            double inverse_matrix[6];
            invert_matrix(&matrix[6 * i], inverse_matrix);
            matrix[6 * i] = inverse_matrix[0];
            matrix[6 * i + 1] = inverse_matrix[1];
            matrix[6 * i + 2] = inverse_matrix[2];
            matrix[6 * i + 3] = inverse_matrix[3];
            matrix[6 * i + 4] = inverse_matrix[4];
            matrix[6 * i + 5] = inverse_matrix[5];
        }
    }

    cout << "Integration ready. Writing output." << endl;

    ofstream volfile("volumes.txt");
    for(unsigned int i = 0; i < particles.gassize(); i++) {
        double r = (particles.gas(i)->get_position() - mid).norm();
        volfile << r << "\t" << Ws[i].rho() << "\t" << Ws[i].p() << "\n";
    }

    delete solver;
    delete simulation_units;
}

/**
 * @brief Entrance point for the GIZMO program
 *
 * Calls the GIZMO constructor.
 *
 * @param argc Number of command line arguments
 * @param argv Command line arguments
 * @return Exit code
 */
int main(int argc, char** argv) {
    MyMPI_Init(&argc, &argv);

    GIZMO(argc, argv);

    return MyMPI_Finalize();
}
