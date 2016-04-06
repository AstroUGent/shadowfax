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
 * @file ICMaker.cpp
 *
 * @brief Sideprogram to generate initial conditions: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "ICMaker.hpp"
#include "ICGenerator.hpp"
#include "MPIGlobal.hpp"
#include "SnapshotWriterFactory.hpp"
#include "io/UnitSet.hpp"
#include "utilities/HelperFunctions.hpp"
#include "utilities/ParticleVector.hpp"
#include <getopt.h>
#include <mpi.h>
#include <string>
#include <vector>
using namespace std;

/**
  * \brief Main routine to generate initial condition files
  *
  * The (optional) arguments to this program are
\verbatim
./icmakerNd [ncell [mode (cart|rand) [filename [seed [name [type]]]]]]
\endverbatim
  *   - ncell: number of cells to be generated. If a cartesian grid is chosen,
  *     the actual number of cells generated can be lower (see
  *     BlockICGenerator::BlockICGenerator) (default: 10000)
  *   - mode: a cartesian grid (cart) or a regularized random sampled grid
  *     (rand) (default: rand)
  *   - filename: name of the xml file containing necessary information for the
  *     BlockICGenerator (default: "overdensity(3d).xml")
  *   - seed: seed for the random generator used (default: 42), ignored for
  *     parallel runs, where every process has a different seed
  *   - name: name of the initial condition file that is generated (default:
  *     icfile.hdf5)
  *   - type: type of the initial condition file (Shadowfax/Gadget, default:
  *     Gadget)
  *
  * @param argc Number of command line arguments
  * @param argv Array of command line arguments
  */
ICMaker::ICMaker(int argc, char** argv) {
    // suppress output for processes other than the process with rank 0
    if(MPIGlobal::rank) {
        cout.rdbuf(NULL);
    }

    unsigned int ncell = 0;
    ICMode mode = IC_RAND;
    string setup;
    string spec;
    unsigned int seed = 42;
    string output_type = "Gadget";
    string output_name;

    static struct option long_options[] = {
            {"ncell", required_argument, NULL, 'n'},
            {"mode", required_argument, NULL, 'm'},
            {"setup", required_argument, NULL, 's'},
            {"predefined", required_argument, NULL, 'p'},
            {"seed", required_argument, NULL, 'r'},
            {"type", required_argument, NULL, 't'},
            {"filename", required_argument, NULL, 'o'},
            {0, 0, 0, 0}};

    int c;
    // force rescan of the arguments
    optind = 1;
    opterr = 0;
    while((c = getopt_long(argc, argv, ":n:m:s:p:r:t:o:", long_options,
                           NULL)) != -1) {
        switch(c) {
            case 'n':
                ncell = atoi(optarg);
                break;
            case 'm':
                if(string(optarg) == "cart") {
                    mode = IC_CART;
                } else {
                    if(string(optarg) == "rand") {
                        mode = IC_RAND;
                    } else {
                        cerr << "Error! Unknown grid mode: " << optarg << endl;
                        return;
                    }
                }
                break;
            case 's':
                setup = optarg;
                break;
            case 'p':
                spec = optarg;
                break;
            case 'r':
                seed = atoi(optarg);
                break;
            case 't':
                output_type = optarg;
                break;
            case 'o':
                output_name = optarg;
                break;
            case ':':
                cerr << "Error! Missing required argument for option " << optopt
                     << "!" << endl;
                return;
        }
    }

    // command line errors
    if(!setup.size() && !spec.size()) {
        cerr << "Error! No setup argument specified!" << endl;
        cerr << "Either specify a file containing a simulations setup or "
                "choose a predefined setup."
             << endl;
        return;
    }
    if(setup.size() && spec.size()) {
        cerr << "Error! Multiple setup arguments specified!" << endl;
        cerr << "You cannot specify both a setup file and a predefined setup!"
             << endl;
        return;
    }

    // command line warnings
    if(!ncell) {
        cout << "Warning! No number of cells specified. Using default value of "
                "10,000."
             << endl;
        ncell = 10000;
    }
    if(!output_name.size()) {
        cout << "Warning! No output filename specified. Using default value "
             << "\"icfile.hdf5\"." << endl;
        output_name = "icfile.hdf5";
    }

    // allocate MPI communication buffer
    unsigned int maxsize = 1 << 30;  // 1 GB
    unsigned int size = ncell * sizeof(GasParticle) * 100;
    // since size is a 32-bit integer, we have to make sure it is small enough
    // before entering the loop below. If not, we risk overflowing, in which
    // case we enter an endless loop...
    if(size > maxsize) {
        size = maxsize;
    }
    MPIGlobal::sendsize = 1;
    while(MPIGlobal::sendsize < size) {
        MPIGlobal::sendsize <<= 1;
    }
    MPIGlobal::sendsize = std::min(MPIGlobal::sendsize, maxsize);
    MPIGlobal::recvsize = MPIGlobal::sendsize;
    delete[] MPIGlobal::sendbuffer;
    delete[] MPIGlobal::recvbuffer;
    MPIGlobal::sendbuffer = new char[MPIGlobal::sendsize];
    MPIGlobal::recvbuffer = new char[MPIGlobal::recvsize];
    cout << "Assigned "
         << HelperFunctions::human_readable_bytes(MPIGlobal::sendsize)
         << " for MPI communication buffer" << endl;

    ICGenerator* icgen;
    if(spec.size()) {
        icgen = new SpecificICGenerator(ncell, 0, atoi(spec.c_str()), seed,
                                        mode);
    } else {
        icgen = new BlockICGenerator(ncell, mode, seed);
        ((BlockICGenerator*)icgen)->read_xml(setup);
    }
    ParticleVector icpart = icgen->generate();

    Header& header = icpart.get_header();
    header.set_time(0);
    header.set_global_timestep(false);
    UnitSet SI_units;

    SnapshotWriter* writer = SnapshotWriterFactory::generate(
            output_type, output_name, SI_units, SI_units, -1);
    writer->write_snapshot(0, icpart);
    delete writer;
    delete icgen;
}

#ifdef ICMAKER

/**
 * @brief Entrance point for the icmaker program
 *
 * Calls the ICMaker constructor.
 *
 * @param argc Number of command line arguments
 * @param argv Command line arguments
 * @return Exit code
 */
int main(int argc, char** argv) {
    MyMPI_Init(&argc, &argv);

    ICMaker(argc, argv);

    return MyMPI_Finalize();
}
#endif
