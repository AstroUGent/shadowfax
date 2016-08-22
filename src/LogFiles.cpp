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
 * @file LogFiles.cpp
 *
 * @brief User-friendly diagnostic log files: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "LogFiles.hpp"
#include "MPIGlobal.hpp"              // for rank
#include "ProgramLog.hpp"             // for LOGS
#include "RestartFile.hpp"            // for RestartFile
#include "StateVector.hpp"            // for StateVector
#include "TimeLine.hpp"               // for TimeLine
#include "Vec.hpp"                    // for Vec
#include "riemann/RiemannSolver.hpp"  // for RiemannSolver
#include "utilities/DMParticle.hpp"
#include "utilities/GasParticle.hpp"
#include "utilities/HelperFunctions.hpp"  // for human_readable_counter, etc
#include "utilities/ParticleVector.hpp"   // for ParticleVector
#include <algorithm>                      // for min
#include <iostream>                       // for cout
#include <stdio.h>                        // for rename
using namespace std;

/**
 * @brief Constructor
 *
 * Open the log files.
 *
 * @param outputdir Directory where the logfiles should be stored
 */
LogFiles::LogFiles(std::string outputdir) {
    if(!MPIGlobal::rank) {
        _timer.start();

        // time log
        _timename = outputdir + string("/timelog.txt");
        _timefile.open(_timename);
        _timefile << "#time\titime\tdt" << endl;

        // energy log
        _energyname = outputdir + string("/energylog.txt");
        _energyfile.open(_energyname);
        _energyfile << "#time\te\tekiny\tu\tepot\tekin\tesum" << endl;

        // statistics log
        _statname = outputdir + string("/statistics.txt");
        _statfile.open(_statname);
        _stat_stepcount = 1;
        _stat_solvercount = 0;

        _timer.stop();

        LOGS("LogFiles opened");
    }
}

LogFiles::~LogFiles() {
    cout << "Spent " << _timer.value() << "s writing log files." << endl;
}

/**
 * @brief Write a line to the log files
 *
 * @param timeline TimeLine of the simulation
 * @param particles Simulation ParticleVector
 * @param solver RiemannSolver used to solve the Riemann problem for the
 * hydrodynamics
 */
void LogFiles::write_to_logs(TimeLine*& timeline, ParticleVector*& particles,
                             RiemannSolver*& solver) {
    _timer.start();
    // gather data

    // energy log
    double totenergy = 0.;
    double ekiny = 0.;
    double e_therm = 0.;
    double e_pot = 0.;
    double e_kin = 0.;
    double e_tot = 0.;
    for(unsigned int i = particles->gassize(); i--;) {
        totenergy += particles->gas(i)->get_Qvec().e();
        ekiny += 0.5 * particles->gas(i)->get_Qvec().m() *
                 particles->gas(i)->get_Wvec().vy() *
                 particles->gas(i)->get_Wvec().vy();
        StateVector Wp = particles->gas(i)->get_Wvec();
        StateVector Qp = particles->gas(i)->get_Qvec();
        double p_e_therm =
                Wp.p() / Wp.rho() / (solver->get_gamma() - 1.) * Qp.m();
        if(!Qp.m() || !Wp.rho()) {
            p_e_therm = 0.;
        }
        double p_e_pot =
                0.5 * Qp.m() * particles->gas(i)->get_gravitational_potential();
        if(!Qp.m()) {
            p_e_pot = 0.;
        }
#if ndim_ == 3
        Vec mom(Qp.px(), Qp.py(), Qp.pz());
#else
        Vec mom(Qp.px(), Qp.py());
#endif
        double p_e_kin = 0.5 * mom.norm2() / Qp.m();
        if(!Qp.m()) {
            p_e_kin = 0.;
        }
        double p_e_tot = p_e_therm + p_e_pot + p_e_kin;
        e_therm += p_e_therm;
        e_pot += p_e_pot;
        e_kin += p_e_kin;
        e_tot += p_e_tot;
    }
    for(unsigned int i = particles->dmsize(); i--;) {
        ekiny += 0.5 * particles->dm(i)->get_mass() *
                 particles->dm(i)->get_velocity().y() *
                 particles->dm(i)->get_velocity().y();
        double p_e_pot = 0.5 * particles->dm(i)->get_mass() *
                         particles->dm(i)->get_gravitational_potential();
        double p_e_kin = 0.5 * particles->dm(i)->get_mass() *
                         particles->dm(i)->get_velocity().norm2();
        double p_e_tot = p_e_pot + p_e_kin;
        e_pot += p_e_pot;
        e_kin += p_e_kin;
        e_tot += p_e_tot;
    }
    double totenergy_glob, ekiny_glob, e_therm_glob, e_pot_glob, e_kin_glob,
            e_tot_glob;
    // better use Reduce, since we only need to gather information on process 0
    MyMPI_Reduce(&totenergy, &totenergy_glob, 1, MPI_DOUBLE, MPI_SUM, 0);
    MyMPI_Reduce(&ekiny, &ekiny_glob, 1, MPI_DOUBLE, MPI_SUM, 0);
    MyMPI_Reduce(&e_therm, &e_therm_glob, 1, MPI_DOUBLE, MPI_SUM, 0);
    MyMPI_Reduce(&e_pot, &e_pot_glob, 1, MPI_DOUBLE, MPI_SUM, 0);
    MyMPI_Reduce(&e_kin, &e_kin_glob, 1, MPI_DOUBLE, MPI_SUM, 0);
    MyMPI_Reduce(&e_tot, &e_tot_glob, 1, MPI_DOUBLE, MPI_SUM, 0);
    totenergy = totenergy_glob;
    ekiny = ekiny_glob;
    e_therm = e_therm_glob;
    e_pot = e_pot_glob;
    e_kin = e_kin_glob;
    e_tot = e_tot_glob;

    // stat log
    unsigned int numactive = particles->get_numactive();
    unsigned int glob_numactive;
    MyMPI_Reduce(&numactive, &glob_numactive, 1, MPI_UNSIGNED, MPI_SUM, 0);
    numactive = glob_numactive;

    unsigned long solvercount = solver->get_neval() - _stat_solvercount;
    _stat_solvercount = solver->get_neval();
    unsigned long glob_solvercount;
    MyMPI_Reduce(&solvercount, &glob_solvercount, 1, MPI_UNSIGNED_LONG, MPI_SUM,
                 0);
    solvercount = glob_solvercount;

    // only rank 0 writes log files
    if(MPIGlobal::rank) {
        return;
    }

    // time log
    _timefile << timeline->get_time() << "\t"
              << HelperFunctions::human_readable_long(
                         timeline->get_integertime())
              << "\t"
              << timeline->get_realtime_interval(timeline->get_timestep())
              << endl;

    // energy log
    _energyfile << timeline->get_time() << "\t" << totenergy << "\t" << ekiny;
    _energyfile << "\t";
    _energyfile << e_therm << "\t" << e_pot << "\t" << e_kin << "\t" << e_tot;
    // force flush
    _energyfile << endl;

    // stat log
    _statfile << "Step " << _stat_stepcount << "\n";
    _stat_stepcount++;
    _statfile << "System time: " << timeline->get_time() << ", System step: "
              << timeline->get_realtime_interval(timeline->get_timestep())
              << "\n";
    _statfile << HelperFunctions::human_readable_counter(numactive)
              << " active particles\n";
    _statfile << HelperFunctions::human_readable_counter(solvercount)
              << " Riemann solver evaluations\n";
    _statfile << endl;

    _timer.stop();
    LOGS("Written to LogFiles");
}

/**
 * @brief Copy the input file to the output file in a single read/write
 * operation
 *
 * This version of the method uses a single buffer with the size of the block to
 * read. If this size is large, this might cause memory problems. It is advised
 * to use copy_file_buffer in this case.
 *
 * @param ifile Input file
 * @param ofile Output file
 * @param pos Size of the block to copy (from the start of the file)
 */
void LogFiles::copy_file_single(ifstream& ifile, ofstream& ofile, long pos) {
    char* contents = new char[pos];
    ifile.read(contents, pos);
    ofile.write(contents, pos);
    delete[] contents;
}

/**
 * @brief Copy the input file to the output file in a buffered read/write
 * operation
 *
 * This version of the method uses a buffer with a fixed size of 1 MB to read.
 * If the block being read is larger, the file is copied in multiple steps.
 *
 * @param ifile Input file
 * @param ofile Output file
 * @param pos Size of the block to copy (from the start of the file)
 */
void LogFiles::copy_file_buffer(ifstream& ifile, ofstream& ofile, long pos) {
    unsigned int bufsize = 1 << 20;
    char* contents = new char[bufsize];
    unsigned int numread = 0;
    while(numread * bufsize < pos) {
        unsigned int blocksize = pos - numread * bufsize;
        blocksize = std::min(bufsize, blocksize);
        ifile.read(contents, blocksize);
        ofile.write(contents, blocksize);
        numread++;
    }
    delete[] contents;
}

/**
 * @brief Copy part of the given file into the new given file
 *
 * @param iname Name of the input file
 * @param oname Name of the (new) output file
 * @param pos Size of the block (from the start of the file) that will be copied
 */
void LogFiles::copy_file(const string& iname, const string& oname, long pos) {
    ifstream ifile(iname);
    ofstream ofile(oname);
    // if the file is larger than 1 MB, we don't read and write all contents
    // at once
    if(pos > (1 << 20)) {
        copy_file_buffer(ifile, ofile, pos);
    } else {
        copy_file_single(ifile, ofile, pos);
    }
}

/**
 * @brief Reset the file with the given name so that all contents after the
 * given position are removed from it
 *
 * This allows us to reopen a log file as if no data was written after the
 * restart dump was made.
 *
 * @param name Name of the file
 * @param pos Position in the file
 */
void LogFiles::reset_file(string name, long pos) {
    string tempname = name + ".temp";
    copy_file(name, tempname, pos);
    rename(tempname.c_str(), name.c_str());
}

/**
 * @brief Restart constructor. Initialize the log files from the given
 * RestartFile
 *
 * We read the name of the files and open the streams in append mode. We then
 * reset every stream to its old position, meaning we start writing again from
 * the point where the file was when the restart dump was made. This does not
 * erase the remainder of the stream however!
 *
 * @param rfile RestartFile to read from
 */
LogFiles::LogFiles(RestartFile& rfile) : _timer(rfile) {
    if(!MPIGlobal::rank) {
        _timer.start();

        // time log
        rfile.read(_timename);
        long timepos;
        rfile.read(timepos);
        reset_file(_timename, timepos);
        _timefile.open(_timename, ios_base::app);

        // energy log
        rfile.read(_energyname);
        long energypos;
        rfile.read(energypos);
        reset_file(_energyname, energypos);
        _energyfile.open(_energyname, ios_base::app);

        // stat log
        rfile.read(_statname);
        long statpos;
        rfile.read(statpos);
        reset_file(_statname, statpos);
        _statfile.open(_statname, ios_base::app);
        rfile.read(_stat_stepcount);
        rfile.read(_stat_solvercount);

        _timer.stop();

        LOGS("Restarted LogFiles");
    }
}

/**
 * @brief Dump the log files to the given RestartFile
 *
 * We only write the names of the files, the streams themselves are not
 * affected.
 *
 * @param rfile RestartFile to write to
 */
void LogFiles::dump(RestartFile& rfile) {
    _timer.dump(rfile);
    if(!MPIGlobal::rank) {
        // time log
        rfile.write(_timename);
        long timepos = _timefile.tellp();
        rfile.write(timepos);

        // energy log
        rfile.write(_energyname);
        long energypos = _energyfile.tellp();
        rfile.write(energypos);

        // stat log
        rfile.write(_statname);
        long statpos = _statfile.tellp();
        rfile.write(statpos);
        rfile.write(_stat_stepcount);
        rfile.write(_stat_solvercount);

        LOGS("LogFiles dumped");
    }
}
