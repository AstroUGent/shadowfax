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
 * @file RestartFile.cpp
 *
 * @brief Restart file support: implementation
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#include "RestartFile.hpp"
#include "../inc/ShadowfaxVersion.hpp"
#include "Error.hpp"
#include "MPIGlobal.hpp"
#include "MPIMethods.hpp"
#include <cstring>
#include <dirent.h>
#include <sstream>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/utsname.h>
using namespace std;

/**
 * @brief Get a zero padded string containing the given digit, with a maximal
 * capacity of limit-1
 *
 * @param digit Digit that needs to be zero padded
 * @param limit Maximal number +1 that should be representable, this sets the
 * number of digits in the string
 * @return String containing the digit with extra zeros added
 */
string RestartFile::get_padded_digit(int digit, int limit) {
    stringstream name;
    // if limit is 1, the procedure below will set the width to zero, which
    // results in strange behaviour
    // in this case, we can just add the digit to the stringstream without
    // padding
    if(limit > 1) {
        int digits = 0;
        int limitcopy = limit;
        int check = 1;
        while(limitcopy) {
            limitcopy /= 10;
            check *= 10;
            digits++;
        }
        check /= 10;
        if(check == limit) {
            digits--;
        }
        name.fill('0');
        name.width(digits);
    }
    name << digit;
    return name.str();
}

/**
 * @brief Get a filename containing the given prefix and the given rank, zero
 * padded to represent up to size-1 different ranks with the same number of
 * digits
 *
 * @param prefix Prefix of the filename
 * @param rank Rank to add to the filename
 * @param size Maximal size of the rank that should be representable with the
 * same number of digits
 * @param suffix Suffix to add to the end of the filename
 * @return Filename
 */
string RestartFile::get_filename(string prefix, int rank, int size,
                                 string suffix) {
    return prefix + string(".") + get_padded_digit(rank, size) + suffix;
}

/**
 * @brief Get a folder name containing the given prefix and the given rank, zero
 * padded to represent up to size-1 different ranks with the same number of
 * digits
 *
 * @param prefix Prefix of the folder name
 * @param rank Rank to add to the folder name
 * @param size Maximal size of the rank that should be representable with the
 * same number of digits
 * @param suffix Suffix to add to the end of the folder name
 * @return Folder name
 */
string RestartFile::get_foldername(string prefix, int rank, int size,
                                   string suffix) {
    return prefix + string("_") + get_padded_digit(rank, size) + suffix;
}

/**
 * @brief Dump constructor
 *
 * The process with local rank 0 makes a backup of existing restart files and
 * then creates an ASCII-file restart.noderank.txt which contains some system
 * information and information about the restart files.
 *
 * Every process then opens a binary file restart.rank.dat for writing.
 *
 * @param outputdir Folder where the restart files will be written. The general
 * file is in this folder, the actual restart files are in a subfolder of this
 * folder
 * @param simtime Current time of the simulation, written to the general file
 */
RestartFile::RestartFile(std::string outputdir, double simtime) {
    if(!MPIGlobal::local_rank) {
        string restartfile =
                get_filename(outputdir + string("/restart"),
                             MPIGlobal::noderank, MPIGlobal::nodesize, ".txt");
        string restartfileback = get_filename(outputdir + string("/restart"),
                                              MPIGlobal::noderank,
                                              MPIGlobal::nodesize, ".txt.back");
        string restartfolder =
                get_foldername(outputdir + string("/restart_files"),
                               MPIGlobal::noderank, MPIGlobal::nodesize);
        string restartfolderback = get_foldername(
                outputdir + string("/restart_files"), MPIGlobal::noderank,
                MPIGlobal::nodesize, "_back");

        // backup existing restart file
        ifstream file(restartfile.c_str());
        // check if the file exists
        if(file) {
            rename(restartfile.c_str(), restartfileback.c_str());

            struct stat info;
            // check if the folder exists
            if(stat(restartfolder.c_str(), &info) != 0) {
                cerr << "Error! Restart file found, but folder is missing!"
                     << endl;
                my_exit();
            }

            // check if there is already a backed up folder
            // if so: delete its contents, since otherwise rename won't work
            DIR* dp;
            dp = opendir(restartfolderback.c_str());
            if(dp) {
                // we have to delete every file in the folder
                // this is accomplished using the ancient code below
                struct dirent* dirp;
                while((dirp = readdir(dp))) {
                    string filename(dirp->d_name);
                    if(filename != "." && filename != "..") {
                        filename = restartfolderback + string("/") + filename;
                        remove(filename.c_str());
                    }
                }
                closedir(dp);
            }
            // rename the folder. If the folder already exists, it is
            // deleted first. Since this only works for empty folders,
            // rename also only works for empty folders.
            rename(restartfolder.c_str(), restartfolderback.c_str());
        }

        // create restart folder
        mkdir(restartfolder.c_str(), S_IRWXU);
        // create restart file
        ofstream rfile(restartfile.c_str());
        rfile << "# Restart file #\n\n";

        rfile << "## System info ##\n";
        // retrieve system information
        struct utsname osinfo;
        uname(&osinfo);
        rfile << "System name: " << osinfo.sysname << " " << osinfo.machine
              << " (" << osinfo.release << ")\n";
        rfile << "Kernel version: " << osinfo.version << "\n";
        rfile << "Node name: " << osinfo.nodename << "\n";
        rfile << "Number of local processes: " << MPIGlobal::local_size << "\n";
        rfile << "Total number of processes: " << MPIGlobal::size << "\n\n";

        rfile << "## Code info ##\n";
        rfile << "Version number: " << GIT_BUILD_STRING << "\n";
        rfile << "Code compiled on " << COMPILE_DATE << ", " << COMPILE_TIME
              << "\n\n";

        rfile << "## Restart info ##\n";
        timeval time;
        // localtime returns a pointer which apparently does not need to be
        // freed
        // (doing so gives an error and is not allowed, see
        // http://stackoverflow.com/questions/8694365/how-is-the-result-struct-
        // of-localtime-allocated-in-c)
        struct tm* date;
        // get the local time in seconds since the epoch (with microsecond
        // precision)
        gettimeofday(&time, NULL);
        // convert it to a more extended format including a date
        date = localtime(&time.tv_sec);
        // convert the date to a string format
        // 22: 3 (month) + 2 (day) + 4 (year) + 2 (hour) + 2 (minutes)
        //      + 2 (seconds) + 3 spaces + 1 comma + 2 colons + '\0'
        char buffer[22];
        strftime(buffer, 22, "%b %e %Y, %H:%M:%S", date);
        rfile << "Restart directory: " << restartfolder << "\n";
        rfile << "Restart file written on " << buffer << "\n";
        rfile << "Simulation time: " << simtime << "\n";
    }
    // make sure the folder exists before trying to write new restart files
    MyMPI_Barrier();

    string restartfolder =
            get_foldername(outputdir + string("/restart_files"),
                           MPIGlobal::noderank, MPIGlobal::nodesize);
    string restartname =
            get_filename(restartfolder + string("/restart"),
                         MPIGlobal::local_rank, MPIGlobal::local_size, ".dat");
    _ofile.open(restartname.c_str(), ios::out | ios::binary);
}

/**
 * @brief Restart constructor
 *
 * The process with rank 0 opens the general restart file with given name and
 * reads in some system information. This information is then checked with
 * information of the current system to see whether a restart is possible.
 * The general file also contains the name of the restart folder, which is read
 * in by the process with rank 0 and then communicated to the other processes.
 *
 * All processes then open the binary file restart.rank.dat in this folder for
 * reading.
 *
 * @param filename
 */
RestartFile::RestartFile(std::string filename) {
    char filenamebuffer[100];
    if(!MPIGlobal::local_rank) {
        string restartfilename = get_filename(filename, MPIGlobal::noderank,
                                              MPIGlobal::nodesize, ".txt");
        ifstream rfile(restartfilename.c_str());
        if(!rfile.is_open()) {
            cerr << "Cannot open file \"" << restartfilename << "\"!" << endl;
            my_exit();
        }

        // read general overview
        string line;
        // ignore first 3 lines
        getline(rfile, line);
        getline(rfile, line);
        getline(rfile, line);

        // retrieve system information
        cout << "## System info ##" << endl;
        char sysnamebuf[100];
        char machinebuf[100];
        char releasebuf[100];
        char nodenamebuf[100];
        string sysname;
        string machine;
        string release;
        string version;
        string nodename;
        int mpisize;
        int local_size;

        getline(rfile, line);
        sscanf(line.c_str(), "System name: %s %s (%s", sysnamebuf, machinebuf,
               releasebuf);
        // remove the trailing ) from syskernel
        releasebuf[strlen(releasebuf) - 1] = '\0';
        sysname = string(sysnamebuf);
        machine = string(machinebuf);
        release = string(releasebuf);
        cout << "System name: " << sysname << " " << machine << " (" << release
             << ")" << endl;

        getline(rfile, line);
        version = line.substr(16);
        cout << "Kernel version: " << version << endl;

        getline(rfile, line);
        sscanf(line.c_str(), "Node name: %s", nodenamebuf);
        nodename = string(nodenamebuf);
        cout << "Node name: " << nodename << endl;

        getline(rfile, line);
        sscanf(line.c_str(), "Number of local processes: %i", &local_size);
        cout << "Number of local processes: " << local_size << endl;

        getline(rfile, line);
        sscanf(line.c_str(), "Total number of processes: %i", &mpisize);
        cout << "Total number of processes: " << mpisize << endl;
        cout << endl;

        // check if system is compatible
        struct utsname osinfo;
        uname(&osinfo);
        bool check = string(osinfo.sysname) == sysname;
        check &= string(osinfo.machine) == machine;
        check &= string(osinfo.release) == release;
        check &= string(osinfo.version) == version;
        check &= string(osinfo.nodename) == nodename;
        check &= MPIGlobal::local_size == local_size;
        check &= MPIGlobal::size == mpisize;
        if(!check) {
            cerr << "Cannot restart, systems are incompatible!" << endl;
            my_exit();
        }

        // skip two lines
        getline(rfile, line);
        getline(rfile, line);

        // retrieve code information
        cout << "## Code info ##" << endl;
        string codeversion;
        getline(rfile, line);
        codeversion = line.substr(16);
        cout << "Version number: " << codeversion << "\n";
        string compiletime;
        getline(rfile, line);
        compiletime = line.substr(17);
        cout << "Code compiled on " << compiletime << endl;
        cout << endl;

        // check if version is compatible
        if(codeversion != GIT_BUILD_STRING) {
            cerr << "Cannot restart, code versions are incompatible!" << endl;
            my_exit();
        }

        // check if code time is compatible
        stringstream compilationstring;
        // not sure if this date and time will be the same as the ones in the
        // method below... (there might possibly be a difference of a few
        // milliseconds)
        compilationstring << COMPILE_DATE << ", " << COMPILE_TIME;
        if(compilationstring.str() != compiletime) {
            cout << "Warning: code compilation time does not match with the "
                    "code used in the restart file!"
                 << endl;
            cout << "This code: " << compilationstring.str() << endl;
            cout << "Restart file: " << compiletime << endl;
            cout << endl;
        }

        // skip two lines
        getline(rfile, line);
        getline(rfile, line);

        // retrieve restart information
        cout << "## Restart info ##" << endl;

        getline(rfile, line);
        sscanf(line.c_str(), "Restart directory: %s", filenamebuffer);
        cout << "Restart directory: " << filenamebuffer << endl;

        getline(rfile, line);
        cout << line << endl;

        getline(rfile, line);
        cout << line << endl;
        cout << endl;
    }
    // wait for process 0 and broadcast restart directory
    MyMPI_Bcast(filenamebuffer, 100, MPI_CHAR, 0, MPIGlobal::nodecomm);

    // everything is ok, we start reading the restart file
    string restartfolder(filenamebuffer);
    string restartname =
            get_filename(restartfolder + string("/restart"),
                         MPIGlobal::local_rank, MPIGlobal::local_size, ".dat");
    cout << "Reading " << restartname << endl;
    _ifile.open(restartname.c_str(), ios::in | ios::binary);
}
