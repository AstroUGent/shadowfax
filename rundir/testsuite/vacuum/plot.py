################################################################################
# This file is part of Shadowfax
# Copyright (C) 2015 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
#
# Shadowfax is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Shadowfax is distributed in the hope that it will be useful,
# but WITOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with Shadowfax. If not, see <http://www.gnu.org/licenses/>.
################################################################################

import h5py
import numpy as np
import pylab as pl
import glob
import vacuum

for f in sorted(glob.glob("2d/vacuum2d*.hdf5")):
    print "processing", f[:-5]
    file = h5py.File(f, "r")
    coords = np.array(file["/PartType0/Coordinates"])
    rho = np.array(file["/PartType0/Density"])
    vs = np.array(file["/PartType0/Velocities"])
    P = np.array(file["/PartType0/InternalEnergy"])
    P = 2.*rho*P/3.
    t = file["/Header"].attrs["Time"]
    
    xsol, rhosol, usol, Psol = vacuum.vacuum([1., 0., 1.], [0., 0., 0.], t)
    
    fig, ax = pl.subplots(2, 2)
    ax[0][0].plot(xsol, rhosol, "r-")
    ax[0][0].plot(coords[:,0], rho, "k.")
    ax[0][1].plot(xsol, usol, "r-")
    ax[0][1].plot(coords[:,0], vs[:,0], "k.")
    ax[1][0].plot(xsol, Psol, "r-")
    ax[1][0].plot(coords[:,0], P, "k.")
    pl.savefig("{name}.png".format(name=f[:-5]))
    pl.close()

for f in sorted(glob.glob("3d/vacuum3d*.hdf5")):
    print "processing", f[:-5]
    file = h5py.File(f, "r")
    coords = np.array(file["/PartType0/Coordinates"])
    rho = np.array(file["/PartType0/Density"])
    vs = np.array(file["/PartType0/Velocities"])
    P = np.array(file["/PartType0/InternalEnergy"])
    P = 2.*rho*P/3.
    t = file["/Header"].attrs["Time"]
    
    xsol, rhosol, usol, Psol = vacuum.vacuum([1., -4., 1.], [1., 4., 1.], t)
    
    fig, ax = pl.subplots(2, 2)
    ax[0][0].plot(xsol, rhosol, "r-")
    ax[0][0].plot(coords[:,0], rho, "k.")
    ax[0][1].plot(xsol, usol, "r-")
    ax[0][1].plot(coords[:,0], vs[:,0], "k.")
    ax[1][0].plot(xsol, Psol, "r-")
    ax[1][0].plot(coords[:,0], P, "k.")
    pl.savefig("{name}.png".format(name=f[:-5]))
    pl.close()
