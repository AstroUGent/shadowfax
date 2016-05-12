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

import pylab as pl
import numpy as np
import h5py
import glob

nsnap = 101
gamma = 5./3.

time = np.zeros(nsnap)
total_energy = np.zeros(nsnap)
kinetic_energy = np.zeros(nsnap)
thermal_energy = np.zeros(nsnap)
potential_energy = np.zeros(nsnap)
for fname in sorted(glob.glob("3d/evrard3d*.hdf5")):
    print "processing", fname
    f = h5py.File(fname, "r")
    density = np.array(f["/PartType0/Density"])
    x = np.array(f["/PartType0/Coordinates"])

    radius = np.zeros(len(x))
    for j in range(len(x)):
        radius[j] = np.sqrt((x[j,0]-2.)**2 + (x[j,1]-2.)**2 + (x[j,2]-2.)**2)
  
    exact = np.loadtxt("exact3d/profile3d{nr}.txt".format(nr = fname[-8:-5]))
    pl.loglog(exact[:,0], exact[:,1], "r-")
    pl.loglog(radius, density, "k.")
    pl.xlim(0.01, 0.8)
    pl.savefig("{fname}.png".format(fname = fname[:-5]))
    pl.close()
