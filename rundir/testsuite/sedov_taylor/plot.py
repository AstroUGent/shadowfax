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

import scipy as sp
import pylab as pl
import numpy as np
import scipy.interpolate as interpol
import h5py
import subprocess
import glob
import sedov

for fname in sorted(glob.glob("2d/sedov_taylor2d*.hdf5")):
    f = h5py.File(fname, "r")
    rho = np.array(f["/PartType0/Density"])
    x = np.array(f["/PartType0/Coordinates"])
    t = f["/Header"].attrs["Time"][0]

    radius = np.zeros(len(x))
    velocity = np.zeros(len(x))
    for j in range(len(x)):
        radius[j] = np.sqrt((x[j,0]-0.5)**2 + (x[j,1]-0.5)**2)

    xs, rhos = sedov.get_analytical_solution(1., 5./3., 2, t)
    pl.plot(xs, rhos, "r-")
    pl.plot(radius, rho, "k.")

    pl.savefig("{fname}.png".format(fname = fname[:-5]))
    pl.close()

for fname in sorted(glob.glob("3d/sedov_taylor3d*.hdf5")):
    f = h5py.File(fname, "r")
    rho = np.array(f["/PartType0/Density"])
    x = np.array(f["/PartType0/Coordinates"])
    t = f["/Header"].attrs["Time"][0]

    radius = np.zeros(len(x))
    velocity = np.zeros(len(x))
    for j in range(len(x)):
        radius[j] = np.sqrt((x[j,0]-0.5)**2 + (x[j,1]-0.5)**2 + (x[j,2]-0.5)**2)

    xs, rhos = sedov.get_analytical_solution(1., 5./3., 3, t)
    pl.plot(xs, rhos, "r-")
    pl.plot(radius, rho, "k.")

    pl.savefig("{fname}.png".format(fname = fname[:-5]))
    pl.close()
