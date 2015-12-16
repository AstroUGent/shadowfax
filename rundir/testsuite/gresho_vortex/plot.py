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
import h5py
import glob

def get_analytical():
    xs = [0., 0.2, 0.4, 0.6]
    ys = [0., 1., 0., 0.]
    return xs, ys

for fname in sorted(glob.glob("2d/gresho_vortex2d*.hdf5")):
    f = h5py.File(fname, "r")
    v = np.array(f["/PartType0/Velocities"])
    x = np.array(f["/PartType0/Coordinates"])

    radius = np.zeros(len(x))
    velocity = np.zeros(len(x))
    for j in range(len(x)):
        radius[j] = np.sqrt((x[j,0]-0.5)**2 + (x[j,1]-0.5)**2)
        velocity[j] = np.sqrt(v[j,0]**2 + v[j,1]**2)

    xs, vs = get_analytical()
    pl.plot(xs, vs, "r-")
    pl.plot(radius, velocity, "k.")

    pl.savefig("{fname}.png".format(fname = fname[:-5]))
    pl.close()

for fname in sorted(glob.glob("3d/gresho_vortex3d*.hdf5")):
    f = h5py.File(fname, "r")
    v = np.array(f["/PartType0/Velocities"])
    x = np.array(f["/PartType0/Coordinates"])

    radius = np.zeros(len(x))
    velocity = np.zeros(len(x))
    for j in range(len(x)):
        radius[j] = np.sqrt((x[j,0]-0.5)**2 + (x[j,1]-0.5)**2)
        velocity[j] = np.sqrt(v[j,0]**2 + v[j,1]**2)

    xs, vs = get_analytical()
    pl.plot(xs, vs, "r-")
    pl.plot(radius, velocity, "k.")

    pl.savefig("{fname}.png".format(fname = fname[:-5]))
    pl.close()
