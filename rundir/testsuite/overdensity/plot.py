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

for i in range(11):
  f = h5py.File("2d/overdensity2d{0:03d}.hdf5".format(i), "r")
  density = np.array(f["/PartType0/Density"])
  coords = np.array(f["/PartType0/Coordinates"])
  
  radius = np.zeros(len(density))
  for j in range(len(density)):
    radius[j] = np.sqrt((coords[j,0]-0.5)**2 + (coords[j,1]-0.5)**2)
    
  exact = sp.loadtxt("exact2d/csnap{0:03d}.dat".format(i))
  pl.plot(exact[:,0], exact[:,1], "r-")
  pl.plot(radius, density, "k.")
  
  pl.savefig("overdensity2d{0:03d}.png".format(i))
  pl.close()
  
for i in range(11):
  f = h5py.File("3d/overdensity3d{0:03d}.hdf5".format(i), "r")
  density = np.array(f["/PartType0/Density"])
  coords = np.array(f["/PartType0/Coordinates"])
  
  radius = np.zeros(len(density))
  for j in range(len(density)):
    radius[j] = np.sqrt((coords[j,0]-0.5)**2 + (coords[j,1]-0.5)**2 + (coords[j,2]-0.5)**2)
    
  exact = sp.loadtxt("exact3d/csnap{0:03d}.dat".format(i))
  pl.plot(exact[:,0], exact[:,1], "r-")
  pl.plot(radius, density, "k.")
  
  pl.savefig("overdensity3d{0:03d}.png".format(i))
  pl.close()
