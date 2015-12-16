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

class Particle:
    def __init__(self, x, y, z, mass):
        self.pos = [x, y, z]
        self.m = mass

def plummer(r):
    return 3.*1000./4./np.pi*(1.+r**2)**(-5./2.)

for fname in sorted(glob.glob("3d/nbody3d*.hdf5")):
    f = h5py.File(fname, "r")
    x = np.array(f["/PartType1/Coordinates"])
    m = np.array(f["/PartType1/Masses"])

    particles = []
    for j in range(len(x)):
        particles.append(Particle(x[j,0], x[j,1], x[j,2], m[j]))

    # calculate center of mass
    com = [0.,0.,0]
    totmass = 0.
    for p in particles:
        for j in range(3):
            com[j] += p.m*p.pos[j]
        totmass += p.m
    for j in range(3):
        com[j] /= totmass

    # calculate radii
    radius = np.zeros(len(particles))
    for i in range(len(particles)):
        p = particles[i]
        radius[i] = np.sqrt((p.pos[0]-com[0])**2 + (p.pos[1]-com[1])**2 + (p.pos[2]-com[2])**2)

    rmin = min(radius)
    rmax = max(radius)
    dr = 0.01*(rmax - rmin)
    radii = np.arange(rmin, rmax + dr, dr)

    # analytical solution
    dens = []
    for r in radii:
        dens.append(plummer(r))
    pl.semilogy(radii, dens)

    # binned simulation result
    densities = np.zeros(len(radii))
    for p in particles:
        r = np.sqrt((p.pos[0]-com[0])**2 + (p.pos[1]-com[1])**2 + (p.pos[2]-com[2])**2)
        index = 0
        while index+1 < len(radii) and r > radii[index+1]:
            index += 1
        densities[index] += p.m
    for j in range(len(radii)):
        densities[j] /= (4.*np.pi/3.*(3.*radii[j]*dr**2 + 3.*radii[j]**2*dr + dr**3))
    pl.semilogy(radii, densities, "o")

    pl.savefig("{fname}.png".format(fname = fname[:-5]))
    pl.close()
