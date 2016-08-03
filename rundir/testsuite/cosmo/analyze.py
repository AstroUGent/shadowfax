#! /usr/bin/python

################################################################################
# This file is part of Shadowfax
# Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

import numpy as np
import h5py
import time

# snapshot to read
filename = "3d/cosmo3d010.hdf5"
# expected size of largest group
maxgroup_expected = 8
# expected number of particles in groups
numgroup_expected = 204

# Particle. Holds a position and a reference to the group the particle is in
# (if any)
class Particle:
    def __init__(self, pos):
        self._pos = pos
        self._group = None

# function used to calculate the Cartesian distance between two particles
# we use a dirty global variable to wrap around the periodic boundaries
def distance(part1, part2):
    global box
    p1 = part1._pos
    p2 = part2._pos
    r = 0.
    for i in range(3):
        dx = p1[i] - p2[i]
        if dx < -0.5*box[i]:
            dx += box[i]
        if dx > 0.5*box[i]:
            dx -= box[i]
        r += dx**2
    return np.sqrt(r)

# group of particles that are close enough together
# stores a list of member particles and a flag used to determine if the group
# is already in the final group list we construct
class Group:
    def __init__(self, first_member):
        self._members = [first_member]
        self._in_list = False

    def add_member(self, member):
        self._members.append(member)

    # merge this group with another group
    def join(self, group):
        for member in group._members:
            if not member._group == self:
                self.add_member(member)
                member._group = self

# group 2 particles together
# each particle is at most part of one group
# if none of the particles is already part of a group, we create a new group
# for the two particles
# if one of them is, we simply add the other particle to that group
# if both are already part of a group, we need to check if this group is the
# same. If not, we merge both groups.
def group(p1, p2):
    if p1._group:
        if p2._group:
            if not p1._group == p2._group:
                p1._group.join(p2._group)
        else:
            p2._group = p1._group
            p1._group.add_member(p2)
    else:
        if p2._group:
            p1._group = p2._group
            p2._group.add_member(p1)
        else:
            p1._group = Group(p1)
            p2._group = p1._group
            p1._group.add_member(p2)

# Block structure used to speed up neighbour searching
# the entire simulation box is subdivided into a 3D matrix of blocks with equal
# sizes. Every block contains a list of the particles that are inside that
# particular block. We make sure the typical size of a block is significantly
# larger than our search radius, so that we only need to perform our search
# algorithm on individual blocks and on pairs of neighbouring blocks. This
# technique is very similar to what is used in SWIFT.
class Block:
    def __init__(self):
        self._particles = []

    def add(self, p):
        self._particles.append(p)

    # interact this block with another block (which can be self in case of a
    # self-interaction)
    def interact(self, block):
        if block == self:
            # self-interaction
            for i1 in range(len(self._particles)):
                p1 = self._particles[i1]
                for i2 in range(i1+1, len(self._particles)):
                    p2 = self._particles[i2]
                    if distance(p1, p2) < 0.2*ravg:
                        group(p1, p2)
        else:
            # neighbour interaction
            for i1 in range(len(self._particles)):
                p1 = self._particles[i1]
                for i2 in range(len(block._particles)):
                    p2 = block._particles[i2]
                    if distance(p1, p2) < 0.2*ravg:
                        group(p1, p2)

# open the last snapshot and read the coordinates
file = h5py.File(filename, "r")
x = np.array(file["/PartType1/Coordinates"])

# read the boxsize into the global variable
box = file["/Header"].attrs["BoxSize"]

# create and fill the blocks
# since our search radius is 0.2*boxlength/16, 16 blocks is both safe and fast
numblocks = 16

blocks = []
for i in range(numblocks):
    blocksrow = []
    for j in range(numblocks):
        blockscol = []
        for k in range(numblocks):
            blockscol.append(Block())
        blocksrow.append(blockscol)
    blocks.append(blocksrow)

for ix in range(len(x)):
    i = int(x[ix][0]*numblocks/box[0])
    j = int(x[ix][1]*numblocks/box[1])
    k = int(x[ix][2]*numblocks/box[2])
    blocks[i][j][k].add(Particle(x[ix]))

# the average distance between particles is set by the volume of the box and the
# number of particles
ravg = (box[0]*box[1]*box[2]/len(x))**(1./3.)

# time everything
# on my system, a direct computation takes ~50s. The block computation only
# takes 0.4s
start_time = time.time()

# interact the blocks
# each block has one self interaction, and 7 interactions with neighbour blocks
# a block typically only holds a single particle, so although this code looks
# heavy, it is actually very fast
# note that we only interact each pair of blocks once. This is why we do not
# interact with e.g. blocks[i-1][j][k]
for i in range(len(blocks)):
    # get the (periodically) next block in the first direction
    iplus = (i+1)%len(blocks)
    for j in range(len(blocks[i])):
        # get the (periodically) next block in the second direction
        jplus = (j+1)%len(blocks[i])
        for k in range(len(blocks[i][j])):
            # get the (periodically) next block in the third direction
            kplus = (k+1)%len(blocks[i][j])
            blocks[i][j][k].interact(blocks[i][j][k])
            
            blocks[i][j][k].interact(blocks[iplus][j][k])
            blocks[i][j][k].interact(blocks[i][jplus][k])
            blocks[i][j][k].interact(blocks[i][j][kplus])
            
            blocks[i][j][k].interact(blocks[iplus][jplus][k])
            blocks[i][j][k].interact(blocks[iplus][j][kplus])
            blocks[i][j][k].interact(blocks[i][jplus][kplus])
            
            blocks[i][j][k].interact(blocks[iplus][jplus][kplus])

# construct the final list of groups by looping over all particles and adding
# each group once
groups = []
for i in range(len(blocks)):
    for j in range(len(blocks[i])):
        for k in range(len(blocks[i][j])):
          for p in blocks[i][j][k]._particles:
              if p._group:
                  if not p._group._in_list:
                      groups.append(p._group)
                      p._group._in_list = True

# get the largest group
if len(groups) > 0:
    maxgroup = groups[0]
    totpart = 0
    for group in groups:
        totpart += len(group._members)
        if len(group._members) > len(maxgroup._members):
            maxgroup = group
    maxgroup_members = len(maxgroup._members)
else:
    maxgroup_members = 0
    totpart = 0

end_time = time.time()
print (end_time-start_time), "seconds"

print "Maximum group size:", maxgroup_members
print "Expected maximum group size:", maxgroup_expected
if maxgroup_expected == maxgroup_members:
    print "\033[92mPASS\033[0m"
else:
    print "\033[91mFAIL\033[0m"
print "Total number of particles in groups:", totpart
print "Expected number of particles in groups:", numgroup_expected
if numgroup_expected == totpart:
    print "\033[92mPASS\033[0m"
else:
    print "\033[91mFAIL\033[0m"
