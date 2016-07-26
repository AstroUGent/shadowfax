# -*- coding: utf-8 -*-

################################################################################
# This file is part of Shadowfax
# Copyright (C) 2016 Yorick Van Den Bossche (yorick.vandenbossche@ugent.be)
#                    Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

"""
Created on Wed Feb 10 16:12:03 2016

@author: Yorick Van Den Bossche (yorick.vandenbossche@ugent.be)
"""
import numpy as np

Fe_vals = np.array([-99., -0.5, 0.5])
Mg_vals = np.array([0., 0.47, -0.55])
kperm = 8248.98310351/(4./(1.+3.*0.76))
rho = 1.e-19
def cool(T, nH):
    return 150.*kperm*T*(10.**nH)

def getrho(nH):
    return nH

for Fe in Fe_vals:
    for Mg in Mg_vals:
        for z in np.arange(0., 12., 6.):
              for nH in np.arange(-24., -16., 4.):
                file = open("coolingtables/RadLoss_{Fe}_{Mg}_{z}_{nH}.rates"
                            .format(Fe = Fe, Mg = Mg, z = z, nH = 10.**nH), "w")
                file.write("{rho}\n".format(rho = getrho(10.**nH)))
                file.write("bla\n")
                for T in np.logspace(1., 5., num = 10):
                    file.write("{T}\t{cool}\n".format(T = T,
                                                      cool = cool(T, nH)))
                file.close()
