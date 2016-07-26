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

def f(x1, x2, x3, x4, x5):
    return 10.**(x1+x2+x3+x4+x5)

def g(x1, x2, x3):
    return x3

for i in xrange(0, 3):
    for j in xrange(0, 3):
        for k in xrange(0, 3):
             for l in xrange(0, 3):
                file = open("coolingtables/RadLoss_{i}_{j}_{k}_{g}.rates"
                            .format(i = i, j = j, k = k, g = g(i, j, l)), "w")
                file.write("{l}\n".format(l = l))
                file.write("bla\n")
                for m in xrange(0, 3):
                    file.write("{m}\t{f}\n".format(m = m, f = f(i, j, k, l, m)))
                file.close()
