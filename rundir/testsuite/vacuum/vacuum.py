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

import numpy as np
import pylab as pl

gamma = 5./3.

def get_solution(WL, WR, x, t):
    if t > 0.:
        u = x/t
    else:
        if x > 0.:
            u = 1.e9
        else:
            u = -1.e9
    uL = WL[1]
    uR = WR[1]
    if WL[0] > 0.:
        aL = np.sqrt(gamma*WL[2]/WL[0])
    else:
        aL = 0.
    if WR[0] > 0.:
        aR = np.sqrt(gamma*WR[2]/WR[0])
    else:
        aR = 0.
    
    if WR[0] == 0.:
        # vacuum right state
        if u <= uL - aL:
            solution = WL
        else:
            SstarL = uL + 2.*aL/(gamma-1.)
            if u >= SstarL:
                solution = np.zeros(3)
            else:
                solution = np.zeros(3)
                fac = 2./(gamma+1.) + (gamma-1.)*(uL-u)/(gamma+1.)/aL
                solution[0] = WL[0]*fac**(2./(gamma-1.))
                solution[1] = 2.*(aL + 0.5*(gamma-1.)*uL + u)/(gamma+1.)
                solution[2] = WL[2]*fac**(2.*gamma/(gamma-1.))
    else:
        if WL[0] == 0.:
            # vacuum left state
            if u >= uR + aR:
                solution = WR
            else:
                SstarR = uR - 2.*aR/(gamma-1.)
                if u <= SstarR:
                    solution = np.zeros(3)
                else:
                    solution = np.zeros(3)
                    fac = 2./(gamma+1.) - (gamma-1.)*(uR-u)/(gamma+1.)/aR
                    solution[0] = WR[0]*fac**(2./(gamma-1.))
                    solution[1] = 2.*(0.5*(gamma-1.)*uR - aR + u)/(gamma+1.)
                    solution[2] = WR[2]*fac**(2.*gamma/(gamma-1.))
        else:
            # vacuum generation
            SstarL = uL + 2.*aL/(gamma-1.)
            SstarR = uR - 2.*aR/(gamma-1.)
            if u <= SstarL:
                if u <= uL - aL:
                    solution = WL
                else:
                    solution = np.zeros(3)
                    fac = 2./(gamma+1.) + (gamma-1.)*(uL-u)/(gamma+1.)/aL
                    solution[0] = WL[0]*fac**(2./(gamma-1.))
                    solution[1] = 2.*(aL + 0.5*(gamma-1.)*uL + u)/(gamma+1.)
                    solution[2] = WL[2]*fac**(2.*gamma/(gamma-1.))
            else:
                if u >= SstarR:
                    if u >= uR + aR:
                        solution = WR
                    else:
                        solution = np.zeros(3)
                        fac = 2./(gamma+1.) - (gamma-1.)*(uR-u)/(gamma+1.)/aR
                        solution[0] = WR[0]*fac**(2./(gamma-1.))
                        solution[1] = 2.*(0.5*(gamma-1.)*uR - aR + u)/(gamma+1.)
                        solution[2] = WR[2]*fac**(2.*gamma/(gamma-1.))
                else:
                    solution = np.zeros(3)
    return solution

def vacuum(WL, WR, t):
    WL = np.array(WL)
    WR = np.array(WR)
    x = np.arange(0.005, 1., 0.01)-0.5
    rho = np.zeros(len(x))
    u = np.zeros(len(x))
    P = np.zeros(len(x))
    for i in range(len(x)):
        Wstar = get_solution(WL, WR, x[i], t)
        rho[i] = Wstar[0]
        u[i] = Wstar[1]
        P[i] = Wstar[2]

    return x+0.5, rho, u, P

WL = np.array([1., -4., 1.])
WR = np.array([1., 4., 1.])

x = np.arange(0.005, 1., 0.01)-0.5
rho = np.zeros(len(x))
count = 0
for t in np.arange(0., 0.1, 0.01):
    for i in range(len(x)):
        rho[i] = get_solution(WL, WR, x[i], t)[0]
    pl.plot(x, rho)
    pl.savefig("snap{0:03d}.png".format(count))
    pl.close()
    count += 1
