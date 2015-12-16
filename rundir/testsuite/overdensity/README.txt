= The overdensity test =

The overdensity test consists of cubic box which contains a spherical
overdensity. The sperical region has radius 0.25, center (0.5,0.5[,0.5]) and
following hydrodynamical values:
rho1 = 1., u1 = 0., v1 = 0.[, w1 = 0.], p1 = 1.
The outside region has following hydrodynamical values:
rho2 = 0.125, u2 = 0., v2 = 0.[, w2 = 0.], p2 = 0.1

The solution consists of an inward travelling rarefaction wave and an outward
travelling shockwave, separated by a contact discontinuity.

This problem is a good basic test for the Riemann solver, since it involves
all three wave patterns present in the general solution of the Riemann problem.

It also makes a good converence and scaling test, since the results are easily
compared to an "exact" solution (obtained by a high resolution 1D equivalent
solver) and good convergence is reached with limited computational resources.

To run the basic problem, run the command:
bash run.sh

Results can be plotted using the provided python script:
python plot.py

The initial condition was based upon Test 1 on page 129 of Toro, E.F. 2009,
Riemann Solvers and Numerical Methods for Fluid Dynamics, 3rd edn. (Springer-
Verlag), e-ISBN 978-3-540-49834-6.
The two- and three-dimensional case are described in chapter 17 of the same
book.

-- Written by Bert Vandenbroucke, Thursday 3 December 2015 --
