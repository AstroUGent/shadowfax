= The Sedov-Taylor shockwave =

This test problem consists of a box filled with a homogeneous fluid at zero or
very low temperature. At time zero, a large amount of energy is inserted into
the center (0.5, 0.5[, 0.5]) of the box. This causes a strong shock wave to
emerge and spread radially outward.

The initial hydrodynamical values in the box are
rho = 1., u = 0., v = 0.[, w = 0.], p = 1.e-6
and an energy E0 = 1. is inserted into the center.

Taylor (1950) and Sedov (1959) showed that the resulting system is self-similar
and that the position of the shock wave as a function of time is given by
rs = (E/rho)^0.2 t^0.4 for spherical symmetry and
rs = (E/rho)^0.25 t^0.5 for cylindrical symmetry
where E = E0/alpha, with alpha some constant which can be calculated.

The setup as given here is based on the problem presented in Springel and
Hernquist (2002).

This problem is very demanding for the hydro solver and is also a good test
for the vacuum Riemann solver. The two dimensional case currently only works
well with cartesian initial conditions. The three dimensional case can handle
random setups.

To run the problem, just run:
bash run.sh

Results can be plotted together with the analytical solution using
python plot.py

== References ==

Sedov, L. 1977, Similitude et dimensions en mecanique, 7th edn. (Editions Mir)

Springel, V., & Hernquist, L. 2002, MNRAS, 333, 649

Taylor, G. 1950, Royal Society of London Proceedings A, 201, 159

-- Written by Bert Vandenbroucke, Tuesday 12 August 2014 --
