Revisions:
==========
Jan 1 2008, Original Version

Contact Person:
===============
Federico Halpern, feh2@lehigh.edu

General comments:
=================
The code contained in this package is a version of the Weiland model
current as of Nov 2006. The computation of diffusion coefficients
dates from Oct 2007. The code was written in free-form style by 
F. Halpern, based upon the original F77 code by J.Weiland.

The code has been built to compute diffusivities for a single point.

This package may be revised in the future to output the drift-wave
growth rates, and for numerical stability and speed. Also, there 
may be issues in the computation of toroidal momentum diffusivities.

However, the thermal diffusivities have been satisfactorily 
compared with the original version of the code.

Contents of package:
====================
w20_README.txt
w20data.f90
w20main.f90
w20wguess.f90
w20disp.f90
w20diff.f90
r8tomsqz.f
aitken.f90

Input parameters:
=================


Integer i_print	 Switch for text output. Set this variable to 9 to obtain
                 printout of the dispersion relation at each iteration

Real*8  z_te	 Electron temperature in keV
Real*8  z_ne	 Electron density in m^-3
Real*8  z_vtor   Toroidal velocity in m/s
Real*8  z_vpol   Poloidal velocity in m/s
Real*8  z_btor   Toroidal magnetic field in T
Real*8  z_rmaj   Local major radius in m
Real*8  z_eps    Aspect ratio a/R
Real*8  z_aimp   Impurity mass in a.m.u
Real*8  z_ahyd   Hydrogenic mass in a.m.u
Real*8  z_zimp   Impurity charge number Z
Real*8  z_gte    Normalized electron temperature gradient (R/Te)(dTe/dr)
                 Derivative is taken respect to minor radius r (not rho)
Real*8  z_gti    Normalized hydrogenic ion temperature gradient (R/Th)(dTh/dr)
Real*8  z_gtz    Normalized impurity ion temperature gradient (R/Tz)(dTz/dr)
Real*8  z_gne    Normalized electron density gradient (R/ne)(dne/dr)
Real*8  z_gni    Normalized hydrogenic ion density gradient (R/nh)(dnh/dr)
Real*8  z_gnz    Normalized impurity ion density gradient (R/nz)(dnz/dr)
Real*8  z_gvt    Normalized toroidal velocity gradient (R/vphi)(dvphi/dr)
Real*8  z_gvp    Normalized poloidal velocity gradient (R/vtht)(dvtht/dr)
Real*8  z_kyrho  ( k_y * rho_s ), should be 0.316 as used by Weiland
Real*8  z_tauh   Ratio of hydrogenic ion temperature to electron temperature Ti/Te
Real*8  z_tauz   Ratio of impurity ion temperature to electron temperature Tz/Te
Real*8  z_nz_ne  Ratio of impurity ion density to electron density nz/ne
Real*8  z_ns_ne  Ratio of fast ion density to electron density ns/ne
Real*8  z_fte    Fraction of trapped electrons
Real*8  z_q      Magnetic safety factor q
Real*8  z_shear  Magnetic shear (r/q)(dq/dr)
Real*8  z_kappa  Local elongation kappa
Real*8  z_wexb   ExB shearing rates in 1/s

Output
======

Real*8, dimension(6) diffs  Diffusion coefficients in m^2/s
Real*8, dimension(6) vconv  Convective velocities in m/s 

Transport channels:
diffs(1): Ion thermal diffusivity
diffs(2): Hydrogenic particle diffusivity
diffs(3): Electron thermal diffusivity
diffs(4): Impurity particle diffusivity
diffs(5): Toroidal momentum diffusivity
diffs(6): Poloidal momentum diffusivity

WARNING:
========
The convective velocities are still in experimental stage as far
as coding goes. Their use is highly discouraged, especially for
the momentum transport, as they should be included in diffs(5)
and diffs(6) already.