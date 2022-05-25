# multiscale_solver
Fluid kinetic multiscale solver

This C++ code is a coupled version of the Lattice-Boltzmann (LB) RD3Q41 model and Direct Simulation Monte Carlo (DSMC) fluid solvers.

Some portions of the code are SIMD vectorized and require avx support.

Physical domain layout: The code is written only for the channel flow geometry with the DSMC region restricted to tunable equal
widths from either walls. The LB solver governs the rest of the domain.
