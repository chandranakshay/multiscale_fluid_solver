# multiscale_solver
Fluid kinetic multiscale solver

This C++ code is a coupled version of the Lattice-Boltzmann (LB) RD3Q41 model and Direct Simulation Monte Carlo (DSMC) fluid solvers.

Some portions of the code are SIMD vectorized and require avx support.

Physical domain layout: The code is written only for the channel flow geometry (default setup: gravity driven flow) with the DSMC region restricted to tunable equal
widths from either walls. The LB solver governs the rest of the domain.

Running a simulation: requires you to setup nondimensional parameters, number of DSMC cells and equivalent LB nodes in the domain, number of particles per cell, 
the buffer layer distance, and perlin-noise parameters (if required) in ./DSMC_Solver/main.C.
