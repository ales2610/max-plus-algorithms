# max-plus-algorithms

This repository contains MATLAB implementations of algorithms for the analysis of max-plus linear systems, with a focus on spectral properties and transient behaviour.

Max-plus algebra replaces classical addition and multiplication by:

Addition: a ⊕ b = max(a,b)

Multiplication: a ⊗ b = a + b

Such models naturally arise in discrete-event systems, scheduling, transportation networks, and timed Petri nets.

Files in this repository
howard.m

Implementation of Howard’s policy iteration algorithm for computing a generalised max-plus eigenmode.

Given a weighted max-plus system with multiple durations, the algorithm computes:

the max-plus eigenvalue (cycle mean),

the associated eigenvector,

a critical node belonging to a critical cycle.

This corresponds to solving a mean-payoff / max-plus spectral problem.

lea.m

Monte Carlo algorithm for estimating transient deviation bounds in max-plus linear systems.

The method:

applies random perturbations to a max-plus matrix,

compares two trajectories with a time shift,

estimates lower and upper bounds on their deviation,

provides 95% confidence intervals.

This approach is related to Lyapunov-type and coupling-time analyses for stochastic max-plus systems.

example_lea.m

Example script illustrating how to use lea.m.

The script:

defines a max-plus system matrix,

sets algorithm parameters,

runs the Monte Carlo estimation,

returns deviation bounds and confidence intervals.

This file serves as a minimal working example.

Requirements

MATLAB (no external toolboxes required)

References

F. Baccelli, G. Cohen, G. J. Olsder, J.-P. Quadrat, Synchronization and Linearity, Wiley.

Mean-payoff games and max-plus spectral theory literature.
