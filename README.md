# max-plus-algorithms

This repository contains MATLAB implementations for the analysis of  
**max-plus linear systems**, with a focus on spectral properties and transient
behavior.

Max-plus algebra replaces classical addition and multiplication by:

- **Addition:**  
  a ⊕ b = max(a, b)

- **Multiplication:**  
  a ⊗ b = a + b

Such models arise naturally in **discrete-event systems**, **scheduling**,
**transportation networks**, and **timed Petri nets**.

---

## Repository Contents

### `howard.m`

Implementation of **Howard’s policy iteration algorithm** for computing a
**generalised max-plus eigenmode**.

Given a heigher order max-plus system, the algorithm computes:

- the **max-plus eigenvalue** (cycle mean),
- the associated **eigenvector**,
- a **critical node** belonging to a critical cycle.

---

### `lea.m`

Monte Carlo algorithm for estimating **Lyapunov exponent** of
max-plus linear systems.

The algorithm:

- applies randomness to a max-plus matrix,
- propagates two trajectories with a fixed time shift that represent lower and upper bound of Lyapunov exponent,
- estimates their mean values to get bounds for Lyapunov exponent,
- computes **95% confidence intervals** from Monte Carlo samples.

In this algorithm we assume that the delays added to the random matrix A with some probability and are exponetialy distributed. This can easily changed to any other distribution.

---

### `example_lea.m`

Example script demonstrating the use of `lea.m` for a simple timed events graph visible in the picture below. The net has four transitions and eight places. We can describe this net with a recursive equation of heigher order, but this can be further described as a equation of first order with a single matrix A, that you can found in `example_lea.m`.

The script:

1. Defines a max-plus system matrix
2. Sets algorithm parameters
3. Runs the Monte Carlo estimation
4. Returns deviation bounds and confidence intervals

![Timed event graph](docs/example_net.svg)

---

## Requirements

- MATLAB  
- No external toolboxes required

---

## References

- B. Heidergott, G. J. Olsder, J. van der Woude,  
  *Max Plus at Work*, Princeton University Press, 2006.

- R. M. P. Goverde, B. Heidergott, G. Merlet,  
  *A coupling approach to estimating the Lyapunov exponent of stochastic  
  max-plus linear systems*.

---

## License

This code is provided for research and educational purposes.
