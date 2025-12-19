# max-plus-algorithms

This repository contains MATLAB implementations for the analysis of  
**max-plus linear systems**, with a focus on spectral properties and transient
behaviour.

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

Given a max-plus system with multiple durations, the algorithm computes:

- the **max-plus eigenvalue** (cycle mean),
- the associated **eigenvector**,
- a **critical node** belonging to a critical cycle.

This corresponds to solving a **mean-payoff / max-plus spectral problem**.

---

### `lea.m`

Monte Carlo algorithm for estimating **transient deviation bounds** in
max-plus linear systems.

The algorithm:

- applies random perturbations to a max-plus matrix,
- propagates two trajectories with a fixed time shift,
- estimates **lower and upper bounds** on their deviation,
- computes **95% confidence intervals** from Monte Carlo samples.

This method is closely related to **Lyapunov exponent estimation** and
**coupling-time analysis** for stochastic max-plus systems.

---

### `example_lea.m`

Example script demonstrating the use of `lea.m`.

The script:

1. Defines a max-plus system matrix
2. Sets algorithm parameters
3. Runs the Monte Carlo estimation
4. Returns deviation bounds and confidence intervals

This file serves as a **minimal working example**.

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

