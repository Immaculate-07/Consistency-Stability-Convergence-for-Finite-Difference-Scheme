# Finite Difference Method in MATLAB

## Overview

This repository contains a MATLAB implementation of four finite difference schemes for solving the **Heat Equation**:

$$
\frac{\partial u}{\partial t} = \alpha \frac{\partial^2 u}{\partial x^2}
$$

where:
-  $u(x,t)$  is the temperature at position \( x \) and time \( t \),
-  $\alpha$ is the thermal diffusivity constant.

The following four finite difference schemes are implemented to solve this equation:

- Crank-Nicolson Scheme:     $$\frac{{U_i^{n + 1} - U_i^n}}{\tau } = \frac{1}{{2{h^2}}}\left[ {U_{i + 1}^{n + 1} - 2U_i^{n + 1} + U_{i - 1}^{n + 1} + U_{i + 1}^n - 2U_i^n + U_{i - 1}^n} \right]$$
 
- Forward Euler Scheme:    $$\frac{{U_i^{n + 1} - U_i^n}}{\tau } = \frac{1}{{{h^2}}}\left( {U_{i - 1}^n - 2U_i^n + U_{i + 1}^n} \right)$$
- Central Difference for Time Derivative:     $$\frac{{U_i^{n + 1} - U_i^{n-1}}}{{2\tau }} = \frac{{U_{i+1}^n - 2U_i^n + U_{i - 1}^n}}{{{h^2}}}$$
- Stabilized Central Time Scheme:    $$\frac{{U_i^{n + 1} - U_i^{n - 1}}}{{2\tau }} = \frac{1}{{{h^2}}}\left( {U_{i + 1}^n - U_i^{n + 1} - U_i^{n - 1} + U_{i - 1}^n} \right)$$

The main objective of this project is to implement these schemes to solve the **1D Heat Equation** with zero Dirichlet boundary conditions and an initial condition of your choice, over the interval  (0,1).

$ Getting Started 
1. Clone the repository:
   git clone https://github.com/Immaculate-07/Consitency-Stability-Convergence-for-Finite-Difference-Scheme
   cd Consitency-Stability-Convergence-for-Finite-Difference-Scheme

## Problem

Consider the following Partial Differential Equation (PDE) for the **Heat Equation**:

$$
u_t = u_{xx}
$$

with the boundary and initial conditions:

$$
u(0, t) = u(1, t) = 0 \quad \text{for} \quad t > 0
$$

and the initial condition:

$$
u(x, 0) = 
\begin{cases} 
\sin(\pi x) & \text{for} \quad x \in [0,1] \\
0 & \text{otherwise}
\end{cases}
$$

where:
- $u(x, t)$  represents the temperature distribution at position \( x \) and time \( t \),
- $u_t$  is the partial derivative of \( u \) with respect to time,
- $u_{xx}$ is the second spatial derivative of \( u \).

The exact solution to this problem is:

$$
u(x,t) = \exp(-\pi^2 t) \sin(\pi x)
$$



## Documentation
The project is accompanied by a detailed report which contains the Analysis of Stability and Consistency.

