# Finite Difference Method in MATLAB

## Overview

This repository contains a MATLAB implementation of three finite difference schemes for solving the **Heat Equation**:

$$
\frac{\partial u}{\partial t} = \alpha \frac{\partial^2 u}{\partial x^2}
$$

where:
-  $u(x,t)$  is the temperature at position \( x \) and time \( t \),
-  $\alpha$ is the thermal diffusivity constant.

The following three finite difference schemes are implemented to solve this equation:

- Crank-Nicolson Scheme:     $$\frac{{U_i^{n + 1} - U_i^n}}{\tau } = \frac{1}{{2{h^2}}}\left[ {U_{i + 1}^{n + 1} - 2U_i^{n + 1} + U_{i - 1}^{n + 1} + U_{i + 1}^n - 2U_i^n + U_{i - 1}^n} \right]$$
 
- Forward Euler Scheme:    $$\frac{{U_i^{n + 1} - U_i^n}}{\tau } = \frac{1}{{{h^2}}}\left( {U_{i - 1}^n - 2U_i^n + U_{i + 1}^n} \right)$$

- Stabilized Central Time Scheme:    $$\frac{{U_i^{n + 1} - U_i^{n - 1}}}{{2\tau }} = \frac{1}{{{h^2}}}\left( {U_{i + 1}^n - U_i^{n + 1} - U_i^{n - 1} + U_{i - 1}^n} \right)$$

The main objective of this project is to implement these schemes to solve the **1D Heat Equation** with zero Dirichlet boundary conditions and an initial condition of your choice, over the interval  (0,1).

## Prerequisites
* MATLAB R2015 or later

## Project Structure
1. Crank Nicholson Scheme Implementation/:
   - CNS.m: Crank-Nicholson function which generates
     * A Mesh Plot to the solution to the heat equation
     * A plot of $U_{Approx}$ and $U_{exact}$ at the final time
   - MainCNS.m: A script that contains the implementation of CNS for different values of N which generates
     * Results for each N
     * Convergence plot.
 
   


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


## Results

Here are the results from applying the finite difference schemes to the heat equation:

### Result 1: Crank-Nicolson Scheme


| <img src="src/cn1.png" alt="Image 1" style="width:75%;"/> | <img src="src/cn2.png" alt="Image 2" style="width:75%;"/> |
|-----------------------------------------------------------|-----------------------------------------------------------|
| Crank-Nicolson solution with dt=0.0125                    | Plot of Exact and Approximated solution at T = 0.5        |


| <img src="src/cn3.png" alt="Image 3" style="width:50%;"/> |
|-----------------------------------------------------------|
| Convergence Plot for Crank-Nicolson Solution              |


### Result 2: Forward Euler Scheme

| <img src="src/fe1.png" alt="Image 4" style="width:75%;"/> | <img src="src/fe2.png" alt="Image 5" style="width:75%;"/> |
|-----------------------------------------------------------|-----------------------------------------------------------|
| Forward Euler solution with dt=0.0025                     | Plot of Exact and Approximated solution at T = 0.5        |


| <img src="src/fe3.png" alt="Image 6" style="width:50%;"/> |
|-----------------------------------------------------------|
| Convergence Plot for Forward Euler                        |

### Result 3: Stabilized Centered Time

| <img src="src/st1.png" alt="Image 7" style="width:75%;"/> | <img src="src/st2.png" alt="Image 8" style="width:75%;"/> |
|-----------------------------------------------------------|-----------------------------------------------------------|
| Stabilaized Central Time Solution of Heat Equation with dt=0.0053  | Plot of Exact and Approximated solution at T = 0.5 |


| <img src="src/st3.png" alt="Image 9" style="width:50%;"/> |
|-----------------------------------------------------------|
| Convergence Plot for Stabilized Central Time             |



## Documentation
The project is accompanied by a detailed report which contains the Analysis of Stability and Consistency.
![File](docs/ConsistencyStabilityConvergence.pdf)



