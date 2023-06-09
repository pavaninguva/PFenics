# Classic Cahn-Hilliard Equation

## Derivation

The classic Cahn-Hilliard equation is
$$
\frac{\partial c}{\partial t} = \nabla \cdot (M \nabla \mu),
$$
where $c$ is the phase field variable, $M$ is the mobility coefficient, and $\mu$ is regarded as the chemical potential. The CH equation is derived by considering the total free energy of the system,
$$
F(c) = \int_{V} f(c) + \frac{\kappa}{2}(\nabla c)^{2} \ dV.
$$
As is with the AC equation, $\mu$ is the variational derivative of $F$,
$$
\mu = f'(c) - \kappa \nabla^{2}c.
$$
We thus have,
$$
\frac{\partial c}{\partial t} = \nabla \cdot \left( M \nabla \left(f'(c) - \kappa \nabla^{2}c \right) \right)
$$
Assuming a constant mobility and $\kappa$ (which is the case for the classic CH equation), we get,
$$
\frac{\partial c}{\partial t} = \nabla^{2} \left(f'(c)\right) - \kappa \nabla^{4} c.
$$
Typically, $f(c) = \frac{1}{4} (c^{2} -1)^{2}$ is specified. As the CH equation is a 4th order PDE, this makes it difficult to naively implement it in FEniCS, hence it is typically solved as a set of coupled 2nd order PDEs,
$$
\frac{\partial c}{\partial t} = \nabla \cdot \left( M \nabla \mu \right), \\
%
\mu = f'(c) - \kappa \nabla^{2}c.
$$

## Relevant Features

### Mass Conservation

By formulation, the Cahn-Hilliard equation is conservative, i.e., it is based on the conservation equation,
$$
\frac{\partial c}{\partial t} + \nabla \cdot \mathbf{J} = g,
$$
where $\mathbf{J}$ is the flux and $g$ is a source/sink term. In the case of the classic Cahn-Hilliard equation, no source/sink terms are considered. Correspondingly with the use of the Gauss divergence theorem, 
$$
\frac{d}{dt}\int_{V} c \ dV = \int_{V} \frac{\partial c}{\partial t} \ dV = -\int_{V} \nabla \cdot \mathbf{J} \ dV   = -\int_{S} \mathbf{J}\cdot \hat{\mathbf{n}} \ dS.
$$
With appropriate boundary conditions e.g., no-flux i.e., $\mathbf{J}\cdot \hat{\mathbf{n}} = 0$, we obtain the property that total mass is conserved,
$$
\frac{d}{dt} \int_{V} c \ dV = 0
$$

### Energy Conservation

The total energy of the system decreases and follows the same arguments made in `Allen_Cahn_2D`. 

## Simulation Setup

### Case 1:  7 Circles

The first case considered assumes a constant mobility with $M = 1$. A simulation domain of $V\in [0, 2\pi]^{2}$ is specified. The initial conditions are specified with 7 circles,
$$
c(x,y,t=0) = -1 + \sum_{i=1}^{7} f\left( \sqrt{(x-x_{i})^{2} + (y-y_{i})^{2}} -r_{i} \right) \\
%
f(s) = 
\begin{cases}
2\exp{\frac{-\kappa}{s^{2}}} & s < 0 \\
%
0 & \text{otherwise}
\end{cases}
$$

### Case 2: Energy Decay

We track the evolution of the energy on the domain $V \in [0,2\pi]^{2}$ with $\kappa = 0.01$ and the initial conditions
$$
c(x,y,0) = 0.0 + \epsilon,
$$

where $\epsilon$ is random noise that is added to the homogeneous profile to initiate demixing. 

## Weak Form

Recall the Cahn-Hilliard equation
$$
\frac{\partial c}{\partial t} = \nabla \cdot \left( M \nabla \mu \right), \\
%
\mu = f'(c) - \kappa \nabla^{2}c.
$$
Using a backwards Euler scheme for the time derivative and treating all terms implicitly,
$$
\frac{c^{n+1} - c^{n}}{\Delta t} - \nabla \cdot (M \nabla \mu^{n+1}) = 0, \\
%
\mu^{n+1} - f'(c^{n+1}) + \kappa \nabla^{2}c^{n+1} = 0.
$$
Introducing test functions $q$ and $v$, 
$$
\int_{V} c^{n+1} q - c^{n}q + M\Delta t\nabla\mu^{n+1} \cdot \nabla q \ dV = 0 \\
%
\int_{V} \mu^{n+1}v - f'(c^{n+1})v - \kappa \nabla c^{n+1} \cdot \nabla v \ dV = 0
$$

## References

1. Wu, H., 2021. A Review on the Cahn-Hilliard Equation: Classical Results and Recent Advances in Dynamic Boundary Conditions. *arXiv preprint arXiv:2112.13812*.
2. Church, J.M., Guo, Z., Jimack, P.K., Madzvamuse, A., Promislow, K., Wetton, B., Wise, S.M. and Yang, F., 2019. High accuracy benchmark problems for Allen-Cahn and Cahn-Hilliard dynamics. *Communications in computational physics*, *26*(4).