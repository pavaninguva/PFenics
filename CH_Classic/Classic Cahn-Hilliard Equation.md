# Classic Cahn-Hilliard Equation

## Derivation

The classic Cahn-Hilliard equation is
$$
\frac{\partial c}{\partial t} = \nabla \cdot (M \nabla \mu),
$$
where $c$ is the phase field variable, $M$ is the mobility coefficient, and $\mu$ is regarded as the chemical potential. 

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
With appropriate boundary conditions e.g., no-flux i.e., $\mathbf{J}\cdot \hat{\mathbf{n}} = 0$, we obtain the property that total mass is conserved, i.e., 
$$
\frac{d}{dt} \int_{V} c \ dV = 0
$$

### Energy Conservation

### Interfacial Tension

## Weak Form

## References

1. Wu, H., 2021. A Review on the Cahn-Hilliard Equation: Classical Results and Recent Advances in Dynamic Boundary Conditions. *arXiv preprint arXiv:2112.13812*.