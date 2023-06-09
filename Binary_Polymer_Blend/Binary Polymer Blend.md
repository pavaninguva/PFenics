# Binary Polymer Blend

## Derivation

This model is based on the classic Cahn-Hilliard equation and makes a few modifications to the original equation to adapt it for describing the demixing of a binary polymer blend. In a binary polymer blend, there are two species present (We can refer to them as species A and B for convenience). The order parameter for this system is the volume fraction which is conserved. By tracking the volume fraction of either $A$ or $B$, we can infer the volume fraction of the other species simply by recognizing,
$$
c_{A} + c_{B} = 1
$$
Therefore, only a single Cahn-Hilliard equation is necessary to track the evolution of the system. Recall the Cahn-Hilliard equation
$$
\frac{\partial c}{\partial t} = \nabla \cdot (M \nabla \mu),
$$
where $c$ is now the volume fraction of species $A$. The free energy of the system, $F$ is given as
$$
F(c) = \int_{V} f(c) + \frac{\kappa}{2}(\nabla c)^{2} \ dV.
$$
 To better represent the thermodynamics of the polymer blend, we use the Flory-Huggins free energy for $f(c)$,
$$
f(c) = \frac{c \ln{c}}{N} + \frac{(1-c) \ln{(1-c)}}{N} + \chi c(1-c),
$$
where $N$ is the degree of polymerization, and $\chi$ is the Flory-Huggins interaction parameter. Often, $N$ is taken to be the chain length. The use of the Flory-Huggins equation (or similar equations with a logarithmic term) in the Cahn-Hilliard equation is not uncommon in the literature and is often referred to as a logarithmic potential.

When evaluated, the CH equation is
$$
\frac{\partial c}{\partial t} = \nabla \cdot (M \nabla \mu), \\
%
\mu = \frac{1}{N} \left( \ln{\frac{c}{1-c}} + N\chi -2N\chi c \right) - \kappa\nabla^{2}c.
$$
The use of the FH equation imposes certain requirements on the mobility $M$ and gradient energy parameter $\kappa$. It can be shown that for a binary polymer blend, 
$$
\kappa = \underbrace{\frac{1}{3}(R_{G,1}^{2}+R_{G,2}^{2})\chi}_{\text{Enthalpic Contribution}} +   \underbrace{\frac{R_{G,1}^{3}}{3} \left(  \frac{1}{N_{1}c} +  \frac{1}{N_{2}(1-c)} \right)}_{\text{Entropic Contribution}},
$$
where the subscript $1,2$ corresponds to species A or B respectively and $R_{G,i}$ is the radius of gyration for species $i$. In most cases, the entropic contribution is neglected which is justifiable in the limit of large $N$. In addition, for a symmetric polymer blend, we can assume $N_{1} = N_{2}$ and $R_{G,1} = R_{G,2}$, which simplifies $\kappa$ to
$$
\kappa = \frac{2R_{G}^{2}}{3} \chi.
$$
$M$ is given as
$$
M = Dc(1-c),
$$
where $D$ is the diffusion coefficient. To understand this form of $M$, in the limit of an ideal solution where $\chi = 0$ and $N = 1$, 
$$
\mu = \ln{\frac{c}{1-c}}.
$$
Using the chain rule,
$$
\nabla \mu = \frac{\partial \mu}{\partial c}\nabla c = \frac{1}{c(1-c)}\nabla c.
$$
Therefore, 
$$
\frac{\partial c}{\partial t} = \nabla \cdot \left( Dc(1-c) \frac{1}{c(1-c)}\nabla c \right) = \nabla \cdot (D \nabla c).
$$
We see that we recover the classical diffusion equation when $M$ is specified accordingly. This form of $M$ has also been explored in the literature and is referred to as a degenerate mobility. 

### Defining The Model

Putting things together, we have
$$
\frac{\partial c}{\partial t} = \nabla \cdot (Dc(1-c) \nabla \mu) \\
%
\mu =\frac{1}{N} \left( \ln{\frac{c}{1-c}} + N\chi -2N\chi c \right) - \frac{2 R_{G}^{2}\chi}{3}\nabla^{2}c.
$$
Introducing the following scalings
$$
\mathbf{x} = R_{G}\tilde{\mathbf{x}}, \quad t = \frac{NR_{G}^{2}}{D}\tilde{t}.
$$
The model equations become,
$$
\tilde{\mu} = \frac{1}{N} \left( \ln{\frac{c}{1-c}} + N\chi -2N\chi c \right) - \frac{2\chi}{3}\tilde{\nabla}^{2}c.
$$

$$
N\tilde{\mu} = \ln{\frac{c}{1-c}} + N\chi - 2N\chi x - \frac{2}{3}N\chi \tilde{\nabla}^{2} c.
$$

$$
\frac{D}{NR_{G}^{2}}\frac{\partial c}{\partial \tilde{t}}  = \frac{D}{R_{G}^{2}}\tilde{\nabla} \cdot (c(1-c) \tilde{\nabla} \tilde{\mu}) \\
%
\frac{\partial c}{\partial \tilde{t}} = \tilde{\nabla} \cdot (c(1-c) \nabla (N\tilde{\mu}))
$$

Treating $N\chi$ as a single parameter and dropping the $\sim$, we have
$$
\frac{\partial c}{\partial t} = \nabla \cdot (M \nabla \mu), \\
\mu = f'(c) - \kappa \nabla^{2}c,
$$
where
$$
M = c(1-c), \quad f(c) = c\ln{c} + (1-c)\ln{(1-c)} + N\chi c(1-c), \quad \kappa = \frac{2}{3}N\chi.
$$


## References

1. Inguva, P.K., Walker, P.J., Yew, H.W., Zhu, K., Haslam, A.J. and Matar, O.K., 2021. Continuum-scale modelling of polymer blends using the Cahn–Hilliard equation: transport and thermodynamics. *Soft Matter*, *17*(23), pp.5645-5665.
2. Nauman, E.B. and He, D.Q., 1994. Morphology predictions for ternary polymer blends undergoing spinodal decomposition. *Polymer*, *35*(11), pp.2243-2255.
3. Ariyapadi, M.V. and Nauman, E.B., 1990. Gradient energy parameters for polymer–polymer–solvent systems and their application to spinodal decomposition in true ternary systems. *Journal of Polymer Science Part B: Polymer Physics*, *28*(12), pp.2395-2409.
4. Vasishtha, N. and Nauman, E.B., 1994. Hydrodynamic effects in the phase separation of binary polymer mixtures. *Chemical Engineering Communications*, *129*(1), pp.29-39.