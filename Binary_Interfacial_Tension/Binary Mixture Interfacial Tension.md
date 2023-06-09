# Binary Mixture Interfacial Tension

## Derivation

Recall the Allen-Cahn and Cahn-Hilliard equations

Allen-Cahn:
$$
\frac{\partial c}{\partial t} = - M\mu, \\
\mu = f'(c) - \kappa \nabla^{2}c.
$$
Cahn-Hilliard:
$$
\frac{\partial c}{\partial t} = \nabla \cdot (M \nabla \mu), \\
\mu = f'(c) - \kappa \nabla^{2}c.
$$
An important result arising from this work (which incidentally is related to the work of van der Waals) is that of the square-gradient theory which can be used to compute the interfacial tension in a mixture. Assuming a planar interface (i.e., 1-D where we only have a single interfacial region in the domain), we can compute the interfacial tension of the system,
$$
\sigma = \rho_{m}RT \int_{-\infty}^{\infty} \kappa \left( \frac{\partial c}{\partial x} \right)^{2} dx,
$$
where $\rho_{m}$ is the monomer molar density, $R$ is the ideal gas constant, and $T$ is the temperature. The above expression can be adapted to work for non-dimensionalized units, 
$$
\frac{\sigma}{L_{0} \rho_{m}RT} = \int_{-\infty}^{\infty} \tilde{\kappa} \left( \frac{\partial c}{\partial \tilde{x}} \right)^{2} d\tilde{x},
$$
where $L_{0}$ is the length scale and as considered in the Binary polymer blend case is $R_{G}$. We can call the LHS $\tilde{\sigma}$ to indicate a dimensionless interfacial tension for clarity. 

We consider four cases:

### 1. Allen-Cahn with Quartic Potential and Constant Mobility

Recall the Allen-Cahn equation
$$
\frac{\partial c}{\partial t} = - M\mu, \\
\mu = f'(c) - \kappa \nabla^{2}c.
$$
Specifying a quartic potential and $M = 1$, 
$$
f(c) = 0.25(1-c^{2})^{2},
$$

$$
f'(c) = -(c-c^{3}).
$$

We thus have 
$$
\frac{\partial c}{\partial t} = \kappa\nabla^{2}c +   c-c^{3}.
$$
At equilibrium, in 1-D,
$$
\frac{d^{2}c}{dx^{2}} =  \frac{1}{\kappa} (c^{3} - c).
$$
This equation permits the solution 
$$
c(x) = \tanh{\left(\frac{x}{\sqrt{2\kappa}}\right)}.
$$
Check if the above equation is permitted as a solution:

LHS:
$$
\frac{d^{2}c}{dx^{2}} = -\frac{1}{\kappa} \left( \tanh{\left( \frac{x}{\sqrt{2\kappa}} \right)} \text{sech}^{2} \left( \frac{x}{\sqrt{2\kappa}}\right) \right)
$$

$$
\tanh{\left( \frac{x}{\sqrt{2\kappa}} \right)} \text{sech}^{2} \left( \frac{x}{\sqrt{2\kappa}}\right) = \frac{\sinh{\left( \frac{x}{\sqrt{2\kappa}}\right)}}{\cosh^{3}{\left(\frac{x}{\sqrt{2\kappa}} \right)}}
$$

$$
\Rightarrow \frac{d^{2}c}{dx^{2}} = \frac{-1}{\kappa} \frac{\sinh{\left( \frac{x}{\sqrt{2\kappa}}\right)}}{\cosh^{3}{\left(\frac{x}{\sqrt{2\kappa}} \right)}}
$$

RHS:
$$
c^{3} - c = \tanh^{3}{\left( \frac{x}{\sqrt{2\kappa}} \right)} - \tanh{\left( \frac{x}{\sqrt{2\kappa}} \right)} \\
%
= \frac{\sinh^{3}{\left( \frac{x}{\sqrt{2\kappa}} \right)}}{\cosh^{3}{\left(\frac{x}{\sqrt{2\kappa}} \right)}} - \frac{\sinh{\left( \frac{x}{\sqrt{2\kappa}} \right)}}{\cosh{\left(\frac{x}{\sqrt{2\kappa}} \right)}} \\
%
= \frac{\sinh^{3}{\left(\frac{x}{\sqrt{2\kappa}} \right) - \sinh{\left(\frac{x}{\sqrt{2\kappa}}\right)}} \cosh^{2}{\left(\frac{x}{\sqrt{2\kappa}} \right)}}{\cosh^{3}{\left(\frac{x}{\sqrt{2\kappa}} \right)}} \\
%
= \frac{\sinh^{3}{\left(\frac{x}{\sqrt{2\kappa}} \right) - \sinh{\left(\frac{x}{\sqrt{2\kappa}} \right) \left( 1 + \sinh^{2}{\left(\frac{x}{\sqrt{2\kappa}} \right)} \right)}}}{\cosh^{3}{\left(\frac{x}{\sqrt{2\kappa}} \right)}} \\
%
\frac{-\sinh{\left(\frac{x}{\sqrt{2\kappa}} \right)}}{\cosh^{3}{\left(\frac{x}{\sqrt{2\kappa}} \right)}}
$$

$$
\Rightarrow \frac{1}{\kappa} (c^{3} - c) = \frac{-1}{\kappa}\frac{\sinh{\left( \frac{x}{\sqrt{2\kappa}}\right)}}{\cosh^{3}{\left(\frac{x}{\sqrt{2\kappa}} \right)}}
$$

$$
\therefore LHS = RHS
$$

To obtain an analytical expression for $\sigma$, there are a few more steps. 
$$
\frac{dc}{dx} = \frac{1}{\sqrt{2\kappa}} \text{sech}^{2}{\left( \frac{x}{\sqrt{2\kappa}} \right)}
$$

$$
\left( \frac{dc}{dx} \right)^{2} = \frac{1}{2\kappa}\text{sech}^{4}{\left( \frac{x}{\sqrt{2\kappa}} \right)}
$$

$$
\int_{-\infty}^{\infty} \kappa \left( \frac{dc}{dx} \right)^{2} dx =  \frac{1}{2}\int_{-\infty}^{\infty}\text{sech}^{4}{\left( \frac{x}{\sqrt{2\kappa}} \right)} dx \\
%
= \frac{1}{2}\left[\sqrt{2\kappa}\tanh{\left( \frac{x}{\sqrt{2\kappa}}\right)} - \frac{\sqrt{2\kappa}}{3} \tanh^{3}{\left(\frac{x}{\sqrt{2\kappa}}\right)} \right]^{\infty}_{-\infty}\\
%
= \frac{1}{2} \left[ \left( \sqrt{2\kappa} - \frac{1}{3}\sqrt{2\kappa} \right) - \left( -\sqrt{2\kappa} + \frac{1}{3}\sqrt{2\kappa} \right) \right] \\
%
= \frac{1}{2} \left( \frac{4}{3}\sqrt{2\kappa} \right) = \frac{2\sqrt{2}}{3}\sqrt{\kappa}
$$

$$
\therefore \sigma = \frac{2\sqrt{2}}{3}\sqrt{\kappa}
$$

### 2. Allen-Cahn with Logarithmic Potential and Constant Mobility

At equilibrium, in 1-D, the model equation is 
$$
\frac{1}{\kappa}\frac{d^{2}c}{dx^{2}} = \ln{\left( \frac{x}{1-x} \right)} - 2 N\chi c + N\chi,
$$
where $\kappa = \frac{2}{3}N\chi$ as previously discussed when using a F-H type free energy. 

The solution to this equation is very sensitive to the initial guess. For a given value of $N\chi$, it is possible to find the equilibrium composition of each phase using a function like `fsolve`. This is equivalent to finding the minima of the potential. Assuming we found the two minima $c_{0}$ and $c_{1}$, the initial guess is specified as a modified Boltzmann sigmoid curve
$$
c^{*}(x) = \frac{c_{0} - c_{1}}{1 + \exp(N\chi x)} + c_{1}
$$


## Weak Form Representation

### Equilibrium Allen-Cahn equation

Recall the Allen-Cahn equation at equilibrium,
$$
\kappa \frac{d^{2}c}{dx^{2}} + c - c^{3} =  0
$$
Introduce test function $v$ and integrate by parts



## References

1. Wazwaz, A.M., 2007. The tanh–coth method for solitons and kink solutions for nonlinear parabolic equations. *Applied Mathematics and Computation*, *188*(2), pp.1467-1475.
2. Yang, J., Li, Y., Lee, C., Choi, Y. and Kim, J., 2023. Fast evolution numerical method for the Allen–Cahn equation. *Journal of King Saud University-Science*, *35*(1), p.102430.
3. Chen, L.Q. and Zhao, Y., 2022. From classical thermodynamics to phase-field method. *Progress in Materials Science*, *124*, p.100868.
4. Inguva, P.K., Walker, P.J., Yew, H.W., Zhu, K., Haslam, A.J. and Matar, O.K., 2021. Continuum-scale modelling of polymer blends using the Cahn–Hilliard equation: transport and thermodynamics. *Soft Matter*, *17*(23), pp.5645-5665.
5. Nauman, E.B. and He, D.Q., 2001. Nonlinear diffusion and phase separation. *Chemical Engineering Science*, *56*(6), pp.1999-2018.

