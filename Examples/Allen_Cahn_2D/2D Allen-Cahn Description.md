# 2D Allen-Cahn Equation

## Derivation

The classic Allen-Cahn equation is given by
$$
\frac{\partial c}{\partial t} = -M \mu,
$$
where $c$ is the phase field variable, $M$ is the mobility, and $\mu$ is the chemical potential. To obtain $\mu$, first consider the total free energy of the system $F$,
$$
F(c) = \int_{V} f(c) + \frac{\kappa}{2}(\nabla c)^{2} dV.
$$
$\mu$ is defined as the variational derivative of $F$,
$$
\mu = \frac{\delta F}{\delta c} = \frac{\partial F}{\partial c} - \nabla \frac{\partial F}{\partial \nabla c} = f'(c) - \kappa\nabla^{2}c.
$$
We thus have
$$
\frac{\partial c}{\partial t} = M(\kappa\nabla^{2}c - f'(c)).
$$
Assuming a constant mobility, the Allen-Cahn equation is
$$
\frac{\partial c}{\partial t} = \kappa\nabla^{2}c - f'(c).
$$
In most cases in the literature, a quartic polynomial for $f(c)$ is specified,
$$
f(c) = \frac{1}{4}(c^{2}-1)^{2},
$$
which gives us the final form of the Allen-Cahn equation we shall consider
$$
\frac{\partial c}{\partial t} = \kappa \nabla^{2}c - c^{3} + c
$$

## Relevant Features

### Mass Conservation

The first property of note is that the total mass is not conserved. Taking the time derivative of the volume integral of $c$ over the whole domain,
$$
\frac{d}{dt}\int_{V} c \ dV = \int_{V} \frac{\partial c}{\partial t} \ dV \\
= \int_{V} \kappa \nabla^{2}c - f'(c) \ dV
$$
We can split up the integral to handle the gradient term separately. Recognizing that $\nabla^{2}c = \nabla \cdot \nabla c$ and using the Divergence theorem, 
$$
\int_{V} \kappa \nabla^{2}c \ dV = \kappa\int_{V} \nabla \cdot \nabla c \ dV \\
%
= \kappa \int_{S} (\nabla c \cdot \hat{\mathbf{n}}) \ dV
$$
Using a suitable boundary condition e.g., no flux, the above integral is simply zero. We thus get,
$$
\frac{d}{dt}\int_{V} c \ dV  = - \int_{V} f'(c) \ dV
$$
We can see that mass conservation is not guaranteed even with a no-flux boundary condition on the domain. 

### Energy Conservation

The second property of note is that the energy of the system will decrease over time. 

Taking the time derivative of $F$,
$$
\frac{d}{dt} F = \frac{d}{dt}\int_{V} f(c) + \frac{\kappa}{2}(\nabla c)^{2} dV \\
%
= \int_{V}\frac{d}{dt}f(c) + \frac{d}{dt} \left( \frac{\kappa}{2}(\nabla c)^{2} \right) dV \\
%
= \int_{V} \frac{\partial c}{\partial t} f'(c)  + \kappa (\nabla c \cdot \frac{d}{dt} (\nabla c)) \ dV \\
%
= \int_{V} \frac{\partial c}{\partial t} f'(c)  + \kappa \left(\nabla c \cdot \nabla \left( \frac{\partial c}{\partial t} \right) \right) dV.
$$
To proceed further, we need to simplify the dot product. For convenience,  $\frac{\partial c}{\partial t} \equiv c_{t}$. Hence the complicated dot product is
$$
\nabla c \cdot \nabla c_{t}
$$
Using a vector calculus identity,
$$
\nabla c \cdot \nabla c_{t} = \nabla \cdot (c_{t} \nabla c) - c_{t}\nabla^{2} c
$$
Recall that the expression is in an integral,
$$
\int_{V} \nabla c \cdot \nabla c_{t} \ dV = \\
%
\int_{V} \nabla \cdot (c_{t} \nabla c) - c_{t}\nabla^{2} c \ dV = \\
\int_{V} \nabla \cdot (c_{t} \nabla c) \ dV - \int_{V} c_{t}\nabla^{2} c \ dV
$$
We can address the first term using the Divergence theorem,
$$
\int_{V} \nabla \cdot (c_{t} \nabla c) \ dV = \int_{S} c_{t}\nabla c \cdot \hat{\mathbf{n}} \ dS
$$
We recognize that $c_{t}$ is a scalar which allows us to use the scalar multiplication property of the dot product,
$$
\int_{S} c_{t}\nabla c \cdot \hat{\mathbf{n}} \ dS = \int_{S} c_{t} (\nabla c \cdot \hat{\mathbf{n}}) \ dS
$$
At this point, we can employ our boundary condition which could be a no-flux boundary condition i.e., ($\hat{\mathbf{n}} \cdot \nabla c = 0$) to conclude that the above integral is simply zero. A more general argument can be made about the validity of the above conclusion, but it is a bit more complicated. We thus get
$$
\frac{d}{dt} F = \int_{V} \frac{\partial c}{\partial t} f'(c) - \kappa \frac{\partial c}{\partial t} \nabla^{2}c \ dV \\
%
= \int_{V} \frac{\partial c}{\partial t} \left ( f'(c) - \kappa \nabla^{2}c \right) dV \\
%
= \int_{V} \frac{\partial c}{\partial t} \left ( -\frac{\partial c}{\partial t} \right) dV \\
%
= -\int_{V}\left (\frac{\partial c}{\partial t} \right)^{2} \leq 0.
$$

### Mean Curvature Flow

The last feature we shall note is that for a curved surface, we can track the evolution of the mean curvature analytically. This means for a benchmarking case of a circle / sphere, we can track the evolution of its radius $r$,
$$
r(t) = \sqrt{r_{0}^{2} + 2(1-d)t},
$$
where $r_{0}$ is the initial radius of the circle, $d$ is the dimension of the system (i.e., 2 for a 2D simulation), and $t$ is the time. 

## Weak Form

To implement the Allen-Cahn equation in FEniCS, we need to obtain the weak form of the equation. 

Recall the Allen-Cahn equation,
$$
\frac{\partial c}{\partial t} = \kappa\nabla^{2}c - f'(c),
$$
with the initial conditions,
$$
c(\mathbf{x},0) = c_{0}(\mathbf{x})
$$
We first discretize the time derivative with a backwards difference scheme,
$$
\frac{\partial c}{\partial t} \approx \frac{c^{n+1} - c^{n}}{\Delta t}.
$$
Treating all terms implicitly,
$$
\frac{c^{n+1} - c^{n}}{\Delta t} = \kappa \nabla^{2}c^{n+1} - f'(c^{n+1}).
$$
Reordering to move all terms to the LHS, 
$$
c^{n+1} - c^{n} -\Delta t \kappa \nabla^{2}c^{n+1} + \Delta t f^{n+1} = 0
$$
We seek to transform the problem into the standard form
$$
a(c^{n+1},v) + L_{n+1}(v) = 0,
$$
where $v$ is the finite element test function that is introduced. The LHS contains all the terms with the variable being solved for (i.e., the terms that are treated implicitly), and the LHS contains all the terms that are treated explicitly / are known (e.g., $c^{n}$). Correspondingly, we have
$$
a(c^{n+1},v) = \int_{V} c^{n+1} v + \kappa \Delta t \nabla c^{n+1} \cdot \nabla v + \Delta t f^{n+1}v \ dV,
$$
and
$$
L_{n+1}(v) = \int_{V} -c^{n} v \ dV.
$$
We are thus in a position to put implement the Allen-Cahn equation into FEniCS. 

## References

1. Kim, Y., Ryu, G. and Choi, Y., 2021. Fast and accurate numerical solution of Allen–Cahn equation. *Mathematical Problems in Engineering*, *2021*, pp.1-12.
2. Church, J.M., Guo, Z., Jimack, P.K., Madzvamuse, A., Promislow, K., Wetton, B., Wise, S.M. and Yang, F., 2019. High accuracy benchmark problems for Allen-Cahn and Cahn-Hilliard dynamics. *Communications in computational physics*, *26*(4).
3. Li, Y., Lee, H.G., Jeong, D. and Kim, J., 2010. An unconditionally stable hybrid numerical method for solving the Allen–Cahn equation. *Computers & Mathematics with Applications*, *60*(6), pp.1591-1606.