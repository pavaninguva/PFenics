# Planar Interface

We consider two cases.

## Allen-Cahn with Constant Mobility and Quartic Potential

Recall the Allen-Cahn equation
$$
\frac{\partial c}{\partial t} = - M\mu, \\
\mu = f'(c) - \kappa \nabla^{2}c.
$$
Specifying $M=1$ and $f(c) = 0.25(1-c^{2})^2$, 
$$
\frac{\partial c}{\partial t} = \kappa \nabla^{2}c + c - c^{3}
$$
where in 1-D
$$
\frac{\partial c}{\partial t} = \kappa \frac{\partial^{2}c}{\partial x^{2}} + c - c^{3}.
$$
We seek traveling wave solutions of the form
$$
c(x,t) = \frac{1}{2} \left(1 - \tanh{\left( \frac{x - ut}{2\sqrt{2\kappa}} \right)} \right),
$$
where $u$ is the speed of the traveling wave. Considering each term in the 1-D equation

$\frac{\partial c}{\partial t}$:
$$
\frac{\partial c}{\partial t} = -\frac{1}{2}\frac{\partial}{\partial t} \left( \tanh{\left( \frac{x - ut}{2\sqrt{2\kappa}} \right)} \right) \\
%
= -\frac{1}{2} \frac{-u}{2\sqrt{2\kappa}} \text{sech}^{2}{\left( \frac{x - ut}{2\sqrt{2\kappa }} \right)} \\
%
= \frac{u}{4\sqrt{2\kappa}}\text{sech}^{2}{\left( \frac{x - ut}{2\sqrt{2\kappa }} \right)}
$$
$\kappa\frac{\partial^{2}c}{\partial x^{2}}$:
$$
\frac{\partial c}{\partial x} = \frac{-1}{4\sqrt{2\kappa}}\text{sech}^{2}{\left( \frac{x - ut}{2\sqrt{2\kappa }} \right)}
$$

$$
\kappa \frac{\partial^{2}c}{\partial x^{2}} = \frac{1}{8}\tanh{\left( \frac{x - ut}{2\sqrt{2\kappa}} \right)}\text{sech}^{2}{\left( \frac{x - ut}{2\sqrt{2\kappa }} \right)}
$$

$c$:
$$
c = \frac{1}{2}\left(1-\tanh{\left( \frac{x - ut}{2\sqrt{2\kappa}} \right)}\right)
$$
$c^{3}$:
$$
c^{3} = \frac{1}{8} \left( 1 - \tanh^{3}\left( \frac{x - ut}{2\sqrt{2\kappa}}\right) + 3\tanh^{2}{\left( \frac{x - ut}{2\sqrt{2\kappa}} \right)} - 3\tanh{\left( \frac{x - ut}{2\sqrt{2\kappa}} \right)} \right)
$$


Simplifying each term by dividing by $\text{sech}^{2}$,

$\frac{\partial c}{\partial t}$: 
$$
\frac{u}{4\sqrt{2\kappa}}
$$
$\kappa \frac{\partial^{2}c}{\partial x^{2}}$:
$$
\frac{1}{16} \tanh{ \left( \frac{x - ut}{2\sqrt{2\kappa}} \right)}
$$
$c^{3}$:
$$
\frac{1}{8} \cosh^{2}{\left( \frac{x - ut}{2\sqrt{2\kappa}} \right)} - \frac{1}{8} \tanh{\left( \frac{x - ut}{2\sqrt{2\kappa}} \right)}\sinh^{2}{\left( \frac{x - ut}{2\sqrt{2\kappa}} \right)} + \frac{3}{8} \sinh^{2}{\left( \frac{x - ut}{2\sqrt{2\kappa}} \right)} \\
%
-\frac{3}{8} \sinh{\left( \frac{x - ut}{2\sqrt{2\kappa}} \right)} \cosh{\left( \frac{x - ut}{2\sqrt{2\kappa}} \right)}
$$

$$
\frac{1}{8} \left( 1 + 4\sinh^{2}{\left( \frac{x - ut}{2\sqrt{2\kappa}} \right) } \right) - \frac{3}{16} \sinh{\left( 2 \frac{x - ut}{2\sqrt{2\kappa}} \right)} \\ 
 - \frac{1}{8} \tanh{\left( \frac{x - ut}{2\sqrt{2\kappa}} \right)}\sinh^{2}{\left( \frac{x - ut}{2\sqrt{2\kappa}} \right)}
$$

$c$:
$$
\frac{1}{2}\cosh^{2}{\left( \frac{x - ut}{2\sqrt{2\kappa}} \right)} - \frac{1}{2}\sinh{\left( \frac{x - ut}{2\sqrt{2\kappa}} \right)}\cosh{\left( \frac{x - ut}{2\sqrt{2\kappa}} \right)} 
$$

$$
\frac{1}{2}\cosh^{2}{\left( \frac{x - ut}{2\sqrt{2\kappa}} \right)} - \frac{1}{4} \sinh{\left( 2\frac{x - ut}{2\sqrt{2\kappa}} \right)}
$$

$c-c^{3}$:
$$
\frac{1}{2}\cosh^{2}{\left( \frac{x - ut}{2\sqrt{2\kappa}} \right)} - \frac{1}{4} \sinh{\left( 2\frac{x - ut}{2\sqrt{2\kappa}} \right)} - \frac{1}{8} - \frac{1}{2} \sinh^{2}{\left( \frac{x - ut}{2\sqrt{2\kappa}} \right)} + \frac{3}{16}\sinh{\left( 2\frac{x - ut}{2\sqrt{2\kappa}} \right)} \\ 
%
+ \frac{1}{8} \tanh{\left( \frac{x - ut}{2\sqrt{2\kappa}} \right)}\sinh^{2}{\left( \frac{x - ut}{2\sqrt{2\kappa}} \right)}
$$

$$
\frac{1}{2}-\frac{1}{8}-\frac{1}{16}\sinh{\left( 2\frac{x - ut}{2\sqrt{2\kappa}} \right)}  + \frac{1}{8} \tanh{\left( \frac{x - ut}{2\sqrt{2\kappa}} \right)} \left( \cosh^{2}{\left( \frac{x - ut}{2\sqrt{2\kappa}} \right)} - 1 \right)
$$

$$
\frac{3}{8}-\frac{1}{16}\sinh{\left( 2\frac{x - ut}{2\sqrt{2\kappa}} \right)} - \frac{1}{8} \tanh{\left( \frac{x - ut}{2\sqrt{2\kappa}} \right)} \\
%
+ \frac{1}{8}\sinh{\left( \frac{x - ut}{2\sqrt{2\kappa}} \right)}\cosh{\left( \frac{x - ut}{2\sqrt{2\kappa}} \right)}
$$

$$
-\frac{1}{8} -\frac{1}{16}\sinh{\left( 2\frac{x - ut}{2\sqrt{2\kappa}} \right)} - \frac{1}{8} \tanh{\left( \frac{x - ut}{2\sqrt{2\kappa}} \right)} + \frac{1}{16}\sinh{\left( 2\frac{x - ut}{2\sqrt{2\kappa}} \right)}
$$

$$
\frac{3}{8} -  \frac{1}{8} \tanh{\left( \frac{x - ut}{2\sqrt{2\kappa}} \right)}
$$

$\kappa \frac{\partial^{2}c}{\partial x^{2}} + c -c^{3}$:
$$
\frac{3}{8} -  \frac{1}{8} \tanh{\left( \frac{x - ut}{2\sqrt{2\kappa}} \right)} + \frac{1}{8} \tanh{\left( \frac{x - ut}{2\sqrt{2\kappa}} \right)}
$$

$$
\kappa \frac{\partial^{2}c}{\partial x^{2}} + c -c^{3} = \frac{3}{8}
$$
Putting things together
$$
\frac{u}{4\sqrt{2\kappa}} = \frac{3}{8} \quad \Rightarrow u = \frac{3\sqrt{\kappa}}{\sqrt{2}} 
$$
The traveling wave has the above velocity and satisfies the Allen-Cahn equation. 

We use the initial condition
$$
c(t,x=0) = \frac{1}{2} \left(1 - \tanh{\left( \frac{x}{2\sqrt{2\kappa}} \right)} \right)
$$
The simulation is simulated to t = 5. 



## References 

1. Choi, J.W., Lee, H.G., Jeong, D. and Kim, J., 2009. An unconditionally gradient stable numerical method for solving the Allen–Cahn equation. *Physica A: Statistical Mechanics and its Applications*, *388*(9), pp.1791-1803.



