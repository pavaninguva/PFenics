# Finite Differencing For CH With Log Potential

## Equation

Consider the Cahn-Hilliard equation with a Log potential
$$
\begin{align}
\frac{\partial \phi}{\partial t} &= \nabla \cdot (M \nabla \mu), \\
\mu &= \frac{\partial f}{\partial \phi} - \kappa \nabla^{2}\phi.
\end{align}
$$
We have $f$ which the the homogeneous free energy of mixing given by
$$
f = \phi \ln{\phi} + (1-\phi)\ln{(1-\phi)} + \chi \phi (1-\phi),
$$
and $M$, which is the mobility given by 
$$
M = \phi (1-\phi).
$$
The gradient energy parameter is given by 
$$
\kappa = \frac{2}{3}\chi.
$$
We have Neumann BCs for both $\mu$ and $\phi$ at both ends of the domain, 
$$
\frac{\partial \phi}{\partial x}\bigg|_{x=0} = \frac{\partial \phi}{\partial x}\bigg|_{x=L} = \frac{\partial \mu}{\partial x}\bigg|_{x=0} = \frac{\partial \mu}{\partial x}\bigg|_{x=L} = 0.
$$


## Applying Finite Differences

Casting the equations in 1D,
$$
\begin{align}
\frac{\partial \phi}{\partial t} &= \frac{\partial}{\partial x} \left( M \frac{\partial \mu}{\partial x} \right), \\
%
\mu &= \frac{\partial f}{\partial \phi} - \kappa \frac{\partial^{2}\phi}{\partial x^{2}}. 
\end{align}
$$
As a first iteration, we treat all terms implicitly with a backwards Euler approximation for time stepping. Using the $i$ superscript to denote time and the $j$ subscript to denote space,
$$
\frac{\phi^{i+1}_{j} - \phi^{i}_{j}}{\Delta t} = \frac{q^{i+1}_{j+ 1/2} - q_{j-1/2}^{i+1}}{\Delta x},
$$
where
$$
q^{i+1}_{j+1/2} = M_{j+1/2}^{i+1}\frac{\mu_{j+1} - \mu_{j}}{\Delta x}, \\
q^{i+1}_{j-1/2} = M_{j-1/2}^{i+1}\frac{\mu_{j} - \mu_{j-1}}{\Delta x}.
$$
Note that the above is still a central difference scheme about the half way point. 

Combining things
$$
\phi_{j}^{i+1} - \phi_{j}^{i} = \frac{\Delta t}{(\Delta x)^{2}} \left( M^{i+1}_{j+1/2}( \mu_{j+1}^{i+1} -\mu_{j}^{i+1}) - M_{j-1/2}^{i+1}(\mu_{j}^{i+1} -  \mu_{j-1}^{i+1}) \right)
$$
There are multiple choices for $M_{j+1/2}^{i+1}$ :
$$
M_{j+1/2}^{i+1} = \frac{1}{2}\left( M(\phi_{j+1}^{i+1}) + M(\phi_{j}^{i+1}) \right), 
\\
M_{j+1/2}^{i+1} = M\left( \frac{\phi_{j+1}^{i+1} + \phi_{j}^{i+1}}{2} \right),
\\
M_{j+1/2}^{i+1} = \frac{2M(\phi_{j+1}^{i+1})M(\phi_{j}^{i+1})}{M(\phi_{j+1}^{i+1}) + M(\phi_{j}^{i+1})}.
$$
And for the chemical potential,
$$
\mu_{j}^{i+1} = \frac{\partial f}{\partial \phi}_{j}^{i+1} - \kappa \frac{\phi^{i+1}_{j+1} - 2\phi_{j}^{i+1} + \phi_{j-1}^{i+1}}{(\Delta x)^{2}}
$$
To handle the boundary conditions, which are set at $j = 0$, and $j=N$, we have
$$
\frac{\partial \phi}{\partial x} \bigg|_{x=0} \approx \frac{\phi_{1}^{i+1} - \phi_{-1}^{i+1}}{2\Delta x} = 0 \quad \Rightarrow \phi_{-1}^{i+1} = \phi_{1}^{i+1}
$$

$$
\frac{\partial \mu}{\partial x}\bigg|_{x=0} \approx \frac{\mu_{1}^{i+1} - \mu_{-1}^{i+1}}{2\Delta x} = 0 \quad \Rightarrow \mu_{-1}^{i+1} = \mu_{1}^{i+1}
$$

where $\phi_{-1}$ , $\mu_{-1}$ are the ghost node. We can now eliminate the ghost nodes in the discretized equations at the boundaries
$$
\mu_{0}^{i+1} = \frac{\partial f}{\partial \phi}\bigg|_{0}^{i+1} - 2\kappa \frac{\phi_{1}^{i+1} - \phi_{0}^{i+1}}{(\Delta x)^{2}}
$$
For $\phi$, we note that at $j=0$, $M_{j+1/2}^{i+1} = M_{j-1/2}^{i+1}$, which gives us
$$
\phi_{0}^{i+1} - \phi_{0}^{i} = \frac{2\Delta t}{(\Delta x)^{2}} \left( M_{1/2}^{i+1} (\mu_{1}^{i+1} - \mu_{0}^{i+1}) \right)
$$
Similarly, for $j=N$, 
$$
\frac{\partial \phi}{\partial x} \bigg|_{x=L} \approx \frac{\phi_{N+1}^{i+1} - \phi_{N-1}^{i+1}}{2\Delta x} = 0 \quad \Rightarrow \phi_{N-1}^{i+1} = \phi_{N+1}^{i+1}
$$

$$
\frac{\partial \mu}{\partial x}\bigg|_{x=L} \approx \frac{\mu_{N+1}^{i+1} - \mu_{N-1}^{i+1}}{2\Delta x} = 0 \quad \Rightarrow \mu_{N-1}^{i+1} = \mu_{N+1}^{i+1}
$$

Which gives us
$$
\mu_{N}^{i+1} = \frac{\partial f}{\partial \phi}\bigg|_{N}^{i+1} -2\kappa \frac{\phi_{N-1}^{i+1}- \phi_{N}^{i+1}}{(\Delta x)^{2}}
$$

and, for $\phi$, we similarly note that at $j=N$, $M_{N+1/2}^{i+1} = M_{N-1/2}^{i+1}$,  
$$
\phi_{N}^{i+1} - \phi_{N}^{i} = \frac{2\Delta t}{(\Delta x)^{2}}M_{N-1/2}^{i+1}\left( \mu_{N-1}^{i+1} - \mu_{N}^{i+1} \right)
$$
Putting things together such that we are only solving for $\phi$, i.e., treating things as the fourth order PDE directly, we get
$$
\begin{pmatrix}
\phi_{0}^{i+1} - \phi_{0}^{i} - \frac{2\Delta t}{(\Delta x)^{2}}\left( M_{1/2}^{i+1} \left( \mu_{1}^{i+1} - \mu_{0}^{i+1} \right) \right) \\
%
\vdots\\
%
\phi_{j}^{i+1} - \phi_{j}^{i} - \frac{\Delta t}{(\Delta x)^{2}} \left( M^{i+1}_{j+1/2}( \mu_{j+1}^{i+1} -\mu_{j}^{i+1}) - M_{j-1/2}^{i+1}(\mu_{j}^{i+1} -  \mu_{j-1}^{i+1}) \right) \\
%
\vdots \\
% 
\phi_{N}^{i+1} - \phi_{N}^{i} - \frac{2\Delta t}{(\Delta x)^{2}}M_{N-1/2}^{i+1}\left( \mu_{N-1}^{i+1} - \mu_{N}^{i+1} \right)
\end{pmatrix} = 0
$$

## Trying Method of Lines 

Again starting from the same equations
$$
\begin{align}
\frac{\partial \phi}{\partial t} &= \nabla \cdot (M \nabla \mu), \\
\mu &= \frac{\partial f}{\partial \phi} - \kappa \nabla^{2}\phi.
\end{align}
$$
But with the method of lines, we only employ a spatial discretization, while keeping time continuous. 

For the interior points, we have,
$$
\frac{d \phi_{j}}{\partial t} = \frac{1}{(\Delta x)^{2}} \left( M_{j+1/2}( \mu_{j+1} -\mu_{j}) - M_{j-1/2}(\mu_{j} -  \mu_{j-1}) \right),
$$
with $\mu_{j}$ as an auxiliary variable
$$
\mu_{j} = \frac{\partial f}{\partial \phi}_{j} - \frac{\kappa}{(\Delta x)^{2}} \left(\phi_{j+1} - 2\phi_{j} + \phi_{j-1}\right)
$$
For the left BC ($j = 0, 1$)

For the ghost node, we similarly note that $\phi_{-1} = \phi_{1}$,  at $j=0$, $M_{1/2} = M_{-1/2}$, which gives us
$$
\frac{d\phi_{0}}{dt} = \frac{2}{(\Delta x)^{2}} \left( M_{1/2}(\mu_{1} - \mu_{0}) \right),
$$

$$
\frac{d\phi_{1}}{dt} = \frac{1}{(\Delta x)^{2}}(M_{1+1/2}(\mu_{2} - \mu_{1}) - M_{1/2}(\mu_{1}-\mu_{0}))
$$

with 
$$
\mu_{0} = \frac{\partial f}{\partial \phi}_{0} - \frac{2\kappa}{(\Delta x)^{2}} (\phi_{1} - \phi_{0})
$$

$$
\mu_{1} = \frac{\partial f}{\partial \phi}_{1} -\frac{\kappa}{(\Delta x)^{2}}\left( \phi_{2} - 2\phi_{1} + \phi_{0} \right)
$$

$$
\mu_{2} = \frac{\partial f}{\partial \phi}_{2} - \frac{\kappa}{(\Delta x)^{2}}(\phi_{3} - 2\phi_{2} + \phi_{1})
$$

And for the right BC ($j = N-1, N$)

We similarly note that $\phi_{N+1} = \phi_{N-1}$ and $M_{N+1/2} = M_{N-1/2}$, which gives
$$
\frac{d\phi_{N}}{dt} = \frac{2}{(\Delta x)^{2}} \left(M_{N-1/2}(\mu_{N-1} - \mu_{N}) \right),
$$

$$
\frac{d\phi_{N-1}}{dt} = \frac{1}{(\Delta x)^{2}}\left(M_{N-1/2}(\mu_{N} - \mu_{N-1}) - M_{N-1-1/2}(\mu_{N-1} - \mu_{N-2}) \right)
$$

where
$$
\mu_{N} = \frac{\partial f}{\partial \phi}_{N} - \frac{2\kappa}{(\Delta x)^{2}}(\phi_{N-1} - \phi_{N}),
$$

$$
\mu_{N-1} = \frac{\partial f}{\partial \phi}_{N-1} - \frac{\kappa}{(\Delta x)^{2}}\left(\phi_{N} - 2\phi_{N-1} + \phi_{N-2} \right)
$$

