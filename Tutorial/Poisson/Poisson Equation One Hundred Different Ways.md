# Poisson Equation: One Hundred Different Ways

## Linear Poisson Equation with Pure Neumann Boundary Conditions

Consider the PDE
$$
-\nabla^{2} u = f, \quad \mathbf{x} \in V, \\
\nabla u \cdot n = g, \quad \mathbf{x} \in \partial V.
$$
The solution is determined up to a constant $c$. Physically, imagine a plate / disc with no-flux boundary conditions i.e., the heat/mass cannot flow out of the domain, hence the temperature / concentration is only defined to a constant. Hence another constraint on the system is necessary. For this example, we shall specify
$$
\int_{V} u \ dV = 0.
$$
Incorporating this constraint directly into the code is not so straightforward. We shall proceed by constructing the weak form using a Lagrange multiplier. 

We can define the "energy" functional $\mathcal{J}$, as
$$
\mathcal{J}(u) = \frac{1}{2}\int_{V} |\nabla u|^{2} \ dV - \int_{V}fu \ dV - \int_{\partial V} gu \ d\partial V
$$
The form of the function $\mathcal{J}$ is makes sense when considering the variational derivative for the domain: (i.e., $V$),
$$
\frac{\delta \mathcal{J}}{\delta u} = -\nabla \cdot \nabla u -f 
$$
The solution to the PDE will minimize the energy since it satisfies $\frac{\delta \mathcal{J}}{\delta u} = 0$. We can append the constraint with a Lagrange multiplier since it can be seen as an equality constraint.
$$
\mathcal{L}(u,\lambda) = \mathcal{J}(u) + \lambda \int_{V} u \ dV.
$$
Taking the Gateaux derivative with respect to $\delta u$ and $\delta \lambda$, 
$$
\delta \mathcal{L} = \int_{V} \nabla u \cdot \nabla \delta u \ dV - \int_{V}f \delta u \ dV - \int_{\partial V} g \delta u \  d\partial V + \delta \lambda \int_{V} u \  dV + \lambda \int_{V} \delta u \ dV
$$
Replacing $\delta u$ with our test function $v$ and $\delta \lambda$ with another test function $w$, we get a weak form that can be declared in FEniCS



## Linear Poisson Equation with Integral Term

Consider the PDE
$$
- \nabla^{2}u + \int_{V} u \ dV = f(\mathbf{x}), \quad \mathbf{x} \in V \\
 u = 0, \quad \mathbf{x} \in \partial V
$$
This is a Poisson equation with an integral term and zero Dirichlet boundary conditions on the domain.  We note this PDE is linear. Consider the homogeneous case where $f(\mathbf{x}) = 0$, 
$$
\nabla^{2} u = \int_{V} u \ dV
$$
Suppose $u_{1}$ and $u_{2}$ are solutions to the PDE given above, then $u_{1} + u_{2}$ must be a solution if the PDE is linear. This can be show by direct substitution,
$$
\begin{align}
\text{LHS} &=\nabla^{2}(u_{1} + u_{2}) \\
%
&= \nabla^{2}u_{1} + \nabla^{2}u_{2} \\
&= \int_{V} u_{1} \ dV + \int_{V} u_{2} \ dV \\
&= \int_{V} (u_{1} + u_{2}) \ dV = \text{RHS} \\
\end{align}
$$
This result is not surprising as both the Laplacian and the integral operators are linear. Directly declaring the integero-partial differential equation into FEniCS is not so simple by virtue of the integral term. We can instead conceptualize the solution method as follows: 

1. Neglecting the integral term, the PDE is

$$
-\nabla^{2}u  = f
$$

Expressing the equation in the weak form, we get
$$
\int_V \nabla u \cdot \nabla v \ dV - \int_{V} fv \ dV = 0
$$
The solution 