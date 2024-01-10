# Perturbation Expansion Solution of Allen--Cahn Equation

Consider the Allen-Cahn equstion
$$
\frac{\partial f}{\partial t} = \kappa \nabla^{2}f + f - f^{3}, \quad f(t=0) = f^{0}(\mathbf{x})
$$
where $f$ is the phase-field variable and $\kappa$ is the gradient term. To develop a perturbation expansion solution for the Allen-Cahn equation, we introduce an $\epsilon$ term to the cubic nonlinearity to give
$$
\frac{\partial f}{\partial t} = \kappa \nabla^{2}f + f - \epsilon f^{3}.
$$
We consider a perturbation expansion solution of the form
$$
f(x,y,t) = f_{0} + \epsilon f_{1} + \epsilon^{2}f_{2} + \dots
$$
Introducing the perturbation expansion
$$
\frac{\partial}{\partial t} \left( f_{0} + \epsilon f_{1} + \epsilon^{2}f_{2} + \dots \right) = \\
%
\kappa \nabla^{2}\left( f_{0} + \epsilon f_{1} + \epsilon^{2}f_{2} + \dots \right) 
$$


