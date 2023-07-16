# Anisotropic Dendritic Solidification

## Model Formulation

We are modeling the solidification of a liquid material by cooling. The Allen-Cahn equation is coupled to a thermal transport equation since we are modeling the effect of cooling on solidification. 

The heat transfer equation is
$$
\frac{\partial \Delta T}{\partial t} = D_{T} \nabla^{2}\left(\Delta T \right) + \frac{\partial c}{\partial t}
$$
and the phase field equation is 
$$
\tau_{c}\frac{\partial c}{\partial t} = \nabla \cdot \left(  \mathbf{D} \nabla c \right) + c(1-c)m(c, \Delta T),
$$
where 
$$
m(c,\Delta T) = c - \frac{1}{2} - \frac{k_{1}}{\pi} \arctan{\left( k_{2} \Delta T \right)},
$$

$$
\mathbf{D} = \alpha^{2}(1 + k_{3}\beta) \pmatrix{1+ k_{3}\beta & -k_{3} \frac{\partial \beta}{\partial \psi} \\ k_{3}\frac{\partial \beta}{\partial \psi} & 1 + k_{3}\beta},
$$

$$
\beta = \frac{1- \Phi^{2}}{1 + \Phi^{2}},
$$

$$
\Phi = \tan{\left( \frac{N}{2} \psi \right)},
$$

$$
\psi = k_{4} + \arctan{\left( \frac{ \frac{\partial c}{\partial y} }{\frac{\partial c}{\partial x}} \right)}
$$

The following nominal parameters are used

| Parameter  | Value                |
| ---------- | -------------------- |
| $\tau_{c}$ | $3$$\times$$10^{-4}$ |
| $D_{T}$    | 2.25                 |
| $k_{1}$    |                      |
|            |                      |

