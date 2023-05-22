# Nucleation Benchmark

## Derivation

This nucleation benchmark is based on the work of Wu et al., 2021 and considers nucleation using a single Allen-Cahn type phase field model. 

The total free energy of the system is given by 
$$
F(c) = \int_{V} \frac{\kappa}{2} (\nabla c)^{2} + w c^{2}(1-c)^{2} - \Delta f \left( c^{3} (10 - 15c + 6c^{2}) \right) \ dV,
$$
where $\kappa$ is the gradient energy parameter, $w$ is a parameter that controls the depth of the double well, $\Delta f$ is the driving force for solidification and is positive below the melting point, and the 5th order polynomial is a monotonically increasing function to interpolate the bulk free energies of the solid and the liquid. In the literature, this term is expressed in a few ways (see 1-2),
$$
p(c)f_{S} + (1-p(c))f_{L},
$$

$$
\Delta f p(c),
$$

where $p(c)$ is the 5th order polynomial in equation (1), and $f_{S}$ and $f_{L}$ are the free energies of the solid and liquid bulk phases respectively. The evolution of $c$ is governed by the Allen-Cahn equation,
$$
\frac{\partial c}{\partial t} = -M \frac{\delta F}{\delta c} = M \left( \kappa \nabla^{2}c - w\frac{dg}{dc} + \Delta f \frac{dp}{dc}  \right),
$$
where $g = c^{2}(1-c)^{2}$. 

## Relevant Features

To set up the simulation, we need to explore some theoretical features of classical nucleation theory and the Allen-Cahn equation. 

It can be shown that the interfacial tension $\gamma$ can be analytically derived and is expressed in terms of $\kappa$ and $w$,
$$
\gamma = \frac{\sqrt{\kappa w}}{3\sqrt{2}}.
$$
For a 2D system, we can write down the free energy of a circular solid nucleus, 
$$
\Delta G (r) = \underbrace{2\pi r \gamma}_{\text{Penalty for forming interface}} - \underbrace{\pi r^{2} \Delta f}_{\text{Driving force for forming bulk} }
$$
To identify the critical radius $r^{*}$, we set $\frac{d\Delta G}{dr} = 0$,
$$
2\pi \gamma - 2\pi r^{*}\Delta f = 0,
$$

$$
r^{*} = \frac{\gamma}{\Delta f} = \frac{\sqrt{\kappa w}}{3\sqrt{2} \Delta f}.
$$

Depending on the initial size of the nucleus $r_{0}$, a few things can happen,

1. $r_{0} = r^{*}$: The nucleus will remain static
2. $r_{0} > r^{*}$: The nucleus will grow
3. $r_{0} < r^{*}$: The nucleus will shrink and disappear

Correspondingly, we would also like to track the evolution of the solid fraction $Y(t)$ which has a theoretical relation,
$$
Y(t) = 1 - \exp(-Kt^{n}),
$$
where $K$ is a constant depending on the nucleation and growth rates and $n$ is related to the spatial dimensions of the problem. We can plot,
$$
\ln{\left(1 - Y(t)\right)} = -Kt^{n},
$$
to analyze the results. 

## Simulation Parameters

Recall the original equation,
$$
\frac{\partial c}{\partial t} = M \left( \kappa \nabla^{2}c - w\frac{dg}{dc} + \Delta f \frac{dp}{dc}  \right).
$$
The equation can be non-dimensionalized by introducing the scalings
$$
\mathbf{x} = \sqrt{\frac{\kappa}{w}} \tilde{\mathbf{x}}, \quad t = \frac{1}{Mw} \tilde{t},
$$
to give,
$$
Mw \frac{\partial c}{\partial \tilde{t}} = M\kappa \frac{w}{\kappa}\tilde{\nabla}^{2}c - Mw\frac{dg}{dc} + M\Delta f \frac{dp}{dc},
$$

$$
\frac{\partial c}{\partial \tilde{t}} = \tilde{\nabla}^{2}c - \frac{dg}{dc} + \frac{\Delta f}{w} \frac{dp}{dc}.
$$

By non-dimensionalizing the model, we have reduced the parameter space down to a single parameter $\frac{\Delta f}{w}$. $r^{*}$ can be similarly non-dimensionalized,
$$
\sqrt{\frac{\kappa}{w}} \tilde{r}^{*} = \frac{\sqrt{\kappa w}}{3\sqrt{2} \Delta f}
$$

$$
\tilde{r}^{*} = \frac{1}{3\sqrt{2}}\frac{w}{\Delta f}
$$

For the simulation, we shall specify $\frac{\Delta f}{w} = \frac{\sqrt{2}}{30}$, which gives $\tilde{r}^{*} = 5$. A total of three simulation cases are considered: $\tilde{r}_{0} = \{ \tilde{r}^{*}, 0.99\tilde{r}^{*}, 1.01\tilde{r}^{*} \}$. For the rest of the note, we will drop the ~ above the symbols, but note that we are still working in dimensionless variables. 

The simulation uses the following parameters:

1. $t \in [0, 100]$
2. nx = ny = 100
3. $V \in [-50,50]^{2}$ 
4. $c(t=0) = \frac{1}{2} \left( 1 - \tanh{\frac{\sqrt{x^{2 + y^{2}}} -r_{0}}{\sqrt{2}}} \right)$

## Weak Form

Recall the model equation,
$$
\frac{\partial c}{\partial t} = \nabla^{2}c - \frac{dg}{dc} + \frac{\Delta f}{w} \frac{dp}{dc}.
$$
Discretizing the time derivative with a backwards difference scheme and treating all terms implicitly,
$$
\frac{c^{n+1} -c^{n}}{\Delta t} = \nabla^{2}c^{n+1} - g'(c^{n+1}) + \frac{\Delta f}{w}p'(c^{n+1}) 
$$
Introducing a test function $v$ and integrating,
$$
F = \int_{V} c^{n+1}v - c^{n}v - \Delta t \nabla c^{n+1} \cdot \nabla v + \Delta t g'(c^{n+1})v  - \frac{\Delta f}{w}p'(c^{n+1})v \ dV = 0
$$
We can implement the above into FEniCS. 

## Post Processing

The post processing pipeline in Paraview to extract the solid fraction:

1. Apply the `ExtractBlock` filter to convert the partial dataset into normal cell data
2. Apply a `Clip` filter with a `Scalar` specified for the `Clip Type`. The value set should be 0.5. Do not select the `invert` option. 
3. Apply an `IntegrateVariables` filter. This gives the area of the solid
4. Apply a `Calculator` filter to divide the area by the total volume 
5. Apply the `PlotDataOverTime` filter to obtain the evolution profile. 

To export an animation, click `File`, `Save Animation` and save the `.png` accordingly. Use the following line in the terminal to crop,

```bash
mogrify -trim *.png
```

## References

1. Wu, W., Montiel, D., Guyer, J.E., Voorhees, P.W., Warren, J.A., Wheeler, D., Gránásy, L., Pusztai, T. and Heinonen, O.G., 2021. Phase field benchmark problems for nucleation. *Computational Materials Science*, *193*, p.110371.
2. Takaki, T., 2014. Phase-field modeling and simulations of dendrite growth. *ISIJ international*, *54*(2), pp.437-444.

