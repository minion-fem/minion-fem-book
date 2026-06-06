---
title: "J2 Plasticity"
---

# J2 Plasticity with Return Mapping

## Yield criterion

The J2 (von Mises) yield function with combined isotropic and kinematic hardening is:

$$f(\boldsymbol{\sigma}, \boldsymbol{\beta}, \alpha) = \|\boldsymbol{\xi}\| - \sqrt{\tfrac{2}{3}}\,K(\alpha) \leq 0$$

where $\boldsymbol{\xi} = \text{dev}(\boldsymbol{\sigma}) - \boldsymbol{\beta}$ is the overstress relative to the back-stress $\boldsymbol{\beta}$, and $K(\alpha) = \sigma_{y0} + K'\alpha$ is the isotropic yield radius.

The Voigt norm $\|\boldsymbol{\xi}\|$ uses the inner product $\boldsymbol{\xi}:\boldsymbol{\xi} = \xi_{ij}\xi_{ij}$, which in engineering-shear Voigt form is $\xi_{xx}^2 + \xi_{yy}^2 + \xi_{zz}^2 + \tfrac{1}{2}(\gamma_{xy}^2 + \gamma_{yz}^2 + \gamma_{xz}^2)$.

## Elastic predictor / plastic corrector

Given strain increment $\Delta\boldsymbol{\varepsilon}$ and old state $(\boldsymbol{\sigma}^n, \boldsymbol{\varepsilon}^p_n, \boldsymbol{\beta}^n, \alpha_n)$:

**1. Trial stress** (freeze plastic flow)

$$\boldsymbol{\sigma}^\text{tr} = \mathbf{C}^e : (\boldsymbol{\varepsilon}^n + \Delta\boldsymbol{\varepsilon} - \boldsymbol{\varepsilon}^p_n)$$

$$\boldsymbol{\xi}^\text{tr} = \text{dev}(\boldsymbol{\sigma}^\text{tr}) - \boldsymbol{\beta}^n, \qquad \xi^\text{tr}_\text{eq} = \|\boldsymbol{\xi}^\text{tr}\|$$

**2. Yield check** — if $\xi^\text{tr}_\text{eq} \leq \sqrt{\tfrac{2}{3}}\,K(\alpha_n)$, the step is elastic: accept the trial state.

**3. Plastic corrector** — project back to the yield surface by solving for the consistency parameter $\Delta\gamma$ (Newton iteration on a scalar residual):

$$g(\Delta\gamma) = \xi^\text{tr}_\text{eq} - 2G\Delta\gamma - \sqrt{\tfrac{2}{3}}\bigl[K(\alpha_{n+1}) + H(\alpha_{n+1}) - H(\alpha_n)\bigr] = 0$$

$$\alpha_{n+1} = \alpha_n + \sqrt{\tfrac{2}{3}}\Delta\gamma, \qquad \frac{dg}{d(\Delta\gamma)} = -2G\!\left(1 + \frac{K'+H'}{3G}\right)$$

**4. State update**

$$\mathbf{n}_{n+1} = \frac{\boldsymbol{\xi}^\text{tr}}{\xi^\text{tr}_\text{eq}}, \qquad \Delta\boldsymbol{\varepsilon}^p = \Delta\gamma\,\mathbf{n}_{n+1}$$

$$\boldsymbol{\varepsilon}^p_{n+1} = \boldsymbol{\varepsilon}^p_n + \Delta\boldsymbol{\varepsilon}^p, \qquad \boldsymbol{\beta}_{n+1} = \boldsymbol{\beta}_n + \sqrt{\tfrac{2}{3}}\Delta H\,\mathbf{n}_{n+1}$$

$$\boldsymbol{\sigma}_{n+1} = \text{hydro}(\boldsymbol{\sigma}^\text{tr}) + \text{dev}(\boldsymbol{\sigma}^\text{tr}) - 2G\Delta\gamma\,\mathbf{n}_{n+1}$$

## Consistent algorithmic tangent

The consistent tangent (required for quadratic Newton convergence) is:

$$\mathbf{C}^{ep} = \kappa\,\mathbf{1}\otimes\mathbf{1} + 2G\theta\!\left(\mathbf{I} - \tfrac{1}{3}\mathbf{1}\otimes\mathbf{1}\right) - 2G\bar\theta\,\mathbf{n}\otimes\mathbf{n}$$

where:

$$\theta = 1 - \frac{2G\Delta\gamma}{\xi^\text{tr}_\text{eq}}, \qquad \bar\theta = \frac{1}{1+\frac{K'+H'}{3G}} - (1-\theta)$$

The $\mathbf{n}\otimes\mathbf{n}$ correction makes the tangent depend on the direction of plastic flow; it reduces to the elastic tangent when $\Delta\gamma \to 0$ (continuous transition at the elastic–plastic boundary).
