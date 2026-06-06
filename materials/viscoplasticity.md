---
title: "Anand Viscoplasticity"
---

# Anand Viscoplasticity

The Anand model is a rate-dependent constitutive law for metals and solder alloys at elevated temperature. Unlike J2 plasticity — which has a sharp yield surface — viscoplastic flow is allowed at any stress level; the rate depends on how far the stress exceeds the current deformation resistance. This makes it suitable for creep, stress relaxation, and thermomechanical fatigue.

A single scalar internal variable $s$ (the **deformation resistance**, Pa) captures the entire inelastic history. There is no yield surface and no elastic/plastic split; every increment potentially involves viscoplastic flow.

---

## Constitutive equations

### Flow rule

Viscoplastic flow is directed along the deviatoric stress, with a scalar rate governed by an Arrhenius-sinh law:

$$\dot{\boldsymbol{\varepsilon}}^{vp} = \sqrt{\tfrac{3}{2}}\;\dot{\bar{\varepsilon}}^{vp}\;\frac{\boldsymbol{\sigma}'}{\lVert\boldsymbol{\sigma}'\rVert}$$

$$\dot{\bar{\varepsilon}}^{vp} = A\,\exp\!\left(-\frac{Q}{k_B T}\right) \left[\sinh\!\left(\xi\,\frac{\lVert\boldsymbol{\sigma}'\rVert}{s}\right)\right]^{1/m} \equiv f(\lVert\boldsymbol{\sigma}'\rVert, s, T)$$

The sinh term captures both the low-stress (linear creep) and high-stress (power-law creep) regimes in a single expression. At low stress $\sinh(x) \approx x$, recovering a linear-viscous response; at high stress $\sinh(x) \approx \frac{1}{2}e^x$, recovering power-law creep.

| Parameter | Unit | Meaning |
|---|---|---|
| $A$ | 1/s | Pre-exponential rate constant |
| $Q$ | eV | Activation energy |
| $k_B$ | eV/K | Boltzmann constant |
| $T$ | K | Absolute temperature |
| $\xi$ | — | Stress multiplier (sinh argument scaling) |
| $m$ | — | Strain-rate sensitivity exponent |

### Hardening rule

The deformation resistance $s$ evolves toward a strain-rate- and temperature-dependent saturation value $s^*$:

$$\dot{s} = h_0\,\left|\,1 - \frac{s}{s^*}\right|^{a}\;\operatorname{sgn}\!\left(1 - \frac{s}{s^*}\right)\;\dot{\bar{\varepsilon}}^{vp}$$

$$s^* = \tilde{s}\left(\frac{\dot{\bar{\varepsilon}}^{vp}}{A_{\mathrm{eff}}}\right)^n, \quad A_{\mathrm{eff}} = A\exp\!\left(-\frac{Q}{k_BT}\right)$$

When $s < s^*$ the material hardens; when $s > s^*$ it softens. The exponents $a$ and $n$ control how quickly $s$ converges to $s^*$.

---

## Stress update algorithm

### Collinearity and scalar reduction

The full 6-component tensor update reduces to a 2-scalar system. The key observation is that the updated deviatoric stress is **collinear** with the trial deviatoric stress:

$$\boldsymbol{\sigma}'_{\mathrm{new}} = \boldsymbol{\sigma}'_{\mathrm{tr}} - 3G\,\Delta\bar{\varepsilon}^{vp}\,\mathbf{n}, \quad \mathbf{n} = \sqrt{\tfrac{3}{2}}\frac{\boldsymbol{\sigma}'_{\mathrm{tr}}}{\lVert\boldsymbol{\sigma}'_{\mathrm{tr}}\rVert}$$

This is the same **radial-return** geometry as J2 plasticity, but without a yield surface — $\Delta\bar{\varepsilon}^{vp} = f\,\Delta t$ is always nonzero. Taking norms reduces the update to a 2×2 nonlinear system in the equivalent stress $\bar\sigma = \lVert\boldsymbol{\sigma}'\rVert$ and resistance $s$:

$$\boxed{\begin{aligned}
R_1(\bar\sigma, s) &= \bar\sigma - \bar\sigma_{\mathrm{tr}} + 3G\,f(\bar\sigma, s)\,\Delta t = 0 \\
R_2(\bar\sigma, s) &= s - s_{\mathrm{old}} - g(\bar\sigma, s)\,\Delta t = 0
\end{aligned}}$$

### Newton-Raphson iteration

The 2×2 system is solved by Newton's method with analytical Jacobian entries:

$$J_{11} = 1 + 3G\,\Delta t\,\frac{\partial f}{\partial\bar\sigma}, \quad J_{12} = 3G\,\Delta t\,\frac{\partial f}{\partial s}$$

$$J_{21} = -\left(\frac{n\xi a}{m(s^*-s)}\coth\!\left(\frac{\xi\bar\sigma}{s}\right)g + \frac{\partial f}{\partial\bar\sigma}\,h\right)\Delta t, \quad J_{22} = 1 - \left(h\frac{\partial f}{\partial s} - \frac{a}{s^*-s}hf\right)\Delta t$$

where the partial derivatives of $f$ are:

$$\frac{\partial f}{\partial\bar\sigma} = \frac{\xi}{ms}\coth\!\left(\frac{\xi\bar\sigma}{s}\right)f, \qquad \frac{\partial f}{\partial s} = -\frac{\xi\bar\sigma}{ms^2}\coth\!\left(\frac{\xi\bar\sigma}{s}\right)f$$

### Cutback strategy

After each Newton correction, the algorithm checks for:

1. Negative equivalent stress $(\bar\sigma + \delta\bar\sigma < 0)$
2. Non-positive deformation resistance $(s + \delta s \leq 0)$
3. Unrealistically large stress-to-resistance ratio

If any condition is triggered the step is halved, up to 20 times. This protects convergence near the initial state or after large load jumps.

---

## Consistent algorithmic tangent

To achieve quadratic Newton convergence at the structural level, the material tangent $\mathbf{C}^{ep} = \partial\boldsymbol{\sigma}/\partial\boldsymbol{\varepsilon}$ must be **consistent** with the discrete stress update — not the continuum rate-form tangent.

Differentiating the converged update with respect to $\Delta\boldsymbol{\varepsilon}$ and applying the implicit-function theorem to the Newton residuals gives:

$$\frac{\partial\Delta\bar{\varepsilon}^{vp}}{\partial\boldsymbol{\varepsilon}} = \beta\,\hat{\mathbf{n}}, \qquad \beta = \frac{2}{3}\!\left(1 - \frac{J_{22}}{\det\mathbf{J}}\right)$$

The final tangent has the isotropic Mandel form:

$$\boxed{\mathbf{C}^{ep} = K\,\mathbf{I}\otimes\mathbf{I} + 2G_{\mathrm{eff}}\,\mathbb{I}_4^{\mathrm{dev}} + \alpha\,\hat{\mathbf{n}}\otimes\hat{\mathbf{n}}}$$

$$G_{\mathrm{eff}} = G\!\left(1 - \frac{3G\,\Delta\bar{\varepsilon}^{vp}}{\bar\sigma_{\mathrm{tr}}}\right), \qquad \alpha = \frac{6G^2\,\Delta\bar{\varepsilon}^{vp}}{\bar\sigma_{\mathrm{tr}}} - 2G\!\left(1 - \frac{J_{22}}{\det\mathbf{J}}\right)$$

In the elastic limit $(\Delta\bar{\varepsilon}^{vp}\to 0)$: $G_{\mathrm{eff}}\to G$, $\alpha\to 0$, recovering the isotropic elastic tangent continuously.

---

## Comparison with J2 plasticity

| Feature | J2 Plasticity | Anand Viscoplasticity |
|---|---|---|
| Internal variables | $\boldsymbol{\varepsilon}^p$, $\alpha$, $\boldsymbol{\beta}$ | $\bar{\varepsilon}^{vp}$, $s$ |
| Yield surface | Sharp — elastic/plastic split | None — always viscoplastic |
| Flow rate | Instantaneous at yield | $\dot{\bar{\varepsilon}}^{vp} = A_{\mathrm{eff}}[\sinh(\xi\sigma/s)]^{1/m}$ |
| Temperature dependence | None (in Minion's J2) | Arrhenius factor $\exp(-Q/k_BT)$ |
| Newton system | Scalar $\Delta\gamma$ equation | 2×2 system in $(\bar\sigma, s)$ |
| Typical use | Structural metals, rate-independent | Solder alloys, polymers, creep |
