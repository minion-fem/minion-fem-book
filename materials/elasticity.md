---
title: "Linear Elasticity"
---

# Isotropic Linear Elasticity

The stress–strain relation in Voigt notation is $\boldsymbol{\sigma} = \mathbf{C}^e \boldsymbol{\varepsilon}$, where the 6×6 material Jacobian $\mathbf{C}^e$ is the standard Lamé matrix parameterised by Young's modulus $E$ and Poisson's ratio $\nu$:

$$C^e_{ijkl} = \lambda\,\delta_{ij}\delta_{kl} + \mu(\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk})$$

$$\lambda = \frac{E\nu}{(1+\nu)(1-2\nu)}, \qquad \mu = \frac{E}{2(1+\nu)}$$

This is the base elastic model; it also serves as the elastic predictor inside the inelastic routines.
