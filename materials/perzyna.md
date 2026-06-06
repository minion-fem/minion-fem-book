---
title: "Perzyna Viscoplasticity"
---

# Perzyna Viscoplasticity

For rate-dependent materials, the Perzyna overstress model allows plastic flow outside the rate-independent yield surface. The flow rate is controlled by how far the stress exceeds the yield surface:

$$\dot{\boldsymbol{\varepsilon}}^{vp} = \frac{1}{\eta}\left\langle \frac{f(\boldsymbol{\sigma})}{\sigma_{y0}} \right\rangle^m \mathbf{n}$$

where $\langle \cdot \rangle$ are Macaulay brackets (zero when negative), $\eta$ is the viscosity, $m$ is the rate-sensitivity exponent, and $\mathbf{n} = \partial f / \partial \boldsymbol{\sigma}$ is the flow direction.

Unlike rate-independent plasticity, there is no consistency condition: the yield surface can be violated, and the overstress $f(\boldsymbol{\sigma}) > 0$ drives the flow rate. In the limit $\eta \to 0$ the rate-independent J2 model is recovered.

For a more complete rate-dependent model with Arrhenius-sinh kinetics, temperature dependence, and evolving deformation resistance, see [Anand Viscoplasticity](viscoplasticity.md).
