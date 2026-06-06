---
title: "Material Models"
slug: material-models
---

# Material Models

Constitutive models relate stress to strain — they encode the physical behaviour of a material. This section covers models ranging from simple linear elasticity to rate-dependent viscoplasticity, all formulated with consistent algorithmic tangents for quadratic Newton convergence.

- [Linear Elasticity](elasticity.md) — Lamé formulation, serves as elastic predictor in all inelastic models
- [J2 Plasticity](j2_plasticity.md) — von Mises yield criterion, return mapping, isotropic and kinematic hardening
- [Perzyna Viscoplasticity](perzyna.md) — overstress model, rate-dependent flow outside the yield surface
- [Anand Viscoplasticity](viscoplasticity.md) — Arrhenius-sinh flow rule, temperature-dependent hardening, consistent tangent
