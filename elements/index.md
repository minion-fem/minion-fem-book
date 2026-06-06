---
title: "Finite Element"
slug: finite-element
---

# Finite Element

The finite element method discretises a continuous PDE into a system of algebraic equations by dividing the domain into elements and approximating the solution with piecewise polynomial basis functions. This section covers element formulations from first principles — shape functions, numerical integration, and the modifications needed for large deformation and nearly incompressible materials.

- [Element Formulation](formulation.md) — isoparametric mapping, shape functions, strain-displacement matrix **B**
- [Geometric Nonlinearity](nlgeom.md) — Updated Lagrangian, Hughes-Winget stress rotation, geometric stiffness
- [Numerical Integration](integration.md) — Gauss quadrature rules for triangles, quads, tetrahedra, hexahedra
- [Volumetric Locking](locking.md) — why linear hex locks for incompressible materials, and the B-bar fix
