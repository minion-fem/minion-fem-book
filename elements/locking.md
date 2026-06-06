---
title: "Volumetric Locking and B-bar"
---

# Volumetric Locking and the B-bar Method

Linear hexahedral elements (C3D8) suffer from **volumetric locking** when the material is nearly incompressible ($\nu \to 0.5$). This page explains why, and how Minion's B-bar formulation eliminates it.

---

## What is volumetric locking?

Consider a single linear hex element under pure shear. The standard displacement interpolation $\mathbf{u} = \mathbf{N}\mathbf{d}$ with trilinear shape functions can represent the correct displacement field exactly. Now apply a uniform pressure. The correct response is pure volumetric compression — zero distortion.

For a nearly incompressible material, the volumetric strain $\varepsilon_v = \varepsilon_{xx} + \varepsilon_{yy} + \varepsilon_{zz}$ is energetically very costly (bulk modulus $\kappa = E/[3(1-2\nu)] \to \infty$ as $\nu \to 0.5$). The element displacement field is then **over-constrained**: it cannot simultaneously satisfy the traction boundary conditions and the near-incompressibility constraint with the 8 nodal degrees of freedom available. The element becomes artificially stiff — it **locks**.

The symptom is a displacement result that is orders of magnitude too small and a pressure distribution that oscillates wildly between integration points.

### Why C3D4 does not have this problem

Tetrahedral elements have a constant strain field (one integration point). The volumetric constraint is enforced at the element level as a single scalar equation, leaving enough freedom for the element to deform correctly. For meshes with C3D8 elements and $\nu > 0.45$, locking is often severe enough to invalidate results.

---

## The B-bar method

The B-bar approach, due to Hughes (1980), replaces the volumetric part of the strain-displacement matrix $\mathbf{B}$ with an **element-averaged** version. This allows each integration point to "see" the same volumetric strain as the element as a whole, decoupling the volumetric and deviatoric responses.

### Decomposition of B

Split the B matrix into volumetric and deviatoric parts:

$$\mathbf{B} = \mathbf{B}^{\mathrm{dev}} + \mathbf{B}^{\mathrm{vol}}$$

The volumetric part at integration point $g$ for node $a$ is:

$$B^{\mathrm{vol}}_{ij,a} = \frac{1}{3}\bar{\delta}_{ij}\,N_{a,k}\delta_{kj} = \frac{1}{3}\delta_{ij}\,(\nabla N_a \cdot \mathbf{e}_j)$$

In Voigt notation the first three rows (normal strains) of $\mathbf{B}_a$ contain the volumetric coupling; the last three rows (shear strains) do not.

### Averaged volumetric gradient

Define the element-averaged spatial gradient:

$$\bar{\nabla} N_a = \frac{1}{V^e} \int_{V^e} \nabla N_a \, dV \approx \frac{\sum_g w_g \det(\mathbf{J}_g)\,\nabla N_a|_g}{\sum_g w_g \det(\mathbf{J}_g)}$$

The B-bar matrix at integration point $g$ replaces the local $\nabla N_a$ in the normal-strain rows with $\bar\nabla N_a$:

$$\bar{B}^{(g)}_{a} = \begin{bmatrix}
\tfrac{2}{3}N_{a,x} + \tfrac{1}{3}\bar{N}_{a,x} & -\tfrac{1}{3}N_{a,y} + \tfrac{1}{3}\bar{N}_{a,y} & -\tfrac{1}{3}N_{a,z} + \tfrac{1}{3}\bar{N}_{a,z} \\
-\tfrac{1}{3}N_{a,x} + \tfrac{1}{3}\bar{N}_{a,x} & \tfrac{2}{3}N_{a,y} + \tfrac{1}{3}\bar{N}_{a,y} & -\tfrac{1}{3}N_{a,z} + \tfrac{1}{3}\bar{N}_{a,z} \\
-\tfrac{1}{3}N_{a,x} + \tfrac{1}{3}\bar{N}_{a,x} & -\tfrac{1}{3}N_{a,y} + \tfrac{1}{3}\bar{N}_{a,y} & \tfrac{2}{3}N_{a,z} + \tfrac{1}{3}\bar{N}_{a,z} \\
N_{a,y} & N_{a,x} & 0 \\
0 & N_{a,z} & N_{a,y} \\
N_{a,z} & 0 & N_{a,x}
\end{bmatrix}$$

where $N_{a,x} = \partial N_a/\partial x|_g$ is the local (integration-point) gradient and $\bar{N}_{a,x} = \partial \bar N_a/\partial x$ is the averaged gradient. Shear rows are unchanged.

The element stiffness is then assembled with $\bar{\mathbf{B}}$ in place of $\mathbf{B}$:

$$\mathbf{K}^e = \sum_g w_g \det(\mathbf{J}_g)\; \bar{\mathbf{B}}_g^T\,\mathbf{C}\,\bar{\mathbf{B}}_g$$

---

## Implementation

The averaged gradient $\bar{\nabla}N_a$ is precomputed once per element and cached (`compute_averaged_dndx()`). The B-bar matrix is then assembled per integration point using the formula above (`_update_bee_matrix()`). For linear-elastic problems a vectorised batch version processes all section elements simultaneously with NumPy einsum.

Under **geometric nonlinearity** (NLGEOM), the averaged gradient is recomputed at each Newton iteration because the current configuration changes, and B-bar is applied on top of spatial gradients evaluated at the deformed configuration.

---

## Effect on results

For compressible materials ($\nu < 0.4$), standard and B-bar elements give nearly identical results. The difference becomes large near $\nu = 0.5$:

| $\nu$ | Standard C3D8 | B-bar C3D8 |
|---|---|---|
| 0.3 | Correct | Correct |
| 0.45 | Slightly stiff | Correct |
| 0.499 | Severely locked | Correct |
| 0.4999 | Fails (near-singular $\mathbf{K}$) | Correct |

Using a linear hex without B-bar for rubber-like or metallic plasticity problems ($\nu \approx 0.5$ in the plastic regime) produces incorrect results. The C3D8 implementation here always uses B-bar.
