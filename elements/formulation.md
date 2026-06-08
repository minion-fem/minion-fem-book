---
title: "Element Formulation"
---

# Element Formulation

All standard isoparametric elements share the same structure: physical coordinates $\mathbf{x}$ are interpolated from nodal coordinates $\mathbf{x}_I$ using the same shape functions $N_I(\boldsymbol{\xi})$ that interpolate the displacement field:

$$\mathbf{x}(\boldsymbol{\xi}) = \sum_I N_I(\boldsymbol{\xi})\,\mathbf{x}_I, \qquad \mathbf{u}(\boldsymbol{\xi}) = \sum_I N_I(\boldsymbol{\xi})\,\mathbf{u}_I$$

The Jacobian of the isoparametric mapping and the spatial shape-function gradients are:

$$\mathbf{J} = \frac{\partial N}{\partial \boldsymbol{\xi}}\,\mathbf{x}, \qquad \frac{\partial N}{\partial \mathbf{x}} = \mathbf{J}^{-1}\frac{\partial N}{\partial \boldsymbol{\xi}}, \qquad dV = \det(\mathbf{J})\,w$$

where $w$ is the integration weight.

---

## 3D continuum elements

### C3D4 — 4-node linear tetrahedron

Natural coordinates $(\xi_1, \xi_2, \xi_3)$ are volume coordinates. Shape functions are:

$$N_1 = \xi_1,\quad N_2 = \xi_2,\quad N_3 = \xi_3,\quad N_4 = 1 - \xi_1 - \xi_2 - \xi_3$$

The gradients are constant over the element (constant-strain tetrahedron). One integration point is used (degree-1 exact quadrature over $T_3$).

**Use when:** meshing complex geometry automatically. Cheaper but stiffer than C3D8; requires finer meshes for bending-dominated problems.

### C3D8 — 8-node linear hexahedron

Natural coordinates $(\xi, \eta, \zeta) \in [-1,1]^3$. Shape functions have the trilinear form:

$$N_I(\xi,\eta,\zeta) = \frac{1}{8}(1+\xi_I\xi)(1+\eta_I\eta)(1+\zeta_I\zeta)$$

where $(\xi_I, \eta_I, \zeta_I) \in \{-1,+1\}^3$ are the nodal natural coordinates.

Full $2 \times 2 \times 2$ Gauss quadrature (8 integration points, degree 3 exact) is used for the stiffness matrix. Under large deformation (NLGEOM), the Jacobian is evaluated at the current configuration. The implementation uses the **B-bar method** to prevent volumetric locking for nearly incompressible materials — see [Volumetric Locking](locking.md).

**Use when:** structured meshes or problems where mesh quality can be controlled. Better accuracy per degree of freedom than C3D4 for smooth solutions.

---

## 2D continuum elements

These are plane-strain (CPE) and plane-stress (CPS) formulations. The strain state is reduced: plane strain constrains $\varepsilon_{zz} = \varepsilon_{xz} = \varepsilon_{yz} = 0$; plane stress constrains $\sigma_{zz} = \sigma_{xz} = \sigma_{yz} = 0$.

| Element | Nodes | Integration pts |
|---|---|---|
| CPS3 / CPE3 | 3 (linear triangle) | 1 |
| CPS4 / CPE4 | 4 (bilinear quad) | 4 |

CPS4 is used for 2D topology optimisation.

---

## Beam elements

| Element | Nodes | Formulation |
|---|---|---|
| B31 | 2 | Timoshenko (first-order, 1 integration pt) |
| B32 | 3 | Timoshenko (second-order, 2 integration pts) |

Section properties (area, moments of inertia) are precomputed and stored per element. The element frame rotates with the beam axis; local stiffness is transformed into the global system using direction cosines.

---

## Strain-displacement matrix $\mathbf{B}$

For 3D continuum elements, the per-node block of $\mathbf{B}$ is built from the spatial shape-function gradient $\nabla N_a = (N_{a,x},\, N_{a,y},\, N_{a,z})$:

$$\mathbf{B}_a = \begin{bmatrix}
N_{a,x} & 0 & 0 \\
0 & N_{a,y} & 0 \\
0 & 0 & N_{a,z} \\
N_{a,y} & N_{a,x} & 0 \\
0 & N_{a,z} & N_{a,y} \\
N_{a,z} & 0 & N_{a,x}
\end{bmatrix}$$

Voigt ordering: $[\varepsilon_{xx},\, \varepsilon_{yy},\, \varepsilon_{zz},\, \gamma_{xy},\, \gamma_{yz},\, \gamma_{xz}]$ with engineering shear $\gamma_{ij} = 2\varepsilon_{ij}$.

The element stiffness matrix is then:

$$\mathbf{K}^e = \int_V \mathbf{B}^T \mathbf{C}\, \mathbf{B}\, dV \approx \sum_g w_g \det(\mathbf{J}_g)\, \mathbf{B}_g^T \mathbf{C}_g\, \mathbf{B}_g$$

Under geometric nonlinearity, $\mathbf{B}$ is evaluated at the current deformed configuration and a geometric stiffness term $\mathbf{K}_\text{geo}$ is added.
