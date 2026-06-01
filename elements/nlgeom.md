---
title: "Geometric Nonlinearity (NLGEOM)"
---

# Geometric Nonlinearity in C3D Elements

---

## 1. Governing Equations

### 1.1 Virtual work principle

The equilibrium equation is obtained from the principle of virtual work (where $\delta(\cdot)$ denotes a virtual variation, $\sigma$ the Cauchy stress, and $\mathbf{D}$ the rate-of-deformation tensor — symmetric part of the spatial velocity gradient $\mathbf{L} = \partial\mathbf{v}/\partial\mathbf{x}$):

$$ 0 = \delta\Pi(\mathbf{u}) = \int_V \sigma : \delta\mathbf{D}\,dV  - \int_S \mathbf{t}^T \delta\mathbf{v}\,dS  - \int_V \mathbf{f}^T \delta\mathbf{v}\,dV $$

Equivalently, using Kirchhoff stress $\tau = J\sigma$ (where $J = |dV/dV^0|$ is the volume ratio) and integrating over the reference volume $V^0$:

$$0 = \int_{V^0} \tau : \delta\mathbf{D}\,dV^0  - \int_S \mathbf{t}^T \delta\mathbf{v}\,dS  - \int_V \mathbf{f}^T \delta\mathbf{v}\,dV$$

### 1.2 Linearization (Newton iteration)

Newton's method requires linearizing $\delta\Pi(\mathbf{u}+d\mathbf{u})$ about $\mathbf{u}$, where $d(\cdot)$ denotes the increment within the current Newton iteration:

$$0 \approx \delta\Pi(\mathbf{u})+ \int_{V^0} \bigl(d\tau : \delta\mathbf{D} + \tau : d\delta\mathbf{D}\bigr)\,dV^0 - (\text{load-stiffness terms})$$

The element stiffness contribution comes from the bilinear term. Transforming the integration to the current configuration:

$$\int_{V^0} \bigl(d\tau : \delta\mathbf{D} + \tau : d\delta\mathbf{D}\bigr)\,dV^0 = \int_V \Bigl(\frac{1}{J}d\tau : \delta\mathbf{D} + \sigma : d\delta\mathbf{D}\Bigr)\,dV$$

---

## 2. Material Tangent

### 2.1 Rate-form materials (e.g. elasto-plasticity)

For rate-form materials the material Jacobian $\mathbf{C}$ is defined in the corotational frame. In the global frame, the Jaumann rate of Kirchhoff stress gives (where $\mathbf{W} = \frac{1}{2}(\mathbf{L}-\mathbf{L}^T)$ is the spin tensor, skew-symmetric part of $\mathbf{L}$):

$$d\tau = J\bigl(\mathbf{C} : d\mathbf{D} + d\mathbf{W}\cdot\sigma - \sigma\cdot d\mathbf{W}\bigr)
$$

where

$$\mathbf{C} = \frac{1}{J}\frac{\partial\Delta(J\sigma)}{\partial\Delta\varepsilon}$$

and $\Delta(\cdot)$ denotes the increment over the current load step (accumulated Newton iterations).

### 2.2 Expanding the Jacobian

The Jacobian can be split as:

$$\mathbf{C} = \frac{1}{J}\frac{\partial\Delta(J\sigma)}{\partial\Delta\varepsilon} = \frac{\partial\Delta\sigma}{\partial\Delta\varepsilon}  + \frac{1}{J}\frac{\partial\Delta J}{\partial\Delta\varepsilon}\,\sigma$$

Since $dJ = J(d\varepsilon_{11}+d\varepsilon_{22}+d\varepsilon_{33})$, the second (unsymmetric) term is:

$$ C^{\mathrm{unsym}}_{ijkl} = \delta_{kl}\,\sigma_{ij} $$

This term is skew and is **dropped by default** for solver efficiency. It may be retained when convergence difficulties arise (unsymmetric solver flag).

### 2.3 Unified framework

Both total-form materials (e.g. hyperelasticity) and rate-form materials yield the same assembled stiffness structure.  The only difference is that total-form materials compute stress directly from the deformation gradient rather than through incremental time-integration.

---

## 3. Derivation of the Element Stiffness Matrix

### 3.1 Starting point

Substituting the material tangent into the bilinear form:

$$\int_V \bigl(  \bigl(\mathbf{C}:d\mathbf{D}  + d\mathbf{W}\cdot\sigma - \sigma\cdot d\mathbf{W}\bigr):\delta\mathbf{D}  + \sigma : d\delta\mathbf{D}\Bigr)\,dV$$

We must simplify the spin-related and second-order terms.

### 3.2 Second variation of the strain-rate tensor

Using the identity $\sigma : d\delta\mathbf{W} = 0$ (since $\sigma$ is symmetric and $d\delta\mathbf{W}$ is skew):

$$\sigma : d\delta\mathbf{D}= \sigma : d\delta\mathbf{L}= \sigma : d(\delta\mathbf{F}\cdot\mathbf{F}^{-1})$$

Noting that $d\delta\mathbf{F} = \partial(d\delta\mathbf{x})/\partial\mathbf{X} = 0$(virtual displacements and increments are independent), and using the identity $d\mathbf{F}^{-1} = -\mathbf{F}^{-1}\cdot d\mathbf{F}\cdot\mathbf{F}^{-1}$:

$$ \sigma : d\delta\mathbf{D} = -\sigma : (\delta\mathbf{L}\cdot d\mathbf{L}) $$

### 3.3 Spin-stress commutator terms

The two spin-related terms $(d\mathbf{W}\cdot\sigma - \sigma\cdot d\mathbf{W}):\delta\mathbf{D}$ are simplified using index notation. Since $d\mathbf{W}$ is skew-symmetric:

$$dW_{ij}\sigma_{jk}\delta D_{ik} + \sigma_{ij}dW_{kj}\delta D_{ik}$$

Using the symmetry of $\sigma$ and $\mathbf{D}$, swapping indices $i\leftrightarrow k$ in the first term:

$$= \sigma_{kj}\delta D_{ki}dW_{ij} + \sigma_{ij}\delta D_{ik}dW_{kj}= 2\sigma_{ij}\delta D_{ik}dW_{kj}$$

In tensor form: $\sigma : (2\delta\mathbf{D}\cdot d\mathbf{W})$.

### 3.4 Combining all terms

$$\int_V \bigl(\delta\mathbf{D}:\mathbf{C}:d\mathbf{D}- \sigma:(\delta\mathbf{L}\cdot d\mathbf{L} - 2\delta\mathbf{D}\cdot d\mathbf{W})\bigr)\,dV$$

Substituting $d\mathbf{W} = \frac{1}{2}(d\mathbf{L}-d\mathbf{L}^T) = d\mathbf{L}-d\mathbf{D}$:

$$= \int_V \bigl(\delta\mathbf{D}:\mathbf{C}:d\mathbf{D} - \sigma:(2\delta\mathbf{D}\cdot d\mathbf{D}  + (\delta\mathbf{L}-2\delta\mathbf{D})\cdot d\mathbf{L})\bigr)\,dV$$

Noting that $\delta\mathbf{L}-2\delta\mathbf{D} = -\delta\mathbf{L}^T$:

$$\boxed{\int_V \left(  \delta\mathbf{D}:\mathbf{C}:d\mathbf{D}  - \sigma:\left(2\delta\mathbf{D}\cdot d\mathbf{D}    - \frac{\partial\delta\mathbf{v}}{\partial\mathbf{x}}^T      \cdot\frac{\partial d\mathbf{v}}{\partial\mathbf{x}}\right)\right) dV}$$

This is the **final form** of the element stiffness integral, consisting of:

1. **Material stiffness**: $\delta\mathbf{D}:\mathbf{C}:d\mathbf{D}$
2. **Geometric stiffness**: the remaining $\sigma$-dependent terms

---

## 4. Geometric Stiffness in Component Form

### 4.1 Expansion

Let node indices be $I, J$ and spatial component indices $i, j, k$. Let $g_{Ik} = \partial N_I/\partial x_k$ denote the spatial shape-function gradient.

The geometric-stiffness integrand in component form is:

$$-\frac{1}{2}\sigma_{ij}\Bigl(  \delta v_{Ii}\,\frac{\partial N_I}{\partial x_k}\,\frac{\partial N_J}{\partial x_j}\,dv_{Jk}+ \delta v_{Ik}\,\frac{\partial N_I}{\partial x_i}\,\frac{\partial N_J}{\partial x_k}\,dv_{Jj}+ \delta v_{Ii}\,\frac{\partial N_I}{\partial x_k}\,\frac{\partial N_J}{\partial x_k}\,dv_{Jj}- \delta v_{Ik}\,\frac{\partial N_I}{\partial x_i}\,\frac{\partial N_J}{\partial x_j}\,dv_{Jk}\Bigr)$$

### 4.2 Stiffness matrix entry

Reading off the coefficient of $\delta v_{I\alpha}\,dv_{J\beta}$ gives the geometric tangent contribution at DOF $(I,\alpha)$-$(J,\beta)$:

$$\boxed{K_{\mathrm{geo}}(I,\alpha;\,J,\beta) = -\frac{1}{2}\Bigl[  \sigma_{\alpha\beta}\,(\mathbf{g}_I\cdot\mathbf{g}_J)+ (\boldsymbol{\sigma}_{\alpha:}\cdot\mathbf{g}_J)\,g_{I\beta}+ (\mathbf{g}_I\cdot\boldsymbol{\sigma}_{:\beta})\,g_{J\alpha}- \delta_{\alpha\beta}\,(\mathbf{g}_I\cdot\boldsymbol{\sigma}\,\mathbf{g}_J)\Bigr]\,dV}$$

where:

- $\mathbf{g}_I = \nabla N_I$ (spatial gradient of shape function $I$)
- $\boldsymbol{\sigma}_{\alpha:}$ = row $\alpha$ of the stress tensor
- $\boldsymbol{\sigma}_{:\beta}$ = column $\beta$ of the stress tensor
- $\delta_{\alpha\beta}$ = Kronecker delta

---

## 5. Configurations

The element formulation uses **three distinct configurations** (where $\mathbf{X}$ are reference/undeformed coordinates, $\mathbf{x}$ are current/deformed coordinates, and $\mathbf{u}$ is the displacement): 

| Configuration | Coordinates | Used for |
|---|---|---|
| **Reference** | $\mathbf{X}$ | Small-strain Jacobian |
| **Current** | $\mathbf{X} + \mathbf{u}^{n+1}$ | Assembly $\mathbf{B}$, $\det J$, geometric stiffness |
| **Midpoint** | $\mathbf{X} + \frac{1}{2}(\mathbf{u}^n + \mathbf{u}^{n+1})$ | Constitutive strain increment, Hughes-Winget rotation |

### 5.2 Jacobian and spatial gradients

The Jacobian mapping and spatial shape-function gradients are:

$$\mathbf{J} = \frac{\partial N}{\partial\boldsymbol{\xi}}\,\mathbf{x},\qquad\frac{\partial N}{\partial\mathbf{x}} = \mathbf{J}^{-1}\frac{\partial N}{\partial\boldsymbol{\xi}},\qquad dV = \det(\mathbf{J})\,w $$

Under NLGEOM, the Jacobian is evaluated at the **current configuration**, so $\mathbf{x} = \mathbf{X}+\mathbf{u}^{n+1}$.

### 5.3 Strain-displacement matrix $\mathbf{B}$

For each node $a$ with spatial gradients $N_{a,x}, N_{a,y}, N_{a,z}$, the per-node block is:

$$ \mathbf{B}_a =
\begin{bmatrix}
N_{a,x} & 0 & 0 \\
0 & N_{a,y} & 0 \\
0 & 0 & N_{a,z} \\
N_{a,y} & N_{a,x} & 0 \\
0 & N_{a,z} & N_{a,y} \\
N_{a,z} & 0 & N_{a,x}
\end{bmatrix} $$

Voigt ordering: $[xx, yy, zz, xy, yz, xz]$ with **engineering shear** (i.e. $\gamma_{xy} = \partial u/\partial y + \partial v/\partial x$).

### 5.4 Strain increment

The strain increment is:

$$ \Delta\boldsymbol{\varepsilon} = \mathbf{B}\,\Delta\mathbf{u},\qquad \boldsymbol{\varepsilon}^{n+1} = \boldsymbol{\varepsilon}^n + \Delta\boldsymbol{\varepsilon} $$

Under NLGEOM, the $\mathbf{B}$ matrix used for the constitutive strain path is evaluated at the **midpoint configuration** (see Section 6).

---

## 6. Objective Stress Transport: Hughes-Winget Algorithm

### 6.1 Purpose

Before the constitutive update, the old stress must be rotated into the current configuration to remove rigid-body rotation contamination. The **Hughes-Winget** algorithm provides a second-order accurate, incrementally objective rotation.

### 6.2 Velocity gradient at midpoint

The velocity gradient is computed from the midpoint-configuration spatial gradients:

$$
L_{ij} = \sum_a \Delta u_{ai}\,\frac{\partial N_a}{\partial x_j}\bigg|_{\mathrm{mid}}
$$

where $\Delta u_a = \mathbf{u}_a^{n+1} - \mathbf{u}_a^n$ and the shape function gradients are evaluated at the midpoint configuration.

### 6.3 Spin tensor

$$
\mathbf{W} = \frac{1}{2}(\mathbf{L} - \mathbf{L}^T)
$$

### 6.4 Hughes-Winget rotation matrix

The incremental rotation is obtained via a first-order Pade approximation of the exponential map:

$$
\boxed{
\Delta\mathbf{R} = \left(\mathbf{I} - \frac{1}{2}\mathbf{W}\right)^{-1}
                   \left(\mathbf{I} + \frac{1}{2}\mathbf{W}\right)
}
$$

This is a second-order accurate, unconditionally proper-orthogonal approximation to $\exp(\mathbf{W})$.

### 6.5 Stress rotation

The old Cauchy stress is rotated into the current configuration:

$$
\boldsymbol{\sigma}_{\mathrm{old}}^{\mathrm{rot}}
= \Delta\mathbf{R}\,\boldsymbol{\sigma}_{\mathrm{old}}\,\Delta\mathbf{R}^T
$$

In Voigt form this is implemented as $\boldsymbol{\sigma}' = \mathbf{M}(\Delta\mathbf{R})\,\boldsymbol{\sigma}$, where $\mathbf{M}$ is the $6\times6$ stress-Voigt rotation matrix.

---

## 7. Integration Point Update Algorithm

### 7.1 NLGEOM path (per integration point)

1. Compute spatial gradients and $\det\mathbf{J}$ at the **current** configuration $(\mathbf{x} = \mathbf{X} + \mathbf{u}^{n+1})$; build $\mathbf{B}_{\mathrm{cur}}$.
2. Switch to the **midpoint** configuration $(\mathbf{x} = \mathbf{X} + \tfrac{1}{2}(\mathbf{u}^n+\mathbf{u}^{n+1}))$; compute midpoint spatial gradients.
3. Form the strain increment $\Delta\boldsymbol{\varepsilon} = \mathbf{B}_{\mathrm{mid}}\,\Delta\mathbf{u}$; update $\boldsymbol{\varepsilon}^{n+1} = \boldsymbol{\varepsilon}^n + \Delta\boldsymbol{\varepsilon}$.
4. Compute the Hughes-Winget spin $\mathbf{W}$ from midpoint gradients; form $\Delta\mathbf{R}$.
5. Rotate old stress: $\boldsymbol{\sigma}_{\mathrm{old}}^{\mathrm{rot}} = \Delta\mathbf{R}\,\boldsymbol{\sigma}_{\mathrm{old}}\,\Delta\mathbf{R}^T$.
6. Constitutive update: $(\boldsymbol{\sigma}_{\mathrm{old}}^{\mathrm{rot}},\,\boldsymbol{\varepsilon}^n,\,\boldsymbol{\varepsilon}^{n+1}) \;\to\; (\boldsymbol{\sigma}^{n+1},\,\mathbf{C})$.
7. Assemble material stiffness: $\mathbf{K} \mathrel{+}= \mathbf{B}_{\mathrm{cur}}^T\,\mathbf{C}\,\mathbf{B}_{\mathrm{cur}}\,dV_{\mathrm{cur}}$.
8. Assemble geometric stiffness: $\mathbf{K} \mathrel{+}= \mathbf{K}_{\mathrm{geo}}(\boldsymbol{\sigma}^{n+1},\,\nabla N_{\mathrm{cur}})\,dV_{\mathrm{cur}}$.
9. Assemble internal force: $\mathbf{R} \mathrel{-}= \mathbf{B}_{\mathrm{cur}}^T\,\boldsymbol{\sigma}^{n+1}\,dV_{\mathrm{cur}}$.

### 7.2 Small-strain path

When NLGEOM is disabled:

1. Compute spatial gradients and $\det\mathbf{J}$ at the **reference** configuration $(\mathbf{x} = \mathbf{X})$; build $\mathbf{B}$.
2. Form strain increment $\Delta\boldsymbol{\varepsilon} = \mathbf{B}\,\Delta\mathbf{u}$; update $\boldsymbol{\varepsilon}^{n+1}$.
3. Constitutive update: $(\boldsymbol{\sigma}^n,\,\boldsymbol{\varepsilon}^n,\,\boldsymbol{\varepsilon}^{n+1}) \;\to\; (\boldsymbol{\sigma}^{n+1},\,\mathbf{C})$.
4. $\mathbf{K} \mathrel{+}= \mathbf{B}^T\,\mathbf{C}\,\mathbf{B}\,dV$; $\quad\mathbf{R} \mathrel{-}= \mathbf{B}^T\,\boldsymbol{\sigma}^{n+1}\,dV$.

No stress rotation or geometric stiffness is needed.

---

## 8. Logarithmic (Hencky) Strain for Field Output

### 8.1 Incremental deformation gradient

$$
\mathbf{J}_{\mathrm{old}} = \frac{\partial N}{\partial\boldsymbol{\xi}}\,(\mathbf{X}+\mathbf{u}^n),\qquad
\mathbf{J}_{\mathrm{new}} = \frac{\partial N}{\partial\boldsymbol{\xi}}\,(\mathbf{X}+\mathbf{u}^{n+1})
$$

$$
\mathbf{F} = \mathbf{J}_{\mathrm{old}}^{-1}\,\mathbf{J}_{\mathrm{new}}
$$

### 8.2 Right Cauchy-Green tensor and eigendecomposition

$$
\mathbf{C}_R = \mathbf{F}^T\mathbf{F},\qquad
\mathbf{C}_R \leftarrow \frac{1}{2}(\mathbf{C}_R + \mathbf{C}_R^T)
\quad(\text{symmetrize for numerical stability})
$$

Eigendecomposition: $\mathbf{C}_R\,\mathbf{p}_a = \lambda_a\,\mathbf{p}_a$, with clamping $\lambda_a \leftarrow \max(\lambda_a, 10^{-20})$.

### 8.3 Logarithmic strain

$$
e_a = \frac{1}{2}\ln\lambda_a,\qquad
\mathbf{E}_{\log} = \sum_a e_a\,\mathbf{p}_a\otimes\mathbf{p}_a
$$

Converted to engineering-shear Voigt form:

$$
LE_{xx} = E_{\log,xx},\quad
LE_{xy} = 2E_{\log,xy},\quad\text{etc.}
$$

---

## 9. Summary: Stiffness Matrix Composition

Under NLGEOM, the total element stiffness matrix is:

$$
\boxed{
\mathbf{K} = \underbrace{\mathbf{B}^T\,\mathbf{C}\,\mathbf{B}\,dV}_{\text{material stiffness}}
           + \underbrace{\mathbf{K}_{\mathrm{geo}}(\sigma, \nabla N)\,dV}_{\text{geometric stiffness}}
}
$$

where:

- $\mathbf{B}$ is evaluated at the **current** configuration
- $\mathbf{C}$ is the material tangent from the constitutive update (which
  uses midpoint-configuration strain increments and Hughes-Winget rotated stress)
- $\mathbf{K}_{\mathrm{geo}}$ is the geometric stiffness from Section 4

The internal force (residual) is:

$$
\mathbf{R}_{\mathrm{int}} = -\mathbf{B}^T\,\boldsymbol{\sigma}\,dV
$$
