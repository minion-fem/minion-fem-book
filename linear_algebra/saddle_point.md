---
title: "Saddle Point Systems"
---

# Saddle Point Systems

Constrained structural problems — tied surfaces, multi-point constraints enforced by Lagrange multipliers, contact with Lagrange multipliers — lead to a linear system of the form:

$$\begin{bmatrix} A & B^T \\ B & C \end{bmatrix} \begin{bmatrix} x \\ \lambda \end{bmatrix} = \begin{bmatrix} f \\ g \end{bmatrix}$$

where $A$ is the (symmetric positive semi-definite) stiffness matrix, $B$ encodes the constraint equations, and $\lambda$ are the Lagrange multipliers. In the typical FEM case $C = 0$.

This is a **saddle point system**: the solution $(x_*, \lambda_*)$ is a saddle point of the Lagrangian — a minimum in $x$ and a maximum in $\lambda$. The matrix is indefinite (eigenvalues of both signs), which breaks the assumptions of CG and standard AMG.

:::{figure} images/saddle_point.png
:width: 70%
:align: center

The function $z = x^2 - y^2$ has a saddle point at the origin: it is a minimum along $x$ and a maximum along $y$. The Lagrangian $\mathcal{L}(x,\lambda) = J(x) + \lambda^T(Bx-g)$ has the same geometry — minimised over $x$, maximised over $\lambda$.
:::

---

## Where saddle point systems arise

**Lagrange multiplier constraints** arise whenever equality constraints $Bx = g$ are enforced exactly (as opposed to the penalty method). The constrained minimisation problem

$$\min_x \frac{1}{2}x^T A x - f^T x \quad \text{subject to } Bx = g$$

has first-order conditions

$$Ax + B^T\lambda = f, \qquad Bx = g$$

which assemble directly into the saddle point system above. Examples:

- **Tied surfaces**: displacement continuity across a non-conforming mesh interface
- **Rigid body constraints**: all DOFs of a rigid region move together
- **Incompressibility**: $\text{div}\,u = 0$ enforced via pressure as a Lagrange multiplier (mixed formulation)

---

## Properties of the saddle point matrix

The block matrix $\mathcal{K} = \begin{bmatrix} A & B^T \\ B & 0 \end{bmatrix}$ has the following properties:

- **Indefinite**: eigenvalues of both signs. Standard SPD solvers (Cholesky, CG) cannot be applied directly.
- **Non-singular** if: $A$ is positive definite on $\ker(B)$, and $B$ has full row rank (inf-sup condition).
- **Symmetric**: MINRES (not GMRES) is the appropriate Krylov solver.
- **Condition number large**: scales poorly without a good preconditioner.

---

## Three key techniques

### 1. Augmentation

Add a term $\gamma B^T B$ to the $(1,1)$ block:

$$\begin{bmatrix} A + \gamma B^T B & B^T \\ B & 0 \end{bmatrix} \begin{bmatrix} x \\ \lambda \end{bmatrix} = \begin{bmatrix} f + \gamma B^T g \\ g \end{bmatrix}$$

The modified $(1,1)$ block $A + \gamma B^T B$ is now **symmetric positive definite** (not just semi-definite), even if $A$ is singular. This enables construction of the Schur complement and simplifies preconditioning. The solution is unchanged.

### 2. Regularisation

Add a positive definite matrix $C$ to the $(2,2)$ block:

$$\begin{bmatrix} A & B^T \\ B & -C \end{bmatrix} \begin{bmatrix} x \\ \lambda \end{bmatrix} = \begin{bmatrix} f \\ g \end{bmatrix}$$

Regularisation improves the condition number and numerical stability. It corresponds to a penalised formulation where the constraint is satisfied approximately (with error of order $C$).

### 3. Schur complement

Block elimination of $x$ from the system yields:

$$\begin{bmatrix} A & B^T \\ 0 & -BA^{-1}B^T \end{bmatrix} \begin{bmatrix} x \\ \lambda \end{bmatrix} = \begin{bmatrix} f \\ g - BA^{-1}f \end{bmatrix}$$

The matrix $S = -BA^{-1}B^T$ is the **Schur complement**. Solving the system then proceeds in two steps:

1. Solve $S\lambda = g - BA^{-1}f$ for $\lambda$ (the reduced system on the constraint DOFs)
2. Solve $Ax = f - B^T\lambda$ for $x$

The Schur complement system has dimension equal to the number of constraints (typically much smaller than the full system). This is the basis of FETI and domain decomposition methods.

If $A$ is singular (unconstrained floating body), the Schur complement does not exist. Augmentation ($A \to A + \gamma B^T B$) regularises $A$ to make it invertible before forming $S$.

---

## Augmented Lagrangian iteration

The **augmented Lagrangian method** avoids solving the full coupled saddle point system by iterating between two cheap sub-problems:

$$\mathcal{L}_A(x, \lambda, \gamma) = J(x) + \lambda^T(Bx - g) + \frac{\gamma}{2}\|Bx - g\|^2$$

**Uzawa iteration:**

1. **$x$-update**: solve $(A + \gamma B^T B)\,x_{k+1} = f - B^T \lambda_k + \gamma B^T g$
2. **$\lambda$-update**: $\lambda_{k+1} = \lambda_k + \gamma(B x_{k+1} - g)$

The $x$-update is a standard SPD system (AMG-friendly). The $\lambda$-update is a cheap explicit step with no linear solve.

| Method | Condition number | Constraint satisfaction |
|---|---|---|
| Penalty ($\gamma \to \infty$) | $O(\gamma)$ — explodes | Approximate, error $O(1/\gamma)$ |
| Pure Lagrange multiplier | Indefinite | Exact |
| Augmented Lagrangian | Finite, controlled by $\gamma$ | Exact at convergence |

---

## Visualising the saddle point

The saddle point geometry of the Lagrangian can be reproduced with:

```python
import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(-2, 2, 100)
y = np.linspace(-2, 2, 100)
X, Y = np.meshgrid(x, y)
Z = X**2 - Y**2          # saddle: minimum in x, maximum in y

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='none', alpha=0.8)
fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10, label='Value of z')

ax.scatter([0], [0], [0], color='red', s=100, label='Saddle Point (0,0,0)')
ax.set_title(r'Saddle Point: $z = x^2 - y^2$')
ax.set_xlabel('x (Minimization direction)')
ax.set_ylabel('y (Maximization direction)')
ax.set_zlabel('z')
plt.legend()
plt.tight_layout()
plt.savefig('saddle_point.png', dpi=150)
```
