---
title: "Iterative Solvers"
---

# Iterative Solvers

An iterative solver starts from an initial guess $\mathbf{x}_0$ and produces a sequence of improving approximations $\mathbf{x}_1, \mathbf{x}_2, \ldots$ converging to the solution of $\mathbf{K}\mathbf{x} = \mathbf{b}$. Each iteration costs $O(n_{\text{nnz}})$ — one matrix-vector product — with no fill-in. This makes iterative solvers the only practical option for large 3D problems where direct factorisation is memory-prohibitive.

---

## Krylov subspace methods

The dominant family of iterative methods for sparse linear systems is **Krylov subspace methods**. The idea: at iteration $k$, search for the best solution in the subspace spanned by the first $k$ Krylov vectors

$$\mathcal{K}_k(\mathbf{K}, \mathbf{r}_0) = \operatorname{span}\{\mathbf{r}_0,\, \mathbf{K}\mathbf{r}_0,\, \mathbf{K}^2\mathbf{r}_0,\, \ldots,\, \mathbf{K}^{k-1}\mathbf{r}_0\}$$

where $\mathbf{r}_0 = \mathbf{b} - \mathbf{K}\mathbf{x}_0$ is the initial residual. Theoretically, after $n$ iterations the exact solution is found; in practice convergence is declared much earlier once $\|\mathbf{r}_k\|/\|\mathbf{r}_0\| < \varepsilon$.

The Krylov sequence is orthogonalised as it is built (Arnoldi / Lanczos algorithms) to avoid numerical instability — raw Krylov vectors tend toward the dominant eigenvector and become nearly linearly dependent.

### Conjugate Gradient (CG)

For **symmetric positive definite** (SPD) systems, CG minimises the $\mathbf{K}$-norm of the error at each step:

$$\mathbf{x}_{k+1} = \mathbf{x}_k + \alpha_k \mathbf{p}_k, \qquad \alpha_k = \frac{\mathbf{r}_k^T \mathbf{r}_k}{\mathbf{p}_k^T \mathbf{K} \mathbf{p}_k}$$

The search directions $\mathbf{p}_k$ are $\mathbf{K}$-conjugate ($\mathbf{p}_i^T \mathbf{K} \mathbf{p}_j = 0$ for $i \neq j$), which guarantees that each iteration makes non-redundant progress.

FEM stiffness matrices are SPD after proper constraint application, making CG the natural choice.

### GMRES

For **non-symmetric** or **indefinite** systems (contact with Lagrange multipliers, unsymmetric tangent from viscoplasticity), CG does not apply. GMRES minimises the residual over the full Krylov space at each step, at the cost of storing all previous search directions. A restart parameter (typically 30–50) limits memory use at the expense of some convergence efficiency.

---

## Convergence and condition number

The convergence rate of CG depends on the **condition number** $\kappa(\mathbf{K}) = \lambda_{\max} / \lambda_{\min}$:

$$\frac{\|\mathbf{e}_k\|_{\mathbf{K}}}{\|\mathbf{e}_0\|_{\mathbf{K}}} \leq 2\left(\frac{\sqrt{\kappa}-1}{\sqrt{\kappa}+1}\right)^k$$

The number of iterations to reach tolerance $\varepsilon$ scales as $O(\sqrt{\kappa} \log(1/\varepsilon))$.

For FEM stiffness matrices, the condition number depends on the mesh:

$$\kappa(\mathbf{K}) = O(h^{-2})$$

where $h$ is the element size. Halving the mesh size quadruples the condition number and doubles the iteration count. A 3D mesh with $n$ nodes has $h \sim n^{-1/3}$, so $\kappa \sim n^{2/3}$ — unpreconditioned CG on a 1M-DOF problem requires thousands of iterations.

---

## Preconditioning

A **preconditioner** $\mathbf{M} \approx \mathbf{K}^{-1}$ transforms the system to $\mathbf{M}\mathbf{K}\mathbf{x} = \mathbf{M}\mathbf{b}$, reducing the condition number. The ideal preconditioner is $\mathbf{K}^{-1}$ itself (one iteration), but computing it exactly is as expensive as a direct solve. The goal is to approximate $\mathbf{K}^{-1}$ cheaply enough that the total cost (preconditioner apply × iterations) beats a direct solve.

> Iterative methods only outperform direct solvers when paired with a good preconditioner. Otherwise, the direct solver wins on robustness.

### Diagonal (Jacobi) scaling

The cheapest preconditioner: $\mathbf{M} = \text{diag}(\mathbf{K})^{-1}$. Normalises each equation by its diagonal entry. Effective only when the diagonal dominates, which is not generally true for FEM systems. Used as the smoother in GAMG (see [AMG](amg.md)).

### Incomplete LU (ILU)

Perform LU factorisation but drop fill-in entries below a threshold. Produces an approximate factorisation $\tilde{\mathbf{L}}\tilde{\mathbf{U}} \approx \mathbf{K}$ stored in the same sparsity pattern as $\mathbf{K}$ (ILU(0)) or with limited additional fill (ILU(k)). ILU is a good general-purpose preconditioner but does not scale well to millions of DOFs.

### Algebraic Multigrid (AMG)

The most effective preconditioner for large FEM systems. AMG builds a hierarchy of coarser representations of the problem and uses them to damp errors at all length scales simultaneously, achieving near-$O(n)$ total solve cost. See [Algebraic Multigrid](amg.md) for a detailed treatment.
