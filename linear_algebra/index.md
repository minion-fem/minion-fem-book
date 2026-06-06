---
title: "Linear Algebra"
---

# Linear Algebra

Each Newton iteration reduces to solving a sparse linear system $\mathbf{K}\,\Delta\mathbf{u} = -\mathbf{R}$. The stiffness matrix is large, sparse, and symmetric positive definite — its structure and the choice of solver have a decisive impact on performance at scale. This section covers sparse matrix formats, direct factorisation, Krylov iterative methods, and algebraic multigrid preconditioning.

- [Sparse Matrix](sparse_matrix.md) — sparsity structure, storage formats (COO/CSR/CSC/BSR), assembly
- [Direct Solvers](direct.md) — LU/Cholesky factorisation, fill-in, reordering, multifrontal methods
- [Iterative Solvers](iterative.md) — Krylov methods, CG and GMRES, preconditioning, condition number
- [Algebraic Multigrid](amg.md) — GAMG, near null space, convergence data on million-DOF problems
- [Saddle Point Systems](saddle_point.md) — Lagrange multiplier constraints, augmentation, Schur complement, augmented Lagrangian iteration
