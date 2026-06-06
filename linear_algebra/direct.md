---
title: "Direct Solvers"
---

# Direct Solvers

A direct solver computes an exact factorisation of the stiffness matrix $\mathbf{K} = \mathbf{L}\mathbf{U}$ (or $\mathbf{L}\mathbf{D}\mathbf{L}^T$ for symmetric systems) and solves the two triangular systems $\mathbf{L}\mathbf{y} = \mathbf{b}$, $\mathbf{U}\mathbf{x} = \mathbf{y}$. The result is numerically exact up to floating-point rounding, with no convergence parameter to tune.

---

## LU and Cholesky factorisation

For a general sparse matrix, Gaussian elimination produces an LU factorisation. The key complication for sparse matrices is **fill-in**: entries that are zero in $\mathbf{K}$ but become nonzero in $\mathbf{L}$ or $\mathbf{U}$ during elimination. In the worst case (a dense $n \times n$ matrix), fill-in turns $O(n)$ nonzeros into $O(n^2)$ nonzeros in the factors.

For FEM stiffness matrices:
- In **2D**, fill-in grows as $O(n \log n)$ to $O(n^{1.5})$ depending on mesh structure.
- In **3D**, fill-in grows as $O(n^{4/3})$ to $O(n^2)$.

This is the fundamental scaling limit of direct solvers: beyond ~1–5 million DOFs in 3D, memory and compute cost become prohibitive.

For **symmetric positive definite** (SPD) matrices — which FEM stiffness matrices always are (before constraint modification) — the **Cholesky factorisation** $\mathbf{K} = \mathbf{L}\mathbf{L}^T$ requires roughly half the operations and memory of full LU, and is numerically more stable.

---

## Fill-in and variable reordering

Fill-in depends on the order in which variables are eliminated. Reordering the rows and columns of $\mathbf{K}$ before factorisation can dramatically reduce fill-in without changing the solution.

Finding the optimal reordering (minimum fill-in) is NP-hard. Several practical heuristics are widely used:

**Approximate Minimum Degree (AMD)**: greedily eliminates the variable with the fewest connections at each step, estimated before actual elimination. Used by default in CHOLMOD and UMFPACK.

**Nested Dissection (Metis/SCOTCH)**: recursively bisects the mesh graph using a separator, numbers the separator nodes last. This limits how far fill-in can propagate across the separator. Provides near-optimal $O(n^{4/3})$ fill-in for 3D FEM problems and is preferred for large models. MUMPS and Pardiso both recommend Metis as the default reordering.

**Reverse Cuthill-McKee (RCM)**: reduces the matrix bandwidth by ordering nodes in a BFS traversal. Useful for banded solvers but generally worse than AMD/Metis for sparse direct solvers.

The impact of reordering on memory and time is significant:

| Reordering | Fill-in (relative) | Factorisation time |
|---|---|---|
| None (natural ordering) | 1× (baseline, often bad) | 1× |
| RCM | 0.5–0.8× | 0.5–0.8× |
| AMD | 0.3–0.6× | 0.3–0.6× |
| Nested Dissection (Metis) | 0.1–0.4× | 0.1–0.4× |

---

## Supernodal and multifrontal methods

Naive sparse LU processes one column at a time, which makes poor use of modern CPUs' BLAS-3 capabilities (matrix-matrix operations). Two approaches restructure the factorisation to exploit dense subproblems:

**Supernodal method** (SuperLU): groups columns with identical sparsity patterns into *supernodes* and processes each supernode as a dense BLAS-3 operation. Efficient on shared-memory multicore systems.

**Multifrontal method** (MUMPS, UMFPACK): decomposes the elimination tree into *frontal matrices* — dense subproblems assembled from element contributions and partial factorisation results from child nodes. Each frontal matrix is processed as a dense block, enabling BLAS-3 and natural parallelism across the elimination tree. The main commercial FEM codes (Abaqus, Nastran, Simscale) all use multifrontal solvers.

---

## Solver backends

### CHOLMOD (via scikit-sparse)

Cholesky factorisation for SPD matrices. Uses AMD or nested-dissection reordering, supernodal factorisation, and BLAS-3 kernels. Approximately **2× faster** than SuperLU on typical FEM stiffness matrices.

```bash
brew install suite-sparse && pip install scikit-sparse
```

### SuperLU (via scipy)

Supernodal LU factorisation. No positive-definiteness requirement — handles non-symmetric tangent matrices (e.g., viscoplasticity with the unsymmetric $C^{\text{unsym}}$ term). Available as the SciPy default without additional installation.

### MUMPS

Parallel multifrontal solver with MPI support. Strongly recommends Metis reordering. Handles both symmetric and unsymmetric systems. 

## Benchmark test
Benchmark results on standard FEM matrices (see table below):

| Matrix | DOFs | Sym | MUMPS time | SuperLU time | Pardiso time |
|---|---|---|---|---|---|
| af_shell1 | 505k | Y | — | — | fastest |
| thermal2 | 1.2M | Y | — | — | fastest |
| Cube_Coup_dt0 | 2.2M | Y | failed | failed | failed (>memory) |

All three solvers fail on 2M+ DOF problems on a standard workstation (Intel i5, 16 GB RAM), which is why iterative solvers are needed for large 3D problems.

---

## When to use direct solvers

Direct solvers are the right choice when:
- Limited CPU memory
- Multiple right-hand sides need to be solved with the same $\mathbf{K}$ 
- Robustness matters more than scalability (direct solvers always converge)

For larger problems, see [Iterative Solvers](iterative.md).
