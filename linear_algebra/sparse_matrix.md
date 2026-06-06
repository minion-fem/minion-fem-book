---
title: "Sparse matrix"
---

# Sparse matrix

Each Newton iteration assembles a global sparse linear system $\mathbf{K}\,\Delta\mathbf{u} = -\mathbf{R}$ and solves it. This section covers the structure of the FEM stiffness matrix, how it is stored and assembled, and how to solve it efficiently.

---

## Sparsity: where it comes from

The stiffness matrix $\mathbf{K}$ is assembled from element contributions:

$$K_{ij} = \sum_e K^e_{ij}$$

Element $e$ contributes only to DOF pairs $(i,j)$ that both belong to that element. Two DOFs are coupled if and only if they share at least one element — i.e., the two nodes are **topologically adjacent** in the mesh. All other entries are exactly zero.

For a 3D hexahedral mesh, each node connects to at most 26 neighbours (including itself). With 3 DOFs per node, each row of $\mathbf{K}$ has at most $27 \times 3 = 81$ nonzero entries — regardless of how large the mesh is. The total number of nonzeros therefore grows as $O(n)$ with the number of DOFs $n$.

### What the stiffness matrix looks like

For a structured 3D mesh with nodes numbered sequentially, the nonzero pattern has a **banded structure** around the diagonal, with bandwidth proportional to the number of nodes in one cross-section layer. For unstructured meshes, the band is irregular but the $O(n)$ sparsity still holds.

```
●  ·  ·  ●  ●  ·  ·  ·  ·
·  ●  ·  ●  ●  ●  ·  ·  ·
·  ·  ●  ·  ●  ●  ·  ·  ·
●  ●  ·  ●  ·  ·  ●  ●  ·
●  ●  ●  ·  ●  ·  ●  ●  ●
·  ●  ●  ·  ·  ●  ·  ●  ●
·  ·  ·  ●  ●  ·  ●  ·  ·
·  ·  ·  ●  ●  ●  ·  ●  ·
·  ·  ·  ·  ●  ●  ·  ·  ●
```

*Schematic nonzero pattern for a small 3D elasticity problem. Each `●` represents a 3×3 block coupling two nodes.*

For a 3D elasticity problem with $n_{\text{node}}$ nodes, typical statistics are:

| $n_{\text{node}}$ | DOFs ($\times 3$) | Nonzeros (approx.) | Dense storage |
|---|---|---|---|
| 10k | 30k | ~2M | 7.2 GB → 16 MB |
| 100k | 300k | ~20M | 720 GB → 160 MB |
| 1M | 3M | ~200M | 72 TB → 1.6 GB |

Sparse storage reduces memory by 3–4 orders of magnitude at large scale.

### Sparsity pattern construction

The nonzero structure is determined by the mesh topology alone and is fixed throughout the analysis. It is built once at setup:

1. For each element, collect the DOF index list $\mathcal{I}_e$.
2. For every pair $(i, j) \in \mathcal{I}_e \times \mathcal{I}_e$, mark entry $(i,j)$ as a nonzero.
3. Allocate storage; the pattern does not change across Newton iterations or load steps.

---

## Sparse matrix storage formats

Several formats are in common use. The right choice depends on the operation to be performed.

### COO — Coordinate format

Three arrays store each nonzero as a `(row, col, value)` triple:

```
rows   = [0, 0, 1, 2, 2, ...]
cols   = [0, 3, 1, 2, 4, ...]
values = [k00, k03, k11, k22, k24, ...]
```

**Pros:** easy to build incrementally (just append triples).  
**Cons:** no efficient row or column access; duplicate entries must be summed.  
**Used for:** intermediate assembly; most sparse libraries accept COO as input and convert internally.

### CSR — Compressed Sparse Row

Three arrays:
- `values[nnz]` — nonzero values in row-major order
- `col_indices[nnz]` — column index of each value
- `row_pointers[n+1]` — `row_pointers[i]` is the index in `values` where row $i$ starts

```
K = [10  0   0   3 ]
    [ 0  8   2   0 ]
    [ 0  2   5   0 ]
    [ 3  0   0   7 ]

values      = [10, 3,  8, 2,  2, 5,  3, 7]
col_indices = [ 0, 3,  1, 2,  1, 2,  0, 3]
row_pointers= [ 0,     2,     4,     6,    8]
```

Row $i$ spans `values[row_pointers[i] : row_pointers[i+1]]`.

**Pros:** fast row access and matrix-vector products ($\mathbf{K}\mathbf{x}$); standard format for iterative solvers.  
**Cons:** inserting new nonzeros requires rebuilding; column access is slow.  
**Used for:** global stiffness matrix storage, iterative solver input, PETSc.

### CSC — Compressed Sparse Column

Mirror of CSR with column-major storage. Fast column access and sparse triangular solves.  
**Used for:** direct solver factorisation (CHOLMOD, SuperLU internally use CSC or hybrid formats).

### BSR — Block Sparse Row

Like CSR but stores dense $b \times b$ blocks instead of scalars. For 3D elasticity, $b=3$ (one block per node pair). Reduces index overhead by $b^2$ and enables BLAS-3 operations on the dense blocks.  
**Used for:** PETSc's `BAIJ` format; significant performance advantage on modern CPUs with SIMD.

### Format summary

| Format | Build | Row access | Col access | MatVec | Typical use |
|---|---|---|---|---|---|
| COO | Easy | Slow | Slow | Slow | Assembly input |
| CSR | Hard | Fast | Slow | Fast | Iterative solvers |
| CSC | Hard | Slow | Fast | Fast | Direct solvers |
| BSR | Hard | Fast | Slow | Very fast | Block iterative |

---

## Assembly

Each element computes a local stiffness $\mathbf{K}^e \in \mathbb{R}^{n_e \times n_e}$ and scatters it into the global CSR matrix:

$$K[d_i, d_j] \mathrel{+}= K^e_{ij} \quad \text{for all } (i,j) \in \mathcal{I}_e \times \mathcal{I}_e$$

Because the sparsity pattern is pre-built, each scatter-add is a direct array index into `values[]` — no hash lookup. The column-index lookup `(row, col) → values_index` is cached per element DOF pair so that repeated Newton iterations avoid re-searching.

### Constraint application

After assembly, Dirichlet boundary conditions are applied by the **row-zeroing** method: for each fixed DOF $d$, zero the entire row and column, set the diagonal to 1, and set $R[d]$ to the prescribed increment. Multi-point constraints (MPCs) use either penalty augmentation or constraint elimination via the `AxConstraintManager`.

**Penalty method warning:** MPCs enforced by adding a large penalty $\alpha$ to the diagonal inflate the condition number of $\mathbf{K}$ by $\sim \alpha / K_{ii}$. This can degrade iterative solver convergence dramatically — contact problems with penalty-enforced tie constraints have been observed to require 10× more CG iterations than unconstrained problems of the same size.

---

## Solvers

- [Direct solvers](direct.md) — LU/Cholesky factorisation, fill-in, reordering, CHOLMOD vs SuperLU
- [Iterative solvers](iterative.md) — Krylov methods, preconditioning, CG and GMRES
- [Algebraic multigrid](amg.md) — GAMG, near null space, convergence data
