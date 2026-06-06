---
title: "Algebraic Multigrid"
---

# Algebraic Multigrid (AMG)

Algebraic Multigrid (AMG) is the state-of-the-art preconditioner for large sparse linear systems arising from elliptic PDEs. For 3D structural FEM problems it achieves near-$O(n)$ total solve complexity — qualitatively better than any direct solver or simple preconditioner.

---

## Why simple iterative methods stall

Consider solving $\mathbf{K}\mathbf{x} = \mathbf{b}$ with Gauss-Seidel (GS) as the solver. GS is a local smoother: each step updates one unknown using its immediate neighbours. It damps **high-frequency** error components quickly (errors that oscillate between adjacent nodes), but barely touches **low-frequency** (smooth) error components that vary slowly across the mesh.

This can be understood through the **Rayleigh quotient**: for an SPD matrix $\mathbf{K}$, the quotient

$$R(\mathbf{K}, \mathbf{e}) = \frac{\mathbf{e}^T \mathbf{K} \mathbf{e}}{\mathbf{e}^T \mathbf{e}}$$

lies between $\lambda_{\min}$ and $\lambda_{\max}$. High-frequency error components have large Rayleigh quotient (associated with large eigenvalues); GS reduces them efficiently. Low-frequency components have small Rayleigh quotient — they look like rigid-body modes to the smoother and are barely reduced.

After several GS sweeps the remaining error is smooth. A smoother applied to a smooth error on a fine grid is wasteful: the same smooth function, when sampled on a coarser grid, is just as smooth — and the coarse-grid problem is much smaller.

---

## The multigrid idea

Multigrid exploits this by alternating between fine and coarse grids:

```
Fine grid:   smooth the error with GS (kills high-frequency error)
               ↓  restrict residual to coarse grid
Coarse grid: solve (or recurse) to get a correction
               ↑  interpolate correction back to fine grid
Fine grid:   smooth again
```

This **V-cycle** reduces both high- and low-frequency error in $O(n)$ work. The theoretical result: for model elliptic problems, the number of V-cycles to convergence is independent of mesh size $h$ — so halving $h$ (8× more unknowns in 3D) costs only 8× more work, not $8^{4/3}$ as for a direct solver.

---

## Classical AMG vs GAMG

**Classical AMG** constructs the coarse levels purely from the matrix entries: strong connections (large off-diagonal entries relative to the diagonal) determine which unknowns are coarsened and how the interpolation operators are built. No geometric information is needed.

**GAMG** (Geometric-Algebraic Multigrid, the PETSc implementation) uses the same algebraic coarsening machinery but can additionally use **node coordinates** and **near-null space vectors** to guide the construction of better coarse spaces.

### Why AMG struggles with 3D elasticity

For scalar elliptic problems ($-\nabla^2 u = f$), classical AMG works well. For **vector elasticity** with $d$ DOFs per node, the kernel of the stiffness operator is the space of rigid-body motions: $d$ translations and $d(d-1)/2$ rotations ($6$ modes in 3D).

Classical AMG based on scalar matrix entries does not know that $u_x$, $u_y$, $u_z$ at the same node are physically coupled, nor does it recognise rigid-body modes. This leads to:

- Incorrect coarsening: physically coupled DOFs may be coarsened separately
- Missed near-kernel modes: rigid-body-like error survives the coarse-grid correction

The result is slow or non-convergent AMG.

### Near null space and block structure

The fix: explicitly provide the rigid-body modes as the **near null space** of $\mathbf{K}$. For a 3D elasticity problem with nodes at positions $\mathbf{x}_i = (x_i, y_i, z_i)$, the 6 near-null space vectors are the 3 translations and 3 infinitesimal rotations:

$$\mathbf{t}_x = [1,0,0,\; 1,0,0,\; \ldots], \quad \mathbf{t}_y = [0,1,0,\; \ldots], \quad \mathbf{t}_z = [0,0,1,\; \ldots]$$

$$\mathbf{r}_{xy} = [-y_0,x_0,0,\; -y_1,x_1,0,\; \ldots], \quad \text{etc.}$$

GAMG uses these to build the interpolation operator $\mathbf{P}$ such that the coarse space exactly represents rigid-body motions. Equivalently, specifying `block_size=3` (one block per node, 3 DOFs) and providing node coordinates lets GAMG compute the near-null space internally.

### DOF ordering and block structure

GAMG's block-structured coarsening requires node-interleaved DOF ordering: $[u_{x,0}, u_{y,0}, u_{z,0},\; u_{x,1}, u_{y,1}, u_{z,1},\; \ldots]$. Most FEM codes store DOFs field-separated: $[u_{x,0},\ldots,u_{x,n},\; u_{y,0},\ldots]$.

The solver must therefore permute the matrix and right-hand side before calling PETSc and un-permute the solution afterward.

---

## Matrix structure and AMG performance

The sparsity pattern of the global stiffness matrix directly determines how well AMG can coarsen it. The two figures below compare a pure-elasticity matrix against a contact problem matrix from the same multi-body assembly (7.3M DOFs).

::::::{grid} 2

:::::{grid-item}
:::{figure} images/elastic_matrix.png
:width: 100%

**Pure elasticity** — four diagonal blocks, one per component. Each block has a clean, locally banded structure. AMG sees strong connections only between neighbouring nodes, coarsening proceeds level by level, and the near-null space is well-represented at every level.
:::
:::::

:::::{grid-item}
:::{figure} images/contact_matrix.png
:width: 100%

**With tie and contact constraints** — the same four diagonal blocks, plus dense off-diagonal blocks (rows ~4M–5M) introduced by the penalty-MPC coupling. These long-range connections violate the locality assumption of AMG coarsening and corrupt the near-null space detection.
:::
:::::

::::::

The element-level block structure is visible when zooming into a 5000×5000 sub-block of the interior DOFs:

:::{figure} images/sparsity_zoom.png
:width: 60%
:align: center

Zoomed sparsity pattern (rows/columns 100k–105k). Each small dense block corresponds to one hexahedral element's 24×24 stiffness contribution. The repeating pattern reflects the structured mesh in this region.
:::

At global scale, each dot in the coarse view below represents a 1000×1000 block:

:::{figure} images/contact_sparsity_coarse.png
:width: 60%
:align: center

Full 7.3M×7.3M matrix, 1 dot = 1000×1000 block. The main diagonal band is the elastic stiffness; the scattered off-diagonal dots are the contact/tie coupling terms. AMG's strong-connection graph picks up these long-range connections and misinterprets them as fine-scale features.
:::

---

## V-cycle structure in practice

A typical GAMG setup for a 3D FEM problem builds 4–5 levels. At each level, the smoother is **Chebyshev-accelerated Jacobi** (2 sweeps pre and post). The coarsest level is solved exactly with LU.

Example for an 8.5M-DOF elasticity problem:

| Level | DOFs | Nonzeros |
|---|---|---|
| Fine (L4) | 7,267,200 | 596,079,000 |
| L3 | 131,454 | 26,336,862 |
| L2 | 10,623 | 6,502,221 |
| L1 | 237 | 44,739 |
| Coarse (L0) | 15 | 225 |

Grid complexity (ratio of total DOFs across all levels to fine-grid DOFs): 1.02 — the coarse grids add only 2% overhead in DOF count.

---

## Observed convergence data

The following results use PETSc CG + GAMG with `block_size=3` and node coordinates:

| Problem | DOFs | CG iterations | Solve time |
|---|---|---|---|
| Perforated plate, simple BCs | 1,800,000 | 230 | 529 s (serial) |
| Perforated plate, simple BCs | 8,500,000 | 162 | 2798 s (serial) |
| Multi-body assembly, no contact | 7,267,200 | ~170 | 881 s (16 MPI) |
| Multi-body assembly, with contact (penalty MPC) | 7,267,200 | 1947 | 12419 s (serial) |
| Same, MPI 64 processes | 7,267,200 | 1934 | 575 s |

Key observations:

1. **Iteration count decreases as the mesh is refined** (162 vs 230) — the hallmark of an optimal preconditioner. A perfect AMG would give mesh-independent iteration counts.
2. **Contact/penalty MPC causes a 10× increase in iterations** (1947 vs ~170). The off-diagonal penalty blocks visible in the sparsity pattern above corrupt GAMG's coarsening and near-null space construction.
3. **MPI scaling is reasonable**: 64 processes achieve ~22× speedup on the contact problem.

---

## Practical guidelines

**When GAMG works well:**
- Pure elasticity, simple boundary conditions
- `block_size` = spatial dimension (2 or 3), node coordinates provided
- Near null space provided explicitly or derivable from block structure

**When GAMG struggles:**
- Penalty-method MPCs — the long-range off-diagonal entries corrupt near-null space detection (see sparsity figures above)
- Mixed elements (beam + continuum): block size is inconsistent across the mesh
- Highly heterogeneous materials: large modulus contrasts create strong connections that mislead the coarsening algorithm

For penalty-constrained contact problems, consider switching to constraint elimination or a Schur complement approach — see [Saddle Point Systems](saddle_point.md).
