---
title: "Nonlinear Solver"
---

# Nonlinear Solver

Minion uses an implicit incremental-iterative Newton–Raphson scheme, closely following the Abaqus formulation. The same framework handles both linear elastic problems (one Newton iteration per step) and strongly nonlinear problems (large deformation, plasticity).

---

## Incremental-iterative structure

The total load is divided into $N$ increments. Within each increment $n$, Newton iterations drive the residual to zero:

```
for each load increment n = 1 … N:
    guess: u = u_n  (or extrapolate from previous increment)
    for each Newton iteration k = 1 … max_iter:
        assemble K(u^k) and R(u^k)
        solve:  K Δu = −R
        update: u^{k+1} = u^k + Δu
        if converged: break
    u_{n+1} = u^{k+1}
```

All element stiffness and internal-force contributions are assembled into a global CSR sparse system in the inner loop. The linear solve uses either the SciPy direct solver (UMFPACK) or PETSc CG + GAMG for large 3D models.

---

## Abaqus-compatible convergence check

After each Newton iteration, two norms are evaluated:

**Residual check**

$$|R_\text{max}| < \varepsilon \cdot \tilde{q}$$

where $R_\text{max}$ is the maximum nodal residual force and $\tilde{q}$ is the time-averaged reference force magnitude (computed over the current step to make the tolerance problem-independent).

**Displacement correction check**

$$|c_\text{max}| < c_n \cdot |u_\text{max}|$$

where $c_\text{max}$ is the largest nodal displacement correction in the current iteration and $u_\text{max}$ is the largest total displacement increment in the current step.

Default tolerances: $\varepsilon = c_n = 0.005$. Both checks must pass simultaneously. The displacement check is skipped when $u_\text{max} < 10^{-5}$ (nearly zero displacement increments).

---

## Load stepping and cutback

The load increment size $\Delta t$ is controlled automatically:

- If a step converges in fewer than the target iteration count, the next increment is grown (up to a user-specified maximum).
- If Newton fails to converge, the increment is cut back by a factor (default 0.25) and retried.
- If the increment size falls below the minimum, the analysis terminates with a non-convergence message.

---

## Stiffness update strategies

Two convergence strategies are available and selected per-step in the input:

| Strategy | Behaviour |
|---|---|
| `ABAQUS` | Standard full Newton; stiffness reformed every iteration |
| `SYNOPSYS` | Modified Newton with optional stiffness freeze; fewer factorizations |

The Abaqus strategy provides quadratic convergence near the solution. The Synopsys strategy trades convergence rate for cheaper iterations — useful when the factorisation cost dominates (large models, direct solver).

---

## Linear solver

For small to medium models, SciPy's UMFPACK-based direct solver (`scipy.sparse.linalg.spsolve`) is used. For large 3D models (tens of thousands of DOFs), Minion optionally delegates to PETSc with a CG + GAMG (algebraic multigrid) preconditioner. Node coordinates are passed to GAMG for geometry-aware coarsening, which substantially accelerates convergence on 3D elasticity systems compared to purely algebraic approaches.
