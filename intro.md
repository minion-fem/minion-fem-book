---
title: FEM under the hood
---

# FEM under the hood

This site documents the theory and implementation of nonlinear finite element algorithms — element formulations, constitutive models, implicit solvers, sparse linear algebra, and physics-informed neural topology optimisation. The reference implementation is **MinionFem**, a Python FEM solver built from scratch by C. Ching and Y. Kong.

The goal is to show not just what the equations are, but how they are implemented in a real solver — the design choices, the numerical traps, and the connection from theory to code.

---

## Topics covered

- **Element formulation** — isoparametric elements, strain-displacement matrix $\mathbf{B}$, geometric nonlinearity, B-bar locking fix
- **Material models** — linear elasticity, J2 plasticity with return mapping, Anand viscoplasticity
- **Nonlinear solver** — implicit Newton–Raphson, load stepping, Abaqus convergence criteria
- **Linear algebra** — CSR sparse assembly, CHOLMOD / PETSc GAMG solver backends
- **Neural topology optimisation** — FourierTOuNN, SIMP compliance minimisation, unstructured mesh support


