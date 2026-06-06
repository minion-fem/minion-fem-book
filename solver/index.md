---
title: "Nonlinear Solver"
slug: nonlinear-solver
---

# Nonlinear Solver

Nonlinear structural problems — large deformation, plasticity, contact — are solved incrementally. The total load is divided into steps; within each step, Newton–Raphson iteration drives the residual to zero. This section covers the incremental-iterative framework, convergence criteria, and load stepping strategies.

- [Newton–Raphson Method](newton_raphson.md) — incremental-iterative structure, Abaqus convergence check, automatic load stepping
