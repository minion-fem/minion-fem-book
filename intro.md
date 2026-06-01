---
title: Minion FEM
---

# Minion FEM

Minion is a finite element solver I wrote from scratch in Python. It reads Abaqus-compatible `.inp` files, runs nonlinear implicit structural analyses, and writes field output for post-processing. This site documents the theory behind the key algorithms and shows the solver in action.

---

## Capabilities

### Element library

| Element | Type | Integration pts |
|---|---|---|
| C3D4 | 4-node linear tetrahedron | 1 |
| C3D8 | 8-node linear hexahedron | 8 |
| CPS3 / CPS4 | Plane stress triangle / quad | 1 / 4 |
| CPE3 / CPE4 | Plane strain triangle / quad | 1 / 4 |
| B31 / B32 | 2-node / 3-node Timoshenko beam | 1 / 2 |

### Material models

- **Isotropic linear elasticity** — standard $(\lambda, \mu)$ Lamé formulation
- **J2 plasticity** — elastic predictor / return-mapping corrector with isotropic and kinematic hardening; consistent algorithmic tangent
- **Viscoplasticity** — Perzyna-type overstress model

### Solver features

- Implicit Newton–Raphson with automatic load stepping
- Abaqus-compatible convergence check (time-averaged force norm + displacement correction)
- Geometric nonlinearity (NLGEOM) via Updated Lagrangian: Hughes–Winget objective stress transport, midpoint strain increment, material + geometric tangent stiffness
- Logarithmic (Hencky) strain output for large-deformation field results

### Input / output

- Abaqus `.inp` reader (parts, assemblies, solid sections, beam sections, surfaces, steps, loads)
- Custom `.cmd` / dfISE mesh reader
- JSON intermediate format
- VTK field output

### Topology optimisation (NRTO)

Minion includes a neural reparametrisation topology optimiser (NRTO) that uses a Fourier-feature neural network to represent the density field, trained by back-propagating through the FEM compliance objective.

---

## Architecture

The solver is organised into three phases:

```
Input (.inp / .cmd / .json)
        │
        ▼
  ReadInputPhase       ← parse and validate; produce IxFem data objects
        │
        ▼
  PxFem packager       ← assemble AxFem: mesh, sections, materials, steps
        │
        ▼
  AnalysisPhase        ← time-step loop → Newton loop → element update loop
        │
        ▼
  VTK / JSON output
```

Each analysis step drives a Newton–Raphson loop. Element stiffness and internal force contributions are assembled into a global CSR sparse system and solved with a direct (SciPy) or iterative (PETSc CG + GAMG) solver.

---

## Site map

- **Element formulation** — shape functions, strain-displacement matrix $\mathbf{B}$, geometric nonlinearity, numerical integration
- **Material models** — J2 return mapping, consistent tangent, hardening laws
- **Nonlinear solver** — load increments, Newton iterations, Abaqus convergence criteria
- **Topology optimisation** — FourierTOuNN, compliance minimisation, 2D and 3D demos

The **Appendix** covers the sparse linear algebra background: stiffness-matrix sparsity patterns, direct assembly, Richardson iteration, and conjugate gradients.

```{tableofcontents}
```
