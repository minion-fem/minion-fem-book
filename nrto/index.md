---
title: "Neural Topology Optimisation"
---

# Neural Topology Optimisation

Topology optimisation finds the material layout that minimises structural compliance under a volume constraint. Rather than optimising a grid of independent density variables (SIMP), the neural reparametrisation approach represents the density field with a neural network — giving implicit regularisation and generalising naturally to unstructured meshes.

- [Overview](overview.md) — compliance minimisation, SIMP penalisation, training pipeline, FEM interface
- [FourierTOuNN](fourier_tounn.md) — random Fourier feature map, length-scale control, network architecture
- [Results](results.md) — 2D and 3D cantilever beam examples
