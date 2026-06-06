---
title: "Results"
---

# Topology Optimisation Results

---

## 2D Cantilever Beam

**Setup**: rectangular domain ($64 \times 32$ elements), clamped on the left edge, point load at the bottom-right corner. FourierTOuNN with volume fraction $v_f = 0.4$, penalty continuation $p: 1 \to 8$, 200 iterations.

:::{figure} images/design_evolution.gif
:width: 80%
:align: center

Design evolution over 200 iterations. Starting from a uniform density field, the optimiser progressively discovers load-carrying members and sharpens them toward 0–1 as the SIMP penalty exponent increases. Final compliance $J = 72.64$ — a **30.8% reduction** from the uniform-density baseline $J_0 = 105.0$.
:::

---

## 3D Cantilever Beam

**Setup**: 3D rectangular block, clamped on one face, distributed load on the centre of the opposite face. FourierTOuNN, volume fraction $v_f = 0.3$, 50 iterations.

:::{figure} images/design_evolution3d.gif
:width: 80%
:align: center

3D density field evolution visualised in ParaView. The optimiser carves out internal voids while retaining load-carrying arches and flanges. Final compliance $J \approx 2.68$ — a **17.6% reduction** from the uniform-density baseline $J_0 = 3.251$.
:::

---

## Summary

| Case | Iterations | $v_f$ | $J_0$ | $J_{\text{final}}$ | Reduction |
|---|---|---|---|---|---|
| 2D cantilever | 200 | 0.4 | 105.0 | 72.64 | 30.8% |
| 3D cantilever | 50 | 0.3 | 3.251 | 2.68 | 17.6% |
