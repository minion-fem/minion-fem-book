---
title: "FourierTOuNN"
---

# FourierTOuNN

Minion implements two network architectures. The plain **TOuNN** uses coordinate inputs directly; **FourierTOuNN** first maps coordinates through a random Fourier feature layer that enforces a minimum length scale on the resulting density field.

---

## The length-scale problem

A plain MLP with LeakyReLU activations can in principle approximate any function, including ones with features much finer than the desired minimum member size. Without additional control, the optimiser may produce thin ligaments or point-like void pockets that would be impossible to manufacture.

The key observation: a function's minimum feature size is determined by its spatial frequency content. If we **band-limit the input representation** to frequencies below $1/(2 r_{\min})$, the network output is forced to vary on length scales no finer than $r_{\min}$.

---

## Random Fourier feature map

The Fourier feature map replaces the raw coordinate $\mathbf{x} \in \mathbb{R}^d$ with a higher-dimensional embedding:

$$\phi(\mathbf{x}) = \bigl[\sin(2\pi \mathbf{B}\mathbf{x}),\; \cos(2\pi \mathbf{B}\mathbf{x})\bigr] \in \mathbb{R}^{2M}$$

where $\mathbf{B} \in \mathbb{R}^{M \times d}$ is a fixed random frequency matrix drawn at initialisation. Each row of $\mathbf{B}$ is a frequency vector; $M$ is the embedding size (default 64).

The frequency matrix controls length scale: row frequencies near $1/(2r)$ produce features of size $\sim r$. By sampling frequencies from the interval $[1/(2r_{\max}),\; 1/(2r_{\min})]$, we exclude both very coarse (uninformative) and very fine (sub-manufacturing) spatial modes.

### Frequency initialisation

Two sampling strategies are available:

**Gaussian** (default): rows are drawn from a normal distribution centred at the midpoint frequency with standard deviation set so that 99.7% of mass falls within $[1/(2r_{\max}),\; 1/(2r_{\min})]$:

$$\mu = \frac{1}{2}\!\left(\frac{1}{2r_{\min}} + \frac{1}{2r_{\max}}\right), \qquad \sigma = \frac{1}{6}\!\left(\frac{1}{2r_{\min}} - \frac{1}{2r_{\max}}\right)$$

**Uniform**: rows are drawn with random sign and magnitude uniform in $[1/(2r_{\max}),\; 1/(2r_{\min})]$.

The frequency matrix $\mathbf{B}$ is registered as a buffer (non-trainable) — it is fixed throughout training.

---

## Network architecture

Both architectures share the same MLP backbone:

```
Input → [optional Fourier embedding] → FC → LN → Act → … → FC → Sigmoid → ρ ∈ (0,1)
```

| Component | TOuNN | FourierTOuNN |
|---|---|---|
| Input dimension | $d$ (2 or 3) | $2M$ (default 128 for $M=64$) |
| Hidden layers | 5 | 5 |
| Neurons per layer | 128 | 128 |
| Activation | LeakyReLU | LeakyReLU |
| Normalisation | LayerNorm per layer | LayerNorm per layer |
| Output | Sigmoid → $\rho \in (0,1)$ | Sigmoid → $\rho \in (0,1)$ |
| Weight init | Xavier normal | Xavier normal |

**Why LayerNorm?** The Fourier embedding concentrates activations in the range $[-1,1]$; LayerNorm prevents saturation in deeper layers and stabilises training.

**Why LeakyReLU over ReLU?** The density field must be smooth everywhere; dying ReLU units produce discontinuous density patches. LeakyReLU retains a small gradient for negative inputs.

---

## Generalisation to unstructured meshes

Classical density-filter methods require a structured Cartesian grid or explicit mesh neighbour lists. FourierTOuNN works on any mesh: element centroid coordinates $\mathbf{x}_e$ are the only geometric input. The `MinionToFemInterface` normalises them to $[0,1]^d$ using the bounding box of the mesh before passing them to the network, so that `rmin`/`rmax` are in mesh-normalised units.

This means the same network and the same hyperparameters apply to:
- Regular hexahedral meshes (NRTO paper setting)
- Tetrahedral meshes from automatic meshing
- Curved or non-convex domains

---

## Hard volume constraint

The standard volume constraint uses a soft Lagrangian penalty. Optionally, a **hard volume constraint** can be enforced by adding a learned shift $b$ to the pre-sigmoid logits:

$$\rho_e = \sigma(z_e + b), \quad \text{where } b \text{ satisfies } \frac{1}{n}\sum_e \sigma(z_e + b) = v_f$$

The shift $b$ is found by a root-finding step (bisection) after each forward pass. The root-finding operation is differentiable (via `autograd`) so gradients still flow back through $b$ to the network weights.

---

## When to use which architecture

| Situation | Recommendation |
|---|---|
| Default 2D/3D problem | **FourierTOuNN** — better length-scale control |
| Very coarse mesh ($< 500$ elements) | **TOuNN** — Fourier features over-parameterise |
| Need to inspect frequency content | **FourierTOuNN** with `freq_init='gaussian'` |
| Debugging / fastest iteration | **TOuNN** — simpler, fewer hyperparameters |
