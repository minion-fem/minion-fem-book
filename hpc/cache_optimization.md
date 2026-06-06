---
title: Cache-Friendly Data Layout
---

# Cache-Friendly Data Layout

This is about **rearranging data to match how the CPU likes to read memory**. The layout choice — AoS vs. SoA — determines whether the CPU's vector units run at peak throughput or sit idle waiting for memory.

## AoS vs. SoA

**AoS (Array of Structures)** is the intuitive layout. For finite element node coordinates:

```cpp
struct Node {
    double x, y, z;
};
Node nodes[1000000];
```

**SoA (Structure of Arrays)** splits the fields into separate contiguous arrays:

```cpp
struct Nodes {
    double x[1000000];
    double y[1000000];
    double z[1000000];
};
```

## Why Layout Determines Cache Hit Rate

When the CPU fetches data from RAM, it does not fetch a single value — it fetches an entire **cache line** of 64 bytes at once.

**Consider summing all node X coordinates:**

- **AoS:** fetching `nodes[0]` also loads `nodes[0].y` and `nodes[0].z` into the cache line. Those fields are not needed right now, so two-thirds of every cache line is wasted. The CPU must return to RAM much more frequently. This is a **cache miss**.

- **SoA:** fetching `nodes.x[0]` also loads `nodes.x[1]` through `nodes.x[7]` (8 doubles = 64 bytes). All of them are immediately useful. The CPU processes all 8 values before touching RAM again. This is a **cache hit**.

## Why This Matters for FEM

In a FEM solver, the dominant operation is $\mathbf{F} = \mathbf{K} \mathbf{u}$ — large sparse matrix–vector products with millions of entries. Two consequences of layout:

**AVX2 vectorisation.** Modern CPUs have 256-bit AVX2 registers that process 4 `double` values simultaneously with a single instruction.

- With **AoS** (x, y, z, x, y, z, …), adjacent memory positions hold different fields. The CPU cannot fill an AVX2 register with four consecutive x-values because they are not contiguous. Vectorisation is blocked.
- With **SoA** (x, x, x, x, …), the CPU loads four consecutive x-values directly into the AVX2 register. One instruction, four multiplications.

In practice, Intel MKL exploits exactly this — its sparse routines are written against SoA-style contiguous buffers. Passing data in AoS layout forces an implicit copy or scatter/gather that erases much of MKL's advantage.

**Compiler auto-vectorisation.** With `-O3 -mavx2`, the compiler will attempt to vectorise `for` loops automatically. It can only succeed when the loop body touches contiguous memory. SoA exposes that contiguity; AoS hides it.

## Practical Rule

If an array is iterated over a single field at a time (which is true of most FEM assembly and solver loops), prefer SoA. If elements are always accessed as complete units and the array is small, AoS is fine.

## AVX2 in Brief

AVX2 (Advanced Vector Extensions 2) is the SIMD capability present on modern Intel and AMD CPUs. A 256-bit register holds 4 `double` or 8 `float` values. One `_mm256_add_pd` instruction adds four pairs of doubles simultaneously — a theoretical 4× throughput over scalar code.

Three ways to use it in practice:

1. **Compiler auto-vectorisation** (`-mavx2 -O3`): let the compiler rewrite loops. Works well for SoA loops.
2. **Use a tuned library** (MKL, Eigen, OpenBLAS): the recommended path. The library's authors have already written the intrinsics.
3. **Intrinsic functions** (`_mm256_add_pd`, etc.): hand-written SIMD for hot kernels where the compiler cannot figure it out.
