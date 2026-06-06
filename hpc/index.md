---
title: High-Performance Computing
---

# High-Performance Computing

A FEM solver spends most of its wall-clock time in two places: assembling the global stiffness matrix and solving the resulting linear system. Getting that time down is not about clever algorithms alone — it requires understanding how modern hardware actually moves data.

This chapter covers two practical techniques used in MinionFem's development:

- **Cache-friendly data layout** — how AoS vs. SoA choices determine whether the CPU's vector units run at peak throughput or stall waiting for memory.
- **Performance profiling with VTune** — how to measure where time and memory are actually going, rather than guessing.
