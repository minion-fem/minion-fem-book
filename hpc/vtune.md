---
title: Profiling with Intel VTune
---

# Profiling with Intel VTune

Before optimising anything, measure. VTune tells you exactly which functions allocate the most memory and where the CPU is stalling on memory access — including cross-socket NUMA traffic.

## Build a Debug-Symbols Release Binary

VTune needs symbol information to map samples back to function names and source lines. Do **not** use a plain Debug build: the compiler disables inlining and AVX2 vectorisation in Debug mode, which changes the memory access pattern you are trying to study.

Add these flags in `CMakeLists.txt`:

```cmake
add_compile_options(
    -O3
    -g                        # emit debug symbols
    -fno-omit-frame-pointer   # keep the frame pointer for accurate call stacks
)
```

This gives a Release binary with full optimisations but enough metadata for VTune to reconstruct call stacks.

## Memory Consumption Analysis

To trace every `malloc`/`new` call and find which functions allocate the most heap memory:

```bash
vtune -collect memory-consumption \
      -result-dir ./mem_profile_data \
      -- ./your_solver_exe [args]
```

- `-collect memory-consumption` — records every allocation with its call stack.
- `-result-dir` — writes raw data to this directory.

For a quick terminal summary:

```bash
vtune -report summary -result-dir ./mem_profile_data
```

This prints peak heap usage and the top-5 allocating functions.

## Viewing Results in the GUI

CentOS servers typically have no GUI. Transfer the result directory to a Windows machine for the full VTune interface:

```bash
tar -czvf mem_data.tar.gz ./mem_profile_data
# then scp to local machine
```

Open **Intel VTune Profiler** → **Open Result…** → select the unpacked directory. In the **Bottom-up** tab, sort by **Allocation Size** descending to see the biggest offenders.

## Memory Access Analysis (NUMA Diagnosis)

If peak memory looks reasonable but the solver is still slow, the bottleneck may be cross-socket memory traffic. Switch to memory-access mode:

```bash
vtune -collect memory-access \
      -result-dir ./access_profile_data \
      -- ./your_solver_exe [args]
```

Look for **Remote Hop** counts in the results. A large Remote Hop count means threads on one CPU socket are reading data whose physical pages live on the other socket's memory controller. On a two-socket server this can halve effective memory bandwidth.

The fix is either pinning threads to a single NUMA node (`numactl --cpunodebind=0 --membind=0`) or restructuring data initialisation so that each socket's threads touch only the memory they will later compute on (first-touch policy).
