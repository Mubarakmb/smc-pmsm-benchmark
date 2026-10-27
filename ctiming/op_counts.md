# Per-sample operation counts (from `ctrl_bench.c`, identical laws to `sim.py`)

Counts include the shared output stage `(dwr + a_hat*w + N)/chi`, the clip and the anti-windup test.
`Lm` = fractional memory length (2000 in the benchmark).

| Controller | mul | div | add/sub | special functions | persistent states | memory (double) |
|---|---|---|---|---|---|---|
| CSMC  | 3 | 1 | 3 | sgn | 0 | 0 B |
| ISMC  | 6 | 1 | 7 | sgn | 2 (I, e0) | 16 B |
| TSMC  | 6 | 1 | 6 | sgn x2, **pow** | 1 (I) | 8 B |
| ASMC  | 10 | 2 | 9 | sgn x2, fabs x2, fmax x2 | 2 (K, z) | 16 B |
| FOSMC | 2*Lm + 6 | 1 | 2*Lm + 5 | sgn, **pow x2** | Lm (error buffer) | 8*Lm B buffer + 16*Lm B weights (= 48 kB for Lm = 2000; weights can sit in flash) |
| STSMC | 5 | 1 | 4 | sgn x2, **sqrt**, fabs | 1 (v) | 8 B |

Measured execution time (gcc 13.3 -O2, Intel Xeon @ 2.80 GHz, shared cloud VM, 3 repetitions,
20 M calls each; FOSMC 0.4 M calls):

| Controller | ns per control step (3 runs) |
|---|---|
| CSMC  | 7.39 / 4.70 / 6.22 |
| ISMC  | 7.37 / 9.41 / 7.50 |
| TSMC  | 32.53 / 37.50 / 30.63 |
| ASMC  | 16.12 / 18.32 / 16.18 |
| FOSMC | 2698.75 / 3005.23 / 3434.67 |
| STSMC | 6.15 / 7.18 / 7.47 |

Notes
- Timings are relative (general-purpose CPU with FPU and caches, shared cloud VM). Absolute real-time feasibility on an
  embedded target (e.g. a microcontroller without double-precision FPU) must be measured on that
  target; this is stated as a limitation unless processor-in-the-loop timing is added.
- TSMC cost is dominated by `pow()`; FOSMC cost is linear in Lm (2 MAC per sample per operator).
