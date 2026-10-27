# smcbench: a reproducible benchmark of sliding-mode speed control for PMSM drives

This repository contains the complete simulation framework, tuning and evaluation scripts,
raw results, and the scripts that generate every table and figure of the paper

> M. B. Aremu, A. Ajasa, A. Nasir, *A Unified Benchmarking Framework for Sliding-Mode Control
> Strategies in Permanent Magnet Synchronous Motor Speed Regulation* (submitted).

Repository: <https://github.com/Mubarakmb/smc-pmsm-benchmark>

It is the companion of our review of sliding-mode control for PMSM speed control
(A. Ajasa, M. B. Aremu, A. Nasir, *Arabian Journal for Science and Engineering*, 2026,
doi:[10.1007/s13369-026-11494-6](https://doi.org/10.1007/s13369-026-11494-6)), which identified
the lack of open models, common tuning and standard disturbance profiles as the main obstacle to
comparing published controllers. smcbench is our answer to that gap: six sliding-mode controllers
are compared on one data-sheet motor, one drive model, one tuning protocol and one set of tests,
and anyone can add a controller and compare it under exactly the same conditions.

---

## Contents

1. [What is benchmarked](#1-what-is-benchmarked)
2. [Quick start](#2-quick-start)
3. [Reproducing the paper](#3-reproducing-the-paper)
4. [Repository layout](#4-repository-layout)
5. [Plant and drive model](#5-plant-and-drive-model)
6. [Controllers](#6-controllers)
7. [Tuning protocol](#7-tuning-protocol)
8. [Evaluation](#8-evaluation)
9. [Metrics and objective](#9-metrics-and-objective)
10. [Result files](#10-result-files)
11. [Main results](#11-main-results)
12. [Adding your own controller](#12-adding-your-own-controller)
13. [Computational-cost benchmark (C)](#13-computational-cost-benchmark-c)
14. [Tests](#14-tests)
15. [Limitations](#15-limitations)
16. [How to cite and licence](#16-how-to-cite-and-licence)

---

## 1. What is benchmarked

| Abbrev. | Controller | Sliding variable *s* | Gain condition enforced in tuning |
|---|---|---|---|
| CSMC  | conventional SMC | *e* | *k* ≥ 1.05 *D* |
| ISMC  | integral SMC | *e* − *e*(0) + λ∫*e* | *k* ≥ 1.05 *D* |
| TSMC  | nonsingular integral terminal SMC | *e* + ν∫\|*e*\|^α sgn *e* | *k* ≥ 1.05 *D* |
| ASMC  | adaptive SMC (equivalent-control adaptation, Utkin & Poznyak 2013) | *e* | *K*⁺ ≥ 1.05 *D* |
| FOSMC | fractional-order SMC (fractional-PI surface, Grünwald–Letnikov) | *e* + *k*ᵢ D^−α *e* | *k* ≥ 1.05 *D* |
| STSMC | super-twisting algorithm | *e* | *k*₂ ≥ 1.05 *C*, *k*₁ = ρ√(*k*₂+*C*), ρ ≥ 1.01 (Seeber & Horn 2017) |

*e* = ω\* − ω is the speed error; *D* = *T*_rated/*J* = 4834.4 rad/s² bounds the rated-load
disturbance and *C* = 1.519×10⁶ rad/s³ bounds its derivative for a 5 ms load rise.

The five first-order controllers are each evaluated in two forms: with sgn(*s*) (discontinuous)
and with the same boundary layer sat(*s*/φ) (suffix **-BL**). Together with STSMC this gives
**11 variants**. No controller uses the load torque or its derivative.

## 2. Quick start

Requirements: Python 3 (tested with 3.13) and the packages in `requirements.txt` (the versions listed are those
with which all results were produced). A C compiler is needed only for the timing benchmark.

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python -m smcbench.demo                         # 11 tuned variants, nominal scenario
python -m smcbench.demo ISMC-BL,STSMC "E4 reversal"
PROTOCOL=free python -m smcbench.demo           # gains from the unconstrained protocol
```

The demo prints the main metrics and writes `figures/demo_<scenario>.png`. The first call
compiles the simulator with Numba (a few seconds; the result is cached). One 0.6 s simulation
(300 000 RK4 steps) then takes about 0.2 s.

A single simulation from Python:

```python
from smcbench import tune as T, evaluate as E
best = E.best_params()                         # best seed per variant (PROTOCOL, TAG from env)
r, log = E.evaluate('ISMC-BL', best['ISMC-BL']['x'], E.SCEN['E1 low speed'])
print(r['IAE_load'], r['Te_p2p'])              # metrics dict, see Section 9
# log columns: t, w_ref, w, w_meas, iq_ref, iq, id, Te, TL, s, vd, vq, sat, aux
```

### Windows

Install Python 3.13 from python.org (tick "Add python.exe to PATH"), open the repository folder,
and use the helper script from a terminal (Command Prompt or the VS Code terminal):

```bat
run_windows.bat setup     :: once: creates .venv and installs the packages
run_windows.bat check     :: test suite, should end with "45 passed"
run_windows.bat demo      :: simulates the continuous variants and opens the plots
run_windows.bat outputs   :: regenerates all tables and figures
```

Other demos: `.venv\Scripts\activate` and then, for example,
`python -m smcbench.demo ISMC-BL,STSMC "E4 reversal" --show`.

## 3. Reproducing the paper

`reproduce.sh` has three levels:

| Command | What it does | Time (2 cores) |
|---|---|---|
| `./reproduce.sh check` | test suite: re-simulates the tuned controllers and requires the stored metrics to be reproduced (Section 14) | ~1 min |
| `./reproduce.sh outputs` | regenerates all tables (`tables/`) and figures (`figures/`) from the shipped results | ~1 min |
| `./reproduce.sh full` | re-runs everything from scratch into `rerun/` (tuning of both protocols, load-rise study, evaluation, sensitivity studies, tables and figures) | ~5 h (one tuning run of 1920 evaluations takes ~6 min on one core) |

Results are deterministic: the tuning seeds, the Monte Carlo draws and the measurement-noise
sequences are all seeded, so a full re-run on the same software versions reproduces the shipped
files. The test suite and `outputs` were verified in a clean virtual environment; the regenerated
LaTeX tables are byte-identical to those in the paper. A from-scratch re-tuning of ISMC-BL (three
seeds) and its full evaluation in the clean environment reproduced the shipped gains and all
evaluation results exactly.

Mapping from scripts to the paper:

| Script | Produces | Paper |
|---|---|---|
| `smcbench.tune` | `results/tuning_results_*.json` | tuning protocol, tuned gains |
| `smcbench.evaluate` | `results/eval_results_*.json` | test scenarios, parameter variations, Monte Carlo |
| `smcbench.analyse_eval` | Friedman test, mean ranks, 166-weighting rank sensitivity (stdout) | Monte Carlo section |
| `smcbench.sensitivity` | `results/sensitivity_lyap.json` | solver-step convergence, sampling-time sensitivity |
| `smcbench.disturb` | `results/disturb_lyap.json` | periodic torque ripple, Coulomb friction |
| `smcbench.fo_memory_study` | `results/fo_memory_study.txt` | fractional-memory accuracy |
| `smcbench.fig_protocol` | `figures/fig_protocol.pdf`, `tables/tab_protocol.tex`, `results/protocol_compare.*` | constrained vs unconstrained tuning |
| `smcbench.make_tables` | `tables/tab_gains, tab_nominal, tab_scen, tab_pvar, tab_mc` | result tables |
| `smcbench.figures`, `smcbench.figures2` | `figures/fig_*.pdf` | result figures |
| `smcbench.graphical_abstract` | `figures/graphical_abstract.*` | graphical abstract |
| `ctiming/run_timing.sh` | execution time per control step | computational cost |

## 4. Repository layout

```
smcbench/            Python package
  sim.py             Numba simulation core: plant, drive, six control laws
  config.py          motor data, drive defaults, scenario vector
  metrics.py         metric definitions with event-tied windows
  tune.py            tuning protocol (search spaces, objective, differential evolution)
  evaluate.py        test scenarios, parameter variations, Monte Carlo
  analyse_eval.py    statistics (Friedman test, ranks, weight sensitivity)
  sensitivity.py     solver step and sampling time
  disturb.py         periodic ripple and Coulomb friction
  fo_memory_study.py accuracy of the truncated Grünwald–Letnikov operator
  fig_protocol.py    constrained vs unconstrained comparison
  figures.py, figures2.py, make_tables.py, graphical_abstract.py
  report_tuning.py   console summary of a tuning run (diagnostic)
  demo.py            quick start
  paths.py           output folders (override with SMCBENCH_ROOT)
ctiming/             C implementation of the six laws, timing benchmark, operation counts
results/             all raw results of the paper (JSON) and console logs (results/logs)
tables/, figures/    generated outputs
tests/               regression and physics tests
reproduce.sh         one-command reproduction
```

Environment variables used by the scripts:

| Variable | Default | Meaning |
|---|---|---|
| `PROTOCOL` | `lyap` | `lyap` = stability conditions enforced; `free` = unconstrained |
| `TAG` | `_final` | suffix of result files; `_final` = results reported in the paper |
| `LOAD_RISE` | `5e-3` | load-step rise time (s); sets *C* and the STSMC bounds |
| `N_WORKERS` | all CPUs | parallel tuning jobs |
| `SMCBENCH_ROOT` | repository root | where `results/`, `figures/`, `tables/` are read and written |

`smcbench.tune` and `smcbench.evaluate` write to `results/` under `SMCBENCH_ROOT`. Running them
with the default `TAG` in the repository root therefore works on the shipped files: tuning skips
every job that is already present (it is resumable), and evaluation overwrites the evaluation file
with identical values. Use `SMCBENCH_ROOT` (as `reproduce.sh full` does) or another `TAG` to
keep a re-run separate.

## 5. Plant and drive model

**Motor.** Siemens SIMOTICS S-1FK7063-5AF71, a surface-mounted PMSM (*L_d* = *L_q*). All values
come from the manufacturer data sheet:

| Quantity | Value | Note |
|---|---|---|
| pole pairs *p_n* | 4 | 8 poles |
| rated speed / torque / current | 3000 rpm / 7.3 N·m / 5.6 A rms | |
| maximum torque / current | 35 N·m / 28 A rms | current limit 28√2 = 39.6 A peak |
| *R_s* (20 °C) | 0.65 Ω | |
| *L_s* | 7.7 mH | rotating-field inductance |
| *ψ_f* | 0.1706 Wb | from *k_E* = 87.5 V/1000 rpm (20 °C); *k_T* = 1.37 N·m/A gives 0.1615 Wb hot |
| *J* | 15.1 kg·cm² | |
| *B* | 0 | not in the data sheet; varied in the robustness tests |
| DC link | 600 V | data-sheet family rating |

The data-sheet electrical and mechanical time constants (11.8 ms, 1.56 ms) are reproduced from
*R_s*, *L_s*, *J* and *k_T*, which confirms that *R_s* and *L_s* are per-phase values.

**Drive (identical for every controller).**
- Cascade: speed controller → *i_q*\* (limited to ±39.6 A) → PI current loops with *i_d*\* = 0.
- Current PI by pole–zero cancellation, bandwidth 400 Hz, with decoupling, back-EMF
  feedforward and back-calculation anti-windup.
- Voltage limit |*v_dq*| ≤ *V_dc*/√3 = 346 V (d-axis priority), average-value inverter.
- One-sample computational delay; speed and current loops sampled at *T_s* = 100 µs (ZOH).
- Speed measured by a 20-bit encoder with backward differencing (resolution 0.060 rad/s at
  100 µs); optional additive noise.
- One conditional-integration anti-windup rule for every speed controller: integral-type states
  are frozen while *i_q*\* is saturated and the error drives further into saturation.
- Continuous plant integrated with fixed-step RK4, *h* = 2 µs. All metrics change by at most
  0.4 % between 5, 2 and 1 µs and not at all between 2 and 1 µs.
- Load steps are raised-cosine with rise time *t_r* (default 5 ms), so d*T_L*/d*t* is bounded.
- Optional periodic torque ripple at six times the electrical frequency and smoothed Coulomb
  friction (`Trip`, `Tf` in the motor dictionary).

## 6. Controllers

All six share the output structure

χ̂ *i_q*\* = ω̇\* + â ω + *N*,  χ̂ = 1.5 *p_n* ψ_f / *J*,  â = *B*/*J* (nominal values),

so that, without saturation, ė = *d̃* − *N*, where *d̃* lumps the load and all model mismatch.
The term *N* uses the reaching law −*k* σ(*s*) − *q s* with σ = sgn or sat(·/φ):

| Variant | *N* | Parameters tuned |
|---|---|---|
| CSMC  | *k*σ(*s*) + *q s* | *k*, *q* |
| ISMC  | λ*e* + *k*σ(*s*) + *q s* | λ, *k*, *q* |
| TSMC  | ν\|*e*\|^α sgn *e* + *k*σ(*s*) + *q s* (no \|*e*\|^(α−1) term: nonsingular by construction) | ν, α, *k*, *q* |
| ASMC  | *K*(t)σ(*s*) + *q s*; τż + *z* = σ(*s*), *K̇* = γ₀*K* sgn(\|*z*\| − ᾱ) + projection on [μ, *K*⁺] | *q*, *K*⁺, γ₀, τ, ᾱ |
| FOSMC | *k*ᵢ D^(1−α)*e* + *k*σ(*s*) + *q s*; GL operators with memory *L_m* = 2000 | *k*ᵢ, α, *k*, *q* |
| STSMC | *k*₁\|*s*\|^½ sgn *s* + *v*, *v̇* = *k*₂ sgn *s* | ρ, *k*₂ (constrained) or *k*₁, *k*₂ (unconstrained) |

-BL variants add φ. In ASMC the lower bound and the initial gain are μ = *K*(0) = 0.01 *K*⁺ and the
projection gain is *M* = 10³. The implementation of all laws is in `smcbench/sim.py`
(function `simulate`); `smcbench/tune.py:build` maps the named parameters to the gain vector.

## 7. Tuning protocol

Identical for every variant (`smcbench/tune.py`):

- **Optimiser:** SciPy differential evolution, best/1/bin (SciPy's default strategy), mutation
  dithered in [0.5, 1], crossover 0.7, deferred updating, Sobol initial population of 32
  (SciPy rounds a Sobol population up to a power of two), 60 generations, no polishing.
- **Budget:** exactly 1920 objective evaluations per run, for every variant regardless of its
  number of parameters; three seeds (0, 1, 2) per variant; the best seed is evaluated. (A run
  stops earlier only if its whole population collapses to one point; this happened once, for
  unconstrained ISMC, seed 0, after 1568 evaluations.)
- **Search space:** log-uniform for gains, linear for exponents; switching-gain ranges are scaled
  by *D* (*C* for STSMC) so the same rule applies to every controller.
- **Training scenario:** ramp to rated speed at 5000 rad/s², rated-load step at 0.3 s (5 ms rise),
  run length 0.6 s, 20-bit encoder, all limits active. It is not used for reporting.
- **Protocols:** `lyap` enforces the gain conditions of Section 1 as lower bounds;
  `free` lowers the switching-gain bounds to 0.1 *D* and searches *k*₁ ∈ [1, 10⁴],
  *k*₂ ∈ [10², 10⁸] for STSMC.

Note for STSMC: the stored parameter vector `x` means (ρ, *k*₂) under `lyap` and (*k*₁, *k*₂)
under `free`. Pass `protocol='lyap'|'free'` to `tune.run`/`evaluate.evaluate`, or set
`PROTOCOL`, when simulating stored STSMC gains. The `params` field of each result stores the
physical gains, including *k*₁, in both cases.

## 8. Evaluation

`smcbench/evaluate.py`, all with the best tuned gains and the nominal model inside the controller:

| Set | Content |
|---|---|
| E0 | nominal: rated speed, rated-load step (training conditions, reported for reference) |
| E1 | 10 % rated speed, 50 % load step |
| E2 | 50 % rated speed, 100 % load step |
| E3 | rated speed, 25 / 50 / 150 % load step |
| E4 | speed reversal +rated → −rated at 0.3 s |
| E5 | quasi-step reference (10⁶ rad/s²) with rated load: current saturation and windup |
| P  | one-at-a-time plant changes: *R_s* ×1.39 (100 K hot) and ×0.7; *L* ×0.8, ×1.2; ψ ×0.947 (data-sheet hot/cold ratio) and ×0.9; *J* ×0.5, ×2, ×3; *B* = 0.003 N·m·s |
| MC | 200 random draws, identical for all variants: *R_s* ×U(0.8, 1.5), *L* ×U(0.8, 1.2), ψ ×U(0.9, 1.0), *J* ×U(0.5, 3), *B* ∈ U(0, 0.003), speed U(0.1, 1)×rated, load U(0, 1.5)×rated, speed-noise std 0.02 rad/s |

Further studies: sampling time 50/100/200 µs and solver step 10/5/2/1 µs
(`sensitivity.py`); load-rise time 5/20/50 ms with re-tuning of ISMC-BL, TSMC-BL and STSMC
(`reproduce.sh full`, step 2); periodic torque ripple of 1/3/5 % rated torque at 2 % and 10 %
rated speed, and Coulomb friction of 2/5 % rated torque in a slow reversal (`disturb.py`). The
ripple and friction amplitudes are assumptions (not data-sheet values) and are therefore swept.

## 9. Metrics and objective

Windows are tied to scenario events, so they are identical for all controllers
(`smcbench/metrics.py`): *t_L* = load-step time, *T* = run length, *W* = 50 ms steady window.

| Key | Definition | Unit |
|---|---|---|
| `IAE_trk`, `ITAE_trk` | ∫\|*e*\|, ∫*t*\|*e*\| over [0, *t_L*) | rad, rad·s |
| `OS_%` | overshoot after the ramp | % |
| `t_set` | settling into ±0.5 % band after the ramp (NaN = never) | s |
| `IAE_load` | ∫\|*e*\| over [*t_L*, *T*) | rad |
| `dw_max` | maximum \|*e*\| after the load step | rad/s |
| `t_rec` | recovery into ±0.5 % band after the load step (NaN = not recovered) | s |
| `ess_pre`, `ess_end` | mean \|*e*\| in the window before the load step / at the end | rad/s |
| `E_iq` | ∫*i_q*² | A²·s |
| `sat_%` | share of samples with saturated *i_q*\* | % |
| `TV_iqref` | total variation of *i_q*\* per second, final window | A/s |
| `HF_iq` | rms of *i_q* above 1 kHz, final window | A |
| `sw_rate` | sign changes of *s* per second, final window | 1/s |
| `w_p2p`, `Te_p2p` | peak-to-peak speed and torque, final window | rad/s, N·m |

Objective (training and evaluation):

𝒥 = IAE_trk/*N*₁ + IAE_load/*N*₂ + TV/*N*₃ + 10 ρ_sat + 100·𝟙_fail,

with *N*₁, *N*₂ = (0.5 % of the reference) × window length, *N*₃ = 2 *I*_rated,pk / 1 ms, TV the
total variation of *i_q*\* per second averaged over both steady windows, ρ_sat the saturated
fraction, and 𝟙_fail = 1 if `ess_pre` or `ess_end` exceeds 5 % of the reference. A simulation
that diverges returns 𝒥 = 1000.

## 10. Result files

| File | Content |
|---|---|
| `tuning_results_{lyap,free}_final.json` | list of 33 runs (11 variants × 3 seeds): `name`, `seed`, `J`, `nfev`, `x` (search vector in physical units), `params` (named gains), `metrics` (training scenario), `secs` |
| `tuning_results_lyap_tr20.json`, `_tr50.json` | load-rise study (ISMC-BL, TSMC-BL, STSMC; *t_r* = 20, 50 ms) |
| `eval_results_{lyap,free}_final.json` | per variant: `scen` (E0–E5 metrics), `pvar` (10 plant changes), `mc` (200 Monte Carlo runs; each with metrics, `J`, `w_ref`, `load`, the plant factors `f_Rs`, `f_Ld`, `f_Lq`, `f_psi`, `f_J`, `f_B`, and `TV_both`) |
| `sensitivity_lyap.json` | per variant: `h` (max relative metric change vs *h* = 1 µs) and `Ts` ([*T_s*, 𝒥, torque ripple]) |
| `disturb_lyap.json` | periodic-ripple cases ([speed p-p, mean \|*e*\|]) and Coulomb-friction cases (reversal IAE) |
| `protocol_compare.{json,txt}` | constrained vs unconstrained summary (gains relative to bounds, ranks, paired ratios, Wilcoxon *p*) |
| `fo_memory_study.txt` | GL operator error vs memory length |
| `logs/` | console output of the tuning runs and analyses |

## 11. Main results

Constrained protocol, 200-run Monte Carlo (Friedman χ² = 1555, *p* < 10⁻³⁰⁰):

| Variant | mean rank (1 = best) | median torque ripple (N·m) |
|---|---|---|
| ISMC-BL  | 1.49 | 0.16 |
| TSMC-BL  | 2.29 | 0.18 |
| FOSMC-BL | 2.46 | 0.26 |
| ASMC-BL  | 4.92 | 3.13 |
| STSMC    | 5.34 | 4.91 |
| ASMC     | 6.76 | 8.72 |
| FOSMC    | 7.28 | 15.4 |
| CSMC-BL  | 7.91 | 0.51 |
| TSMC     | 8.72 | 17.5 |
| ISMC     | 8.80 | 17.5 |
| CSMC     | 10.02 | 14.3 |

- ISMC-BL ranks first under all 166 tested weightings of the objective, under both tuning
  protocols.
- With the reaching condition enforced, the discontinuous variants chatter with torque ripple of
  two to three times the rated torque at 100 µs.
- Without it, the optimiser lowers the switching gain below *D* in ten of eleven variants. The
  discontinuous integral variants then stop chattering because they no longer slide (*s* settles
  at (*D* − *k*)/*q*), while the boundary-layer variants become less robust out of sample
  (CSMC-BL fails in 35 % of the Monte Carlo runs).
- STSMC gives the smallest median speed dip (1.03 rad/s) but about 30 times the torque ripple of
  ISMC-BL; its result depends strongly on the assumed load-rise time.
- FOSMC costs about 480 times CSMC per step, and its truncated operator has 42–62 % error with
  2000 samples of memory.

See the paper for the complete analysis.

## 12. Adding your own controller

1. **Control law.** In `smcbench/sim.py`, add a branch `elif ctrl == 7:` in the speed-loop block
   of `simulate` that computes `s` and `N` from the error `e`, the measured speed `w_meas` and the
   gain vector `G` (up to 12 entries). The common output stage, current limit and saturation flag
   follow automatically. Put any state update after the saturation check and respect `freeze`
   (the common anti-windup rule). Use only `e`, `w_meas`, the reference and the nominal model;
   the load is not available to controllers.
2. **Search space.** In `smcbench/tune.py`, add an entry to `SPECS`,
   `'MYSMC': (7, bl, [(name, lo, hi, log), ...])`, and map the parameters to `G` in `build`.
   Scale switching-gain bounds with `D` (or `C_DOT`) and state its stability condition as a
   lower bound, as for the other controllers.
3. **Run.** `python -m smcbench.tune MYSMC` (both protocols if relevant), then
   `python -m smcbench.evaluate MYSMC`. Add the name to the `ORDER` lists of the analysis and
   figure scripts to include it in the tables and plots.
4. **Report** the same quantities as Section 9, and keep the budget (1920 evaluations, 3 seeds)
   unchanged so that the comparison stays fair.

## 13. Computational-cost benchmark (C)

`ctiming/ctrl_bench.c` re-implements the six laws in C (same equations as `sim.py`, sgn
versions, including output stage, saturation and anti-windup) and times them on fixed
pseudo-random inputs:

```bash
./ctiming/run_timing.sh          # gcc -O2; 3 runs of 2e7 calls (FOSMC 4e5)
```

`ctiming/op_counts.md` lists the exact per-sample operation counts and memory, and the timings
reported in the paper (gcc 13.3 -O2, Intel Xeon @ 2.8 GHz; medians of 3 runs: CSMC 6.2 ns,
ISMC 7.5 ns, TSMC 32.5 ns, ASMC 16.2 ns, FOSMC 3.0 µs, STSMC 7.2 ns). The timings are relative
general-purpose-processor figures; embedded execution times must be measured on the target processor.

## 14. Tests

```bash
python -m pytest -q tests
```

The 45 tests check that:
- re-simulating the stored best gains reproduces the stored training objective for all 11 variants
  under both protocols, and selected test-scenario and Monte Carlo metrics (relative tolerance
  10⁻⁶);
- the boundary-layer equilibrium *s*\* = *D*/(*k*/φ + *q*) holds for CSMC-BL, ISMC-BL, TSMC-BL and FOSMC-BL;
- with unconstrained gains (*k* < *D*), ISMC, TSMC and FOSMC settle at *s* = (*D* − *k*)/*q* without
  changing sign;
- the steady-state torque balances the load, and current and voltage limits are respected;
- the Grünwald–Letnikov weights reduce to the identity and the backward difference for orders 0 and 1.

## 15. Limitations

The benchmark is simulation-based. The drive model includes sampling, delay, current and voltage
limits, encoder quantisation, noise and current-loop dynamics, but not PWM ripple, dead time,
fixed-point arithmetic, magnetic saturation or thermal dynamics. One motor is used, the
disturbance amplitudes of the ripple and friction studies are assumptions, the rankings depend on
the stated objective (their sensitivity to the weights is reported), and the timings are
relative. Details are in the limitations section of the paper.

## 16. How to cite and licence

If you use smcbench, please cite the paper above (full reference to be added on publication) and
the companion review:

```bibtex
@article{ajasa2026sliding,
  title   = {Sliding-Mode Control Strategies for {PMSM} Speed Control: A Comprehensive Review, Taxonomy, and Research Gaps},
  author  = {Ajasa, Abdullah and Aremu, Mubarak Badamasi and Nasir, Ali},
  journal = {Arabian Journal for Science and Engineering},
  year    = {2026},
  doi     = {10.1007/s13369-026-11494-6}
}
```

`CITATION.cff` contains the software citation.

- **Code:** MIT License (`LICENSE`).
- **Results data** (`results/`) and generated tables and figures: Creative Commons Attribution 4.0
  International (CC BY 4.0).

### Key references for the implemented methods

- V. I. Utkin, A. S. Poznyak, Adaptive sliding mode control with application to super-twist
  algorithm: Equivalent control method, *Automatica* 49 (2013) 39–47.
- R. Seeber, M. Horn, Stability proof for a well-established super-twisting parameter setting,
  *Automatica* 84 (2017) 241–243, doi:10.1016/j.automatica.2017.07.002.
- J. A. Moreno, M. Osorio, Strict Lyapunov functions for the super-twisting algorithm,
  *IEEE Trans. Automatic Control* 57 (2012) 1035–1040, doi:10.1109/TAC.2012.2186179.
- Y. Feng, X. Yu, Z. Man, Non-singular terminal sliding mode control of rigid manipulators,
  *Automatica* 38 (2002) 2159–2167, doi:10.1016/S0005-1098(02)00147-4.
- I. Podlubny, *Fractional Differential Equations*, Academic Press, 1999.
- A. Levant, Sliding order and sliding accuracy in sliding mode control, *Int. J. Control* 58
  (1993) 1247–1263.

### Contact

Mubarak Badamasi Aremu, Control and Instrumentation Engineering Department, King Fahd University
of Petroleum & Minerals (KFUPM), Dhahran, Saudi Arabia.
