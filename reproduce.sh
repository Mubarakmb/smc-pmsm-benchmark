#!/usr/bin/env bash
# Reproduce the results of the paper.
#   ./reproduce.sh check     run the test suite (re-simulates tuned controllers, ~1 min)
#   ./reproduce.sh outputs   regenerate all tables and figures from the shipped results (~5 min)
#   ./reproduce.sh full      re-run everything from scratch into ./rerun (several CPU-hours)
set -euo pipefail
cd "$(dirname "$0")"
MODE=${1:-check}

outputs () {
  python -m smcbench.make_tables
  python -m smcbench.figures
  python -m smcbench.figures2
  python -m smcbench.fig_protocol
  python -m smcbench.graphical_abstract
  PROTOCOL=lyap python -m smcbench.analyse_eval
  PROTOCOL=free python -m smcbench.analyse_eval
  python -m smcbench.fo_memory_study
}

case "$MODE" in
  check)
    python -m pytest -q tests ;;
  outputs)
    outputs ;;
  full)
    export SMCBENCH_ROOT="$PWD/rerun" TAG=_final
    mkdir -p "$SMCBENCH_ROOT"
    ALL=CSMC,ISMC,TSMC,ASMC,FOSMC,STSMC,CSMC-BL,ISMC-BL,TSMC-BL,ASMC-BL,FOSMC-BL
    # 1) tuning: 11 variants x 3 seeds x 1920 evaluations, both protocols
    PROTOCOL=lyap python -m smcbench.tune $ALL | tee "$SMCBENCH_ROOT/tune_log_lyap.txt"
    PROTOCOL=free python -m smcbench.tune $ALL | tee "$SMCBENCH_ROOT/tune_log_free.txt"
    # 2) load-rise study: re-tune the three leading continuous designs for t_r = 20 and 50 ms
    PROTOCOL=lyap LOAD_RISE=20e-3 TAG=_tr20 python -m smcbench.tune ISMC-BL,TSMC-BL,STSMC
    PROTOCOL=lyap LOAD_RISE=50e-3 TAG=_tr50 python -m smcbench.tune ISMC-BL,TSMC-BL,STSMC
    # 3) out-of-sample evaluation (test scenarios, parameter variations, 200-run Monte Carlo)
    PROTOCOL=lyap python -m smcbench.evaluate
    PROTOCOL=free python -m smcbench.evaluate
    # 4) sampling-time / solver-step sensitivity and non-constant disturbances
    PROTOCOL=lyap python -m smcbench.sensitivity | tee "$SMCBENCH_ROOT/sens.txt"
    PROTOCOL=lyap python -m smcbench.disturb | tee "$SMCBENCH_ROOT/disturb.txt"
    # 5) tables and figures (written to rerun/tables and rerun/figures)
    outputs ;;
  *)
    echo "usage: $0 [check|outputs|full]"; exit 1 ;;
esac
