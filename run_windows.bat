@echo off
REM Windows helper for smcbench. Run from the repository folder (double-click or from a terminal):
REM   run_windows.bat setup     create .venv and install the packages (once)
REM   run_windows.bat check     run the test suite (re-simulates the tuned controllers)
REM   run_windows.bat demo      simulate the boundary-layer variants and STSMC in the nominal scenario
REM   run_windows.bat outputs   regenerate all tables and figures from the stored results
setlocal
cd /d "%~dp0"
set MODE=%1
if "%MODE%"=="" set MODE=check

if "%MODE%"=="setup" (
    py -3 -m venv .venv || python -m venv .venv
    call .venv\Scripts\activate.bat
    python -m pip install --upgrade pip
    python -m pip install -r requirements.txt
    goto :eof
)

if not exist .venv\Scripts\activate.bat (
    echo The virtual environment does not exist yet. Run:  run_windows.bat setup
    exit /b 1
)
call .venv\Scripts\activate.bat

if "%MODE%"=="check" (
    python -m pytest -q tests
    goto :eof
)
if "%MODE%"=="demo" (
    python -m smcbench.demo CSMC-BL,ISMC-BL,TSMC-BL,ASMC-BL,FOSMC-BL,STSMC "E0 nominal" --show
    goto :eof
)
if "%MODE%"=="outputs" (
    python -m smcbench.make_tables
    python -m smcbench.figures
    python -m smcbench.figures2
    python -m smcbench.fig_protocol
    python -m smcbench.graphical_abstract
    set PROTOCOL=lyap
    python -m smcbench.analyse_eval
    set PROTOCOL=free
    python -m smcbench.analyse_eval
    python -m smcbench.fo_memory_study
    goto :eof
)
echo Unknown option "%MODE%". Use: setup, check, demo or outputs.
