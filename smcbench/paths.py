"""Output locations. By default everything lives under the repository root; set the
environment variable SMCBENCH_ROOT to write elsewhere."""
import os

ROOT = os.environ.get('SMCBENCH_ROOT', os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RESULTS = os.path.join(ROOT, 'results')
FIGS = os.path.join(ROOT, 'figures')
TABLES = os.path.join(ROOT, 'tables')
for _d in (RESULTS, FIGS, TABLES):
    os.makedirs(_d, exist_ok=True)
