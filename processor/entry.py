"""
Container entry point: fixes the thread count (THREADS, else the cgroup CPU quota, via the shared
safe_inputs.cpu_threads) before numpy is imported, then runs main() once.
"""
import os
import sys

from processor.safe_inputs import cpu_threads

raw = os.environ.get("THREADS", "").strip()
n = cpu_threads(int(raw) if raw.isdigit() else 0)
for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[var] = str(n)
os.environ["SC_QC_THREADS"] = str(n)

from processor.main import main  # noqa: E402

sys.exit(main())
