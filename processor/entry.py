"""
Container entry point: fixes the thread count from THREADS (or the cgroup CPU
quota) before numpy is imported, then runs main() once.
"""
import math
import os
import sys


def cgroup_cpus() -> int | None:
    """CPUs from the cgroup quota (v2 cpu.max, else v1 cfs), None when unlimited or unknown."""
    try:
        quota, period = open("/sys/fs/cgroup/cpu.max").read().split()[:2]
        if quota != "max":
            return max(1, math.ceil(int(quota) / int(period)))
        return None
    except (OSError, ValueError):
        pass
    try:
        quota = int(open("/sys/fs/cgroup/cpu/cpu.cfs_quota_us").read())
        period = int(open("/sys/fs/cgroup/cpu/cpu.cfs_period_us").read())
        if quota > 0:
            return max(1, math.ceil(quota / period))
    except (OSError, ValueError):
        pass
    return None


def set_threads() -> int:
    raw = os.environ.get("THREADS", "").strip()
    n = int(raw) if raw.lstrip("-").isdigit() else 0
    source = "THREADS"
    if n <= 0:
        n = cgroup_cpus()
        source = "cgroup quota"
        if n is None:
            n, source = 2, "fallback (no cgroup CPU quota)"
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS",
                "NUMEXPR_NUM_THREADS"):
        os.environ[var] = str(n)
    os.environ["SC_QC_THREADS"] = str(n)
    os.environ["SC_QC_THREADS_SOURCE"] = source
    return n


set_threads()

from processor.main import main  # noqa: E402

sys.exit(main())
