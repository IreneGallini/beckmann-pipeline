"""
R/F prediction and wCNmax-vs-R(N-O) plot for one molecule's PySCF scan series
(backend.pipeline.wcnmax_pyscf.run_scan_series output). The prediction rule
itself is vendored unchanged from beckmann.dft.descriptors/wcnmax_rule -- no
reimplementation, just called on freshly-computed rows instead of ones read
from a benchmark CSV.
"""
import base64
import io

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from beckmann_core.wcnmax_rule import find_wcnmax_minimum, predict_from_wcnmax, resolve_series


def predict_outcome(mol: str, rows: list[dict]) -> tuple[str, dict | None]:
    """('R'|'F', minimum-dict-or-None) -- same two vendored functions the
    main pipeline's own benchmark run uses (beckmann.dft.wcnmax_rule.
    run_wcnmax_benchmark), applied to one molecule's rows directly instead
    of looping over a benchmark corpus."""
    minimum = find_wcnmax_minimum(mol, rows)
    return predict_from_wcnmax(minimum), minimum


def plot_wcnmax_series(mol: str, rows: list[dict]) -> str:
    """wCNmax vs R(N-O), base64-encoded PNG (embedded directly in the API
    response -- no separate static file, no frontend charting library)."""
    by_stage = {r["stage"]: r for r in rows if r["mol"] == mol and r["channel"] == "cn"}
    series = resolve_series(by_stage)
    r_values = [r["R_NO"] for r in series]
    w_values = [r["weight"] for r in series]

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(r_values, w_values, marker="o", color="#2b6cb0")
    ax.set_xlabel("R(N-O) / Å")
    ax.set_ylabel("wCNmax")
    ax.set_title(f"{mol}: wCNmax vs R(N-O)\n(PySCF, AIMNet2 geometry -- no Gaussian/NBO7)")
    fig.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=120)
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("ascii")
