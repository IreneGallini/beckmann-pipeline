"""
R/F prediction for one molecule's PySCF scan series
(backend.pipeline.wcnmax_pyscf.run_scan_series output). The prediction rule
itself is vendored unchanged from beckmann.dft.descriptors/wcnmax_rule -- no
reimplementation, just called on freshly-computed rows instead of ones read
from a benchmark CSV.
"""
from beckmann_core.wcnmax_rule import find_wcnmax_minimum, predict_from_wcnmax


def predict_outcome(mol: str, rows: list[dict]) -> tuple[str, dict | None]:
    """('R'|'F', minimum-dict-or-None) -- same two vendored functions the
    main pipeline's own benchmark run uses (beckmann.dft.wcnmax_rule.
    run_wcnmax_benchmark), applied to one molecule's rows directly instead
    of looping over a benchmark corpus."""
    minimum = find_wcnmax_minimum(mol, rows)
    return predict_from_wcnmax(minimum), minimum
