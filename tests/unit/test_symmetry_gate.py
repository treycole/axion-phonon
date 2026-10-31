"""tests/symmetry/symmetry_lib.failures: passes a clean result, and catches a broken operation, a part over its
threshold, a control that does not fail, and a non-zero Tr F."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "symmetry"))
import symmetry_lib  # noqa: E402


def clean():
    n_ops = 12
    rows = np.zeros((n_ops + 2, 5))
    rows[1:n_ops, 1:] = [0.0, 4e-4, 1.5e-3, 2e-4]          # internal, cross, external, total (the 10-31 values)
    rows[n_ops:, 4] = [0.4, 0.02]                            # controls fail
    return dict(n_ops=n_ops, names=np.array([f"op{i}" for i in range(n_ops)] + ["C4z", "-E"]),
                dft_ueV=np.full((n_ops, 2), 0.06), generic=rows.copy(), tube=rows.copy(),
                trace_generic=1e-6, trace_tube=1e-6)


def test_clean_passes():
    assert symmetry_lib.failures(clean()) == []


def test_each_defect_is_caught():
    cases = {
        "total": lambda r: r["generic"].__setitem__((5, 4), 0.01),
        "internal": lambda r: r["tube"].__setitem__((3, 1), 1e-3),
        "bands": lambda r: r["generic"].__setitem__((2, 0), 0.5),
        "DFT": lambda r: r["dft_ueV"].__setitem__((4, 0), 30.0),
        "does not fail": lambda r: r["generic"].__setitem__((12, 4), 1e-3),
        "Tr F": lambda r: r.__setitem__("trace_tube", 1e-2),
    }
    for word, break_it in cases.items():
        r = clean()
        break_it(r)
        bad = symmetry_lib.failures(r)
        assert any(word in b for b in bad), (word, bad)
