"""Regression: abnumber installed but ANARCI backend missing must FALL BACK.

pip-installing abnumber without the conda-only `anarci` backend makes
Chain() raise ModuleNotFoundError at first use. find_cdrs.py must degrade to
the regex anchors instead of aborting the run (it used to abort).

Run: python3 -I tests/test_anarci_fallback.py
"""
import sys
import types
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

VH = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
      "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS")


def _install_broken_abnumber():
    mod = types.ModuleType("abnumber")

    class Chain:  # noqa: D401 - test double
        def __init__(self, *a, **k):
            raise ModuleNotFoundError("No module named 'anarci'")

    mod.Chain = Chain
    sys.modules["abnumber"] = mod


def main() -> None:
    _install_broken_abnumber()
    import find_cdrs

    cdrs, engine, resolved, note = find_cdrs.annotate(VH, "VH", "kabat")
    assert engine == "anchor-regex", f"expected fallback, got {engine}"
    assert cdrs[2] == "WGGDGFYAMDY", cdrs
    assert "anarci" in note.lower(), note
    assert "conda" in note.lower(), "note should tell the user how to fix it"

    # And with no abnumber at all: fallback, but no noisy note.
    sys.modules.pop("abnumber", None)
    sys.modules["abnumber"] = None  # force ImportError on import
    import importlib
    importlib.reload(find_cdrs)
    cdrs2, engine2, _, note2 = find_cdrs.annotate(VH, "VH", "kabat")
    assert engine2 == "anchor-regex" and cdrs2[2] == "WGGDGFYAMDY"
    assert note2 == "", f"unexpected note when abnumber absent: {note2!r}"
    print("ok: ANARCI-missing falls back with guidance; absent abnumber is silent")


if __name__ == "__main__":
    main()
