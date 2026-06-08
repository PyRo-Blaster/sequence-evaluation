"""Shared helpers for the sequence-evaluation scripts.

Stdlib-only by design so it can be imported from any of the uv single-file
scripts without adding to their dependency sets. When a script is run with
``uv run scripts/<name>.py`` the script directory is on ``sys.path``, so a
plain ``import _common`` resolves.
"""

from __future__ import annotations

import json
import sys

# The 20 standard amino acids.
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
# Ambiguous / non-standard codes that ProtParam cannot weigh.
AMBIGUOUS_AA = set("BJOUXZ")


class SequenceError(ValueError):
    """Raised when an input sequence cannot be used for analysis."""


def clean_sequence(seq: str, *, name: str = "sequence", allow_stop: bool = True) -> str:
    """Normalise an amino-acid string and validate its alphabet.

    Uppercases, strips whitespace/newlines, removes a single trailing stop
    (``*``). Raises :class:`SequenceError` on ambiguous or unknown residues so
    callers fail loudly instead of crashing deep inside a dependency.
    """
    if seq is None:
        raise SequenceError(f"{name}: no sequence provided")

    # Drop all whitespace, including internal (FASTA line wrapping pasted in).
    cleaned = "".join(seq.split()).upper()
    if allow_stop:
        cleaned = cleaned.rstrip("*")

    if not cleaned:
        raise SequenceError(f"{name}: empty after cleaning")

    bad = sorted({c for c in cleaned if c not in STANDARD_AA})
    if bad:
        ambiguous = [c for c in bad if c in AMBIGUOUS_AA]
        unknown = [c for c in bad if c not in AMBIGUOUS_AA]
        parts = []
        if ambiguous:
            parts.append(f"ambiguous code(s) {','.join(ambiguous)}")
        if unknown:
            parts.append(f"invalid character(s) {','.join(unknown)}")
        raise SequenceError(
            f"{name}: contains {' and '.join(parts)}; "
            "MW/pI and most analyses require the 20 standard amino acids"
        )
    return cleaned


def emit_json(payload: object) -> None:
    """Print a JSON document to stdout (used for ``--json`` modes)."""
    json.dump(payload, sys.stdout, indent=2, default=str)
    sys.stdout.write("\n")
