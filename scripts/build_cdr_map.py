"""Build the ``{chain: [{name, seq}]}`` CDR map consumed by the other scripts.

``find_cdrs.py --json`` emits one document per variable domain
(``{name, domain_type, cdr1, cdr2, cdr3, ...}``), but ``analyze_interfaces.py``,
``scan_surface_patches.py`` and ``scan_liabilities.py`` all consume a map keyed
by chain. Converting by hand is error-prone and a wrong chain ID silently
produces "no paratope residue mapped to a CDR", so do it here instead.

CDR names are derived from the domain type: VH/VHH -> CDR-H1/2/3,
VL -> CDR-L1/2/3. Use ``--prefix`` to override per entry.

Usage:
    uv run scripts/build_cdr_map.py A=vl.json B=vh.json > cdrs.json
    uv run scripts/build_cdr_map.py A=vhh.json --prefix A=N   # CDR-N1/2/3

# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from _common import emit_json

LETTER = {"VH": "H", "VHH": "H", "VL": "L"}


def main() -> None:
    p = argparse.ArgumentParser(
        description="Assemble a {chain: [{name, seq}]} CDR map from find_cdrs.py JSON"
    )
    p.add_argument("pairs", nargs="+", metavar="CHAIN=FILE",
                   help="Chain ID and the find_cdrs.py --json output for the domain on it")
    p.add_argument("--prefix", action="append", default=[], metavar="CHAIN=LETTER",
                   help="Override the CDR letter for a chain (default from domain type)")
    args = p.parse_args()

    overrides = {}
    for item in args.prefix:
        if "=" not in item:
            p.error(f"--prefix expects CHAIN=LETTER, got {item!r}")
        cid, letter = item.split("=", 1)
        overrides[cid] = letter

    out: dict[str, list[dict]] = {}
    for item in args.pairs:
        if "=" not in item:
            p.error(f"expected CHAIN=FILE, got {item!r}")
        cid, path = item.split("=", 1)
        try:
            d = json.loads(Path(path).read_text())
        except Exception as e:
            print(f"[!] could not read {path}: {e}", file=sys.stderr)
            sys.exit(1)
        missing = [k for k in ("cdr1", "cdr2", "cdr3") if k not in d]
        if missing:
            print(f"[!] {path} is not find_cdrs.py --json output (missing {', '.join(missing)})",
                  file=sys.stderr)
            sys.exit(1)
        letter = overrides.get(cid) or LETTER.get(d.get("domain_type", ""), "X")
        entries = []
        for i in (1, 2, 3):
            seq = d[f"cdr{i}"]
            if not seq:
                print(f"[!] {path}: CDR{i} is empty; skipping", file=sys.stderr)
                continue
            entries.append({"name": f"CDR-{letter}{i}", "seq": seq})
        if cid in out:
            print(f"[!] chain {cid} given more than once", file=sys.stderr)
            sys.exit(1)
        out[cid] = entries

    emit_json(out)


if __name__ == "__main__":
    main()
