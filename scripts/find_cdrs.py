"""Identify CDR regions in antibody variable domain sequences.

Two engines:

1. **abnumber/ANARCI** (preferred) — true Kabat/IMGT/Chothia numbering. Used
   automatically when importable. Install on the fly with:
       uv run --with abnumber scripts/find_cdrs.py ...
   (ANARCI also requires HMMER; abnumber bundles a usable build on most
   platforms.)
2. **Anchor regex** (fallback) — conserved framework motifs, Kabat-compatible
   boundaries, no external dependency. CDR3 N-terminal boundary may be +/-1 for
   some sequences.

Supports VH, VHH (nanobody), VL-kappa and VL-lambda. Chain type is
auto-detected from FR1 motifs when not given.

Usage:
    uv run scripts/find_cdrs.py --vh  EVQLVES... --name VH_B
    uv run scripts/find_cdrs.py --vhh EVQLVES... --name VHH_A
    uv run scripts/find_cdrs.py --vl  DIQMTQ...  --name VL_C
    uv run scripts/find_cdrs.py --seq EVQLVES... --name X --json   # auto-detect

# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""

import argparse
import re
import sys

from _common import SequenceError, clean_sequence, emit_json

VH_FR1 = ("EVQLVE", "QVQLVE", "QVQLQE", "QVQLQQ", "EVQLLE", "QVQLVQ")
VL_FR1 = ("DIQMTQ", "EIVLTQ", "DIVMTQ", "DIVLTQ", "EIVMTQ", "QSVLTQ", "QSALTQ", "SYELTQ", "DIQLTQ")


def detect_chain_type(seq: str) -> str:
    """Return 'VH' or 'VL' from FR1 motifs (VHH reported as VH)."""
    head = seq[:8]
    if any(head.startswith(m) for m in VL_FR1):
        return "VL"
    if any(head.startswith(m) for m in VH_FR1):
        return "VH"
    # Heuristic fallback: kappa/lambda usually carry the FGGGT/FGQGT J motif.
    if re.search(r"[FW]G[GQ]G[TK]", seq[-20:]):
        return "VL"
    return "VH"


def find_vh_cdrs(seq: str) -> tuple[str, str, str]:
    """Return (CDR1, CDR2, CDR3) for a VH or VHH sequence (anchor method)."""
    cys1 = seq.index("C", 15)
    trp_fr2 = seq.index("W", cys1 + 1)
    cdr1 = seq[cys1 + 4 : trp_fr2]

    fr2_end = trp_fr2 + 15
    m = re.search(r"[RK]FT[IV][SA]", seq[fr2_end:])
    if not m:
        raise ValueError("could not find RFTIS anchor for CDR2 boundary")
    cdr2_end = fr2_end + m.start()
    cdr2 = seq[fr2_end:cdr2_end]

    # CDR-H3: bounded by the conserved Cys92 and the J-region WGxGT motif.
    # Anchor on the Cys itself rather than a literal YYC, so FR3 variants such
    # as VYFCAR / VYHCAR / FYCAR are handled. Kabat places two framework
    # residues (93-94, typically AR/AK/TR/SR) between Cys92 and CDR-H3 start.
    j_matches = list(re.finditer(r"WG.G[TAS]", seq))
    if not j_matches:
        raise ValueError("could not find WGxGT J-region anchor for CDR3 end")
    cdr3_end = j_matches[-1].start()

    cys92 = seq.rfind("C", cdr2_end, cdr3_end)
    if cys92 == -1:
        raise ValueError("could not find conserved Cys92 before the J-region")
    cdr3_start = cys92 + 3
    if cdr3_start > cdr3_end:
        raise ValueError("implausible CDR3 boundaries (Cys92 too close to J-region)")
    cdr3 = seq[cdr3_start:cdr3_end]
    return cdr1, cdr2, cdr3


def find_vl_cdrs(seq: str) -> tuple[str, str, str]:
    """Return (CDR1, CDR2, CDR3) for a VL-kappa or VL-lambda sequence (anchor method)."""
    cys1 = seq.index("C", 15)
    wyq_match = re.search(r"W[YF][QK]", seq[cys1 + 1:])
    if not wyq_match:
        raise ValueError("could not find WYQ/WFQ anchor for VL CDR1 end")
    cdr1_end = cys1 + 1 + wyq_match.start()
    cdr1 = seq[cys1 + 1 : cdr1_end]

    kpg_match = re.search(r"KPG[QK][SA][PL]", seq[cdr1_end:])
    if not kpg_match:
        raise ValueError("could not find KPG anchor for VL CDR2 start")
    fr2_anchor_end = cdr1_end + kpg_match.end()
    iy_match = re.search(r"[ILV]Y", seq[fr2_anchor_end:])
    if not iy_match:
        raise ValueError("could not find xY anchor before CDR2")
    cdr2_start = fr2_anchor_end + iy_match.end()
    cdr2 = seq[cdr2_start : cdr2_start + 7]

    # CDR-L3: bounded by the conserved Cys88 and the J-region FGxGT motif.
    # Anchored on the Cys itself (not a literal YYC) so YFC/FYC variants work.
    # Kabat CDR-L3 starts immediately after Cys88 (no intervening framework).
    j_matches = list(re.finditer(r"[FW]G.G[TKQ]", seq))
    if not j_matches:
        raise ValueError("could not find FGxGT J-region anchor for VL CDR3 end")
    cdr3_end = j_matches[-1].start()

    cys88 = seq.rfind("C", cdr2_start, cdr3_end)
    if cys88 == -1:
        raise ValueError("could not find conserved Cys88 before the J-region")
    cdr3_start = cys88 + 1
    if cdr3_start > cdr3_end:
        raise ValueError("implausible CDR3 boundaries (Cys88 too close to J-region)")
    cdr3 = seq[cdr3_start:cdr3_end]
    return cdr1, cdr2, cdr3


def cdrs_via_abnumber(seq: str, scheme: str):
    """Try abnumber/ANARCI.

    Returns ((cdr1, cdr2, cdr3), chain_type) on success, or (None, reason) when
    the engine is unusable so the caller can fall back to the regex anchors.

    NOTE: importing `abnumber` is NOT sufficient. abnumber is a thin wrapper and
    calls `anarci.anarci` lazily at Chain() time; the original ANARCI is not on
    PyPI (conda-only, needs HMMER), so `pip install abnumber` alone yields a
    ModuleNotFoundError at first use. That must degrade to the fallback, not
    abort the run.
    """
    try:
        from abnumber import Chain  # type: ignore
    except Exception:
        return None, "abnumber not installed"
    try:
        chain = Chain(seq, scheme=scheme)
    except (ImportError, ModuleNotFoundError) as e:
        # abnumber present but its ANARCI backend is missing.
        return None, (f"abnumber installed but backend unavailable ({e}); "
                      "install ANARCI via conda: conda install -c bioconda anarci")
    except Exception as e:
        # ANARCI ran but could not number this sequence; the regex anchors may
        # still succeed, so fall back rather than abort.
        return None, f"abnumber could not number the sequence ({e})"
    ctype = "VL" if chain.chain_type in ("K", "L") else "VH"
    return (chain.cdr1_seq, chain.cdr2_seq, chain.cdr3_seq), ctype


def annotate(seq: str, domain_type: str, scheme: str):
    """Return (cdrs, engine, resolved_type, note)."""
    cdrs, info = cdrs_via_abnumber(seq, scheme)
    if cdrs is not None:
        return cdrs, f"abnumber/ANARCI ({scheme})", info, ""
    note = "" if info == "abnumber not installed" else f"ANARCI unavailable: {info}"
    if domain_type == "VL":
        return find_vl_cdrs(seq), "anchor-regex", "VL", note
    return find_vh_cdrs(seq), "anchor-regex", domain_type, note


def main() -> None:
    parser = argparse.ArgumentParser(description="Identify CDRs in antibody variable domains")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--vh", metavar="SEQ", help="VH amino acid sequence")
    group.add_argument("--vhh", metavar="SEQ", help="VHH (nanobody) sequence")
    group.add_argument("--vl", metavar="SEQ", help="VL (kappa or lambda) sequence")
    group.add_argument("--seq", metavar="SEQ", help="Sequence with auto chain-type detection")
    parser.add_argument("--name", default="Query", help="Domain name for display")
    parser.add_argument("--scheme", default="kabat", choices=["kabat", "imgt", "chothia"],
                        help="Numbering scheme for the abnumber engine (default: kabat)")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    args = parser.parse_args()

    raw = args.vh or args.vhh or args.vl or args.seq
    try:
        seq = clean_sequence(raw, name=args.name)
    except SequenceError as e:
        print(f"[!] {e}", file=sys.stderr)
        sys.exit(1)

    if args.vhh:
        domain_type = "VHH"
    elif args.vh:
        domain_type = "VH"
    elif args.vl:
        domain_type = "VL"
    else:
        domain_type = detect_chain_type(seq)

    try:
        cdrs, engine, resolved, note = annotate(
            seq, "VL" if domain_type == "VL" else domain_type, args.scheme)
    except (ValueError, IndexError) as e:
        print(f"[!] CDR detection failed for {args.name}: {e}", file=sys.stderr)
        sys.exit(1)

    # Preserve VHH label if the user asserted it.
    display_type = "VHH" if domain_type == "VHH" else resolved
    cdr1, cdr2, cdr3 = cdrs

    if args.json:
        emit_json({
            "name": args.name, "domain_type": display_type, "engine": engine,
            "engine_note": note, "length": len(seq),
            "cdr1": cdr1, "cdr2": cdr2, "cdr3": cdr3,
            "cdr1_len": len(cdr1), "cdr2_len": len(cdr2), "cdr3_len": len(cdr3),
        })
        return

    print(f"{args.name} ({display_type}), length {len(seq)} aa  [engine: {engine}]")
    if note:
        print(f"  note: {note}")
    print(f"  CDR1 ({len(cdr1):>2} aa): {cdr1}")
    print(f"  CDR2 ({len(cdr2):>2} aa): {cdr2}")
    print(f"  CDR3 ({len(cdr3):>2} aa): {cdr3}")


if __name__ == "__main__":
    main()
