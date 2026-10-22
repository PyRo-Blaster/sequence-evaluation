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

from pathlib import Path

# uv does not put the script's directory on sys.path (uv >=0.13), so make the
# sibling _common module importable regardless of how this file is invoked.
sys.path.insert(0, str(Path(__file__).resolve().parent))

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

    # FR3 ends at the conserved Cys92 (Kabat). The residue before it is an
    # aromatic in essentially every human VH gene, but it is Tyr in most and Phe
    # in a substantial minority (YYC vs YFC) — anchoring on `YY[CA]` silently
    # failed on every YFC heavy chain. Anchor on the Cys and accept either.
    cys2_match = re.search(r"[YFH]C", seq[cdr2_end:])
    if not cys2_match:
        # Last resort: the first cysteine in FR3 is Cys92 in a normal domain.
        idx = seq.find("C", cdr2_end)
        if idx == -1:
            raise ValueError("could not find the conserved FR3 cysteine for CDR3 start")
        cys_end = idx + 1 - cdr2_end
    else:
        cys_end = cys2_match.end()
    # CDR-H3 begins three residues after that Cys (Cys, then the two FR3
    # residues such as A-R), so skip two past the end of the anchor match.
    cdr3_start = cdr2_end + cys_end + 2

    wg_match = re.search(r"WG[QR]G", seq[cdr3_start:])
    if not wg_match:
        raise ValueError("could not find WGQG anchor for CDR3 end")
    cdr3_end = cdr3_start + wg_match.start()
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

    # Same conserved-Cys logic as the heavy chain: accept YYC or YFC.
    cys2_match = re.search(r"[YFH]C", seq[cdr2_start:])
    if not cys2_match:
        raise ValueError("could not find the conserved FR3 cysteine for VL CDR3 start")
    cdr3_start = cdr2_start + cys2_match.end()
    fgggt = re.search(r"[FW]G.G[TK]", seq[cdr3_start:])
    if not fgggt:
        raise ValueError("could not find FGxGT anchor for VL CDR3 end")
    cdr3 = seq[cdr3_start : cdr3_start + fgggt.start()]
    return cdr1, cdr2, cdr3


def cdrs_via_abnumber(seq: str, scheme: str):
    """Try abnumber/ANARCI. Returns (cdrs, chain_type) or None if unavailable."""
    try:
        from abnumber import Chain  # type: ignore
    except Exception:
        return None
    try:
        chain = Chain(seq, scheme=scheme)
    except Exception as e:  # not a valid antibody domain for ANARCI
        raise ValueError(f"abnumber could not parse sequence: {e}")
    ctype = "VL" if chain.chain_type in ("K", "L") else "VH"
    return (chain.cdr1_seq, chain.cdr2_seq, chain.cdr3_seq), ctype


def annotate(seq: str, domain_type: str, scheme: str):
    """Return (cdrs, engine, resolved_type).

    The anchor-regex engine is a fallback for *both* abnumber being unavailable
    and abnumber failing to parse the domain. Previously a parse failure aborted
    the run, which made installing abnumber strictly worse for unusual domains.
    """
    degraded = ""
    try:
        ab = cdrs_via_abnumber(seq, scheme)
    except ValueError as e:
        ab = None
        degraded = f" [abnumber present but could not parse: {e}]"
    if ab is not None:
        cdrs, detected = ab
        return cdrs, f"abnumber/ANARCI ({scheme})", detected
    if domain_type == "VL":
        return find_vl_cdrs(seq), "anchor-regex" + degraded, "VL"
    return find_vh_cdrs(seq), "anchor-regex" + degraded, domain_type


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
        cdrs, engine, resolved = annotate(seq, "VL" if domain_type == "VL" else domain_type, args.scheme)
    except (ValueError, IndexError) as e:
        print(f"[!] CDR detection failed for {args.name}: {e}", file=sys.stderr)
        sys.exit(1)

    # Preserve VHH label if the user asserted it.
    display_type = "VHH" if domain_type == "VHH" else resolved
    cdr1, cdr2, cdr3 = cdrs

    if args.json:
        emit_json({
            "name": args.name, "domain_type": display_type, "engine": engine,
            "length": len(seq),
            "cdr1": cdr1, "cdr2": cdr2, "cdr3": cdr3,
            "cdr1_len": len(cdr1), "cdr2_len": len(cdr2), "cdr3_len": len(cdr3),
        })
        return

    print(f"{args.name} ({display_type}), length {len(seq)} aa  [engine: {engine}]")
    print(f"  CDR1 ({len(cdr1):>2} aa): {cdr1}")
    print(f"  CDR2 ({len(cdr2):>2} aa): {cdr2}")
    print(f"  CDR3 ({len(cdr3):>2} aa): {cdr3}")


if __name__ == "__main__":
    main()
