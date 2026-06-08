"""Identify CDR regions in antibody variable domain sequences.

Supports VH, VHH (nanobody), VL-kappa, and VL-lambda using conserved
framework anchor residues (Kabat/Chothia hybrid approach).

Usage:
    uv run scripts/find_cdrs.py --vh EVQLVES... --name VH_B
    uv run scripts/find_cdrs.py --vhh EVQLVES... --name VHH_A
    uv run scripts/find_cdrs.py --vl DIQMTQ... --name VL_C

# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""

import argparse
import re
import sys


def find_vh_cdrs(seq: str) -> tuple[str, str, str]:
    """Return (CDR1, CDR2, CDR3) for a VH or VHH sequence."""
    # FR1 ends at the conserved Cys. FR2 starts at the conserved Trp.
    cys1 = seq.index("C", 15)  # first conserved Cys ~position 22
    trp_fr2 = seq.index("W", cys1 + 1)  # conserved Trp starting FR2

    cdr1 = seq[cys1 + 4 : trp_fr2]

    # FR2 is 15 residues: W.RQ / W.LQ / W.VR / W.FR + 1 extra position
    fr2_end = trp_fr2 + 15

    # CDR2 ends just before the conserved RFTIS / RFTIA sequence
    m = re.search(r"[RK]FT[IV][SA]", seq[fr2_end:])
    if not m:
        raise ValueError("Could not find RFTIS anchor for CDR2 boundary")
    cdr2_end = fr2_end + m.start()
    cdr2 = seq[fr2_end:cdr2_end]

    # CDR3: between second conserved Cys (YYC) and WGQG / WGQGT / WGRQ
    cys2_match = re.search(r"YY[CA]", seq[cdr2_end:])
    if not cys2_match:
        raise ValueError("Could not find YYC anchor for CDR3 start")
    # CDR3 starts 1 residue after the conserved Cys (the residue right after C
    # is the last FR3 position in most Kabat-aligned sequences).
    cdr3_start = cdr2_end + cys2_match.end() + 1

    wg_match = re.search(r"WG[QR]G", seq[cdr3_start:])
    if not wg_match:
        raise ValueError("Could not find WGQG anchor for CDR3 end")
    cdr3_end = cdr3_start + wg_match.start()
    cdr3 = seq[cdr3_start:cdr3_end]

    return cdr1, cdr2, cdr3


def find_vl_cdrs(seq: str) -> tuple[str, str, str]:
    """Return (CDR1, CDR2, CDR3) for a VL-kappa or VL-lambda sequence."""
    # CDR1: between first conserved Cys and WYQ / WFQ / WYL
    cys1 = seq.index("C", 15)
    wyq_match = re.search(r"W[YF][QK]", seq[cys1 + 1:])
    if not wyq_match:
        raise ValueError("Could not find WYQ/WFQ anchor for VL CDR1 end")
    cdr1_end = cys1 + 1 + wyq_match.start()
    cdr1 = seq[cys1 + 1 : cdr1_end]

    # CDR2: exactly 7 residues after the conserved IY/LY/VY in FR2
    # FR2 contains KPG.SP / KPGQSP / KPGKAP, then CDR2 is 7 aa
    kpg_match = re.search(r"KPG[QK][SA][PL]", seq[cdr1_end:])
    if not kpg_match:
        raise ValueError("Could not find KPG anchor for VL CDR2 start")
    fr2_anchor_end = cdr1_end + kpg_match.end()
    # The conserved xY pair (IY, VY, LY) ends FR2; CDR2 starts right after.
    iy_match = re.search(r"[ILV]Y", seq[fr2_anchor_end:])
    if not iy_match:
        raise ValueError("Could not find xY anchor before CDR2")
    cdr2_start = fr2_anchor_end + iy_match.end()
    cdr2 = seq[cdr2_start : cdr2_start + 7]

    # CDR3: between second conserved Cys (YYC motif) and FGGGT
    cys2_match = re.search(r"YYC", seq[cdr2_start:])
    if not cys2_match:
        raise ValueError("Could not find YYC anchor for VL CDR3 start")
    cdr3_start = cdr2_start + cys2_match.end()
    fgggt = re.search(r"[FW]GGG[TK]", seq[cdr3_start:])
    if not fgggt:
        raise ValueError("Could not find FGGGT anchor for VL CDR3 end")
    cdr3 = seq[cdr3_start : cdr3_start + fgggt.start()]

    return cdr1, cdr2, cdr3


def main() -> None:
    parser = argparse.ArgumentParser(description="Identify CDRs in antibody variable domains")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--vh", metavar="SEQ", help="VH or VHH amino acid sequence")
    group.add_argument("--vhh", metavar="SEQ", help="VHH (nanobody) sequence (same logic as VH)")
    group.add_argument("--vl", metavar="SEQ", help="VL (kappa or lambda) amino acid sequence")
    parser.add_argument("--name", default="Query", help="Domain name for display")
    args = parser.parse_args()

    try:
        if args.vh or args.vhh:
            seq = (args.vh or args.vhh).strip().upper()
            cdr1, cdr2, cdr3 = find_vh_cdrs(seq)
            domain_type = "VHH" if args.vhh else "VH"
        else:
            seq = args.vl.strip().upper()
            cdr1, cdr2, cdr3 = find_vl_cdrs(seq)
            domain_type = "VL"

        print(f"{args.name} ({domain_type}), length {len(seq)} aa")
        print(f"  CDR1 ({len(cdr1):>2} aa): {cdr1}")
        print(f"  CDR2 ({len(cdr2):>2} aa): {cdr2}")
        print(f"  CDR3 ({len(cdr3):>2} aa): {cdr3}")

    except (ValueError, IndexError) as e:
        print(f"[!] CDR detection failed for {args.name}: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
