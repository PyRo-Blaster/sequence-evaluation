"""Calculate physical properties (MW, pI) for protein chains from a FASTA file.

Usage:
    uv run scripts/analyze_properties.py <fasta_file> [--signal-peptide MGWSCI...]
                                                       [--complex]
                                                       [--disulfide-bonds N]

# /// script
# requires-python = ">=3.10"
# dependencies = ["biopython>=1.84"]
# ///
"""

import argparse
import sys
from typing import NamedTuple

from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis


class ChainProps(NamedTuple):
    name: str
    precursor_len: int
    mature_len: int
    mature_mw: float
    mature_pi: float


def calc_chain(name: str, seq: str, signal_peptide: str) -> ChainProps:
    precursor_len = len(seq)
    if signal_peptide and seq.startswith(signal_peptide):
        mature = seq[len(signal_peptide):]
    else:
        mature = seq
    pa = ProteinAnalysis(mature)
    return ChainProps(
        name=name,
        precursor_len=precursor_len,
        mature_len=len(mature),
        mature_mw=pa.molecular_weight(),
        mature_pi=pa.isoelectric_point(),
    )


def complex_pi(mature_seqs: list[str]) -> float:
    def net_charge(ph: float) -> float:
        return sum(ProteinAnalysis(s).charge_at_pH(ph) for s in mature_seqs)

    lo, hi = 0.0, 14.0
    for _ in range(100):
        mid = (lo + hi) / 2.0
        if net_charge(mid) > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def main() -> None:
    parser = argparse.ArgumentParser(description="Protein chain physical properties")
    parser.add_argument("fasta", help="Path to input FASTA file")
    parser.add_argument(
        "--signal-peptide",
        default="",
        help="Signal peptide sequence to strip (default: none)",
    )
    parser.add_argument(
        "--complex",
        action="store_true",
        help="Calculate assembled complex MW and pI",
    )
    parser.add_argument(
        "--disulfide-bonds",
        type=int,
        default=0,
        help="Number of disulfide bonds for MW correction (-2.016 Da each)",
    )
    args = parser.parse_args()

    records = list(SeqIO.parse(args.fasta, "fasta"))
    if not records:
        print(f"[!] No sequences found in {args.fasta}", file=sys.stderr)
        sys.exit(1)

    print(f"{'Chain':<30} {'Precursor':>10} {'Mature':>8} {'MW (Da)':>14} {'pI':>6}")
    print("-" * 72)

    chains: list[ChainProps] = []
    mature_seqs: list[str] = []

    for rec in records:
        seq = str(rec.seq)
        props = calc_chain(rec.id, seq, args.signal_peptide)
        chains.append(props)
        mature_seqs.append(seq[len(args.signal_peptide):] if args.signal_peptide and seq.startswith(args.signal_peptide) else seq)
        print(f"{props.name:<30} {props.precursor_len:>10} {props.mature_len:>8} {props.mature_mw:>14,.2f} {props.mature_pi:>6.2f}")

    if args.complex and len(chains) > 1:
        total_mw = sum(p.mature_mw for p in chains)
        ds_correction = args.disulfide_bonds * 2.016
        corrected_mw = total_mw - ds_correction
        pi = complex_pi(mature_seqs)

        print()
        print("=== Assembled Complex ===")
        print(f"Total MW (sum):          {total_mw:>14,.2f} Da")
        if args.disulfide_bonds:
            print(f"Disulfide correction:    {-ds_correction:>14,.2f} Da ({args.disulfide_bonds} bonds)")
            print(f"Corrected MW:            {corrected_mw:>14,.2f} Da")
        print(f"Complex pI:              {pi:>14.2f}")


if __name__ == "__main__":
    main()
