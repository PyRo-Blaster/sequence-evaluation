"""Calculate physical properties (MW, pI, A280, developability indices) for
protein chains from a FASTA file.

Usage:
    uv run scripts/analyze_properties.py <fasta_file> [--signal-peptide MGWSCI...]
                                                       [--complex]
                                                       [--disulfide-bonds N]
                                                       [--json]

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

from pathlib import Path

# uv does not put the script's directory on sys.path (uv >=0.13), so make the
# sibling _common module importable regardless of how this file is invoked.
sys.path.insert(0, str(Path(__file__).resolve().parent))

from _common import SequenceError, clean_sequence, emit_json

# Molar extinction coefficients at 280 nm (M^-1 cm^-1), Pace et al. 1995.
EXT_TRP = 5500
EXT_TYR = 1490
EXT_CYSTINE = 125  # per disulfide bond (pair of Cys)


class ChainProps(NamedTuple):
    name: str
    precursor_len: int
    mature_len: int
    mature_mw: float
    mature_pi: float
    ext_reduced: int
    ext_oxidized: int
    a280_reduced: float
    a280_oxidized: float
    instability: float
    aliphatic_index: float
    gravy: float
    n_cys: int


def aliphatic_index(seq: str) -> float:
    """Ikai (1980) aliphatic index: relative volume from Ala/Val/Ile/Leu."""
    n = len(seq)
    if n == 0:
        return 0.0
    mol_pct = {aa: 100.0 * seq.count(aa) / n for aa in "AVIL"}
    return mol_pct["A"] + 2.9 * mol_pct["V"] + 3.9 * (mol_pct["I"] + mol_pct["L"])


def calc_chain(name: str, seq: str, signal_peptide: str) -> ChainProps:
    precursor_len = len(seq)
    if signal_peptide and seq.startswith(signal_peptide):
        mature = seq[len(signal_peptide):]
    else:
        mature = seq

    pa = ProteinAnalysis(mature)
    mw = pa.molecular_weight()

    n_trp = mature.count("W")
    n_tyr = mature.count("Y")
    n_cys = mature.count("C")
    ext_reduced = n_trp * EXT_TRP + n_tyr * EXT_TYR
    # ExPASy convention: assume all cysteines pair into cystines.
    ext_oxidized = ext_reduced + (n_cys // 2) * EXT_CYSTINE

    return ChainProps(
        name=name,
        precursor_len=precursor_len,
        mature_len=len(mature),
        mature_mw=mw,
        mature_pi=pa.isoelectric_point(),
        ext_reduced=ext_reduced,
        ext_oxidized=ext_oxidized,
        a280_reduced=ext_reduced / mw if mw else 0.0,
        a280_oxidized=ext_oxidized / mw if mw else 0.0,
        instability=pa.instability_index(),
        aliphatic_index=aliphatic_index(mature),
        gravy=pa.gravy(),
        n_cys=n_cys,
    )


def complex_pi(mature_seqs: list[str], copies: list[int] | None = None) -> float:
    """Solve for the pH at which the assembled complex carries zero net charge.

    ``copies`` weights each chain by its stoichiometry; without it every chain
    counts once, which is wrong for any homomeric assembly (an IgG is 2 HC + 2 LC).
    """
    weights = copies or [1] * len(mature_seqs)

    def net_charge(ph: float) -> float:
        return sum(n * ProteinAnalysis(s).charge_at_pH(ph)
                   for s, n in zip(mature_seqs, weights))

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
    parser.add_argument(
        "--stoichiometry",
        nargs="*",
        default=[],
        metavar="NAME=N",
        help="Copies of each chain in the assembled complex, e.g. HC=2 LC=2. "
             "Without this every chain counts once, which halves the mass of an IgG.",
    )
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    args = parser.parse_args()

    signal = clean_sequence(args.signal_peptide, name="signal peptide") if args.signal_peptide else ""

    records = list(SeqIO.parse(args.fasta, "fasta"))
    if not records:
        print(f"[!] No sequences found in {args.fasta}", file=sys.stderr)
        sys.exit(1)

    chains: list[ChainProps] = []
    mature_seqs: list[str] = []

    for rec in records:
        try:
            seq = clean_sequence(str(rec.seq), name=rec.id)
        except SequenceError as e:
            print(f"[!] {e}", file=sys.stderr)
            sys.exit(1)
        props = calc_chain(rec.id, seq, signal)
        chains.append(props)
        mature_seqs.append(seq[len(signal):] if signal and seq.startswith(signal) else seq)

    stoich: dict[str, int] = {}
    for item in args.stoichiometry:
        if "=" not in item:
            parser.error(f"--stoichiometry expects NAME=N, got {item!r}")
        key, _, val = item.partition("=")
        try:
            stoich[key] = int(val)
        except ValueError:
            parser.error(f"--stoichiometry copy number must be an integer, got {val!r}")
    unknown = [k for k in stoich if k not in {c.name for c in chains}]
    if unknown:
        print(f"[!] --stoichiometry names not in the FASTA: {', '.join(unknown)}", file=sys.stderr)
        sys.exit(1)

    complex_data = None
    if args.complex and len(chains) > 1:
        copies = [stoich.get(p.name, 1) for p in chains]
        total_mw = sum(n * p.mature_mw for p, n in zip(chains, copies))
        ds_correction = args.disulfide_bonds * 2.016
        complex_data = {
            "total_mw": total_mw,
            "stoichiometry": ({p.name: n for p, n in zip(chains, copies)} if stoich else None),
            "disulfide_bonds": args.disulfide_bonds,
            "disulfide_correction": -ds_correction,
            "corrected_mw": total_mw - ds_correction,
            "pi": complex_pi(mature_seqs, copies),
        }

    if args.json:
        emit_json({
            "chains": [c._asdict() for c in chains],
            "complex": complex_data,
        })
        return

    print(f"{'Chain':<24} {'Prec':>5} {'Mat':>5} {'MW (Da)':>13} {'pI':>5} "
          f"{'A280(ox)':>9} {'Instab':>7} {'Aliph':>6} {'GRAVY':>7}")
    print("-" * 92)
    for p in chains:
        print(f"{p.name:<24} {p.precursor_len:>5} {p.mature_len:>5} "
              f"{p.mature_mw:>13,.2f} {p.mature_pi:>5.2f} {p.a280_oxidized:>9.3f} "
              f"{p.instability:>7.1f} {p.aliphatic_index:>6.1f} {p.gravy:>7.3f}")
    print()
    print("A280(ox) = absorbance of a 1 g/L solution at 280 nm assuming all Cys "
          "form cystines (ExPASy convention).")
    print("Instability >40 suggests the protein may be unstable in a test tube.")

    if complex_data:
        print()
        print("=== Assembled Complex ===")
        if complex_data.get("stoichiometry"):
            print("Stoichiometry:           "
                  + ", ".join(f"{k}x{v}" for k, v in complex_data["stoichiometry"].items()))
        else:
            print("Stoichiometry:           1 copy of each chain (assumed; pass --stoichiometry)")
        print(f"Total MW (sum):          {complex_data['total_mw']:>14,.2f} Da")
        if args.disulfide_bonds:
            print(f"Disulfide correction:    {complex_data['disulfide_correction']:>14,.2f} Da "
                  f"({args.disulfide_bonds} bonds)")
            print(f"Corrected MW:            {complex_data['corrected_mw']:>14,.2f} Da")
        print(f"Complex pI:              {complex_data['pi']:>14.2f}")


if __name__ == "__main__":
    main()
