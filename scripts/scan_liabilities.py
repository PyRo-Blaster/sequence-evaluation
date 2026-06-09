"""Scan protein / antibody sequences for common developability liabilities.

Flags sequence-level chemical and post-translational liabilities that matter for
manufacturability and stability of therapeutic proteins: N-glycosylation
sequons, Asn deamidation, Asp isomerization, Met/Trp oxidation, free cysteines,
N-terminal pyroglutamate, and acid-labile Asp-Pro bonds.

Positions are 1-based on the supplied (mature) sequence.

Usage:
    uv run scripts/scan_liabilities.py <fasta_file> [--signal-peptide MGW...] [--json]
    uv run scripts/scan_liabilities.py --seq EVQLVES... --name VH_A [--json]

# /// script
# requires-python = ">=3.10"
# dependencies = ["biopython>=1.84"]
# ///
"""

import argparse
import re
import sys

from _common import SequenceError, clean_sequence, emit_json

# Each rule: (id, compiled regex, severity, human note). Regexes use a capture
# group or lookahead so overlapping motifs are still reported. The flagged span
# is taken from the start of the match.
RULES = [
    ("N-glycosylation sequon", re.compile(r"N[^P][ST]"), "high",
     "N-X-S/T (X!=P): consensus N-linked glycosylation site"),
    ("Deamidation (NG)", re.compile(r"N(?=G)"), "high",
     "Asn-Gly: fastest Asn deamidation motif"),
    ("Deamidation (NS/NT/NN/NH)", re.compile(r"N(?=[STNH])"), "medium",
     "Asn followed by S/T/N/H: moderate deamidation risk"),
    ("Isomerization (DG/DS/DD/DH/DT)", re.compile(r"D(?=[GSDHT])"), "medium",
     "Asp isomerization to iso-Asp"),
    ("Fragmentation (DP)", re.compile(r"D(?=P)"), "medium",
     "Asp-Pro: acid-labile peptide bond, low-pH fragmentation"),
    ("Met oxidation", re.compile(r"M"), "low",
     "Methionine: oxidation-prone (severity depends on solvent exposure)"),
    ("Trp oxidation", re.compile(r"W"), "low",
     "Tryptophan: oxidation-prone (severity depends on solvent exposure)"),
    ("Unpaired cysteine", re.compile(r"C"), "info",
     "Cysteine position (only flagged in summary when the count is odd)"),
]


def scan(seq: str) -> dict:
    findings: list[dict] = []
    for label, rx, severity, note in RULES:
        if label == "Unpaired cysteine":
            continue  # handled separately via parity
        for m in rx.finditer(seq):
            start = m.start()
            findings.append({
                "liability": label,
                "severity": severity,
                "position": start + 1,  # 1-based
                "motif": seq[start:start + 3] if label.startswith(("N-gly", "Deam", "Isom", "Frag")) else seq[start],
                "note": note,
            })

    cys_positions = [i + 1 for i, c in enumerate(seq) if c == "C"]
    cys_unpaired = len(cys_positions) % 2 == 1

    # N-terminal pyroglutamate from N-terminal Gln/Glu.
    nterm_pyroglu = seq[0] in "QE"

    findings.sort(key=lambda f: (f["position"], f["liability"]))
    return {
        "length": len(seq),
        "findings": findings,
        "cysteines": cys_positions,
        "cysteine_count": len(cys_positions),
        "odd_cysteine_count": cys_unpaired,
        "nterm_pyroglutamate": nterm_pyroglu,
        "nterm_residue": seq[0],
    }


SEVERITY_ORDER = {"high": 0, "medium": 1, "low": 2, "info": 3}


def report(name: str, result: dict) -> None:
    print(f"\n=== {name} ({result['length']} aa) ===")
    counts: dict[str, int] = {}
    for f in result["findings"]:
        counts[f["severity"]] = counts.get(f["severity"], 0) + 1
    summary = ", ".join(f"{counts[s]} {s}" for s in ("high", "medium", "low") if s in counts)
    print(f"Liabilities: {summary or 'none in scanned categories'}")

    if result["odd_cysteine_count"]:
        print(f"[!] ODD cysteine count ({result['cysteine_count']}) at positions "
              f"{result['cysteines']}: likely an unpaired/free thiol.")
    else:
        print(f"Cysteines: {result['cysteine_count']} (even — consistent with full pairing)")

    if result["nterm_pyroglutamate"]:
        print(f"[!] N-terminal {result['nterm_residue']}: prone to pyroglutamate formation.")

    if result["findings"]:
        print(f"\n  {'Pos':>5}  {'Sev':<6}  {'Motif':<6}  Liability")
        print(f"  {'-----':>5}  {'------':<6}  {'------':<6}  ---------")
        for f in sorted(result["findings"], key=lambda x: (SEVERITY_ORDER[x["severity"]], x["position"])):
            print(f"  {f['position']:>5}  {f['severity']:<6}  {f['motif']:<6}  {f['liability']}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Developability liability scan")
    parser.add_argument("fasta", nargs="?", help="Path to input FASTA file")
    parser.add_argument("--seq", help="Single amino-acid sequence (alternative to FASTA)")
    parser.add_argument("--name", default="Query", help="Name when using --seq")
    parser.add_argument("--signal-peptide", default="", help="Signal peptide to strip first")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    args = parser.parse_args()

    signal = clean_sequence(args.signal_peptide, name="signal peptide") if args.signal_peptide else ""

    entries: list[tuple[str, str]] = []
    try:
        if args.seq:
            entries.append((args.name, clean_sequence(args.seq, name=args.name)))
        elif args.fasta:
            from Bio import SeqIO
            records = list(SeqIO.parse(args.fasta, "fasta"))
            if not records:
                print(f"[!] No sequences found in {args.fasta}", file=sys.stderr)
                sys.exit(1)
            for rec in records:
                entries.append((rec.id, clean_sequence(str(rec.seq), name=rec.id)))
        else:
            parser.error("provide a FASTA path or --seq")
    except SequenceError as e:
        print(f"[!] {e}", file=sys.stderr)
        sys.exit(1)

    results = {}
    for name, seq in entries:
        if signal and seq.startswith(signal):
            seq = seq[len(signal):]
        results[name] = scan(seq)

    if args.json:
        emit_json(results)
        return
    for name, result in results.items():
        report(name, result)


if __name__ == "__main__":
    main()
