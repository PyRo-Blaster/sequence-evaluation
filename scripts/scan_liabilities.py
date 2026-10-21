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
import json
import re
import sys

from pathlib import Path

# uv does not put the script's directory on sys.path (uv >=0.13), so make the
# sibling _common module importable regardless of how this file is invoked.
sys.path.insert(0, str(Path(__file__).resolve().parent))

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


# Motif span actually described by each rule, for the reported `motif` string.
# N-X-S/T is three residues; NG, N[STNH], D[GSDHT] and DP are two; the oxidation
# rules are single residues.
MOTIF_SPAN = {
    "N-glycosylation sequon": 3,
    "Deamidation (NG)": 2,
    "Deamidation (NS/NT/NN/NH)": 2,
    "Isomerization (DG/DS/DD/DH/DT)": 2,
    "Fragmentation (DP)": 2,
}

SEVERITY_SCALE = ["low", "medium", "high"]
OXIDATION = ("Met oxidation", "Trp oxidation")


def _escalate(severity: str, steps: int = 1) -> str:
    """Move a severity up the scale, saturating at 'high'."""
    if severity not in SEVERITY_SCALE:
        return severity
    i = SEVERITY_SCALE.index(severity)
    return SEVERITY_SCALE[min(i + steps, len(SEVERITY_SCALE) - 1)]


def locate_cdrs(seq: str, cdrs: list[dict]) -> dict[int, str]:
    """Map 1-based sequence position -> CDR name by locating each CDR substring."""
    pos_to_cdr: dict[int, str] = {}
    # CDRs occur in order along the chain, so advance a cursor rather than
    # restarting find() from 0: a short CDR whose sequence also occurs earlier
    # (or is repeated) would otherwise be mapped onto the wrong residues.
    cursor = 0
    for cdr in cdrs or []:
        cdr_seq = (cdr.get("seq") or "").upper()
        name = cdr.get("name") or "CDR"
        if not cdr_seq:
            continue
        i = seq.find(cdr_seq, cursor)
        if i == -1:
            # Out of order or absent: fall back to a global search so a
            # mis-ordered map still maps something instead of dropping it.
            i = seq.find(cdr_seq)
            if i == -1:
                continue
        for offset in range(len(cdr_seq)):
            pos_to_cdr[i + offset + 1] = name
        cursor = i + len(cdr_seq)
    return pos_to_cdr


# The conserved IgG Fc N-glycosylation sequon (EU N297) sits in EEQ[YF]NST.
# It is expected and must not be engineered out unless aglycosylation is the
# goal, so it is reported separately from unwanted (e.g. Fab) sequons and is
# never escalated by CDR proximity.
FC_N297_CONTEXT = re.compile(r"Q[YF]NST")


def scan(seq: str, cdrs: list[dict] | None = None,
         contact_cdrs: set[str] | None = None) -> dict:
    """Scan one chain. ``cdrs`` escalates CDR-localised hits; ``contact_cdrs``
    names the CDRs carrying a paratope contact (CDR-level, not residue-level)."""
    pos_to_cdr = locate_cdrs(seq, cdrs) if cdrs else {}
    contact_cdrs = contact_cdrs or set()

    findings: list[dict] = []
    for label, rx, severity, note in RULES:
        if label == "Unpaired cysteine":
            continue  # handled separately via parity
        span = MOTIF_SPAN.get(label, 1)
        for m in rx.finditer(seq):
            start = m.start()
            position = start + 1  # 1-based
            cdr = pos_to_cdr.get(position)
            final = severity
            is_n297 = (label == "N-glycosylation sequon"
                       and FC_N297_CONTEXT.match(seq, max(start - 2, 0)) is not None)
            if is_n297:
                # Conserved Fc glycan: expected, not a liability to remove, and
                # exempt from CDR escalation.
                final = "expected"
                note = ("conserved Fc N-glycosylation site (EU N297) — expected; "
                        "removing it aglycosylates the Fc and ablates FcgR binding")
            elif cdr:
                # references/developability.md: CDR-localised liabilities are the
                # engineering priority, so they move one step up the scale.
                final = _escalate(severity)
                if label in OXIDATION and cdr in contact_cdrs:
                    # An oxidation hotspot on a contact-bearing CDR gets a
                    # second step: oxidation there can erode potency directly.
                    final = _escalate(final)
            findings.append({
                "liability": label,
                "severity": final,
                "base_severity": severity,
                "escalated": final != severity,
                "position": position,
                "in_cdr": cdr is not None,
                "cdr": cdr,
                "cdr_has_contact": bool(cdr and cdr in contact_cdrs),
                "fc_n297": is_n297,
                "motif": seq[start:start + span],
                "note": note,
            })

    cys_positions = [i + 1 for i, c in enumerate(seq) if c == "C"]
    # NOTE: an odd per-chain cysteine count is NORMAL for antibody chains,
    # because interchain disulfides pair with cysteines on *other* chains
    # (H-L and H-H). Parity is therefore only meaningful across the assembled
    # molecule; see assess_cysteines() for the whole-input judgement.
    cys_odd = len(cys_positions) % 2 == 1

    # N-terminal pyroglutamate from N-terminal Gln/Glu.
    nterm_pyroglu = seq[0] in "QE"

    findings.sort(key=lambda f: (f["position"], f["liability"]))
    return {
        "length": len(seq),
        "cdr_annotated": bool(pos_to_cdr),
        "cdrs_located": sorted(set(pos_to_cdr.values())),
        "contact_cdrs": sorted(contact_cdrs),
        "findings": findings,
        "cysteines": cys_positions,
        "cysteine_count": len(cys_positions),
        "odd_cysteine_count": cys_odd,
        "nterm_pyroglutamate": nterm_pyroglu,
        "nterm_residue": seq[0],
    }


SEVERITY_ORDER = {"high": 0, "medium": 1, "low": 2, "expected": 3, "info": 4}



def report(name: str, result: dict) -> None:
    print(f"\n=== {name} ({result['length']} aa) ===")
    counts: dict[str, int] = {}
    for f in result["findings"]:
        counts[f["severity"]] = counts.get(f["severity"], 0) + 1
    summary = ", ".join(f"{counts[s]} {s}" for s in ("high", "medium", "low") if s in counts)
    print(f"Liabilities: {summary or 'none in scanned categories'}"
          + (f"  (+{counts['expected']} expected Fc glycan)" if "expected" in counts else ""))

    parity = "odd" if result["odd_cysteine_count"] else "even"
    print(f"Cysteines: {result['cysteine_count']} ({parity} for this chain; "
          "interchain partners may be on other chains)")

    if result["nterm_pyroglutamate"]:
        print(f"[!] N-terminal {result['nterm_residue']}: prone to pyroglutamate formation.")

    if result.get("cdr_annotated"):
        located = ", ".join(result.get("cdrs_located") or [])
        print(f"CDRs located in this chain: {located or 'none'}")
        if result.get("contact_cdrs"):
            print(f"Contact-bearing CDRs (from paratope): {', '.join(result['contact_cdrs'])}")
    else:
        print("[note] no CDR map applied — severities are unescalated defaults. "
              "Pass --cdr-json to prioritise CDR-localised liabilities.")

    if result["findings"]:
        print(f"\n  {'Pos':>5}  {'Sev':<6}  {'Motif':<6}  {'CDR':<8}  Liability")
        print(f"  {'-----':>5}  {'------':<6}  {'------':<6}  {'--------':<8}  ---------")
        for f in sorted(result["findings"], key=lambda x: (SEVERITY_ORDER[x["severity"]], x["position"])):
            sev = f["severity"] + ("*" if f.get("escalated") else "")
            print(f"  {f['position']:>5}  {sev:<6}  {f['motif']:<6}  "
                  f"{(f.get('cdr') or '-'):<8}  {f['liability']}")
        if any(f.get("escalated") for f in result["findings"]):
            print("  * escalated from the default severity because the residue is in a CDR.")


def main() -> None:
    parser = argparse.ArgumentParser(description="Developability liability scan")
    parser.add_argument("fasta", nargs="?", help="Path to input FASTA file")
    parser.add_argument("--seq", help="Single amino-acid sequence (alternative to FASTA)")
    parser.add_argument("--name", default="Query", help="Name when using --seq")
    parser.add_argument("--signal-peptide", default="", help="Signal peptide to strip first")
    parser.add_argument("--cdr-json", help=(
        "JSON of {name: [{name, seq}]} from build_cdr_map.py. Keys must match the "
        "FASTA record ids (or --name). CDR-localised liabilities are escalated one "
        "severity step."))
    parser.add_argument("--paratope-json", nargs="*", default=[], metavar="FILE", help=(
        "analyze_interfaces.py --json output(s). Used at CDR level only: CDRs "
        "carrying a paratope contact give oxidation hotspots a second escalation "
        "step. Residue-level mapping is not attempted because structure residue "
        "numbers and 1-based mature-sequence positions are different coordinate "
        "systems."))
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

    cdr_map: dict[str, list[dict]] = {}
    if args.cdr_json:
        try:
            cdr_map = json.loads(Path(args.cdr_json).read_text())
        except Exception as e:
            print(f"[!] could not read {args.cdr_json}: {e}", file=sys.stderr)
            sys.exit(1)

    # Paratope JSONs contribute the set of CDR names that carry a contact.
    contact_cdrs: set[str] = set()
    for path in args.paratope_json:
        try:
            iface = json.loads(Path(path).read_text())
        except Exception as e:
            print(f"[!] could not read {path}: {e}", file=sys.stderr)
            sys.exit(1)
        for res in iface.get("paratope", []):
            region = res.get("region")
            if region and region not in ("FR", "unassigned"):
                contact_cdrs.add(region)

    results = {}
    for name, seq in entries:
        if signal and seq.startswith(signal):
            seq = seq[len(signal):]
        results[name] = scan(seq, cdr_map.get(name), contact_cdrs)

    # Per-chain cysteine parity only means "free thiol" for a single-chain
    # protein. In an assembled antibody the inter-chain bonds (HC-LC, and the
    # hinge HC-HC pair) consume cysteines across chains, so a normal IgG heavy
    # chain has an odd count by construction. Record the assembly context so a
    # report does not read that as an unpaired thiol.
    multi = len(results) > 1
    assembly_total = sum(res["cysteine_count"] for res in results.values())
    for res in results.values():
        res["multi_chain_input"] = multi
        res["assembly_cysteine_total"] = assembly_total
        res["assembly_cysteine_parity"] = "even" if assembly_total % 2 == 0 else "odd"

    unmatched = [k for k in cdr_map if k not in results]
    if unmatched:
        print(f"[!] --cdr-json key(s) {', '.join(unmatched)} match no scanned chain "
              f"(scanned: {', '.join(results)}); those CDRs were not applied.",
              file=sys.stderr)

    if args.json:
        emit_json(results)
        return
    for name, result in results.items():
        report(name, result)


if __name__ == "__main__":
    main()
