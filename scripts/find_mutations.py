"""Map mutations in an IgG heavy-chain constant region vs WT IGHG1*01 (EU numbering).

The constant region is located by the conserved ASTKGPSVF (CH1 start) anchor and
then *globally aligned* to the WT reference, so insertions/deletions in the query
no longer shift every downstream EU number (the failure mode of naive positional
comparison). Recognised engineering mutations are annotated automatically, and a
high divergence count warns when the chain is likely a non-IgG1 isotype/allotype.

Usage:
    uv run scripts/find_mutations.py <full_heavy_chain_sequence> [--label HC1] [--json]

# /// script
# requires-python = ">=3.10"
# dependencies = ["biopython>=1.84"]
# ///
"""

import argparse
import sys

from Bio.Align import PairwiseAligner

from _common import SequenceError, clean_sequence, emit_json

# WT human IgG1 constant region (IGHG1*01): CH1 + Hinge + CH2 + CH3
WT_IgG1_CONST = (
    "ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSS"
    "VVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTL"
    "MISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEY"
    "KCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQP"
    "ENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVLHEALHNHYTQKSLSLSPGK"
)

# EU numbering: CH1 118-215, Hinge 216-230, CH2 231-340, CH3 341-447
EU_NUMBERS: list[int] = (
    list(range(118, 216))   # CH1
    + list(range(216, 231)) # Hinge
    + list(range(231, 341)) # CH2
    + list(range(341, 448)) # CH3
)
assert len(EU_NUMBERS) == len(WT_IgG1_CONST), "EU numbering must match reference length"


def domain_for(eu: int) -> str:
    if 118 <= eu <= 215:
        return "CH1"
    if 216 <= eu <= 230:
        return "Hinge"
    if 231 <= eu <= 340:
        return "CH2"
    return "CH3"


# Per-mutation functional annotation (EU numbering).
MUTATION_NOTES = {
    "L234A": "LALA effector silencing",
    "L235A": "LALA effector silencing",
    "P329G": "completes LALA-PG (abolishes C1q/CDC)",
    "L234F": "FLES silencing",
    "L235E": "FLES silencing",
    "P331S": "FLES silencing",
    "M428L": "LS half-life extension",
    "N434S": "LS / partial half-life extension",
    "T250Q": "QM half-life extension",
    "M252Y": "YTE half-life extension",
    "S254T": "YTE half-life extension",
    "T256E": "YTE half-life extension",
    "T366W": "Knob (knobs-into-holes)",
    "T366S": "Hole (knobs-into-holes)",
    "L368A": "Hole (knobs-into-holes)",
    "Y407V": "Hole (knobs-into-holes)",
    "G236A": "GASDALIE enhanced effector",
    "S239D": "GASDALIE enhanced effector",
    "A330L": "GASDALIE enhanced effector",
    "I332E": "GASDALIE enhanced effector",
    "S267E": "SELF enhanced binding",
    "L328F": "SELF enhanced binding",
    "D356E": "EEM allotype",
    "L358M": "EEM allotype",
}

# Named variants: report present / partial based on the full mutation set.
NAMED_VARIANTS = {
    "LALA": {"L234A", "L235A"},
    "LALA-PG": {"L234A", "L235A", "P329G"},
    "YTE": {"M252Y", "S254T", "T256E"},
    "LS": {"M428L", "N434S"},
    "KiH knob": {"T366W"},
    "KiH hole": {"T366S", "L368A", "Y407V"},
    "GASDALIE": {"G236A", "S239D", "A330L", "I332E"},
    "SELF": {"S267E", "L328F"},
}


def _aligner() -> PairwiseAligner:
    a = PairwiseAligner()
    a.mode = "global"
    a.open_gap_score = -10
    a.extend_gap_score = -0.5
    a.match_score = 2
    a.mismatch_score = -1
    # Don't penalise gaps at either terminus (e.g. clipped K447 / leading slop).
    a.end_gap_score = 0
    return a


def find_mutations(chain_seq: str, label: str) -> dict:
    anchor = "ASTKGPSVF"
    pos = chain_seq.find(anchor)
    if pos == -1:
        raise SequenceError(
            f"{label}: could not find constant-region anchor (ASTKGPSVF). "
            "Pass a full heavy chain including CH1."
        )

    const_seq = chain_seq[pos:]
    aln = _aligner().align(WT_IgG1_CONST, const_seq)[0]
    # aln.aligned -> blocks of (ref_start,ref_end),(qry_start,qry_end)
    ref_aln, qry_aln = aln[0], aln[1]  # gapped strings, same length

    substitutions: list[dict] = []
    deletions: list[dict] = []
    insertions: list[dict] = []

    ref_idx = 0  # index into WT_IgG1_CONST
    last_eu = EU_NUMBERS[0]
    for r, q in zip(ref_aln, qry_aln):
        if r != "-" and q != "-":
            eu = EU_NUMBERS[ref_idx]
            last_eu = eu
            if r != q:
                mut = f"{r}{eu}{q}"
                substitutions.append({
                    "eu": eu, "wt": r, "mut": q, "domain": domain_for(eu),
                    "mutation": mut, "note": MUTATION_NOTES.get(mut, ""),
                })
            ref_idx += 1
        elif q == "-":  # residue present in WT, missing in query -> deletion
            eu = EU_NUMBERS[ref_idx]
            last_eu = eu
            # A clean C-terminal truncation (e.g. des-K447) is expected; record it.
            deletions.append({"eu": eu, "wt": r, "domain": domain_for(eu)})
            ref_idx += 1
        else:  # r == "-" : extra residue in query -> insertion after last_eu
            insertions.append({"after_eu": last_eu, "residue": q})

    mut_set = {s["mutation"] for s in substitutions}
    variants_present = []
    variants_partial = []
    for name, members in NAMED_VARIANTS.items():
        if members <= mut_set:
            variants_present.append(name)
        elif members & mut_set:
            have = sorted(members & mut_set)
            variants_partial.append({"name": name, "have": have,
                                     "missing": sorted(members - mut_set)})

    # Drop pure C-terminal truncation deletions from the "warning" divergence
    # count; only interior deletions are unusual.
    interior_del = [d for d in deletions if d["eu"] < EU_NUMBERS[-1]]
    divergence = len(substitutions) + len(interior_del) + len(insertions)
    likely_non_igg1 = divergence > 20

    return {
        "label": label,
        "substitutions": substitutions,
        "deletions": deletions,
        "insertions": insertions,
        "variants_present": variants_present,
        "variants_partial": variants_partial,
        "divergence": divergence,
        "likely_non_igg1": likely_non_igg1,
    }


def report(result: dict) -> None:
    label = result["label"]
    subs = result["substitutions"]
    print(f"\n{label} — {len(subs)} substitution(s) vs WT IGHG1*01:")

    if result["likely_non_igg1"]:
        print("  [!] High divergence from IgG1 — this chain may be IgG2/IgG4 or a "
              "different allotype. Only an IgG1 reference is shipped, so individual "
              "calls below may be artifacts of isotype, not true engineering.")

    if not subs:
        print("  (none — wild-type IgG1 constant region)")
    else:
        print(f"  {'EU #':>6}  {'Domain':<7}  {'Mut':<8}  Annotation")
        print(f"  {'------':>6}  {'-------':<7}  {'--------':<8}  ----------")
        for s in subs:
            print(f"  {s['eu']:>6}  {s['domain']:<7}  {s['mutation']:<8}  {s['note']}")

    for d in result["deletions"]:
        tag = "C-terminal clip" if d["eu"] == EU_NUMBERS[-1] else "interior deletion"
        print(f"  - deletion: {d['wt']}{d['eu']} ({d['domain']}, {tag})")
    for ins in result["insertions"]:
        print(f"  - insertion: +{ins['residue']} after EU {ins['after_eu']}")

    if result["variants_present"]:
        print(f"\n  Recognised variants present: {', '.join(result['variants_present'])}")
    for vp in result["variants_partial"]:
        print(f"  Partial {vp['name']}: have {vp['have']}, missing {vp['missing']}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Map EU-numbered mutations in an IgG heavy-chain constant region"
    )
    parser.add_argument("sequence", help="Full heavy chain amino acid sequence")
    parser.add_argument("--label", default="Heavy Chain", help="Label for display")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    args = parser.parse_args()

    try:
        seq = clean_sequence(args.sequence, name=args.label)
        result = find_mutations(seq, args.label)
    except SequenceError as e:
        print(f"[!] {e}", file=sys.stderr)
        sys.exit(1)

    if args.json:
        emit_json(result)
    else:
        report(result)


if __name__ == "__main__":
    main()
