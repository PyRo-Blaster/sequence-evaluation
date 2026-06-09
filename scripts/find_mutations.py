"""Map mutations in an IgG heavy-chain constant region vs the wild-type reference
of its isotype, in EU numbering.

The constant region is located by the conserved ASTKGPSVF (CH1 start) anchor. The
isotype (IgG1/IgG2/IgG4) is auto-detected by alignment, the query is *globally
aligned* to that isotype's WT reference (so indels do not shift downstream EU
numbers), and recognised engineering mutations are annotated automatically.

EU numbering for every reference is derived by aligning it to the verified IgG1
reference, which is reliable across CH1/CH2/CH3 and the conserved CPxCP hinge core
(including the IgG4 S228 site). The exact EU numbers of the non-conserved
N-terminal hinge residues of IgG2/IgG4 are approximate.

Reference data note: only the IgG1 reference is independently validated (it
reproduces an approved IgG1 therapeutic). The IgG2/IgG4 references are
reconstructed from canonical isotype differences. For definitive / regulatory
work pass an authoritative sequence via --reference-fasta.

Usage:
    uv run scripts/find_mutations.py <full_heavy_chain_sequence> [--label HC1] [--json]
    uv run scripts/find_mutations.py <seq> --isotype IgG4
    uv run scripts/find_mutations.py <seq> --reference-fasta my_wt.fasta

# /// script
# requires-python = ">=3.10"
# dependencies = ["biopython>=1.84"]
# ///
"""

import argparse
import sys

from Bio.Align import PairwiseAligner

from _common import SequenceError, clean_sequence, emit_json

# WT human IgG1 constant region (IGHG1*01): CH1 + Hinge + CH2 + CH3.
# Validated: VH + this sequence reproduces the trastuzumab heavy chain.
WT_IgG1_CONST = (
    "ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSS"
    "VVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTL"
    "MISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEY"
    "KCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQP"
    "ENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVLHEALHNHYTQKSLSLSPGK"
)

# WT human IgG2 (IGHG2*01) and IgG4 (IGHG4*01) constant regions, reconstructed
# from canonical isotype differences (shorter hinges; CH1/CH2/CH3 substitutions).
# Best-effort — verify with --reference-fasta for definitive work.
WT_IgG2_CONST = (
    "ASTKGPSVFPLAPCSRSTSESTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSS"
    "VVTVPSSNFGTQTYTCNVDHKPSNTKVDKTVERKCCVECPPCPAPPVAGPSVFLFPPKPKDTLMISR"
    "TPEVTCVVVDVSHEDPEVQFNWYVDGVEVHNAKTKPREEQFNSTFRVVSVLTVVHQDWLNGKEYKCK"
    "VSNKGLPAPIEKTISKTKGQPREPQVYTLPPSREEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPEN"
    "NYKTTPPMLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK"
)
WT_IgG4_CONST = (
    "ASTKGPSVFPLAPCSRSTSESTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSS"
    "VVTVPSSSLGTKTYTCNVDHKPSNTKVDKRVESKYGPPCPSCPAPEFLGGPSVFLFPPKPKDTLMIS"
    "RTPEVTCVVVDVSQEDPEVQFNWYVDGVEVHNAKTKPREEQFNSTYRVVSVLTVLHQDWLNGKEYKC"
    "KVSNKGLPSSIEKTISKAKGQPREPQVYTLPPSQEEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPE"
    "NNYKTTPPVLDSDGSFFLYSRLTVDKSRWQEGNVFSCSVMHEALHNHYTQKSLSLSLGK"
)

BUILTIN_REFERENCES = {
    "IgG1": WT_IgG1_CONST,
    "IgG2": WT_IgG2_CONST,
    "IgG4": WT_IgG4_CONST,
}

# EU numbering for the IgG1 reference (CH1 118-215, Hinge 216-230, CH2 231-340,
# CH3 341-447). Other isotypes inherit EU numbers via alignment to IgG1.
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
    "S228P": "IgG4 hinge stabilization (prevents Fab-arm exchange)",
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
    "IgG4 S228P": {"S228P"},
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


def eu_array_for(ref_seq: str) -> list[int | None]:
    """EU number per residue of ``ref_seq``, derived by alignment to IgG1.

    IgG1 returns its canonical numbering directly. For other references each
    residue inherits the EU number of the IgG1 residue it aligns to; residues
    aligning to an IgG1 gap (insertions vs IgG1) get None.
    """
    if ref_seq == WT_IgG1_CONST:
        return list(EU_NUMBERS)
    aln = _aligner().align(WT_IgG1_CONST, ref_seq)[0]
    igg1_aln, ref_aln = aln[0], aln[1]
    out: list[int | None] = []
    igg1_idx = 0
    for g, r in zip(igg1_aln, ref_aln):
        if g != "-" and r != "-":
            out.append(EU_NUMBERS[igg1_idx])
            igg1_idx += 1
        elif g != "-":  # gap in ref -> consume an IgG1 position
            igg1_idx += 1
        else:  # residue in ref aligned to IgG1 gap (rare for IgG2/4)
            out.append(None)
    return out


def detect_isotype(const_seq: str, references: dict[str, str]) -> tuple[str, dict[str, float]]:
    """Return (best isotype name, {name: alignment score})."""
    aligner = _aligner()
    scores = {name: float(aligner.score(seq, const_seq)) for name, seq in references.items()}
    best = max(scores, key=scores.get)
    return best, scores


def find_mutations(
    chain_seq: str,
    label: str,
    isotype: str = "auto",
    references: dict[str, str] | None = None,
) -> dict:
    references = references or BUILTIN_REFERENCES
    anchor = "ASTKGPSVF"
    pos = chain_seq.find(anchor)
    if pos == -1:
        raise SequenceError(
            f"{label}: could not find constant-region anchor (ASTKGPSVF). "
            "Pass a full heavy chain including CH1."
        )
    const_seq = chain_seq[pos:]

    scores = None
    if isotype == "auto":
        isotype, scores = detect_isotype(const_seq, references)
    elif isotype not in references:
        raise SequenceError(f"{label}: unknown isotype '{isotype}' "
                            f"(have {', '.join(references)})")
    ref_seq = references[isotype]
    ref_eu = eu_array_for(ref_seq)

    aln = _aligner().align(ref_seq, const_seq)[0]
    ref_aln, qry_aln = aln[0], aln[1]  # gapped strings, same length

    substitutions: list[dict] = []
    deletions: list[dict] = []
    insertions: list[dict] = []

    ref_idx = 0
    last_eu = ref_eu[0] or EU_NUMBERS[0]
    for r, q in zip(ref_aln, qry_aln):
        if r != "-" and q != "-":
            eu = ref_eu[ref_idx]
            if eu is not None:
                last_eu = eu
                if r != q:
                    mut = f"{r}{eu}{q}"
                    substitutions.append({
                        "eu": eu, "wt": r, "mut": q, "domain": domain_for(eu),
                        "mutation": mut, "note": MUTATION_NOTES.get(mut, ""),
                    })
            ref_idx += 1
        elif q == "-":  # residue present in WT, missing in query -> deletion
            eu = ref_eu[ref_idx]
            if eu is not None:
                last_eu = eu
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

    eu_last = next((e for e in reversed(ref_eu) if e is not None), EU_NUMBERS[-1])
    interior_del = [d for d in deletions if d["eu"] < eu_last]
    divergence = len(substitutions) + len(interior_del) + len(insertions)
    # After choosing the right isotype, lots of leftover divergence suggests an
    # unusual allotype or a reference mismatch.
    high_divergence = divergence > 20

    return {
        "label": label,
        "isotype": isotype,
        "isotype_scores": scores,
        "substitutions": substitutions,
        "deletions": deletions,
        "insertions": insertions,
        "variants_present": variants_present,
        "variants_partial": variants_partial,
        "divergence": divergence,
        "high_divergence": high_divergence,
    }


def report(result: dict) -> None:
    label = result["label"]
    subs = result["substitutions"]
    iso = result["isotype"]
    print(f"\n{label} — isotype {iso}; {len(subs)} substitution(s) vs WT {iso}:")
    if result.get("isotype_scores"):
        ranked = sorted(result["isotype_scores"].items(), key=lambda x: -x[1])
        print("  isotype match scores: " + ", ".join(f"{n} {s:.0f}" for n, s in ranked))

    if result["high_divergence"]:
        print("  [!] High residual divergence even against the best isotype — may be "
              "an unusual allotype or a reference mismatch. Consider --reference-fasta.")

    if not subs:
        print(f"  (none — wild-type {iso} constant region)")
    else:
        print(f"  {'EU #':>6}  {'Domain':<7}  {'Mut':<8}  Annotation")
        print(f"  {'------':>6}  {'-------':<7}  {'--------':<8}  ----------")
        for s in subs:
            print(f"  {s['eu']:>6}  {s['domain']:<7}  {s['mutation']:<8}  {s['note']}")

    eu_last = EU_NUMBERS[-1]
    for d in result["deletions"]:
        tag = "C-terminal clip" if d["eu"] == eu_last else "interior deletion"
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
    parser.add_argument("--isotype", default="auto",
                        choices=["auto", "IgG1", "IgG2", "IgG4"],
                        help="Reference isotype (default: auto-detect)")
    parser.add_argument("--reference-fasta",
                        help="Custom WT constant-region FASTA to compare against "
                             "(overrides built-in references; recommended for "
                             "definitive IgG2/IgG4 work)")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    args = parser.parse_args()

    references = dict(BUILTIN_REFERENCES)
    isotype = args.isotype
    if args.reference_fasta:
        from Bio import SeqIO
        rec = next(iter(SeqIO.parse(args.reference_fasta, "fasta")), None)
        if rec is None:
            print(f"[!] no sequence in {args.reference_fasta}", file=sys.stderr)
            sys.exit(1)
        try:
            custom = clean_sequence(str(rec.seq), name="reference")
        except SequenceError as e:
            print(f"[!] {e}", file=sys.stderr)
            sys.exit(1)
        references = {"custom": custom}
        isotype = "custom"

    try:
        seq = clean_sequence(args.sequence, name=args.label)
        result = find_mutations(seq, args.label, isotype=isotype, references=references)
    except SequenceError as e:
        print(f"[!] {e}", file=sys.stderr)
        sys.exit(1)

    if args.json:
        emit_json(result)
    else:
        report(result)


if __name__ == "__main__":
    main()
