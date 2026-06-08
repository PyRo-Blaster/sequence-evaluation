"""Assemble a Markdown sequence-evaluation report from per-stage JSON outputs.

Each analysis script supports ``--json``; redirect those to files and pass them
here so the final report is built deterministically rather than retyped by hand.

Usage:
    uv run scripts/compile_report.py --title "mAb-X evaluation" \
        --properties props.json \
        --cdrs vh.json vl.json \
        --mutations hc1.json hc2.json \
        --liabilities liab.json \
        --interfaces iface_A.json iface_B.json \
        --homology homology_section.md \
        > report.md

All inputs are optional; sections are emitted only for the data provided.

# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""

import argparse
import json
import sys
from pathlib import Path


def load(path: str):
    try:
        return json.loads(Path(path).read_text())
    except Exception as e:
        print(f"[!] could not read {path}: {e}", file=sys.stderr)
        sys.exit(1)


def section_properties(out, data):
    out.append("## Chain Properties\n")
    out.append("| Chain | Precursor | Mature | MW (Da) | pI | A280 (1 g/L, ox) | Instability |")
    out.append("|-------|-----------|--------|---------|----|------------------|-------------|")
    for c in data.get("chains", []):
        out.append(
            f"| {c['name']} | {c['precursor_len']} | {c['mature_len']} | "
            f"{c['mature_mw']:,.2f} | {c['mature_pi']:.2f} | {c['a280_oxidized']:.3f} | "
            f"{c['instability']:.1f} |"
        )
    cx = data.get("complex")
    if cx:
        out.append("")
        out.append(f"**Assembled complex:** MW {cx['corrected_mw']:,.2f} Da "
                   f"({cx['disulfide_bonds']} disulfide bonds), pI {cx['pi']:.2f}")
    out.append("")


def section_cdrs(out, items):
    out.append("## CDR Annotation\n")
    for d in items:
        out.append(f"### {d['name']} ({d['domain_type']}) — engine: {d['engine']}\n")
        out.append("| CDR | Length | Sequence |")
        out.append("|-----|--------|----------|")
        for n in (1, 2, 3):
            out.append(f"| CDR{n} | {d[f'cdr{n}_len']} | `{d[f'cdr{n}']}` |")
        out.append("")


def section_mutations(out, items):
    out.append("## Fc Mutations (EU numbering)\n")
    for d in items:
        subs = d.get("substitutions", [])
        out.append(f"### {d['label']} — {len(subs)} substitution(s) vs WT IGHG1*01\n")
        if d.get("likely_non_igg1"):
            out.append("> ⚠️ High divergence from IgG1 — chain may be a different "
                       "isotype/allotype; calls below may be artifacts.\n")
        if subs:
            out.append("| EU # | Domain | Mutation | Annotation |")
            out.append("|------|--------|----------|------------|")
            for s in subs:
                out.append(f"| {s['eu']} | {s['domain']} | {s['mutation']} | {s['note'] or '—'} |")
        else:
            out.append("Wild-type IgG1 constant region (no substitutions).")
        if d.get("variants_present"):
            out.append(f"\n**Recognised variants:** {', '.join(d['variants_present'])}")
        for ins in d.get("insertions", []):
            out.append(f"- insertion: +{ins['residue']} after EU {ins['after_eu']}")
        for dele in d.get("deletions", []):
            out.append(f"- deletion: {dele['wt']}{dele['eu']} ({dele['domain']})")
        out.append("")


def section_liabilities(out, data):
    out.append("## Developability Liabilities\n")
    for name, res in data.items():
        out.append(f"### {name} ({res['length']} aa)\n")
        flags = []
        if res.get("odd_cysteine_count"):
            flags.append(f"⚠️ odd cysteine count ({res['cysteine_count']}) — likely free thiol")
        if res.get("nterm_pyroglutamate"):
            flags.append(f"⚠️ N-terminal {res['nterm_residue']} (pyroglutamate-prone)")
        if flags:
            out.append("- " + "\n- ".join(flags))
        findings = res.get("findings", [])
        if findings:
            out.append("\n| Pos | Severity | Motif | Liability |")
            out.append("|-----|----------|-------|-----------|")
            order = {"high": 0, "medium": 1, "low": 2, "info": 3}
            for f in sorted(findings, key=lambda x: (order.get(x["severity"], 9), x["position"])):
                out.append(f"| {f['position']} | {f['severity']} | {f['motif']} | {f['liability']} |")
        out.append("")


def section_interfaces(out, items):
    out.append("## Structural Predictions — Interface Analysis\n")
    for d in items:
        ab = ",".join(d["antibody_chains"])
        ag = ",".join(d["antigen_chains"])
        out.append(f"### {Path(d['cif']).name} (Ab {ab} vs Ag {ag}, cutoff {d['cutoff']} Å)\n")
        out.append(f"- Atom-atom contacts: **{d['total_atom_contacts']}**")
        para = d.get("paratope", [])
        epi = d.get("epitope", [])
        para_str = ", ".join(
            f"{p['chain']}:{p['resname']}{p['resnum']}"
            + (f"[{p['region']}]" if p.get("region") and p["region"] != "FR" else "")
            for p in para
        )
        epi_str = ", ".join(f"{e['chain']}:{e['resname']}{e['resnum']}" for e in epi)
        out.append(f"- Paratope ({len(para)}): {para_str or 'None'}")
        out.append(f"- Epitope ({len(epi)}): {epi_str or 'None'}")
        out.append("")


def main() -> None:
    p = argparse.ArgumentParser(description="Compile a Markdown evaluation report from stage JSON")
    p.add_argument("--title", default="Sequence Evaluation Report")
    p.add_argument("--properties")
    p.add_argument("--cdrs", nargs="*", default=[])
    p.add_argument("--mutations", nargs="*", default=[])
    p.add_argument("--liabilities")
    p.add_argument("--interfaces", nargs="*", default=[])
    p.add_argument("--homology", help="Path to a pre-written Markdown homology section")
    args = p.parse_args()

    out: list[str] = [f"# {args.title}\n"]

    if args.properties:
        section_properties(out, load(args.properties))
    if args.cdrs:
        section_cdrs(out, [load(x) for x in args.cdrs])
    if args.mutations:
        section_mutations(out, [load(x) for x in args.mutations])
    if args.liabilities:
        section_liabilities(out, load(args.liabilities))
    if args.homology:
        out.append("## Homology Search\n")
        out.append(Path(args.homology).read_text().rstrip())
        out.append("")
    if args.interfaces:
        section_interfaces(out, [load(x) for x in args.interfaces])

    sys.stdout.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
