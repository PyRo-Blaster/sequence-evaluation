"""Assemble a phased Markdown evaluation report from per-stage JSON outputs.

The report walks the molecule from triage to binder, and every phase emits two
things: **Findings** (the data) and an **Assessment** (what it means). Assessments
are derived deterministically from the JSON by the rules in `references/` — they
are not free prose, and where an input is absent the phase says so rather than
inferring.

Phases
    1  Molecule Triage          what this is, chain composition, architecture
    2  Chain Properties         MW, pI, A280, instability
    3  Domain Annotation        CDRs and the numbering engine behind them
    4  Engineering Audit        Fc mutations, heterodimerization coherence
    5  Sequence Liabilities     CDR-escalated, ranked
    6  Structural Developability surface patches on the selected model
    7  Lineage                  homology search
    8  Verdict                  ranked risk register and overall call

Paratope/epitope interface reporting is deliberately NOT a phase: it depends on
a cofolded antibody-antigen complex, which is a GPU-required job. Run
`analyze_interfaces.py` directly when a complex exists; its output still feeds
`scan_liabilities.py --paratope-json` to escalate oxidation hotspots.

Usage:
    uv run scripts/compile_report.py --title "mAb-X evaluation" \
        --properties props.json \
        --cdrs vh.json vl.json \
        --mutations hc1.json hc2.json \
        --liabilities liab.json \
        --patches patches.json \
        --homology homology_section.md \
        > report.md

All inputs are optional; phases are emitted for the data provided and the
Verdict records which phases were not assessed.

# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""

import argparse
import datetime
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

SEV_ORDER = {"high": 0, "medium": 1, "low": 2, "expected": 3, "info": 4}

# Thresholds used by the assessment rules: the ProtParam conventions already
# cited in SKILL.md.
INSTABILITY_UNSTABLE = 40.0
GRAVY_SOLUBLE_MAX = 0.4


def load(path: str):
    try:
        return json.loads(Path(path).read_text())
    except Exception as e:
        print(f"[!] could not read {path}: {e}", file=sys.stderr)
        sys.exit(1)


def _fmt_list(items, empty="none"):
    items = [str(i) for i in items if i]
    return ", ".join(items) if items else empty


class Report:
    """Accumulates Markdown plus the cross-phase facts the Verdict needs."""

    def __init__(self, title: str):
        self.out: list[str] = []
        self.title = title
        self.phase_no = 0
        # Collected for the Verdict section.
        self.risks: list[dict] = []          # {severity, where, what, action}
        self.not_assessed: list[str] = []
        self.blocking: list[str] = []        # construct-level errors

    # ---------- emission helpers ----------
    def phase(self, name: str, assessed: bool = True):
        self.phase_no += 1
        self.out.append(f"\n## Phase {self.phase_no} — {name}\n")
        if not assessed:
            self.not_assessed.append(f"Phase {self.phase_no} — {name}")

    def findings(self):
        self.out.append("**Findings**\n")

    def assessment(self, lines: list[str]):
        self.out.append("\n**Assessment**\n")
        for line in lines:
            self.out.append(f"- {line}")
        self.out.append("")

    def skipped(self, reason: str):
        self.out.append(f"_Not assessed — {reason}._\n")

    def add_risk(self, severity: str, where: str, what: str, action: str):
        self.risks.append({"severity": severity, "where": where,
                           "what": what, "action": action})


# ---------------------------------------------------------------- phase 1
def phase_triage(r: Report, props, cdrs, muts):
    r.phase("Molecule Triage", assessed=bool(props or cdrs or muts))
    if not (props or cdrs or muts):
        r.skipped("no properties, CDR or mutation input supplied")
        return

    chains = (props or {}).get("chains", [])
    domains = cdrs or []
    heavy = [d for d in domains if d.get("domain_type") in ("VH", "VHH")]
    light = [d for d in domains if d.get("domain_type") == "VL"]
    vhh = [d for d in domains if d.get("domain_type") == "VHH"]
    isotypes = sorted({m.get("isotype") for m in (muts or []) if m.get("isotype")})

    r.findings()
    if chains:
        r.out.append(f"- Chains in input: **{len(chains)}** — {_fmt_list(c['name'] for c in chains)}")
    if domains:
        r.out.append(f"- Variable domains annotated: **{len(domains)}** "
                     f"({len(heavy)} heavy-chain or VHH, {len(light)} light-chain)")
    if muts:
        r.out.append(f"- Constant regions analysed: **{len(muts)}** "
                     f"({_fmt_list(m.get('label') for m in muts)})")
    if isotypes:
        r.out.append(f"- Isotype(s) detected: **{_fmt_list(isotypes)}**")
    r.out.append("")

    # Modality inference from what is actually present.
    a = []
    if vhh:
        modality = f"VHH / nanobody ({len(vhh)} domain(s))"
    elif heavy and light and len(heavy) >= 2:
        modality = f"bispecific or multispecific — {len(heavy)} heavy chains paired with {len(light)} light chain(s)"
    elif heavy and light:
        modality = "monospecific antibody (one VH + one VL)"
    elif heavy or light:
        modality = f"single variable domain ({_fmt_list(d.get('domain_type') for d in domains)})"
    elif muts:
        modality = "Fc-containing construct (constant region only; no variable domain supplied)"
    else:
        modality = "general protein — no antibody markers in the supplied analyses"
    a.append(f"Modality read: **{modality}**.")

    if len(isotypes) > 1:
        a.append(f"Mixed isotypes across chains ({_fmt_list(isotypes)}) — confirm this is "
                 "intended; a mismatched pair changes effector profile and purification.")
        r.add_risk("medium", "architecture", f"mixed isotypes {_fmt_list(isotypes)}",
                   "confirm intended; verify effector and Protein A behaviour")
    if domains and not muts:
        a.append("No constant region supplied, so this is an Fv/scFv-level assessment: "
                 "Fc engineering and isotype are **not assessed**.")
    if chains and not domains:
        a.append("Chains were measured but no variable domain was annotated — "
                 "CDR-dependent phases below are unescalated.")
    r.assessment(a)


# ---------------------------------------------------------------- phase 2
def phase_properties(r: Report, data):
    r.phase("Chain Properties", assessed=bool(data))
    if not data:
        r.skipped("no analyze_properties.py JSON supplied")
        return

    chains = data.get("chains", [])
    r.findings()
    r.out.append("| Chain | Precursor | Mature | MW (Da) | pI | A280 (1 g/L, ox) | Instability | GRAVY |")
    r.out.append("|-------|-----------|--------|---------|----|------------------|-------------|-------|")
    for c in chains:
        r.out.append(
            f"| {c.get('name','?')} | {c.get('precursor_len','?')} | {c.get('mature_len','?')} | "
            f"{c.get('mature_mw',0):,.2f} | {c.get('mature_pi',0):.2f} | "
            f"{c.get('a280_oxidized',0):.3f} | {c.get('instability',0):.1f} | "
            f"{c.get('gravy',0):.3f} |"
        )
    cx = data.get("complex")
    if cx:
        r.out.append("")
        stoich = cx.get("stoichiometry")
        stoich_s = (", stoichiometry " + ", ".join(f"{k}x{v}" for k, v in stoich.items())
                    if stoich else "")
        r.out.append(f"**Assembled complex:** MW {cx.get('corrected_mw',0):,.2f} Da "
                     f"({cx.get('disulfide_bonds',0)} disulfide bonds{stoich_s}), "
                     f"pI {cx.get('pi',0):.2f}")
    r.out.append("")

    a = []
    unstable = [c for c in chains if c.get("instability", 0) > INSTABILITY_UNSTABLE]
    if unstable:
        names = _fmt_list(f"{c['name']} ({c['instability']:.0f})" for c in unstable)
        a.append(f"Instability index above {INSTABILITY_UNSTABLE:.0f} for {names}. "
                 "The index is a sequence heuristic, not a stability measurement — "
                 "treat it as a prompt for a thermal/accelerated study, not a finding.")
        for c in unstable:
            r.add_risk("low", f"chain {c['name']}",
                       f"instability index {c['instability']:.0f} (>40)",
                       "confirm with DSF/accelerated stability rather than acting on the index")
    else:
        a.append(f"All chains below the instability threshold of {INSTABILITY_UNSTABLE:.0f}.")

    greasy = [c for c in chains if c.get("gravy", 0) > GRAVY_SOLUBLE_MAX]
    if greasy:
        a.append(f"GRAVY above {GRAVY_SOLUBLE_MAX} for {_fmt_list(c['name'] for c in greasy)} — "
                 "unusually hydrophobic for a secreted protein; watch expression and aggregation.")
        for c in greasy:
            r.add_risk("medium", f"chain {c['name']}", f"GRAVY {c['gravy']:.2f} (>{GRAVY_SOLUBLE_MAX})",
                       "monitor aggregation and expression titre")

    if cx:
        pi = cx.get("pi", 0)
        if pi and (pi < 6.0 or pi > 9.0):
            a.append(f"Complex pI {pi:.2f} sits outside the usual 6–9 therapeutic window — "
                     "relevant to formulation pH and ion-exchange polishing.")
        if not cx.get("stoichiometry"):
            a.append("Complex MW assumes **one copy of each chain in the FASTA**. For an IgG "
                     "pass `--stoichiometry` (e.g. `HC=2 LC=2`) or the mass is half the real value.")
    else:
        a.append("No assembled-complex figure — rerun with `--complex` for MW and complex pI.")

    a.append("A280 is the 1 g/L oxidized value (ExPASy convention); use it directly for "
             "concentration determination.")
    r.assessment(a)


# ---------------------------------------------------------------- phase 3
def phase_domains(r: Report, items):
    r.phase("Domain Annotation", assessed=bool(items))
    if not items:
        r.skipped("no find_cdrs.py JSON supplied")
        return

    r.findings()
    r.out.append("| Domain | Type | Engine | CDR1 | CDR2 | CDR3 | CDR3 len |")
    r.out.append("|--------|------|--------|------|------|------|----------|")
    for d in items:
        r.out.append(
            f"| {d.get('name','?')} | {d.get('domain_type','?')} | {d.get('engine','?')} | "
            f"`{d.get('cdr1','')}` | `{d.get('cdr2','')}` | `{d.get('cdr3','')}` | "
            f"{d.get('cdr3_len',0)} |"
        )
    r.out.append("")

    a = []
    engines = {d.get("engine", "?") for d in items}
    regex = [d for d in items if "anchor-regex" in (d.get("engine") or "")]
    if regex:
        a.append(f"{len(regex)} of {len(items)} domain(s) annotated with the **anchor-regex "
                 "fallback**, not ANARCI. CDR1 and CDR2 are reliable; the CDR3 N-terminal "
                 "boundary may be ±1. Install `abnumber` for true scheme numbering before "
                 "quoting these boundaries externally.")
        r.add_risk("low", "annotation", "CDR boundaries from the regex fallback",
                   "re-run with abnumber/ANARCI before publishing CDR definitions")
    else:
        a.append(f"Scheme numbering from {_fmt_list(engines)} — boundaries are scheme-exact.")

    long_h3 = [d for d in items if d.get("domain_type") in ("VH", "VHH") and d.get("cdr3_len", 0) >= 20]
    if long_h3:
        a.append(f"Long CDR-H3 in {_fmt_list(d['name'] for d in long_h3)} "
                 f"({_fmt_list(str(d['cdr3_len']) for d in long_h3)} aa) — long H3 loops are the "
                 "least reliable part of any structure prediction and a common source of "
                 "model-dependent patch and paratope calls.")
    short = [d for d in items if d.get("cdr3_len", 99) <= 4]
    if short:
        a.append(f"Very short CDR3 in {_fmt_list(d['name'] for d in short)} — check the domain "
                 "was passed without its constant region and that the anchors resolved.")
        r.add_risk("medium", "annotation", "implausibly short CDR3",
                   "verify the input was a clean variable domain")
    r.assessment(a)


# ---------------------------------------------------------------- phase 4
# Heterodimerization sets that must appear on opposite chains of a bispecific.
KIH_KNOB = {"T366W"}
KIH_HOLE = {"T366S", "L368A", "Y407V"}


def _architecture_check(muts):
    """Verify heterodimerization mutation sets are structurally coherent.

    Returns (findings, errors). A knob on both arms, or a knob with no hole,
    cannot assemble as intended — that is a construct error, not a risk flag.
    """
    findings, errors = [], []
    if len(muts) < 2:
        return findings, errors
    arms = {}
    for m in muts:
        subs = {s.get("mutation") for s in m.get("substitutions", [])}
        label = m.get("label", "?")
        arms[label] = {
            "knob": bool(subs & KIH_KNOB),
            "hole": bool(subs & KIH_HOLE),
            "hole_partial": 0 < len(subs & KIH_HOLE) < len(KIH_HOLE),
        }
    knobs = [k for k, v in arms.items() if v["knob"]]
    holes = [k for k, v in arms.items() if v["hole"]]
    if not knobs and not holes:
        return findings, errors

    findings.append(f"Knob (T366W) on: {_fmt_list(knobs)}")
    findings.append(f"Hole (T366S/L368A/Y407V) on: {_fmt_list(holes)}")
    if len(knobs) > 1:
        errors.append(f"Knob mutation present on {len(knobs)} chains ({_fmt_list(knobs)}) — "
                      "two knobs cannot heterodimerize.")
    if len(holes) > 1:
        errors.append(f"Hole mutation set present on {len(holes)} chains ({_fmt_list(holes)}) — "
                      "two holes cannot heterodimerize.")
    if knobs and not holes:
        errors.append("Knob present with no complementary hole chain — heterodimer cannot form "
                      "as designed.")
    if holes and not knobs:
        errors.append("Hole set present with no complementary knob chain — heterodimer cannot "
                      "form as designed.")
    for label, v in arms.items():
        if v["hole_partial"]:
            errors.append(f"{label} carries an incomplete hole set — all three of "
                          "T366S/L368A/Y407V are required.")
    if not errors and len(knobs) == 1 and len(holes) == 1 and knobs != holes:
        findings.append("Knob and hole are on opposite chains, exactly once each — coherent.")
    return findings, errors


def phase_engineering(r: Report, items):
    r.phase("Engineering Audit", assessed=bool(items))
    if not items:
        r.skipped("no find_mutations.py JSON supplied (no constant region analysed)")
        return

    r.findings()
    for d in items:
        subs = d.get("substitutions", [])
        iso = d.get("isotype", "IgG1")
        r.out.append(f"\n**{d.get('label','?')}** — isotype {iso}; "
                     f"{len(subs)} substitution(s) vs WT {iso}\n")
        if d.get("high_divergence"):
            r.out.append("> High residual divergence even against the best-matching isotype — "
                         "possible unusual allotype or reference mismatch.\n")
        if subs:
            r.out.append("| EU # | Domain | Mutation | Annotation |")
            r.out.append("|------|--------|----------|------------|")
            for s in subs:
                r.out.append(f"| {s.get('eu','?')} | {s.get('domain','?')} | "
                             f"{s.get('mutation','?')} | {s.get('note') or '—'} |")
        else:
            r.out.append(f"Wild-type {iso} constant region (no substitutions).")
        if d.get("variants_present"):
            r.out.append(f"\nRecognised variants: {_fmt_list(d['variants_present'])}")
        for ins in d.get("insertions", []):
            r.out.append(f"- insertion: +{ins.get('residue')} after EU {ins.get('after_eu')}")
        for dele in d.get("deletions", []):
            r.out.append(f"- deletion: {dele.get('wt')}{dele.get('eu')} ({dele.get('domain')})")
    arch_findings, arch_errors = _architecture_check(items)
    if arch_findings:
        r.out.append("\n**Heterodimerization architecture**\n")
        for f in arch_findings:
            r.out.append(f"- {f}")
    r.out.append("")

    a = []
    all_variants = sorted({v for d in items for v in (d.get("variants_present") or [])})
    if all_variants:
        a.append(f"Engineering present: **{_fmt_list(all_variants)}**.")
    wt = [d for d in items if not d.get("substitutions")]
    if wt:
        a.append(f"{_fmt_list(d.get('label') for d in wt)} carries an unmodified constant region.")
    for e in arch_errors:
        a.append(f"**Construct error:** {e}")
        r.blocking.append(e)
        r.add_risk("high", "architecture", e, "correct the construct before expression")
    if arch_findings and not arch_errors:
        a.append("Heterodimerization set is complementary and correctly distributed.")
    if any(d.get("high_divergence") for d in items):
        a.append("At least one chain diverges strongly from every built-in reference — pass an "
                 "authoritative IMGT/UniProt sequence via `--reference-fasta` before relying on "
                 "the EU calls.")
        r.add_risk("medium", "Fc numbering", "high divergence vs built-in reference",
                   "supply --reference-fasta from IMGT/UniProt")
    isos = {d.get("isotype") for d in items}
    if isos - {"IgG1"} and isos:
        a.append("Non-IgG1 isotype present: only the IgG1 reference in this skill is "
                 "independently validated; IgG2/IgG4 references are reconstructed. Treat those "
                 "EU numbers as provisional for regulatory work.")
    r.assessment(a)


# ---------------------------------------------------------------- phase 5
def phase_liabilities(r: Report, data):
    r.phase("Sequence Liabilities", assessed=bool(data))
    if not data:
        r.skipped("no scan_liabilities.py JSON supplied")
        return

    # scan_liabilities.py emits {"chains": ..., "cysteine_assessment": ...};
    # tolerate a flat {chain: result} mapping too rather than silently
    # iterating the wrapper keys and reporting "no hits".
    cys_assessment = data.get("cysteine_assessment") if isinstance(data, dict) else None
    data = data.get("chains", data) if isinstance(data, dict) else data

    any_cdr = any(res.get("cdr_annotated") for res in data.values())
    r.findings()
    for name, res in data.items():
        r.out.append(f"\n**{name}** ({res.get('length','?')} aa)\n")
        flags = []
        if res.get("odd_cysteine_count"):
            if res.get("multi_chain_input"):
                flags.append(
                    f"{res.get('cysteine_count')} cysteines (odd) — expected in an assembled "
                    "multi-chain antibody, where inter-chain disulfides pair across chains")
            else:
                flags.append(f"odd cysteine count ({res.get('cysteine_count')}) — likely free thiol")
        if res.get("nterm_pyroglutamate"):
            flags.append(f"N-terminal {res.get('nterm_residue')} (pyroglutamate-prone)")
        for f in flags:
            r.out.append(f"- {f}")
        findings = res.get("findings", [])
        ranked = sorted(findings, key=lambda x: (SEV_ORDER.get(x["severity"], 9), x["position"]))
        notable = [f for f in ranked if f["severity"] in ("high", "medium")]
        if notable:
            r.out.append("\n| Pos | Severity | Motif | CDR | Liability |")
            r.out.append("|-----|----------|-------|-----|-----------|")
            for f in notable:
                sev = f["severity"] + ("\\*" if f.get("escalated") else "")
                label = f["liability"] + (" — conserved Fc N297, expected"
                                          if f.get("fc_n297") else "")
                r.out.append(f"| {f['position']} | {sev} | `{f['motif']}` | "
                             f"{f.get('cdr') or '—'} | {label} |")
            low_n = len(ranked) - len(notable)
            if low_n:
                r.out.append(f"\n_{low_n} further low-severity hit(s) (unescalated "
                             "Met/Trp oxidation) omitted from this table._")
        elif ranked:
            r.out.append(f"\n{len(ranked)} hit(s), all low severity.")
        else:
            r.out.append("\nNo hits in the scanned categories.")
    r.out.append("")
    if any_cdr:
        r.out.append("\\* escalated one step because the residue lies in a CDR "
                     "(two steps for oxidation on a contact-bearing CDR).\n")

    a = []
    if not any_cdr:
        a.append("**No CDR map was applied**, so every severity here is the unescalated default "
                 "and CDR-localised liabilities are not distinguished from framework ones. "
                 "Rebuild with `build_cdr_map.py` and pass `--cdr-json` — this is the single "
                 "highest-value input to this phase.")
        r.add_risk("medium", "liability scan", "no CDR context applied",
                   "re-run scan_liabilities.py with --cdr-json")
    cdr_high = [(name, f) for name, res in data.items()
                for f in res.get("findings", [])
                if f["severity"] == "high" and f.get("in_cdr")]
    if cdr_high:
        items = _fmt_list(f"{name} {f['cdr']} pos {f['position']} ({f['liability'].split(' (')[0]})"
                          for name, f in cdr_high)
        a.append(f"**{len(cdr_high)} high-severity liability/liabilities inside a CDR**: {items}. "
                 "These are the primary engineering targets — a deamidation or isomerization "
                 "motif in a CDR can erode potency on storage.")
        for name, f in cdr_high:
            r.add_risk("high", f"{name} {f['cdr']} pos {f['position']}",
                       f"{f['liability']} (motif {f['motif']})",
                       "engineer out, or justify with a forced-degradation study")
    glyc = [(name, f) for name, res in data.items()
            for f in res.get("findings", [])
            if f["liability"].startswith("N-glyc")]
    var_glyc = [(n, f) for n, f in glyc if f.get("in_cdr")]
    if var_glyc:
        where = _fmt_list(f"{n} {f.get('cdr')} pos {f['position']}" for n, f in var_glyc)
        a.append(f"N-glycosylation sequon inside a CDR ({where}) — remove unless "
                 "variable-domain glycosylation is intended; it brings heterogeneity and "
                 "can alter binding.")
        for n, f in var_glyc:
            r.add_risk("high", f"{n} {f.get('cdr')} pos {f['position']}",
                       "N-glyc sequon in a CDR",
                       "remove the sequon unless glycosylation is intended")
    elif glyc:
        a.append(f"{len(glyc)} N-glyc sequon(s), none in a CDR. The conserved Fc N297 sequon is "
                 "expected and should not be removed unless aglycosylation is the goal.")
    odd = [n for n, res in data.items() if res.get("odd_cysteine_count")]
    multi = any(res.get("multi_chain_input") for res in data.values())
    parity = next((res.get("assembly_cysteine_parity") for res in data.values()
                   if res.get("assembly_cysteine_parity")), None)
    total = next((res.get("assembly_cysteine_total") for res in data.values()
                  if res.get("assembly_cysteine_total")), None)
    if odd and multi:
        if parity == "even":
            a.append(f"Odd per-chain cysteine counts on {_fmt_list(odd)} are **expected** here: "
                     f"the {total} cysteines across one copy of each chain pair evenly "
                     f"({total // 2} disulfides), and an antibody heavy chain is odd by "
                     "construction because the HC–LC and hinge bonds pair across chains. "
                     "Not a free-thiol finding.")
        else:
            a.append(f"Cysteine count across the assembly is **odd** ({total} over one copy of "
                     "each chain), so at least one thiol cannot pair. Check for a truncated "
                     "chain or an engineered conjugation site.")
            r.add_risk("medium", "assembly", f"odd total cysteine count ({total})",
                       "identify the unpaired thiol; confirm intentional")
    elif odd:
        a.append(f"Odd cysteine count on {_fmt_list(odd)} — almost always an intentional "
                 "conjugation handle or an engineered free thiol. Confirm it is deliberate; "
                 "if not, it is a mispairing and aggregation risk.")
        for n in odd:
            r.add_risk("medium", f"chain {n}", "unpaired cysteine (odd count)",
                       "confirm intentional; otherwise pair or remove")
    pyro = [n for n, res in data.items() if res.get("nterm_pyroglutamate")]
    if pyro:
        a.append(f"N-terminal Gln/Glu on {_fmt_list(pyro)} — pyroglutamate formation is usually "
                 "benign but shows up as a charge variant; expect it in the CEX profile.")
    if not a:
        a.append("No high- or medium-severity liabilities in the scanned categories.")
    r.assessment(a)


# ---------------------------------------------------------------- phase 6
def phase_patches(r: Report, data):
    r.phase("Structural Developability", assessed=bool(data))
    if not data:
        r.skipped("no predicted model supplied, so charged and hydrophobic surface patches "
                  "are **blocked** — they cannot be inferred from sequence")
        return

    r.findings()
    r.out.append(
        f"Method `{data.get('method', 'surface-patch-v1')}` on chains "
        f"{_fmt_list(data.get('chains', []))} of `{Path(str(data.get('cif','?'))).name}`. "
        "This is **not** TAP PPC/PNC/PSH and **not** SAP; do not apply those cutoffs.\n")
    excluded = data.get("chains_excluded_from_sasa") or []
    if excluded:
        r.out.append(f"Chains {_fmt_list(excluded)} were present in the file and excluded before "
                     "the SASA calculation, so exposure is that of the isolated antibody "
                     "chains (apo), not the bound complex.\n")
    if data.get("coarse_model"):
        r.out.append("> Coarse model (missing side chains). These patch calls are not usable.\n")
    for w in data.get("warnings") or []:
        if data.get("coarse_model") and "side-chain" in w:
            continue
        if excluded and "excluded before the SASA" in w:
            continue
        r.out.append(f"> {w}\n")

    patches = data.get("patches") or []
    if patches:
        r.out.append("| Severity | Kind | n | CDR res | Score | Residues |")
        r.out.append("|----------|------|---|---------|-------|----------|")
        for p in patches:
            residues = ", ".join(
                f"{x['chain']}:{x['resname']}{x['resnum']}"
                + (f" [{x['region']}]" if x.get("region") not in ("FR", "unassigned", None) else "")
                for x in p.get("residues", [])
            )
            r.out.append(f"| {p['severity']} | {p['kind']} | {p['n_residues']} | "
                         f"{p['cdr_residues']} | {p['score']:.3f} | {residues} |")
    else:
        r.out.append("No surface patches at or above the minimum size.")
    r.out.append("")

    a = []
    if data.get("coarse_model"):
        a.append("**Model has no side chains** — nothing in this phase is interpretable. "
                 "Re-predict with an all-atom model before assessing surface patches.")
        r.add_risk("medium", "structural developability", "coarse model, patches unusable",
                   "re-predict an all-atom Fv model")
        r.assessment(a)
        return

    high = [p for p in patches if p["severity"] == "high"]
    if high:
        for p in high:
            res = _fmt_list(f"{x['chain']}:{x['resname']}{x['resnum']}" for x in p["residues"])
            a.append(f"**High-severity {p['kind'].replace('_',' ')} patch** "
                     f"({p['n_residues']} residues, {p['cdr_residues']} in a CDR, "
                     f"score {p['score']:.2f}): {res}.")
            r.add_risk("high", "Fv surface", f"{p['kind']} patch touching a CDR",
                       "review against viscosity/aggregation data; consider charge-spreading mutations")
    medium = [p for p in patches if p["severity"] == "medium"]
    cdr_medium = [p for p in medium if p.get("cdr_residues")]
    if cdr_medium:
        for p in cdr_medium:
            res = _fmt_list(f"{x['chain']}:{x['resname']}{x['resnum']}" for x in p["residues"])
            a.append(f"Medium {p['kind'].replace('_',' ')} patch touching a CDR "
                     f"({p['n_residues']} residues, score {p['score']:.2f}): {res}. "
                     "Below the high bar only because it has fewer than three residues.")
            r.add_risk("medium", "Fv surface",
                       f"{p['kind']} patch touching a CDR ({p['n_residues']} residues)",
                       "watch in viscosity and self-interaction screening")
    elif medium:
        a.append(f"{len(medium)} medium-severity patch(es), none touching a CDR — "
                 "framework surface, lower priority.")
    if not data.get("chains_excluded_from_sasa") and len(data.get("chains_in_file") or []) <= len(data.get("chains", [])):
        a.append("Exposure computed on the antibody chains alone, as intended.")
    if not patches:
        a.append("No clustered charged or hydrophobic surface detected at the default "
                 "0.20 exposure and 7.0 Å clustering settings.")
    a.append("Patch membership is read from **one static model**; a predicted CDR-H3 can move "
             "it. Treat these as review flags for this model, not measurements of the molecule.")
    if any("No --cdr-json" in w for w in (data.get("warnings") or [])):
        a.append("No CDR map was passed, so nothing could be flagged high regardless of size — "
                 "pass `--cdr-json` to separate CDR-local from framework patches.")
    r.assessment(a)


# ---------------------------------------------------------------- phase 7
def phase_lineage(r: Report, homology_path):
    r.phase("Lineage", assessed=bool(homology_path))
    if not homology_path:
        r.skipped("no homology section supplied (Stage 3 has no script — write the section "
                  "from the search output and pass it with `--homology`)")
        return
    r.findings()
    r.out.append(Path(homology_path).read_text().rstrip())
    r.out.append("")
    r.assessment([
        "Homology text is passed through as written — this phase is the one place the report "
        "is not machine-derived, so the search engine (MMseqs2 or EBI BLAST) and the hit "
        "metrics must be stated in the section above.",
    ])


# ---------------------------------------------------------------- phase 8
def phase_verdict(r: Report):
    r.phase("Verdict")
    order = {"high": 0, "medium": 1, "low": 2}
    risks = sorted(r.risks, key=lambda x: order.get(x["severity"], 9))
    r.findings()
    if risks:
        r.out.append("| # | Severity | Where | Finding | Suggested action |")
        r.out.append("|---|----------|-------|---------|------------------|")
        for i, risk in enumerate(risks, 1):
            r.out.append(f"| {i} | {risk['severity']} | {risk['where']} | "
                         f"{risk['what']} | {risk['action']} |")
    else:
        r.out.append("No risks were raised by the phases that ran.")
    r.out.append("")

    a = []
    n_high = sum(1 for x in risks if x["severity"] == "high")
    n_med = sum(1 for x in risks if x["severity"] == "medium")
    if r.blocking:
        call = "**STOP — construct error.** The architecture cannot assemble as designed; fix "
        call += "the construct before any further work."
    elif n_high:
        call = (f"**ENGINEER** — {n_high} high-severity item(s) to resolve before committing to "
                "this sequence.")
    elif n_med:
        call = (f"**ENGINEER / ACCEPT WITH MONITORING** — no high-severity findings; {n_med} "
                "medium item(s) worth a decision.")
    else:
        call = "**PROCEED** — nothing above low severity in the phases that ran."
    a.append(call)

    top = [x for x in risks if x["severity"] == "high"][:5]
    if top:
        a.append("Top engineering targets, in order: "
                 + "; ".join(f"{i}) {x['where']} — {x['what']}" for i, x in enumerate(top, 1)) + ".")
    if r.not_assessed:
        a.append("**Not assessed in this run:** " + _fmt_list(r.not_assessed)
                 + ". The verdict covers only the phases above it.")
    a.append("All assessments are rule-derived from the stage JSON by "
             "`compile_report.py`; none are measurements. Nothing here substitutes for "
             "forced-degradation, DSF, SEC or binding data.")
    r.assessment(a)


# ---------------------------------------------------------------- provenance
def provenance(r: Report, args, inputs: dict):
    r.out.append("\n---\n")
    r.out.append("## Provenance\n")
    r.out.append(f"- Generated: {datetime.datetime.now().astimezone().isoformat(timespec='seconds')}")
    r.out.append(f"- Compiler: `compile_report.py` (phase-report-v1)")
    for key, paths in inputs.items():
        if key.startswith("_") or not paths:
            continue
        r.out.append(f"- {key}: {_fmt_list(Path(str(p)).name for p in paths)}")
    notes = []
    props = inputs.get("_props_data")
    if isinstance(props, dict):
        cx = props.get("complex") or {}
        if cx.get("stoichiometry"):
            notes.append("complex stoichiometry "
                         + ", ".join(f"{k}x{v}" for k, v in cx["stoichiometry"].items()))
        elif cx:
            notes.append("complex stoichiometry assumed 1 per chain")
    patches = inputs.get("_patch_data")
    if isinstance(patches, dict) and patches:
        notes.append(f"patch params: exposure ≥ {patches.get('exposure_cutoff')}, "
                     f"cluster ≤ {patches.get('cluster_cutoff_A')} Å, "
                     f"min size {patches.get('min_size')}, "
                     f"His included: {patches.get('include_his')}")
    cdrs = inputs.get("_cdr_data")
    if cdrs:
        engines = sorted({d.get("engine", "?") for d in cdrs})
        notes.append(f"CDR engine(s): {_fmt_list(engines)}")
    muts = inputs.get("_mut_data")
    if muts:
        notes.append(f"Fc reference isotype(s): {_fmt_list(sorted({d.get('isotype') for d in muts if d.get('isotype')}))}")
    for n in notes:
        r.out.append(f"- {n}")
    r.out.append("")


def main() -> None:
    p = argparse.ArgumentParser(
        description="Compile a phased Markdown evaluation report from stage JSON")
    p.add_argument("--title", default="Sequence Evaluation Report")
    p.add_argument("--properties")
    p.add_argument("--cdrs", nargs="*", default=[])
    p.add_argument("--mutations", nargs="*", default=[])
    p.add_argument("--liabilities")
    p.add_argument("--patches", help="JSON from scan_surface_patches.py")
    p.add_argument("--homology", help="Path to a pre-written Markdown homology section")
    args = p.parse_args()

    props = load(args.properties) if args.properties else None
    cdrs = [load(x) for x in args.cdrs]
    muts = [load(x) for x in args.mutations]
    liab = load(args.liabilities) if args.liabilities else None
    patches = load(args.patches) if args.patches else None

    r = Report(args.title)
    r.out.append(f"# {args.title}\n")
    r.out.append("Phased evaluation: each phase reports its findings and an assessment "
                 "derived from them. Phases without input say so rather than inferring.\n")

    phase_triage(r, props, cdrs, muts)
    phase_properties(r, props)
    phase_domains(r, cdrs)
    phase_engineering(r, muts)
    phase_liabilities(r, liab)
    phase_patches(r, patches)
    phase_lineage(r, args.homology)
    phase_verdict(r)
    provenance(r, args, {
        "Properties": [args.properties] if args.properties else [],
        "CDRs": args.cdrs,
        "Mutations": args.mutations,
        "Liabilities": [args.liabilities] if args.liabilities else [],
        "Patches": [args.patches] if args.patches else [],
        "Homology": [args.homology] if args.homology else [],
        "_props_data": props,
        "_patch_data": patches,
        "_cdr_data": cdrs,
        "_mut_data": muts,
    })

    sys.stdout.write("\n".join(r.out) + "\n")


if __name__ == "__main__":
    main()
