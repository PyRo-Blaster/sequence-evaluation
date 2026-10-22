# Trastuzumab Fv — phased evaluation

Phased evaluation: each phase reports its findings and an assessment derived from them. Phases without input say so rather than inferring.


## Phase 1 — Molecule Triage

**Findings**

- Chains in input: **2** — VL_A, VH_B
- Variable domains annotated: **2** (1 heavy-chain or VHH, 1 light-chain)


**Assessment**

- Modality read: **monospecific antibody (one VH + one VL)**.
- No constant region supplied, so this is an Fv/scFv-level assessment: Fc engineering and isotype are **not assessed**.


## Phase 2 — Chain Properties

**Findings**

| Chain | Precursor | Mature | MW (Da) | pI | A280 (1 g/L, ox) | Instability | GRAVY |
|-------|-----------|--------|---------|----|------------------|-------------|-------|
| VL_A | 110 | 110 | 12,052.35 | 9.01 | 1.208 | 51.6 | -0.366 |
| VH_B | 120 | 120 | 13,164.53 | 8.03 | 2.699 | 32.7 | -0.305 |

**Assembled complex:** MW 25,212.84 Da (2 disulfide bonds, stoichiometry VL_Ax1, VH_Bx1), pI 8.63


**Assessment**

- Instability index above 40 for VL_A (52). The index is a sequence heuristic, not a stability measurement — treat it as a prompt for a thermal/accelerated study, not a finding.
- A280 is the 1 g/L oxidized value (ExPASy convention); use it directly for concentration determination.


## Phase 3 — Domain Annotation

**Findings**

| Domain | Type | Engine | CDR1 | CDR2 | CDR3 | CDR3 len |
|--------|------|--------|------|------|------|----------|
| VL_A | VL | anchor-regex | `RASQDVNTAVA` | `SASFLYS` | `QQHYTTPPT` | 9 |
| VH_B | VH | anchor-regex | `GFNIKDTYIH` | `IYPTNGYTRYADSVKG` | `WGGDGFYAMDY` | 11 |


**Assessment**

- 2 of 2 domain(s) annotated with the **anchor-regex fallback**, not ANARCI. CDR1 and CDR2 are reliable; the CDR3 N-terminal boundary may be ±1. Install `abnumber` for true scheme numbering before quoting these boundaries externally.


## Phase 4 — Engineering Audit

_Not assessed — no find_mutations.py JSON supplied (no constant region analysed)._


## Phase 5 — Sequence Liabilities

**Findings**


**VL_A** (110 aa)


| Pos | Severity | Motif | CDR | Liability |
|-----|----------|-------|-----|-----------|
| 30 | high\* | `NT` | CDR-L1 | Deamidation (NS/NT/NN/NH) |

_2 further low-severity hit(s) (unescalated Met/Trp oxidation) omitted from this table._

**VH_B** (120 aa)

- N-terminal E (pyroglutamate-prone)

| Pos | Severity | Motif | CDR | Liability |
|-----|----------|-------|-----|-----------|
| 31 | high\* | `DT` | CDR-H1 | Isomerization (DG/DS/DD/DH/DT) |
| 55 | high | `NG` | CDR-H2 | Deamidation (NG) |
| 62 | high\* | `DS` | CDR-H2 | Isomerization (DG/DS/DD/DH/DT) |
| 99 | high\* | `W` | CDR-H3 | Trp oxidation |
| 102 | high\* | `DG` | CDR-H3 | Isomerization (DG/DS/DD/DH/DT) |
| 107 | high\* | `M` | CDR-H3 | Met oxidation |
| 73 | medium | `DT` | — | Isomerization (DG/DS/DD/DH/DT) |
| 77 | medium | `NT` | — | Deamidation (NS/NT/NN/NH) |
| 84 | medium | `NS` | — | Deamidation (NS/NT/NN/NH) |
| 90 | medium | `DT` | — | Isomerization (DG/DS/DD/DH/DT) |

_4 further low-severity hit(s) (unescalated Met/Trp oxidation) omitted from this table._

\* escalated one step because the residue lies in a CDR (two steps for oxidation on a contact-bearing CDR).


**Assessment**

- **7 high-severity liability/liabilities inside a CDR**: VL_A CDR-L1 pos 30 (Deamidation), VH_B CDR-H1 pos 31 (Isomerization), VH_B CDR-H2 pos 55 (Deamidation), VH_B CDR-H2 pos 62 (Isomerization), VH_B CDR-H3 pos 99 (Trp oxidation), VH_B CDR-H3 pos 102 (Isomerization), VH_B CDR-H3 pos 107 (Met oxidation). These are the primary engineering targets — a deamidation or isomerization motif in a CDR can erode potency on storage.
- N-terminal Gln/Glu on VH_B — pyroglutamate formation is usually benign but shows up as a charge variant; expect it in the CEX profile.


## Phase 6 — Structural Developability

**Findings**

Method `surface-patch-v1` on chains A, B of `1n8z.cif`. This is **not** TAP PPC/PNC/PSH and **not** SAP; do not apply those cutoffs.

Chains C were present in the file and excluded before the SASA calculation, so exposure is that of the isolated antibody chains (apo), not the bound complex.

| Severity | Kind | n | CDR res | Score | Residues |
|----------|------|---|---------|-------|----------|
| medium | negative_charge | 2 | 1 | 1.320 | A:ASP1, B:ASP62 [CDR-H2] |
| medium | hydrophobic | 2 | 1 | 1.399 | A:TYR49, A:PHE53 [CDR-L2] |
| low | positive_charge | 2 | 0 | 0.638 | A:LYS39, A:LYS42 |
| low | positive_charge | 2 | 0 | 0.515 | A:LYS190, A:ARG211 |
| low | negative_charge | 2 | 0 | 0.559 | A:ASP167, A:ASP170 |
| low | hydrophobic | 2 | 0 | 2.032 | B:VAL5, B:ALA23 |
| low | hydrophobic | 2 | 0 | 1.734 | B:LEU11, B:ALA121 |


**Assessment**

- Medium negative charge patch touching a CDR (2 residues, score 1.32): A:ASP1, B:ASP62. Below the high bar only because it has fewer than three residues.
- Medium hydrophobic patch touching a CDR (2 residues, score 1.40): A:TYR49, A:PHE53. Below the high bar only because it has fewer than three residues.
- Patch membership is read from **one static model**; a predicted CDR-H3 can move it. Treat these as review flags for this model, not measurements of the molecule.


## Phase 7 — Lineage

**Findings**

Search run with **MMseqs2** (ColabFold API) against UniRef; one search per variable domain.

| Domain | Top hit | Identity | Coverage | E-value |
|--------|---------|----------|----------|---------|
| VH_B | Trastuzumab heavy chain variable region | 100% | 100% | 1e-78 |
| VL_A | Trastuzumab light chain variable region | 100% | 100% | 1e-72 |

Both domains are exact matches to the humanized 4D5 framework.


**Assessment**

- Homology text is passed through as written — this phase is the one place the report is not machine-derived, so the search engine (MMseqs2 or EBI BLAST) and the hit metrics must be stated in the section above.


## Phase 8 — Verdict

**Findings**

| # | Severity | Where | Finding | Suggested action |
|---|----------|-------|---------|------------------|
| 1 | high | VL_A CDR-L1 pos 30 | Deamidation (NS/NT/NN/NH) (motif NT) | engineer out, or justify with a forced-degradation study |
| 2 | high | VH_B CDR-H1 pos 31 | Isomerization (DG/DS/DD/DH/DT) (motif DT) | engineer out, or justify with a forced-degradation study |
| 3 | high | VH_B CDR-H2 pos 55 | Deamidation (NG) (motif NG) | engineer out, or justify with a forced-degradation study |
| 4 | high | VH_B CDR-H2 pos 62 | Isomerization (DG/DS/DD/DH/DT) (motif DS) | engineer out, or justify with a forced-degradation study |
| 5 | high | VH_B CDR-H3 pos 99 | Trp oxidation (motif W) | engineer out, or justify with a forced-degradation study |
| 6 | high | VH_B CDR-H3 pos 102 | Isomerization (DG/DS/DD/DH/DT) (motif DG) | engineer out, or justify with a forced-degradation study |
| 7 | high | VH_B CDR-H3 pos 107 | Met oxidation (motif M) | engineer out, or justify with a forced-degradation study |
| 8 | medium | Fv surface | negative_charge patch touching a CDR (2 residues) | watch in viscosity and self-interaction screening |
| 9 | medium | Fv surface | hydrophobic patch touching a CDR (2 residues) | watch in viscosity and self-interaction screening |
| 10 | low | chain VL_A | instability index 52 (>40) | confirm with DSF/accelerated stability rather than acting on the index |
| 11 | low | annotation | CDR boundaries from the regex fallback | re-run with abnumber/ANARCI before publishing CDR definitions |


**Assessment**

- **ENGINEER** — 7 high-severity item(s) to resolve before committing to this sequence.
- Top engineering targets, in order: 1) VL_A CDR-L1 pos 30 — Deamidation (NS/NT/NN/NH) (motif NT); 2) VH_B CDR-H1 pos 31 — Isomerization (DG/DS/DD/DH/DT) (motif DT); 3) VH_B CDR-H2 pos 55 — Deamidation (NG) (motif NG); 4) VH_B CDR-H2 pos 62 — Isomerization (DG/DS/DD/DH/DT) (motif DS); 5) VH_B CDR-H3 pos 99 — Trp oxidation (motif W).
- **Not assessed in this run:** Phase 4 — Engineering Audit. The verdict covers only the phases above it.
- All assessments are rule-derived from the stage JSON by `compile_report.py`; none are measurements. Nothing here substitutes for forced-degradation, DSF, SEC or binding data.


---

## Provenance

- Generated: 2026-10-21T11:34:28+08:00
- Compiler: `compile_report.py` (phase-report-v1)
- Properties: props.json
- CDRs: vl.json, vh.json
- Liabilities: liab_full.json
- Patches: patches.json
- Homology: homology.md
- complex stoichiometry VL_Ax1, VH_Bx1
- patch params: exposure ≥ 0.2, cluster ≤ 7.0 Å, min size 2, His included: False
- CDR engine(s): anchor-regex

