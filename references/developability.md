# Developability Liabilities Reference

`scan_liabilities.py` flags sequence-level chemical and post-translational
liabilities. Severity is a default heuristic — solvent exposure and CDR location
(check against `find_cdrs.py` output) ultimately determine real-world risk.

| Liability | Motif | Default severity | Why it matters |
|-----------|-------|------------------|----------------|
| N-glycosylation sequon | `N-X-[ST]`, X≠P | high | N-linked glycosylation; expected at N297 in Fc, **unwanted** in Fab/CDRs (heterogeneity, altered binding) |
| Deamidation (fast) | `NG` | high | Asn→iso-Asp/Asp; charge change, potency loss; NG is the fastest motif |
| Deamidation (moderate) | `N[STNH]` | medium | Slower Asn deamidation |
| Isomerization | `D[GSDHT]` | medium | Asp→iso-Asp; backbone kink, especially in CDRs |
| Fragmentation | `DP` | medium | Acid-labile Asp-Pro bond; low-pH/formulation fragmentation |
| Met oxidation | `M` | low | Oxidation under stress/light; severity depends on exposure (Fc M252/M428 are sensitive) |
| Trp oxidation | `W` | low | Oxidation; high impact if in CDR (often a paratope residue) |
| Unpaired cysteine | odd Cys count | — | Free thiol → mispairing, aggregation, heterogeneity |
| N-terminal pyroglutamate | N-term `Q`/`E` | — | Pyroglutamate formation; usually benign but a charge-variant source |

## Severity escalation (implemented, not advisory)

`scan_liabilities.py` applies the CDR priority rule itself when given the
context. Build the map once with `build_cdr_map.py`, keyed by the FASTA record
ids being scanned:

```bash
uv run scripts/build_cdr_map.py HC=vh.json LC=vl.json > cdrmap.json
uv run scripts/scan_liabilities.py in.fasta --cdr-json cdrmap.json \
    [--paratope-json iface.json] --json
```

| Context | Effect on severity |
|---------|--------------------|
| Position falls inside a CDR | one step up the low→medium→high scale |
| Met/Trp oxidation in a CDR that carries a paratope contact | a second step up |

Each finding carries `severity` (escalated), `base_severity` (the table default),
`escalated`, `in_cdr`, `cdr` and `cdr_has_contact`, so a report can show both the
default and the context-adjusted call. Without `--cdr-json` nothing is escalated
and the scan says so explicitly rather than implying CDR context was considered.

The `--paratope-json` join is **CDR-level, not residue-level**: it marks which
CDRs carry a contact. Structure residue numbers and 1-based mature-sequence
positions are different coordinate systems, and the scan does not pretend to
reconcile them per residue.

## Interpretation guidance

- **CDR-localised liabilities are the priority.** A deamidation/isomerization or
  oxidation hotspot inside a CDR can erode potency on storage. This is now
  applied automatically by the escalation rules above when a CDR map is passed.
- **N-glyc sequons** in the variable domain are worth removing unless glycosylation
  is intended; the conserved Fc N297 sequon is expected and should *not* be removed
  unless aglycosylation is the design goal.
- **Met/Trp oxidation** flags are low severity by default because most surface Met/Trp
  are tolerated; weight them up if the residue is a known paratope contact (see the
  `analyze_interfaces.py` paratope output).
- **Odd cysteine count** almost always indicates an unpaired/engineered cysteine
  (e.g. a site-specific conjugation handle) — confirm it is intentional.

These are *flags for review*, not automatic rejections. Standard references:
Jarasch et al. 2015 (mAbs developability), Lu et al. 2019.

## Structure-based surface patches

Sequence motifs do not show whether charged or hydrophobic residues form a
surface patch. After a predicted Fv/Fab model exists, run
`scripts/scan_surface_patches.py` on the **antibody chains only** (never the
antigen) of the selected best model. Do not invent patches from sequence when
no structure is available.

The script implements `surface-patch-v1`. It is a review flag. It is **not**
Therapeutic Antibody Profiler PPC/PNC/PSH (Raybould et al. 2019) and **not**
spatial aggregation propensity (Chennamsetty et al. 2009). Do not apply TAP
amber/red cutoffs or SAP cutoffs to these scores.

| Patch | Members (default) | Score |
|-------|-------------------|-------|
| Positive charge | exposed Lys, Arg (formal +1). His excluded unless `--include-his` (+0.5) | Σ \|charge\| × rel_sasa |
| Negative charge | exposed Asp, Glu (formal −1) | Σ \|charge\| × rel_sasa |
| Hydrophobic | exposed Ala, Val, Leu, Ile, Met, Phe, Trp, Tyr | Σ weight × rel_sasa |

- Relative SASA = Shrake-Rupley residue SASA (probe 1.4 Å) / Tien et al. 2013
  theoretical maximum. Default surface cutoff is 0.20.
- Clustering uses side-chain centroid distance, default 7.0 Å, same sign or
  same hydrophobic class only. Reported patches have ≥2 residues.
- Hydrophobic weights are Kyte-Doolittle where that scale is positive. Trp/Tyr
  use fixed aromatic weights (1.0 / 0.8) because raw KD marks them hydrophilic.
- `--cdr-json` uses the same `{chain: [{name, seq}]}` shape as
  `analyze_interfaces.py`. Severity: **high** = ≥3 residues and ≥1 in a CDR;
  **medium** = ≥3 residues, or ≥2 with a CDR residue; otherwise **low**.
  Without CDR tags, nothing is called high.
- A model with no side-chain atoms is `coarse_model`. Do not report those calls.
- One static model, especially a predicted CDR-H3, can move patch membership.
  State that the flag is model-dependent. Do not average an ensemble unless
  the user asks.
