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

## Interpretation guidance

- **CDR-localised liabilities are the priority.** A deamidation/isomerization or
  oxidation hotspot inside a CDR can erode potency on storage. Map positions from
  this scan onto the CDR boundaries from `find_cdrs.py`.
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
