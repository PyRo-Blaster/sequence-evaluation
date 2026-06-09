---
name: sequence-evaluation
description: >
  Use when the user provides a protein FASTA file or amino acid sequence and wants
  sequence characterization, MW/pI calculation, antibody CDR identification,
  IgG Fc mutation mapping, homology search, antibody-antigen structure prediction,
  interface contact analysis, or a consolidated sequence evaluation report. Prefer
  for therapeutic antibodies, bispecifics, trispecifics, and nanobodies/VHH.
---

# Sequence Evaluation

End-to-end workflow for characterizing protein sequences, with a dedicated
antibody sub-workflow. Stages can be run independently or as a full pipeline.

All scripts accept `--json` for machine-readable output; pipe those into
`compile_report.py` (see Reporting) for a deterministic final report.

## Compute Resource Pre-flight (conditional — not every run)

**Skip this entirely for sequence-only work** (Stages 1–3, liabilities,
interface geometry — all CPU-light and always local). Run it **only when you are
about to start a GPU/ML stage** (cofolding, Fv structure modeling, language
models) **and only once per session** — cache the result and reuse it; the
hardware does not change between calls.

```bash
uv run scripts/detect_resources.py [--json]
```

It reports GPU (NVIDIA / Apple MPS), CPU, RAM, Modal CLI availability, and
network egress, then recommends a backend per job class:

| Job class | Tools | Routing rule |
|-----------|-------|--------------|
| `cpu_light` | properties, CDRs (regex), mutations, liabilities, interface geometry, ANARCI/abnumber, IgBLAST, BioPhi/OASis | **always local** — no GPU needed |
| `ab_structure` | IgFold, ImmuneBuilder, ABodyBuilder3, NanoBodyBuilder2, TAP | local GPU if present; else local CPU (slower, fine for a few sequences); else remote |
| `plm_embed` | AbLang, AntiBERTy, language-model scores | local GPU if present; CPU acceptable for small batches |
| `cofold_heavy` | Chai-1, Boltz-1/2, AlphaFold3, RFdiffusion/RFantibody, BindCraft | **GPU-required**: local GPU only if VRAM is sufficient (~16 GB+); else **Modal**; else a **web API/server** |

Decision flow for a GPU-required job:

1. **Local GPU with enough VRAM** → run locally.
2. **No/insufficient local GPU, Modal configured** → offload via the `modal` skill.
3. **No GPU, no Modal, network available** → use a hosted server (AlphaFold3
   server, SAbPred webapps, HelixFold web) — see `references/structural-prediction.md`.
4. **None of the above** → report `blocked`; tell the user what to enable rather
   than attempting a job that will OOM or hang.

Honour explicit user overrides (e.g. "force local" / "use Modal") over the
auto-recommendation, and prefer a slower-but-local path over sending sequences to
an external service when the data is sensitive.

## Stage 1: Physical Properties

Compute MW, pI, extinction coefficient / A280, and developability indices
(instability, aliphatic, GRAVY) per chain, plus the assembled complex.

```bash
uv run scripts/analyze_properties.py <input.fasta> \
    [--signal-peptide MGWSCIILFLVATATGVHS] \
    [--complex] \
    [--disulfide-bonds 18] \
    [--json]
```

- `--signal-peptide`: Strip N-terminal signal peptide before calculating mature chain properties
- `--complex`: Sum chains and solve for complex pI (binary search for net charge = 0)
- `--disulfide-bonds N`: Apply −2.016 Da correction per bond to complex MW
- Reports **A280 (1 g/L, oxidized)** for concentration measurement, plus
  instability index (>40 = potentially unstable), aliphatic index, and GRAVY.
- Inputs are validated: ambiguous/non-standard residues (X/B/Z/U…) fail loudly
  rather than crashing inside the MW/pI calculation.

## Stage 2: Antibody Domain Annotation

### CDR Identification

Run `find_cdrs.py` for each variable domain. Pass raw sequences (no signal peptide, no constant region):

```bash
uv run scripts/find_cdrs.py --vh  EVQLVES...  --name "VH_B"
uv run scripts/find_cdrs.py --vhh EVQLVES...  --name "VHH_A"
uv run scripts/find_cdrs.py --vl  DIQMTQ...   --name "VL_C"
uv run scripts/find_cdrs.py --seq EVQLVES...  --name "X"     # auto-detect chain type
```

**Two engines, auto-selected:**

1. **abnumber/ANARCI (preferred)** — true Kabat/IMGT/Chothia numbering. Used
   automatically when importable. Get it on the fly without a global install:
   ```bash
   uv run --with abnumber scripts/find_cdrs.py --vh EVQLVES... --name VH_B --scheme imgt
   ```
2. **Anchor regex (fallback)** — no dependency, Kabat-compatible boundaries. The
   output reports which engine ran. With the fallback, CDR1 and CDR2 are reliable;
   the CDR3 N-terminal boundary may be ±1 for some VH sequences.

`--scheme {kabat,imgt,chothia}` selects the numbering for the abnumber engine.
See `references/antibody-numbering.md` for scheme comparison and VHH-specific notes.

### Developability Liability Scan

Run `scan_liabilities.py` on each chain (or the whole FASTA) to flag
sequence-level chemical/PTM liabilities relevant to manufacturability:

```bash
uv run scripts/scan_liabilities.py <input.fasta> [--signal-peptide MGW...] [--json]
uv run scripts/scan_liabilities.py --seq EVQLVES... --name VH_A
```

Flags N-glycosylation sequons (N-X-S/T), Asn deamidation (NG/NS…), Asp
isomerization (DG/DS…), Asp-Pro fragmentation, Met/Trp oxidation hotspots, odd
(unpaired) cysteine counts, and N-terminal pyroglutamate. Severity-ranked;
positions are 1-based on the mature chain. Cross-reference high-severity hits in
CDRs against the Stage 2 CDR output — liabilities inside a CDR are the highest
priority to engineer out.

### Fc Mutation Mapping

Run `find_mutations.py` with each full heavy chain sequence to identify all
mutations vs the wild-type reference of its isotype, in EU numbering:

```bash
uv run scripts/find_mutations.py "MGWSCIILFLV...ASTKGPSVF..." --label "HC1" [--json]
uv run scripts/find_mutations.py "<seq>" --isotype IgG4              # force isotype
uv run scripts/find_mutations.py "<seq>" --reference-fasta wt.fasta  # custom WT
```

- Auto-detects the `ASTKGPSVF` (CH1 start) anchor, then **globally aligns** the
  constant region to the reference — insertions/deletions (engineered hinges,
  tags, des-K447) no longer cascade into spurious EU-shifted calls.
- **Isotype-aware**: auto-detects IgG1/IgG2/IgG4 by alignment and compares
  against that isotype's WT, so an unmodified IgG4 reports 0 mutations (not ~30).
  Override with `--isotype` if needed.
- EU numbering is reliable across CH1/CH2/CH3 and the conserved CPxCP hinge core
  (including the IgG4 **S228P** site); the exact EU numbers of the non-conserved
  N-terminal IgG2/IgG4 hinge residues are approximate.
- **Auto-annotates** recognised engineering mutations (LALA, LALA-PG, YTE, LS,
  KiH knob/hole, GASDALIE, SELF, IgG4 S228P…) and reports full vs partial variants.
- **Reference data caveat:** only the IgG1 reference is independently validated
  (it reproduces an approved IgG1 therapeutic). The IgG2/IgG4 references are
  reconstructed from canonical isotype differences. For definitive / regulatory
  work, pass an authoritative IMGT/UniProt sequence via `--reference-fasta`.

Consult `references/antibody-numbering.md` for the full mutation table.

## Stage 3: Homology Search

Use the `protein-sequence-similarity-search` skill for MMseqs2 (default, fast)
or EBI BLAST (fallback). Run one search per variable domain, ideally in parallel,
and capture the top hits, coverage, E-value, and percent identity for the final
report.

After searches complete, review the returned summaries. If the hit descriptions
are incomplete, resolve accessions with the `uniprot-database` skill.

## Stage 4: Structural Prediction & Interface Analysis

This is a `cofold_heavy` (GPU-required) job — run the compute pre-flight (once
per session) and route per its recommendation (local GPU / Modal / web server). See
`references/structural-prediction.md` for full guidance. Summary:

1. **Prepare FASTAs** — one FASTA per complex (VHH + antigen, Fab VH + VL + antigen, etc.)
2. **Run the cofolding model on the backend Stage 0 chose** — Chai-1 via Modal
   (the `chai` + `modal` skills) when offloading; locally if a capable GPU is
   present; or a web server (AlphaFold3 server, HelixFold) if no GPU/Modal
3. **Check scores** — select best model by aggregate score (0.2×pTM + 0.8×ipTM); ipTM > 0.5 is confident
4. **Analyze interface** with `analyze_interfaces.py`:

```bash
uv run scripts/analyze_interfaces.py pred.model_idx_0.cif \
    --antibody A B --antigen C --cutoff 5.0 \
    [--cdr-json cdrs.json] [--json]
```

Chain IDs (A, B, C…) correspond to the order of sequences in the input FASTA.

- Uses a **KD-tree** (`NeighborSearch`) — fast on full Fab/IgG complexes.
- Reports **per-residue-pair contacts** (count + minimum distance), not just a
  flat atom list.
- `--cdr-json` accepts a `{chain: [{name, seq}]}` map (assemble it from the
  `find_cdrs.py --json` outputs) to label each paratope residue with its CDR and
  flag framework-mediated contacts.

## Reporting

Run each stage with `--json`, save the outputs, then assemble them
deterministically with `compile_report.py`:

```bash
uv run scripts/analyze_properties.py in.fasta --complex --json > props.json
uv run scripts/find_cdrs.py --vh ... --name VH --json > vh.json
uv run scripts/find_mutations.py "..." --label HC1 --json > hc1.json
uv run scripts/scan_liabilities.py in.fasta --json > liab.json
uv run scripts/analyze_interfaces.py model.cif --antibody A B --antigen C --json > iface.json

uv run scripts/compile_report.py --title "mAb-X evaluation" \
    --properties props.json --cdrs vh.json vl.json \
    --mutations hc1.json --liabilities liab.json \
    --interfaces iface.json --homology homology_section.md > report.md
```

The compiler emits these sections (only for data provided):

1. **Chain Properties** — MW, pI, A280, instability
2. **CDR Annotation** — CDR1/2/3 per domain, with the engine used
3. **Fc Mutations** — EU-numbered table with functional annotation
4. **Developability Liabilities** — severity-ranked per chain
5. **Homology Search** — passed through as a Markdown section you write from the search results
6. **Structural Predictions** — paratope/epitope with CDR tags

Prepend a **Molecular Architecture** overview (chain composition, targets,
binding domains) and present in a concise lab-note style for downstream design
reviews. Homology search results (Stage 3) have no script, so write that section
as Markdown and pass it via `--homology`.

## Testing / Setup Check

`tests/smoke_test.sh` runs every script against the committed trastuzumab
fixtures in `examples/` and asserts key outputs. Use it to confirm `uv` and the
dependencies work in a fresh environment (e.g. from a SessionStart hook):

```bash
bash tests/smoke_test.sh
```

## Detecting Antibody vs General Protein

Run the antibody-specific stages (2 and 4) when the sequence contains:
- Conserved VH/VHH FR1 motifs: `EVQLVES`, `QVQLVES`, `EVQLVE`
- VL FR1 motifs: `DIQMTQ`, `EIVLTQ`, `DIVMTQ`
- IgG constant region anchor: `ASTKGPSVF` (CH1 start)

For general proteins without these markers, run stages 1 and 3 only.
