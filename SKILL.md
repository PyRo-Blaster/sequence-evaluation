---
name: sequence-evaluation
description: >
  Use when the user provides a protein FASTA file or amino acid sequence and wants
  sequence characterization, MW/pI calculation, antibody CDR identification,
  IgG Fc mutation mapping, homology search, antibody-antigen structure prediction,
  interface contact analysis, structure-based charged or hydrophobic surface
  patches from a predicted Fv model, or a consolidated sequence evaluation report.
  Prefer for therapeutic antibodies, bispecifics, trispecifics, and nanobodies/VHH.
---

# Sequence Evaluation

End-to-end workflow for characterizing protein sequences, with a dedicated
antibody sub-workflow. Stages can be run independently or as a full pipeline.

All scripts accept `--json` for machine-readable output; pipe those into
`compile_report.py` (see Reporting) for a deterministic final report.

## Stage 0: Environment Preflight (run first)

Run this **once at the start of a session**, before any analysis or prediction,
to confirm the environment can actually do the work. It checks tooling, Python
packages, and GPU compute in one pass, then recommends a backend per job class.

```bash
uv run scripts/preflight.py                  # fast check
uv run scripts/preflight.py --provision      # also verify uv can fetch core deps
uv run scripts/preflight.py --check-network  # also probe outbound egress
uv run scripts/preflight.py --modal yes      # assert a host-managed Modal target
uv run scripts/preflight.py --json           # machine-readable
```

What it verifies:

1. **Tooling** — `uv` on PATH (required: every script runs through it);
   `nvidia-smi`; Modal availability **and whether it is authenticated**.
2. **Packages** — per-package status and purpose. Core deps (`biopython`,
   `numpy`) are declared inline in each script and auto-provisioned by `uv run`,
   so they report `auto (uv)`; `--provision` proves that end to end. Optional
   (`abnumber`/`anarci`) and backend (`modal`, `torch`, `chai_lab`, `boltz`)
   packages are reported as installed/missing with the exact install command.
3. **Compute** — GPU (NVIDIA via `nvidia-smi`, Apple MPS) with VRAM, CPU cores, RAM.

Exit code is **0 when every required item is satisfied, 1 otherwise**, so it can
gate a pipeline or CI job. It prints an **Action items** list (a remediation
ladder when a job class is blocked) with the exact commands to fix anything missing.

`--modal auto` (the default) looks for the Modal CLI on PATH, which **misses a
host-managed Modal target** configured under *Settings → Compute* — that one has
no local CLI. After enabling Modal there, re-run with `--modal yes` and
`cofold_heavy` routes to `modal` instead of `blocked`.

**It checks local first** and only escalates when local cannot carry the job.
For sequence-only work the result is almost always `cpu_light -> local` and you
can proceed immediately; re-run it only if the environment changes.

Backend routing per job class:

| Job class | Tools | Routing rule |
|-----------|-------|--------------|
| `cpu_light` | properties, CDRs (regex), mutations, sequence liabilities, surface-patch scan on an existing model, interface geometry, ANARCI/abnumber, IgBLAST, BioPhi/OASis | **always local** — no GPU needed |
| `ab_structure` | IgFold, ImmuneBuilder, ABodyBuilder3, NanoBodyBuilder2, TAP | local GPU if present; else local CPU (slower, fine for a few sequences); else remote |
| `plm_embed` | AbLang, AntiBERTy, language-model scores | local GPU if present; CPU acceptable for small batches |
| `cofold_heavy` | Chai-1, Boltz-1/2, AlphaFold3, RFdiffusion/RFantibody, BindCraft | **GPU-required**: local GPU only if VRAM is sufficient (~16 GB+); else **Modal**; else a **web API/server** |

Decision flow for a GPU-required job:

1. **Local GPU with enough VRAM** (and `chai_lab`/`boltz` installed) → run locally.
2. **No/insufficient local GPU, Modal authenticated** → offload with the ready-made
   apps in `modal/` (see **Cloud GPU via Modal** below).
3. **`modal` installed but not authenticated** → preflight reports `modal-setup`;
   run `modal token new`.
4. **No GPU, no Modal, network available** → use a hosted server (AlphaFold3
   server, SAbPred webapps, HelixFold web) — see `references/structural-prediction.md`.
5. **None of the above** → report `blocked`; tell the user what to enable rather
   than attempting a job that will OOM or hang.

When a job class comes back `blocked`, the pre-flight prints a remediation
ladder — Modal (recommended), a local CUDA GPU, a hosted web server, or an
antibody-specific CPU model — each with the concrete action, where to do it, and
whether it **changes which model runs** rather than just where it runs. The last
two do. In `--json` the same list is under `remediation`, with
`blocked_job_classes` alongside.

**Nothing is enabled automatically.** Turning on a paid backend is the user's
decision; the pre-flight reports and stops. An agent cannot enable Modal —
`host.compute.create('modal')` raises `NotFound` until the user does.

Once Modal is the chosen backend, follow
`references/modal-structure-modeling.md`: which bundled env carries Chai-1 vs
Boltz, the GPU tier, the weights volume, the job-time egress each declares, and
the build-then-submit shape.

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
- `--stoichiometry NAME=N ...`: Copies of each chain in the assembly. **Pass this
  for any homomeric assembly** — without it every chain counts once, so an IgG
  FASTA (one HC + one LC) reports ~74 kDa instead of ~148 kDa. Applies to both
  the MW sum and the complex-pI charge balance.
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
   output reports which engine ran. CDR3 is anchored on the **conserved Cys**
   (Cys92 heavy / Cys88 light) and the J-region `WGxGT`/`FGxGT` motif rather than
   a literal `YYC`, so FR3 variants such as `VYFCAR`/`VYHCAR` are handled and the
   Kabat CDR3 start is exact.

`--scheme {kabat,imgt,chothia}` selects the numbering for the abnumber engine.
See `references/antibody-numbering.md` for scheme comparison and VHH-specific notes.

### Developability Liability Scan

Run `scan_liabilities.py` on each chain (or the whole FASTA) to flag
sequence-level chemical/PTM liabilities relevant to manufacturability:

```bash
uv run scripts/scan_liabilities.py <input.fasta> [--signal-peptide MGW...] \
    [--cdr-json cdrmap_seq.json] [--paratope-json iface.json] [--json]
uv run scripts/scan_liabilities.py --seq EVQLVES... --name VH_A
```

Flags N-glycosylation sequons (N-X-S/T), Asn deamidation (NG/NS…), Asp
isomerization (DG/DS…), Asp-Pro fragmentation, Met/Trp oxidation hotspots, odd
(unpaired) cysteine counts, and N-terminal pyroglutamate. Severity-ranked;
positions are 1-based on the mature chain.

**Pass `--cdr-json`.** The scan applies the CDR priority rule itself: a hit
inside a CDR moves one step up the low→medium→high scale, and Met/Trp oxidation
on a CDR that carries a paratope contact (`--paratope-json`) moves a second step.
Each finding carries both `severity` and `base_severity`, so nothing is hidden.
Keys in the CDR map must match the FASTA record ids being scanned — build it with
`build_cdr_map.py HC=vh.json LC=vl.json`. Without the map every severity is the
unescalated default and the output says so; do not do this join by hand.

The `--paratope-json` link is **CDR-level, not residue-level**: structure
residue numbers and 1-based mature-sequence positions are different coordinate
systems and the scan does not pretend to reconcile them per residue.

### Surface patches (needs a structure)

Sequence motifs cannot show a charged patch or a hydrophobic surface patch.
Run this **only after** a predicted Fv/Fab model exists (Stage 4, or a model
the user supplied). If there is no structure, report the sequence liabilities
and say surface patches are blocked — do not infer them from sequence.

Name the antibody chains in `--chains`. Any other chain in the file (the
antigen, a second Fab) is **detached before the SASA calculation**, so exposure
is always that of the isolated antibody chains. Passing a cofolded complex is
therefore safe — but the reported exposure is apo, not the bound state, which is
the correct substrate for a developability question: a patch matters because it
is exposed in the vial, not because it is buried against target.

```bash
uv run scripts/scan_surface_patches.py pred.model_idx_0.cif \
    --chains A B \
    --cdr-json cdrmap_chains.json \
    [--json]
```

`--cdr-json` is `{chain: [{name, seq}]}` — build it with `build_cdr_map.py`
rather than by hand, and use the same chain IDs as the cofolding FASTA. Method
is `surface-patch-v1` (Shrake-Rupley relative SASA, side-chain centroid
clustering). It reports positive-charge, negative-charge, and hydrophobic
patches. It is **not** TAP PPC/PNC/PSH and **not** SAP — do not apply those
cutoffs. Rules, weights, and severity are in `references/developability.md`.

- **high** only when the patch has ≥3 residues and at least one is in a CDR.
- A `coarse_model` warning means side chains are missing. Do not use that call.
- State that the flag is for the selected model, not an experimental structure.

The **conserved Fc glycan (EU N297)** is recognised by its `Q[YF]NST` context and
reported with severity `expected` rather than `high`, and is never escalated by
CDR proximity — it is required for FcγR binding and must not be engineered out
unless aglycosylation is the goal. Only *other* sequons count as liabilities.

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

## Stage 4: Structural Prediction

This is a `cofold_heavy` (GPU-required) job — run the compute pre-flight (once
per session) and route per its recommendation (local GPU / Modal / web server). See
`references/structural-prediction.md` and `references/modal-setup.md` for
full guidance.

### Model selection is by modality, not by taste

| Input | Model | Skill |
|-------|-------|-------|
| VHH / nanobody (single domain) | **Boltz-1** | `boltz` |
| Fab / Fv (VH + VL) | **Chai-1** | `chai` |

Resolve it with the script rather than by hand — it also stages the FASTA in the
chain order every downstream step depends on:

```bash
uv run scripts/select_structure_model.py --cdrs vh.json vl.json \
    --seqs VH=vh.seq VL=vl.seq --out-fasta fv.fasta [--antigen antigen.fasta] [--json]
```

It reports the modality, the chosen model, the skill to load, the chain-ID map
(A = first record, B = second, …) and the job class. A lone VH with no VL is
rejected rather than guessed at — pass `--vhh` to `find_cdrs.py` if it is a
single-domain antibody.

**Predict the apo Fv/VHH unless you specifically need the complex.** Surface
patches — the developability output — want the unbound antibody: a patch matters
because it is exposed in the vial, not because it is buried against target. An
antigen is only needed for interface work, which is no longer a report phase
(see below).

Steps:

1. **Stage the FASTA** with `select_structure_model.py`
2. **Run the chosen model** on the backend the pre-flight picked — Boltz-1 for VHH,
   Chai-1 for Fab — locally with a capable GPU, else via Modal, else a web server
3. **Check scores** — select the best model by aggregate score (0.2×pTM + 0.8×ipTM)
4. **Scan surface patches** on that selected model (Phase 6 of the report)

### Interface analysis (no longer a report phase)

`analyze_interfaces.py` remains available for paratope/epitope work when a
cofolded **antibody–antigen complex** exists, but it is deliberately not a
report phase: it needs a complex prediction, which is a materially heavier job
than the apo Fv the rest of the workflow runs on. Its output still feeds
`scan_liabilities.py --paratope-json` to escalate oxidation hotspots on
contact-bearing CDRs.

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

5. **Scan surface patches** on that same selected model. Pass antibody chains
   only — never the antigen. Command and severity rules are under Developability
   Liability Scan. If a model exists and developability was requested, do not
   stop after the interface table; if no model exists, do not invent patches.

### Cloud GPU via Modal (Chai-1 / Boltz)

When the preflight says `cofold_heavy -> modal`, use the ready-made apps in
`modal/` — no bespoke setup needed:

```bash
pip install modal && modal token new      # one time

modal run modal/modal_chai1.py \
    --input-faa fastas/vhh_antigen.fasta --out-dir predictions/vhh_antigen

modal run modal/modal_boltz.py \
    --input-path fastas/fab_antigen.fasta --out-dir predictions/fab_antigen
```

- Models run on a Modal GPU; mmCIF models and score files are streamed back to
  `--out-dir`, ready for `analyze_interfaces.py`.
- Weights are cached in a persistent Modal Volume, so only the first run pays
  the download.
- GPU is selected with `MODAL_GPU` (default `A100-40GB`; use `A100-80GB`/`H100`
  for large complexes, `L4`/`A10G` to cut cost). Also `MODAL_TIMEOUT`,
  `CHAI_VERSION`, `BOLTZ_SPEC`.
- **Privacy:** Chai-1 is MSA-free, so sequences stay in your Modal workspace.
  Boltz `--use-msa-server` is **off by default** because it sends sequences to
  the public ColabFold MMseqs2 server — keep it off for confidential sequences.

See `references/modal-setup.md` for GPU sizing, cost notes, and troubleshooting.

## Reporting — phased, with an opinion per phase

Run each stage with `--json`, save the outputs, then assemble them
deterministically with `compile_report.py`. Build the CDR map first: it is the
input that lets the liability scan, the patch scan and the interface analysis
all speak about CDRs instead of bare positions.

```bash
uv run scripts/analyze_properties.py in.fasta --complex --stoichiometry HC=2 LC=2 --json > props.json
uv run scripts/find_cdrs.py --vh ... --name VH --json > vh.json
uv run scripts/find_cdrs.py --vl ... --name VL --json > vl.json

# One map, three consumers. Key it by whatever the consumer calls the chain:
uv run scripts/build_cdr_map.py HC=vh.json LC=vl.json > cdrmap_seq.json    # FASTA record ids
uv run scripts/build_cdr_map.py A=vl.json  B=vh.json > cdrmap_chains.json  # model chain ids

uv run scripts/find_mutations.py "..." --label HC1 --json > hc1.json
uv run scripts/scan_liabilities.py in.fasta --cdr-json cdrmap_seq.json \
    [--paratope-json iface.json] --json > liab.json
uv run scripts/select_structure_model.py --cdrs vh.json vl.json \
    --seqs VH=vh.seq VL=vl.seq --out-fasta fv.fasta        # -> Chai-1 (Fab) or Boltz-1 (VHH)
# ... run that model on a GPU backend, then:
uv run scripts/scan_surface_patches.py model.cif --chains A B \
    --cdr-json cdrmap_chains.json --json > patches.json

uv run scripts/compile_report.py --title "mAb-X evaluation" \
    --properties props.json --cdrs vh.json vl.json \
    --mutations hc1.json --liabilities liab.json \
    --patches patches.json \
    --homology homology_section.md > report.md
```

The report walks the molecule from triage to binder. **Every phase emits
`**Findings**` (the data) and `**Assessment**` (what it means).** Phases with no
input print `Not assessed — <reason>` and are listed again in the Verdict, so a
missing input is visible rather than silently absent.

| Phase | Reads | Assessment answers |
|-------|-------|--------------------|
| 1 Molecule Triage | properties, CDRs, mutations | what modality is this, is the chain composition coherent |
| 2 Chain Properties | `analyze_properties.py` | anything outside the formulation/expression envelope |
| 3 Domain Annotation | `find_cdrs.py` | can these CDR boundaries be trusted, which engine produced them |
| 4 Engineering Audit | `find_mutations.py` | what the mutation set does; is the heterodimerization set coherent |
| 5 Sequence Liabilities | `scan_liabilities.py` | which liabilities are in CDRs, what to engineer first |
| 6 Structural Developability | `scan_surface_patches.py` | which patches touch a CDR, how model-dependent the call is |
| 7 Lineage | `--homology` Markdown | origin and closest known molecules |
| 8 Verdict | all of the above | ranked risk register and a STOP / ENGINEER / PROCEED call |

Phase 8 aggregates every risk raised by phases 1–7 into one table (severity,
where, finding, suggested action) and returns an explicit call:

- **STOP** — a construct error was found (e.g. knob without hole); fix before expression.
- **ENGINEER** — at least one high-severity item to resolve.
- **ENGINEER / ACCEPT WITH MONITORING** — only medium items.
- **PROCEED** — nothing above low severity *in the phases that ran*.

A **Provenance** block closes the report: timestamp, compiler version, every
input file, complex stoichiometry, patch parameters, CDR engine and Fc reference
isotype. Phase 1 replaces the hand-written Molecular Architecture overview —
chain composition and modality are now derived, not typed.

**The assessment layer is rule-derived and must stay that way.** The verdict
strings come from thresholds in `references/`; do not edit them into free prose,
and do not add a conclusion the JSON does not support. Homology (Phase 7) is the
one hand-written section — state the engine (MMseqs2 or EBI BLAST) and the hit
metrics there.

## Testing / Setup Check

`tests/smoke_test.sh` runs every script against the committed trastuzumab
fixtures in `examples/` and asserts key outputs. Use it to confirm `uv` and the
dependencies work in a fresh environment (e.g. from a SessionStart hook):

```bash
bash tests/smoke_test.sh
```

Assert report text on a token that markdown emphasis cannot split — `**not** TAP`
does not contain the substring `not TAP`.

## Detecting Antibody vs General Protein

Run the antibody-specific stages (2 and 4) when the sequence contains:
- Conserved VH/VHH FR1 motifs: `EVQLVES`, `QVQLVES`, `EVQLVE`
- VL FR1 motifs: `DIQMTQ`, `EIVLTQ`, `DIVMTQ`
- IgG constant region anchor: `ASTKGPSVF` (CH1 start)

For general proteins without these markers, run stages 1 and 3 only.
