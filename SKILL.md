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

## Stage 1: Physical Properties

Compute MW and pI per chain, and for the assembled complex.

```bash
uv run scripts/analyze_properties.py <input.fasta> \
    [--signal-peptide MGWSCIILFLVATATGVHS] \
    [--complex] \
    [--disulfide-bonds 18]
```

- `--signal-peptide`: Strip N-terminal signal peptide before calculating mature chain properties
- `--complex`: Sum chains and solve for complex pI (binary search for net charge = 0)
- `--disulfide-bonds N`: Apply −2.016 Da correction per bond to complex MW

## Stage 2: Antibody Domain Annotation

### CDR Identification

Run `find_cdrs.py` for each variable domain. Pass raw sequences (no signal peptide, no constant region):

```bash
uv run scripts/find_cdrs.py --vh  EVQLVES...  --name "VH_B"
uv run scripts/find_cdrs.py --vhh EVQLVES...  --name "VHH_A"
uv run scripts/find_cdrs.py --vl  DIQMTQ...   --name "VL_C"
```

Uses conserved anchor residues (not fixed Kabat positions). CDR1 and CDR2 are reliable; CDR3 N-terminal boundary may be ±1 residue for some VH sequences and CDR1 may include/exclude a boundary residue depending on the numbering scheme. For precise Kabat/IMGT numbering, use ANARCI (`pip install abnumber` or `conda install -c bioconda anarci`). See `references/antibody-numbering.md` for scheme comparison and VHH-specific notes.

### Fc Mutation Mapping

Run `find_mutations.py` with each full heavy chain sequence to identify all
mutations vs WT IGHG1*01 in EU numbering:

```bash
uv run scripts/find_mutations.py "MGWSCIILFLV...ASTKGPSVF..." --label "HC1"
```

The script auto-detects the ASTKGPSVF constant region anchor. Consult
`references/antibody-numbering.md` for a table of common therapeutic mutations
(LALA-PG, KiH, LS, YTE, etc.) to annotate functional significance.

## Stage 3: Homology Search

Use the `protein-sequence-similarity-search` skill for MMseqs2 (default, fast)
or EBI BLAST (fallback). Run one search per variable domain, ideally in parallel,
and capture the top hits, coverage, E-value, and percent identity for the final
report.

After searches complete, review the returned summaries. If the hit descriptions
are incomplete, resolve accessions with the `uniprot-database` skill.

## Stage 4: Structural Prediction & Interface Analysis

See `references/structural-prediction.md` for full guidance. Summary:

1. **Prepare FASTAs** — one FASTA per complex (VHH + antigen, Fab VH + VL + antigen, etc.)
2. **Run Chai-1 in parallel** via Modal (use the `chai` skill together with the `modal` skill)
3. **Check scores** — select best model by aggregate score (0.2×pTM + 0.8×ipTM); ipTM > 0.5 is confident
4. **Analyze interface** with `analyze_interfaces.py`:

```bash
uv run scripts/analyze_interfaces.py pred.model_idx_0.cif \
    --antibody A B --antigen C --cutoff 5.0
```

Chain IDs (A, B, C…) correspond to the order of sequences in the input FASTA.

## Reporting

After running all stages, compile results into a Markdown report with these sections:

1. **Molecular Architecture** — chain composition, targets, binding domains
2. **Chain Properties** — table: chain name, precursor/mature length, MW, pI
3. **CDR Annotation** — one sub-section per variable domain; CDR1/2/3 sequences
4. **Fc Mutations** — EU-numbered mutation table with functional annotation
5. **Homology Search** — top 3–5 hits per domain; Q-Cov, E-value, % identity, inferred origin
6. **Structural Predictions** — pTM/ipTM confidence table; interface contact summary

Present the report in a concise lab-note style with clear tables and per-domain
subsections so that the output can be reused in downstream design reviews.

## Detecting Antibody vs General Protein

Run the antibody-specific stages (2 and 4) when the sequence contains:
- Conserved VH/VHH FR1 motifs: `EVQLVES`, `QVQLVES`, `EVQLVE`
- VL FR1 motifs: `DIQMTQ`, `EIVLTQ`, `DIVMTQ`
- IgG constant region anchor: `ASTKGPSVF` (CH1 start)

For general proteins without these markers, run stages 1 and 3 only.
