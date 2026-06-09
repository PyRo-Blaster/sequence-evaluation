# Structural Prediction & Interface Analysis

## Overview

After sequence-level analysis, co-fold antibody-antigen complexes to predict
binding modes, paratopes, and epitopes. The primary tool is **Chai-1** via Modal GPU
(same setup used in the T002 trispecific antibody workflow). AlphaFold3 or
Boltz-1 are alternatives.

## Choosing a backend (run Stage 0 first)

Cofolding is GPU-required. Run `scripts/detect_resources.py` and route the job:

| Local resource | Recommended backend |
|----------------|--------------------|
| GPU with ~16 GB+ VRAM | Run Chai-1 / Boltz locally |
| No/small GPU, Modal CLI configured | Offload to Modal (`chai` + `modal` skills) |
| No GPU, no Modal, network available | Web server: AlphaFold3 server, HelixFold-Multimer |
| None of the above | Report `blocked`; do not attempt locally (will OOM/hang) |

Web servers impose limits (e.g. AlphaFold3 server ~20 jobs/day) and send the
sequence to a third party — avoid for sensitive sequences; prefer Modal or local.
For unbound Fv/VHH modeling (the `ab_structure` class) IgFold/ImmuneBuilder run
acceptably on CPU for a handful of sequences; reserve GPU/Modal for cofolding.

## Chai-1 via Modal

Chai-1 is a multi-chain structure prediction model that accepts multi-sequence
FASTA and produces ranked mmCIF models with per-chain confidence scores.

### Skill Dependencies

- `modal` skill — for GPU server setup
- `chai` skill — for Chai-1-specific guidance and `modal_chai1.py`

### FASTA Preparation

Each co-folding run requires one FASTA with all chains labeled by role:

```fasta
>antibody_VH|A
EVQLVES...
>antibody_VL|B
EIVLTQ...
>antigen|C
MKVLY...
```

Typical pairings for an antibody-antigen analysis:
- **Nanobody (VHH) vs antigen**: 2 chains (VHH + antigen)
- **Fab vs antigen**: 3 chains (VH + VL + antigen)
- **Full IgG vs antigen**: 4+ chains

### Running Predictions

Run all complexes in parallel using `subprocess.Popen`:

```python
import subprocess

runs = [
    ("vhh_antigen_A", "fastas/vhh_antigen_A.fasta", "predictions/vhh_antigen_A"),
    ("fab_antigen_B", "fastas/fab_antigen_B.fasta", "predictions/fab_antigen_B"),
]

processes = []
for name, fasta, outdir in runs:
    cmd = ["modal", "run", "modal_chai1.py",
           "--input-faa", fasta,
           "--out-dir", outdir]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    processes.append((name, p))

for name, p in processes:
    stdout, stderr = p.communicate()
    print(f"--- {name}: exit {p.returncode} ---")
```

### Output Files

Chai-1 writes 5 ranked models per run (indices 0–4):
```
predictions/
└── vhh_antigen_A/
    └── 2606071616/
        ├── pred.model_idx_0.cif   ← best model
        ├── pred.model_idx_1.cif
        ├── ...
        └── scores.json            ← pTM, ipTM per model
```

### Confidence Score Interpretation

| Metric | Meaning | Threshold |
|--------|---------|-----------|
| pTM    | Overall fold confidence | >0.5 confident |
| ipTM   | Interface confidence (binding) | >0.6 high confidence; 0.4–0.6 moderate |
| Aggregate | 0.2×pTM + 0.8×ipTM | Use for ranking models |

- **ipTM < 0.4**: Low confidence — binding mode unreliable
- **ipTM 0.4–0.6**: Moderate — GPCRs and membrane proteins typically fall here (membrane absent)
- **ipTM > 0.6**: High confidence — interface geometry likely accurate

## Interface Analysis with `analyze_interfaces.py`

After prediction, run `scripts/analyze_interfaces.py` to extract paratope/epitope residues.

```bash
uv run scripts/analyze_interfaces.py pred.model_idx_0.cif \
    --antibody A B \
    --antigen C \
    --cutoff 5.0
```

**Chain assignment** must match the FASTA used as input to Chai-1 (chain A = first sequence, B = second, etc.).

### Distance Cutoffs

| Cutoff | Use case |
|--------|----------|
| 4.0 Å  | Strict contacts (direct hydrogen bonds / salt bridges) |
| 5.0 Å  | Standard (recommended; captures VdW contacts) |
| 8.0 Å  | Extended neighborhood (for buried interface residues) |

### Interpreting the Output

- **Paratope residues**: CDR-derived contacts are expected. Non-CDR (framework) contacts indicate an unusual binding mode worth flagging.
- **Epitope residues**: Check if they overlap with known functional sites (receptor binding sites, active sites) — this suggests functional blockage.
- **Contact count**: Typical Fab-protein interfaces have 100–400 atom contacts. <50 may indicate a poor or artifactual prediction.

## AlphaFold3 / Boltz-1 Alternatives

If Chai-1/Modal is unavailable:

| Tool | Access | Notes |
|------|--------|-------|
| AlphaFold3 | alphafold.ebi.ac.uk (web) | 20 job/day limit; no GPU required |
| Boltz-1 | `boltz` skill | Local or Modal; open weights |
| ESMFold | `esm` skill | Monomer only, very fast |

AlphaFold3 produces `.cif` files compatible with `analyze_interfaces.py`. Boltz-1 uses the same FASTA multi-chain input format as Chai-1.
