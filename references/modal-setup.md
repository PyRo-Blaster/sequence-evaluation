# Modal Setup for Cloud GPU Cofolding (Chai-1 / Boltz)

The cofolding models are GPU-required. When `scripts/preflight.py` reports
`cofold_heavy -> modal` (or `modal-setup`), use the ready-made Modal apps in
`modal/`.

## One-time setup

```bash
pip install modal        # or: uv tool install modal
modal token new          # opens a browser to authenticate
```

Verify with the preflight — `modal CLI` should read `yes (authenticated)`:

```bash
uv run scripts/preflight.py
```

Modal's free tier includes monthly credit; beyond that GPU time is billed per
second, so cost scales with the number of complexes and sampling settings.

## Running a prediction

```bash
# Chai-1 (MSA-free, nothing leaves Modal)
modal run modal/modal_chai1.py \
    --input-faa fastas/vhh_antigen.fasta \
    --out-dir predictions/vhh_antigen

# Boltz
modal run modal/modal_boltz.py \
    --input-path fastas/fab_antigen.fasta \
    --out-dir predictions/fab_antigen
```

Outputs (mmCIF models plus score files) are streamed back and written to
`--out-dir` on your machine, ready for `scripts/analyze_interfaces.py`.

## Choosing a GPU

The GPU is read from the `MODAL_GPU` environment variable at import time:

```bash
MODAL_GPU=H100 modal run modal/modal_chai1.py --input-faa in.fasta --out-dir out/
```

| GPU | VRAM | Good for |
|-----|------|----------|
| `L4` | 24 GB | small complexes (VHH + small antigen); cheapest |
| `A10G` | 24 GB | small/medium complexes |
| `A100-40GB` | 40 GB | **default** — Fab + antigen, most antibody work |
| `A100-80GB` | 80 GB | large multi-chain complexes, full IgG |
| `H100` | 80 GB | fastest; large complexes |

Other knobs: `MODAL_TIMEOUT` (seconds, default 3600), `CHAI_VERSION`,
`BOLTZ_SPEC` (e.g. `boltz==2.1.1`).

## Weight caching

Each app mounts a persistent Modal Volume (`chai1-weights` / `boltz-weights`)
and points the model cache at it, so weights download once and later runs start
immediately. Inspect or clear them with:

```bash
modal volume ls chai1-weights
modal volume delete chai1-weights     # forces a fresh download next run
```

## Privacy

- **Chai-1** runs MSA-free using ESM embeddings — sequences stay within your
  Modal workspace.
- **Boltz** `--use-msa-server` is **off by default** here because it sends
  sequences to the public ColabFold MMseqs2 server. For confidential sequences
  leave it off and supply precomputed MSAs via a Boltz YAML input.

Prefer Modal or a local GPU over public web servers (AlphaFold3 server, SAbPred)
whenever the sequence is unpublished.

## Boltz-2 affinity prediction

Boltz-2 can predict binding affinity alongside structure when the YAML input
declares an affinity property. Pass the YAML to `--input-path`; the extra
affinity outputs are returned with the rest of the prediction files.

## Troubleshooting

| Symptom | Fix |
|---------|-----|
| `modal: command not found` | `pip install modal` |
| Auth / 401 errors | `modal token new` |
| CUDA OOM | raise `MODAL_GPU` (e.g. `A100-80GB`, `H100`) or reduce chains |
| Function timeout | raise `MODAL_TIMEOUT` |
| Slow first run | expected — weights are downloading into the Volume |
