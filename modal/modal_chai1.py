"""Run Chai-1 complex structure prediction on a Modal cloud GPU.

Setup (once):
    pip install modal
    modal token new

Run:
    modal run modal/modal_chai1.py \
        --input-faa fastas/vhh_antigen.fasta \
        --out-dir predictions/vhh_antigen

Pick a GPU by setting MODAL_GPU before the call (default A100-40GB):
    MODAL_GPU=H100 modal run modal/modal_chai1.py --input-faa in.fasta --out-dir out/

Model weights are cached in a persistent Modal Volume, so only the first run
pays the download cost. Chai-1 runs MSA-free by default (ESM embeddings), so
sequences are NOT sent to any third-party MSA server.
"""

import os
from pathlib import Path

import modal

# The decorator is evaluated at import time, so the GPU is chosen via env var.
GPU = os.environ.get("MODAL_GPU", "A100-40GB")
TIMEOUT_S = int(os.environ.get("MODAL_TIMEOUT", "3600"))
CHAI_VERSION = os.environ.get("CHAI_VERSION", "0.6.1")

WEIGHTS_DIR = "/weights"

image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("git")
    .pip_install(f"chai_lab=={CHAI_VERSION}")
    # Keep downloaded weights on the Volume instead of the ephemeral container.
    .env({"CHAI_DOWNLOADS_DIR": f"{WEIGHTS_DIR}/chai"})
)

weights = modal.Volume.from_name("chai1-weights", create_if_missing=True)
app = modal.App("sequence-evaluation-chai1", image=image)


@app.function(gpu=GPU, volumes={WEIGHTS_DIR: weights}, timeout=TIMEOUT_S)
def fold(
    fasta_text: str,
    num_trunk_recycles: int = 3,
    num_diffn_timesteps: int = 200,
    seed: int = 42,
) -> dict[str, bytes]:
    """Fold one multi-chain FASTA; return {relative_path: file_bytes}."""
    import shutil

    from chai_lab.chai1 import run_inference

    workdir = Path("/tmp/chai_run")
    if workdir.exists():
        shutil.rmtree(workdir)
    out_dir = workdir / "out"
    out_dir.mkdir(parents=True)

    fasta_path = workdir / "input.fasta"
    fasta_path.write_text(fasta_text)

    run_inference(
        fasta_file=fasta_path,
        output_dir=out_dir,
        num_trunk_recycles=num_trunk_recycles,
        num_diffn_timesteps=num_diffn_timesteps,
        seed=seed,
        device="cuda:0",
        use_esm_embeddings=True,
    )

    # Persist newly downloaded weights for the next run.
    weights.commit()

    return {
        str(p.relative_to(out_dir)): p.read_bytes()
        for p in sorted(out_dir.rglob("*"))
        if p.is_file()
    }


@app.local_entrypoint()
def main(
    input_faa: str,
    out_dir: str,
    num_trunk_recycles: int = 3,
    num_diffn_timesteps: int = 200,
    seed: int = 42,
) -> None:
    fasta_text = Path(input_faa).read_text()
    if not fasta_text.lstrip().startswith(">"):
        raise SystemExit(f"{input_faa} does not look like a FASTA file")

    print(f"[chai-1] submitting {input_faa} to Modal (GPU={GPU}) ...")
    results = fold.remote(
        fasta_text,
        num_trunk_recycles=num_trunk_recycles,
        num_diffn_timesteps=num_diffn_timesteps,
        seed=seed,
    )

    dest = Path(out_dir)
    dest.mkdir(parents=True, exist_ok=True)
    for rel, blob in results.items():
        target = dest / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(blob)
    print(f"[chai-1] wrote {len(results)} file(s) to {dest}")
    for rel in sorted(results):
        if rel.endswith((".cif", ".json", ".npz")):
            print(f"  {rel}")
    print("\nNext: score the models, then run\n"
          f"  uv run scripts/analyze_interfaces.py {dest}/<best>.cif "
          "--antibody A B --antigen C --json")
