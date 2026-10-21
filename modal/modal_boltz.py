"""Run Boltz (Boltz-1 / Boltz-2) structure prediction on a Modal cloud GPU.

Setup (once):
    pip install modal
    modal token new

Run (FASTA or Boltz YAML input):
    modal run modal/modal_boltz.py \
        --input-path fastas/fab_antigen.fasta \
        --out-dir predictions/fab_antigen

    MODAL_GPU=H100 modal run modal/modal_boltz.py \
        --input-path complexes/binder.yaml --out-dir predictions/binder

PRIVACY: --use-msa-server (default OFF here) sends your sequences to the public
ColabFold MMseqs2 server. Leave it off for confidential sequences and supply
precomputed MSAs in the YAML input instead.

Boltz-2 additionally predicts binding affinity when the YAML declares an
affinity property; see references/modal-setup.md.
"""

import os
from pathlib import Path

import modal

GPU = os.environ.get("MODAL_GPU", "A100-40GB")
TIMEOUT_S = int(os.environ.get("MODAL_TIMEOUT", "3600"))
# Pin a release for reproducibility; override to track Boltz-2 updates.
BOLTZ_SPEC = os.environ.get("BOLTZ_SPEC", "boltz==2.1.1")

WEIGHTS_DIR = "/weights"

image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("git", "wget")
    .pip_install(BOLTZ_SPEC)
    .env({"BOLTZ_CACHE": f"{WEIGHTS_DIR}/boltz"})
)

weights = modal.Volume.from_name("boltz-weights", create_if_missing=True)
app = modal.App("sequence-evaluation-boltz", image=image)


@app.function(gpu=GPU, volumes={WEIGHTS_DIR: weights}, timeout=TIMEOUT_S)
def predict(
    input_text: str,
    filename: str,
    use_msa_server: bool = False,
    diffusion_samples: int = 1,
    recycling_steps: int = 3,
) -> dict[str, bytes]:
    """Run `boltz predict`; return {relative_path: file_bytes}."""
    import shutil
    import subprocess

    workdir = Path("/tmp/boltz_run")
    if workdir.exists():
        shutil.rmtree(workdir)
    out_dir = workdir / "out"
    out_dir.mkdir(parents=True)

    input_path = workdir / filename
    input_path.write_text(input_text)

    cmd = [
        "boltz", "predict", str(input_path),
        "--out_dir", str(out_dir),
        "--cache", f"{WEIGHTS_DIR}/boltz",
        "--accelerator", "gpu",
        "--diffusion_samples", str(diffusion_samples),
        "--recycling_steps", str(recycling_steps),
        "--output_format", "mmcif",
    ]
    if use_msa_server:
        cmd.append("--use_msa_server")

    print("[boltz] " + " ".join(cmd), flush=True)
    proc = subprocess.run(cmd, capture_output=True, text=True)
    print(proc.stdout[-4000:], flush=True)
    if proc.returncode != 0:
        print(proc.stderr[-4000:], flush=True)
        raise RuntimeError(f"boltz predict failed (exit {proc.returncode})")

    weights.commit()

    return {
        str(p.relative_to(out_dir)): p.read_bytes()
        for p in sorted(out_dir.rglob("*"))
        if p.is_file()
    }


@app.local_entrypoint()
def main(
    input_path: str,
    out_dir: str,
    use_msa_server: bool = False,
    diffusion_samples: int = 1,
    recycling_steps: int = 3,
) -> None:
    src = Path(input_path)
    if src.suffix.lower() not in {".fasta", ".fa", ".faa", ".yaml", ".yml"}:
        raise SystemExit(f"{input_path}: expected a .fasta/.fa/.faa or .yaml/.yml input")

    if use_msa_server:
        print("[boltz] WARNING: --use-msa-server sends sequences to the public "
              "ColabFold MMseqs2 server. Do not use for confidential sequences.")

    print(f"[boltz] submitting {input_path} to Modal (GPU={GPU}, {BOLTZ_SPEC}) ...")
    results = predict.remote(
        src.read_text(),
        src.name,
        use_msa_server=use_msa_server,
        diffusion_samples=diffusion_samples,
        recycling_steps=recycling_steps,
    )

    dest = Path(out_dir)
    dest.mkdir(parents=True, exist_ok=True)
    for rel, blob in results.items():
        target = dest / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(blob)
    print(f"[boltz] wrote {len(results)} file(s) to {dest}")
    for rel in sorted(results):
        if rel.endswith((".cif", ".json")):
            print(f"  {rel}")
    print("\nNext: pick the top-ranked model, then run\n"
          f"  uv run scripts/analyze_interfaces.py {dest}/<best>.cif "
          "--antibody A B --antigen C --json")
