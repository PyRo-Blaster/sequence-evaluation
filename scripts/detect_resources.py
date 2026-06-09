"""Detect local compute resources and recommend an execution backend per job class.

Run this as a pre-flight before any structure-prediction or ML stage. It reports
GPU (NVIDIA via nvidia-smi, or Apple MPS), CPU, RAM, and whether the Modal CLI is
available, then maps the skill's job classes to a backend:

    local      - run on this machine (GPU if present)
    local-cpu  - run here on CPU (slower; acceptable for small inputs)
    modal      - offload to a Modal GPU (the `modal` skill)
    web-api    - use a hosted server / API (e.g. AlphaFold3 server, SAbPred)

Job classes:
    cpu_light    - properties, CDRs, mutations, liabilities, interface geometry
    ab_structure - Fv/VHH modeling (IgFold, ImmuneBuilder, ABodyBuilder3): GPU-preferred
    plm_embed    - antibody language models (AbLang/AntiBERTy/BioPhi): GPU-preferred, CPU-OK
    cofold_heavy - complex cofolding / diffusion (Chai-1, Boltz, AF3, RFdiffusion): GPU-required

Usage:
    uv run scripts/detect_resources.py [--json] [--no-network-check]

# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""

import argparse
import os
import platform
import shutil
import socket
import subprocess
import sys

# Rough VRAM (GB) guidance. Cofolding a real Fab/IgG complex wants a large card;
# small complexes fit in less. These are routing heuristics, not hard limits.
VRAM_COFOLD_OK = 16
VRAM_AB_STRUCTURE_OK = 6


def detect_nvidia() -> list[dict]:
    exe = shutil.which("nvidia-smi")
    if not exe:
        return []
    try:
        out = subprocess.run(
            [exe, "--query-gpu=name,memory.total", "--format=csv,noheader,nounits"],
            capture_output=True, text=True, timeout=10,
        )
    except (subprocess.SubprocessError, OSError):
        return []
    gpus = []
    for line in out.stdout.strip().splitlines():
        if not line.strip():
            continue
        parts = [p.strip() for p in line.split(",")]
        name = parts[0]
        try:
            vram_gb = round(float(parts[1]) / 1024, 1)
        except (IndexError, ValueError):
            vram_gb = None
        gpus.append({"name": name, "vram_gb": vram_gb, "kind": "cuda"})
    return gpus


def detect_apple_mps() -> list[dict]:
    if platform.system() != "Darwin" or platform.machine() != "arm64":
        return []
    # Apple Silicon: unified memory, MPS backend. Report as a GPU of unknown VRAM.
    return [{"name": f"Apple Silicon ({platform.processor() or 'arm64'})",
             "vram_gb": None, "kind": "mps"}]


def detect_ram_gb() -> float | None:
    try:
        if platform.system() == "Linux":
            with open("/proc/meminfo") as fh:
                for line in fh:
                    if line.startswith("MemTotal:"):
                        return round(int(line.split()[1]) / (1024 ** 2), 1)
        elif platform.system() == "Darwin":
            out = subprocess.run(["sysctl", "-n", "hw.memsize"],
                                 capture_output=True, text=True, timeout=5)
            return round(int(out.stdout.strip()) / (1024 ** 3), 1)
    except (OSError, ValueError, subprocess.SubprocessError):
        pass
    try:  # POSIX fallback
        return round(os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / (1024 ** 3), 1)
    except (ValueError, OSError, AttributeError):
        return None


def check_network() -> bool | None:
    """Best-effort egress check (short timeout). None if it can't be determined."""
    for host, port in (("1.1.1.1", 443), ("8.8.8.8", 53)):
        try:
            with socket.create_connection((host, port), timeout=2):
                return True
        except OSError:
            continue
    return False


def best_gpu(gpus: list[dict]) -> dict | None:
    if not gpus:
        return None
    # Prefer the card with the most known VRAM; unknown VRAM (MPS) sorts last.
    return max(gpus, key=lambda g: (g["vram_gb"] is not None, g["vram_gb"] or 0))


def recommend(gpus, cpu_count, network) -> dict:
    gpu = best_gpu(gpus)
    has_gpu = gpu is not None
    vram = gpu["vram_gb"] if gpu else None
    # Unknown VRAM (e.g. Apple MPS) is treated as "enough for ab_structure,
    # uncertain for cofold" — cofold heavy is routed remote to be safe.
    modal = shutil.which("modal") is not None

    def remote_choice(reason_local_no_gpu):
        if modal:
            return "modal", "no suitable local GPU; Modal CLI present"
        if network:
            return "web-api", "no local GPU and no Modal; use a hosted server/API"
        return "blocked", ("no local GPU, no Modal, no network egress — "
                           "GPU job cannot run; install locally or enable a backend")

    rec = {}
    rec["cpu_light"] = {"backend": "local", "reason": "no GPU needed"}

    if has_gpu and (vram is None or vram >= VRAM_AB_STRUCTURE_OK):
        rec["ab_structure"] = {"backend": "local", "reason": f"GPU available ({gpu['name']})"}
        rec["plm_embed"] = {"backend": "local", "reason": f"GPU available ({gpu['name']})"}
    elif cpu_count and cpu_count >= 4:
        rec["ab_structure"] = {"backend": "local-cpu",
                               "reason": f"no GPU; {cpu_count} cores — OK for a few sequences, slower"}
        rec["plm_embed"] = {"backend": "local-cpu", "reason": "no GPU; CPU fine for small batches"}
    else:
        b, why = remote_choice(True)
        rec["ab_structure"] = {"backend": b, "reason": why}
        rec["plm_embed"] = {"backend": "local-cpu", "reason": "CPU fine for small batches"}

    if has_gpu and vram is not None and vram >= VRAM_COFOLD_OK:
        rec["cofold_heavy"] = {"backend": "local",
                               "reason": f"GPU with {vram} GB VRAM ({gpu['name']})"}
    else:
        reason_detail = (f"local GPU has {vram} GB (<{VRAM_COFOLD_OK} GB recommended)"
                         if (has_gpu and vram is not None) else "no/unsized local GPU")
        b, why = remote_choice(True)
        rec["cofold_heavy"] = {"backend": b, "reason": f"{reason_detail}; {why}"}

    return rec


def main() -> None:
    ap = argparse.ArgumentParser(description="Detect compute resources and recommend backends")
    ap.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    ap.add_argument("--no-network-check", action="store_true",
                    help="Skip the outbound connectivity probe")
    args = ap.parse_args()

    gpus = detect_nvidia() + detect_apple_mps()
    cpu_count = os.cpu_count()
    ram_gb = detect_ram_gb()
    network = None if args.no_network_check else check_network()
    modal = shutil.which("modal") is not None
    rec = recommend(gpus, cpu_count, network)

    report = {
        "platform": platform.platform(),
        "cpu_count": cpu_count,
        "ram_gb": ram_gb,
        "gpus": gpus,
        "modal_cli": modal,
        "network_egress": network,
        "recommendations": rec,
    }

    if args.json:
        import json
        json.dump(report, sys.stdout, indent=2)
        sys.stdout.write("\n")
        return

    print("=== Local compute resources ===")
    print(f"Platform : {report['platform']}")
    print(f"CPU cores: {cpu_count}    RAM: {ram_gb if ram_gb else '?'} GB")
    if gpus:
        for g in gpus:
            vram = f"{g['vram_gb']} GB" if g["vram_gb"] is not None else "unknown VRAM"
            print(f"GPU      : {g['name']} ({g['kind']}, {vram})")
    else:
        print("GPU      : none detected")
    print(f"Modal CLI: {'yes' if modal else 'no'}    "
          f"Network  : {'yes' if network else ('no' if network is False else 'not checked')}")
    print("\n=== Recommended backend per job class ===")
    for job, info in rec.items():
        print(f"  {job:<13} -> {info['backend']:<10}  ({info['reason']})")
    if any(i["backend"] == "blocked" for i in rec.values()):
        print("\n[!] One or more GPU job classes cannot run here. Options: install a "
              "local GPU stack, configure the `modal` skill, or run the relevant "
              "web tool (e.g. AlphaFold3 server, SAbPred).")


if __name__ == "__main__":
    main()
