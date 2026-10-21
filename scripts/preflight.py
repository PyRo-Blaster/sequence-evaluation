"""Environment preflight: verify tooling, Python packages, and GPU compute.

Run this FIRST, before any analysis or prediction stage. It answers three
questions in one pass:

  1. Tooling   - is `uv` present? is the Modal CLI installed and authenticated?
  2. Packages  - are the packages each stage needs available (and can uv fetch
                 the ones the scripts declare inline)?
  3. Compute   - is there a usable GPU, how much VRAM, how many cores / RAM?

It then recommends an execution backend per job class and prints concrete
remediation commands for anything missing.

Usage:
    uv run scripts/preflight.py                 # fast local check
    uv run scripts/preflight.py --json          # machine-readable
    uv run scripts/preflight.py --provision     # actually fetch core deps via uv
    uv run scripts/preflight.py --check-network # probe outbound egress

Exit code is 0 when every REQUIRED item is satisfied, 1 otherwise, so this can
gate a pipeline or a CI job.

# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""

import argparse
import importlib.util
import json
import os
import platform
import shutil
import socket
import subprocess
import sys
from pathlib import Path

# ---------------------------------------------------------------- packages ---
# tier: required  -> the CPU-light stages cannot run without it
#       antibody  -> improves antibody annotation accuracy (ANARCI numbering)
#       backend   -> needed to offload GPU work to Modal
#       local-gpu -> only needed to run cofolding ON THIS MACHINE
PACKAGES = [
    {"module": "Bio",       "pip": "biopython>=1.84", "tier": "required",
     "purpose": "MW/pI, pairwise alignment, mmCIF parsing"},
    {"module": "numpy",     "pip": "numpy>=1.26",     "tier": "required",
     "purpose": "interface distance math"},
    {"module": "abnumber",  "pip": "abnumber",        "tier": "antibody",
     "purpose": "ANARCI Kabat/IMGT/Chothia CDR numbering (else regex fallback); "
                "HEAVY: pulls anarcii+PyTorch, multi-GB first install"},
    {"module": "modal",     "pip": "modal",           "tier": "backend",
     "purpose": "cloud GPU for Chai-1 / Boltz cofolding"},
    {"module": "torch",     "pip": "torch",           "tier": "local-gpu",
     "purpose": "local GPU inference"},
    {"module": "chai_lab",  "pip": "chai_lab",        "tier": "local-gpu",
     "purpose": "run Chai-1 locally instead of on Modal"},
    {"module": "boltz",     "pip": "boltz",           "tier": "local-gpu",
     "purpose": "run Boltz locally instead of on Modal"},
]

# Core deps are declared inline by each script, so `uv run` provisions them
# automatically. --provision proves that actually works in this environment.
CORE_PROVISION = ["biopython>=1.84", "numpy>=1.26"]

VRAM_COFOLD_OK = 16   # GB; a Fab/IgG complex wants a large card
VRAM_AB_STRUCTURE_OK = 6


# ----------------------------------------------------------------- compute ---
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
        try:
            vram = round(float(parts[1]) / 1024, 1)
        except (IndexError, ValueError):
            vram = None
        gpus.append({"name": parts[0], "vram_gb": vram, "kind": "cuda"})
    return gpus


def detect_apple_mps() -> list[dict]:
    if platform.system() != "Darwin" or platform.machine() != "arm64":
        return []
    return [{"name": f"Apple Silicon ({platform.machine()})",
             "vram_gb": None, "kind": "mps"}]


def detect_ram_gb() -> float | None:
    try:
        if platform.system() == "Linux":
            for line in Path("/proc/meminfo").read_text().splitlines():
                if line.startswith("MemTotal:"):
                    return round(int(line.split()[1]) / (1024 ** 2), 1)
        elif platform.system() == "Darwin":
            out = subprocess.run(["sysctl", "-n", "hw.memsize"],
                                 capture_output=True, text=True, timeout=5)
            return round(int(out.stdout.strip()) / (1024 ** 3), 1)
    except (OSError, ValueError, subprocess.SubprocessError):
        pass
    try:
        return round(os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / (1024 ** 3), 1)
    except (ValueError, OSError, AttributeError):
        return None


def best_gpu(gpus: list[dict]) -> dict | None:
    if not gpus:
        return None
    return max(gpus, key=lambda g: (g["vram_gb"] is not None, g["vram_gb"] or 0))


# ----------------------------------------------------------------- tooling ---
def modal_authenticated() -> bool:
    """True if Modal credentials are present (token file or env vars)."""
    if os.environ.get("MODAL_TOKEN_ID") and os.environ.get("MODAL_TOKEN_SECRET"):
        return True
    cfg = Path(os.environ.get("MODAL_CONFIG_PATH", Path.home() / ".modal.toml"))
    try:
        return cfg.is_file() and "token_id" in cfg.read_text()
    except OSError:
        return False


def check_package(entry: dict, uv_present: bool) -> dict:
    """Is the module importable in THIS interpreter? Report version if so.

    Required deps are declared inline by each script, so `uv run` provisions
    them per-script even when they are absent from this interpreter. In that
    case report "auto (uv)" rather than a misleading "missing".
    """
    found, version = False, None
    try:
        spec = importlib.util.find_spec(entry["module"])
        found = spec is not None
    except (ImportError, ValueError):
        found = False
    if found:
        try:
            mod = __import__(entry["module"])
            version = getattr(mod, "__version__", None)
        except Exception:
            version = None
    if found:
        status = "installed"
    elif entry["tier"] == "required" and uv_present:
        status = "auto (uv)"
    else:
        status = "missing"
    return {**entry, "installed": found, "version": version, "status": status}


def provision_core(timeout: int = 600) -> dict:
    """Prove `uv` can actually fetch the inline core deps in this environment."""
    if not shutil.which("uv"):
        return {"ok": False, "detail": "uv not on PATH"}
    cmd = ["uv", "run", "--quiet"]
    for dep in CORE_PROVISION:
        cmd += ["--with", dep]
    cmd += ["python", "-c", "import Bio, numpy; print(Bio.__version__, numpy.__version__)"]
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return {"ok": False, "detail": f"timed out after {timeout}s (slow or blocked network?)"}
    except OSError as e:
        return {"ok": False, "detail": str(e)}
    if out.returncode == 0:
        return {"ok": True, "detail": out.stdout.strip()}
    return {"ok": False, "detail": (out.stderr or out.stdout).strip().splitlines()[-1:] or ["unknown error"]}


def check_network() -> bool | None:
    for host, port in (("1.1.1.1", 443), ("8.8.8.8", 53)):
        try:
            with socket.create_connection((host, port), timeout=2):
                return True
        except OSError:
            continue
    return False


# --------------------------------------------------------------- reasoning ---
def recommend(gpus, cpu_count, pkgs, modal_ready, network) -> dict:
    gpu = best_gpu(gpus)
    vram = gpu["vram_gb"] if gpu else None
    have = {p["module"]: p["installed"] for p in pkgs}
    local_cofold_pkg = have.get("chai_lab") or have.get("boltz")

    def remote(prefix: str) -> tuple[str, str]:
        if modal_ready:
            return "modal", f"{prefix}; Modal CLI installed and authenticated"
        if have.get("modal"):
            return "modal-setup", f"{prefix}; modal installed but not authenticated — run `modal token new`"
        if network:
            return "web-api", f"{prefix}; no Modal — use a hosted server (AlphaFold3 server, SAbPred)"
        if network is None:
            return "needs-backend", (f"{prefix}; no Modal configured and network not checked "
                                     "— re-run with --check-network, or set up Modal")
        return "blocked", f"{prefix}; no GPU, no Modal, no network egress"

    rec = {"cpu_light": {"backend": "local", "reason": "no GPU needed"}}

    if gpu and (vram is None or vram >= VRAM_AB_STRUCTURE_OK):
        rec["ab_structure"] = {"backend": "local", "reason": f"GPU available ({gpu['name']})"}
        rec["plm_embed"] = {"backend": "local", "reason": f"GPU available ({gpu['name']})"}
    elif cpu_count and cpu_count >= 4:
        rec["ab_structure"] = {"backend": "local-cpu",
                               "reason": f"no GPU; {cpu_count} cores — OK for a few sequences, slower"}
        rec["plm_embed"] = {"backend": "local-cpu", "reason": "no GPU; CPU fine for small batches"}
    else:
        b, why = remote("too few cores and no GPU")
        rec["ab_structure"] = {"backend": b, "reason": why}
        rec["plm_embed"] = {"backend": "local-cpu", "reason": "CPU fine for small batches"}

    if gpu and vram is not None and vram >= VRAM_COFOLD_OK and local_cofold_pkg:
        rec["cofold_heavy"] = {"backend": "local",
                               "reason": f"{vram} GB GPU ({gpu['name']}) with chai_lab/boltz installed"}
    elif gpu and vram is not None and vram >= VRAM_COFOLD_OK:
        b, why = remote(f"{vram} GB GPU present but neither chai_lab nor boltz installed locally")
        rec["cofold_heavy"] = {"backend": b, "reason": why}
    else:
        detail = (f"local GPU has {vram} GB (<{VRAM_COFOLD_OK} GB recommended)"
                  if (gpu and vram is not None) else "no/unsized local GPU")
        b, why = remote(detail)
        rec["cofold_heavy"] = {"backend": b, "reason": why}
    return rec


def remediation(blocked: list[str], modal: bool, network) -> list[dict]:
    """Concrete next steps for blocked job classes, best option first.

    Ported from detect_resources.py (PR #3) when the two tools were merged
    into this single preflight.
    """
    if not blocked:
        return []
    return [
        {
            "option": "Modal",
            "recommended": True,
            "action": "Enable Modal, then re-run this preflight with --modal yes",
            "where": "Settings -> Compute -> Modal",
            "note": ("Managed GPU; nothing to install locally. This is the configured "
                     "route for structure modelling: VHH -> Boltz-1, Fab/Fv -> Chai-1 "
                     "(see select_structure_model.py). Ready-made apps live in modal/."),
            "changes_model": False,
            "available_now": bool(modal),
        },
        {
            "option": "Local CUDA GPU",
            "recommended": False,
            "action": f"Install the model stack on a GPU with >= {VRAM_COFOLD_OK} GB VRAM",
            "where": "this machine, or an SSH host registered under Settings -> Compute",
            "note": "Keeps sequences entirely local — the right choice for sensitive material.",
            "changes_model": False,
            "available_now": False,
        },
        {
            "option": "Hosted web server",
            "recommended": False,
            "action": "Submit to AlphaFold3 server / HelixFold by hand",
            "where": "see references/structural-prediction.md",
            "note": ("Sends the sequence to a third party and imposes job quotas. "
                     "Do not use for sensitive sequences."),
            "changes_model": True,
            "available_now": bool(network),
        },
        {
            "option": "Antibody-specific CPU model",
            "recommended": False,
            "action": "Model the Fv with IgFold / ImmuneBuilder on local CPU",
            "where": "ab_structure job class — already routed local",
            "note": ("Unblocks the structure-dependent developability checks without a "
                     "GPU, but it is NOT the configured predictor: it replaces "
                     "Boltz-1/Chai-1 rather than giving them a backend."),
            "changes_model": True,
            "available_now": True,
        },
    ]


def action_items(uv_ok, pkgs, modal_ready, rec) -> list[str]:
    todo = []
    if not uv_ok:
        todo.append("Install uv (required to run every script): "
                    "curl -LsSf https://astral.sh/uv/install.sh | sh")
    for p in pkgs:
        if p["installed"]:
            continue
        if p["tier"] == "required" and p.get("status") != "auto (uv)":
            todo.append(f"REQUIRED {p['module']}: pip install '{p['pip']}' "
                        f"(or install uv so scripts self-provision)")
        elif p["tier"] == "antibody":
            todo.append(
                f"Optional {p['module']}: only needed for exact ANARCI numbering; "
                "the regex fallback handles standard and YFC-style frameworks. "
                "First install pulls anarcii+PyTorch (multi-GB, can take tens of "
                "minutes): uv run --with abnumber scripts/find_cdrs.py ...")
    if rec["cofold_heavy"]["backend"] in ("modal-setup", "blocked", "web-api"):
        if not any(p["module"] == "modal" and p["installed"] for p in pkgs):
            todo.append("Cloud GPU for Chai-1/Boltz: pip install modal  (then `modal token new`)")
        elif not modal_ready:
            todo.append("Modal installed but not authenticated: run `modal token new`")
    return todo


def main() -> None:
    ap = argparse.ArgumentParser(description="Environment preflight: tooling, packages, GPU")
    ap.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    ap.add_argument("--provision", action="store_true",
                    help="Actually fetch the inline core deps via uv to prove it works")
    ap.add_argument("--check-network", action="store_true", help="Probe outbound egress")
    ap.add_argument("--modal", choices=["auto", "yes", "no"], default="auto",
                    help="Whether a Modal backend is available. 'auto' looks for the "
                         "Modal CLI on PATH, which MISSES a host-managed Modal target "
                         "(Settings -> Compute) — pass 'yes' when one is enabled.")
    args = ap.parse_args()

    uv_path = shutil.which("uv")
    gpus = detect_nvidia() + detect_apple_mps()
    cpu_count = os.cpu_count()
    ram_gb = detect_ram_gb()
    pkgs = [check_package(p, uv_path is not None) for p in PACKAGES]
    modal_cli = shutil.which("modal") is not None
    # modal_source records HOW availability was established: an explicit
    # --modal yes is "asserted" even where the CLI also exists, because a
    # host-managed Modal target has no CLI at all.
    if args.modal == "yes":
        modal_ready, modal_source = True, "asserted"
    elif args.modal == "no":
        modal_ready, modal_source = False, "none"
    else:
        modal_ready = modal_cli and modal_authenticated()
        modal_source = "cli" if modal_ready else "none"
    network = check_network() if args.check_network else None
    provision = provision_core() if args.provision else None
    rec = recommend(gpus, cpu_count, pkgs, modal_ready, network)
    blocked = [job for job, info in rec.items()
               if info["backend"] in ("blocked", "needs-backend")]
    fixes = remediation(blocked, modal_ready, network)
    todo = action_items(uv_path is not None, pkgs, modal_ready, rec)

    required_ok = uv_path is not None and (provision["ok"] if provision else True)
    report = {
        "system": {"platform": platform.platform(), "python": platform.python_version(),
                   "cpu_count": cpu_count, "ram_gb": ram_gb},
        "modal_available": modal_ready,
        "modal_source": modal_source,
        "blocked_job_classes": blocked,
        "remediation": fixes,
        "tooling": {"uv": uv_path, "modal_cli": modal_cli, "modal_authenticated": modal_ready,
                    "nvidia_smi": shutil.which("nvidia-smi") is not None},
        "gpus": gpus,
        "packages": pkgs,
        "provision": provision,
        "network_egress": network,
        "recommendations": rec,
        "action_items": todo,
        "required_ok": required_ok,
    }

    if args.json:
        json.dump(report, sys.stdout, indent=2)
        sys.stdout.write("\n")
        sys.exit(0 if required_ok else 1)

    print("=== System ===")
    print(f"Platform : {report['system']['platform']}")
    print(f"Python   : {report['system']['python']}   CPU cores: {cpu_count}   "
          f"RAM: {ram_gb if ram_gb else '?'} GB")

    print("\n=== Tooling ===")
    print(f"uv            : {uv_path or 'NOT FOUND (required)'}")
    print(f"nvidia-smi    : {'yes' if report['tooling']['nvidia_smi'] else 'no'}")
    print(f"modal CLI     : {'yes' if modal_cli else 'no'}"
          f"   available: {'yes' if modal_ready else 'no'} ({modal_source})")
    if args.check_network:
        print(f"network egress: {'yes' if network else 'no'}")

    print("\n=== GPU ===")
    if gpus:
        for g in gpus:
            vram = f"{g['vram_gb']} GB" if g["vram_gb"] is not None else "unknown VRAM"
            print(f"  {g['name']} ({g['kind']}, {vram})")
    else:
        print("  none detected")

    print("\n=== Packages ===")
    print(f"  {'module':<10} {'tier':<10} {'status':<16} purpose")
    print(f"  {'-'*10} {'-'*10} {'-'*16} {'-'*40}")
    for p in pkgs:
        status = f"installed {p['version']}" if p["installed"] and p["version"] else p["status"]
        print(f"  {p['module']:<10} {p['tier']:<10} {status:<16} {p['purpose']}")
    print("  note: 'required' packages are declared inline per script and fetched")
    print("        automatically by `uv run`; --provision verifies that end to end.")

    if provision:
        print(f"\n=== Provision core deps via uv ===\n  "
              f"{'OK — ' + str(provision['detail']) if provision['ok'] else 'FAILED — ' + str(provision['detail'])}")

    print("\n=== Recommended backend per job class ===")
    for job, info in rec.items():
        print(f"  {job:<13} -> {info['backend']:<12} ({info['reason']})")

    if fixes:
        print("\n=== Remediation ladder (blocked: "
              + ", ".join(blocked) + ") ===")
        for f in fixes:
            mark = "*" if f["recommended"] else " "
            avail = "available" if f["available_now"] else "not available yet"
            warn = "  [changes the model]" if f["changes_model"] else ""
            print(f"  {mark} {f['option']:<28} {avail}{warn}")
            print(f"      {f['action']}  ({f['where']})")

    if todo:
        print("\n=== Action items ===")
        for t in todo:
            print(f"  - {t}")
    else:
        print("\nAll checks satisfied.")

    sys.exit(0 if required_ok else 1)


if __name__ == "__main__":
    main()
