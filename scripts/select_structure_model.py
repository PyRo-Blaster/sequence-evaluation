"""Pick the structure-prediction backend for an antibody and stage its input.

Model choice is by modality, not by taste:

    VHH / nanobody (single domain) -> Boltz-1
    Fab / Fv (VH + VL)             -> Chai-1

Both are `cofold_heavy` jobs and require a GPU. This script does not run them;
it resolves the modality, writes the multi-chain FASTA in the right chain order,
and reports which skill to load and whether this machine can execute it. Pair it
with `detect_resources.py` for the backend decision (local GPU / Modal / web).

Chain order in the emitted FASTA is the contract for everything downstream:
chain A is the first record, B the second, and `analyze_interfaces.py` /
`scan_surface_patches.py` must be given those same IDs.

Usage:
    uv run scripts/select_structure_model.py --cdrs vh.json vl.json \
        --seqs VH=vh.seq VL=vl.seq --out-fasta fv.fasta
    uv run scripts/select_structure_model.py --cdrs vhh.json \
        --seqs VHH=vhh.seq --antigen antigen.fasta --out-fasta cplx.fasta --json

# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from _common import SequenceError, clean_sequence, emit_json

# modality -> (model, skill to load, note)
ROUTES = {
    "VHH": ("Boltz-1", "boltz",
            "single-domain input; Boltz-1 is the configured predictor for VHH/nanobody"),
    "FAB": ("Chai-1", "chai",
            "paired VH+VL input; Chai-1 is the configured predictor for Fab/Fv"),
}


def resolve_modality(domains: list[dict]) -> tuple[str, str]:
    """Return (modality_key, human description) from find_cdrs.py outputs."""
    types = [d.get("domain_type") for d in domains]
    vhh = [t for t in types if t == "VHH"]
    vh = [t for t in types if t == "VH"]
    vl = [t for t in types if t == "VL"]

    if vhh and not vl:
        return "VHH", f"{len(vhh)} VHH domain(s), no light chain"
    if vh and vl:
        return "FAB", f"{len(vh)} VH + {len(vl)} VL"
    if vh and not vl:
        raise SystemExit(
            "[!] a lone VH with no VL is not a modality this routing covers. "
            "Pass --vhh to find_cdrs.py if it is a single-domain antibody, or "
            "supply the partner light chain."
        )
    raise SystemExit(f"[!] could not resolve a modality from domain types {types}")


def main() -> None:
    p = argparse.ArgumentParser(
        description="Select the structure-prediction model for an antibody and stage its FASTA")
    p.add_argument("--cdrs", nargs="+", required=True, metavar="FILE",
                   help="find_cdrs.py --json output, one per variable domain")
    p.add_argument("--seqs", nargs="+", required=True, metavar="NAME=FILE",
                   help="Sequence file per chain, in the chain order to emit "
                        "(A, B, ... follow this order)")
    p.add_argument("--antigen", metavar="FASTA",
                   help="Antigen FASTA appended after the antibody chains")
    p.add_argument("--out-fasta", required=True, help="Where to write the staged FASTA")
    p.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    args = p.parse_args()

    domains = []
    for f in args.cdrs:
        try:
            domains.append(json.loads(Path(f).read_text()))
        except Exception as e:
            raise SystemExit(f"[!] could not read {f}: {e}")

    modality, detail = resolve_modality(domains)
    model, skill, note = ROUTES[modality]

    records: list[tuple[str, str]] = []
    for item in args.seqs:
        if "=" not in item:
            p.error(f"--seqs expects NAME=FILE, got {item!r}")
        name, path = item.split("=", 1)
        try:
            raw = Path(path).read_text()
        except Exception as e:
            raise SystemExit(f"[!] could not read {path}: {e}")
        if raw.lstrip().startswith(">"):
            raw = "".join(l for l in raw.splitlines() if not l.startswith(">"))
        try:
            records.append((name, clean_sequence(raw, name=name)))
        except SequenceError as e:
            raise SystemExit(f"[!] {e}")

    antigen_records: list[tuple[str, str]] = []
    if args.antigen:
        name, buf = None, []
        for line in Path(args.antigen).read_text().splitlines():
            if line.startswith(">"):
                if name:
                    antigen_records.append((name, "".join(buf)))
                name, buf = line[1:].split()[0], []
            elif line.strip():
                buf.append(line.strip())
        if name:
            antigen_records.append((name, "".join(buf)))

    all_records = records + antigen_records
    chain_ids = [chr(ord("A") + i) for i in range(len(all_records))]

    out = Path(args.out_fasta)
    with out.open("w") as fh:
        for cid, (name, seq) in zip(chain_ids, all_records):
            fh.write(f">{name}|{cid}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")

    payload = {
        "modality": modality,
        "modality_detail": detail,
        "model": model,
        "skill": skill,
        "note": note,
        "job_class": "cofold_heavy",
        "gpu_required": True,
        "fasta": str(out),
        "chains": [{"chain": cid, "name": name, "length": len(seq)}
                   for cid, (name, seq) in zip(chain_ids, all_records)],
        "antibody_chains": chain_ids[:len(records)],
        "antigen_chains": chain_ids[len(records):],
        "next_steps": [
            f"Run `detect_resources.py` once to choose the backend for cofold_heavy.",
            f"Load the `{skill}` skill and predict {out.name} with {model}.",
            "Then: scan_surface_patches.py on the antibody chains "
            f"({', '.join(chain_ids[:len(records)])}) of the selected model.",
        ],
    }

    if args.json:
        emit_json(payload)
        return

    print(f"Modality : {modality}  ({detail})")
    print(f"Model    : {model}  -> load the `{skill}` skill")
    print(f"Reason   : {note}")
    print(f"Job class: cofold_heavy (GPU required)")
    print(f"FASTA    : {out}")
    for c in payload["chains"]:
        print(f"  chain {c['chain']}: {c['name']} ({c['length']} aa)")
    if payload["antigen_chains"]:
        print(f"  antigen chains: {', '.join(payload['antigen_chains'])}")
    else:
        print("  no antigen supplied — this is an apo Fv/VHH model, which is the "
              "correct substrate for surface-patch developability")
    print("\nNext:")
    for s in payload["next_steps"]:
        print(f"  - {s}")


if __name__ == "__main__":
    main()
