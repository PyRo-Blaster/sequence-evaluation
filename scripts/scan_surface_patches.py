"""Flag charged and hydrophobic surface patches on a predicted structure.

This is a review flag, not a TAP or SAP score. It does not reproduce Raybould
PPC/PNC/PSH or Chennamsetty SAP, and those published cutoffs must not be applied
to these numbers.

Method (``surface-patch-v1``):

1. Residue SASA from Shrake-Rupley (probe 1.4 Å). Relative SASA is residue SASA
   divided by the Tien et al. 2013 theoretical maximum (Gly-X-Gly, ALLOWED).
   A residue is surface-exposed when relative SASA ≥ ``--exposure`` (default 0.20).
2. Positive patch members are exposed Lys/Arg (formal charge +1). Negative patch
   members are exposed Asp/Glu (−1). His is excluded unless ``--include-his``
   (then +0.5). Hydrophobic members are exposed Ala/Val/Leu/Ile/Met/Phe/Trp/Tyr.
3. Same-kind residues are clustered by side-chain centroid distance
   (Cβ and beyond; Cα if the side chain is missing). An edge exists at
   ≤ ``--cluster-cutoff`` Å (default 7.0). A patch is a connected component of
   size ≥ ``--min-size`` (default 2).
4. Score is the sum of per-residue contributions: |formal charge| × rel_sasa
   for charge patches, or a hydrophobicity weight × rel_sasa for hydrophobic
   patches. Weights are Kyte-Doolittle where that scale is positive; Trp/Tyr
   use fixed aromatic weights (1.0 / 0.8) because raw KD marks them hydrophilic.
5. Optional ``--cdr-json`` (same shape as ``analyze_interfaces.py``) tags each
   residue CDR vs framework by exact sequence match. Severity is a heuristic
   for review, not a developability gate:

   - high: ≥3 residues and ≥1 of them in a CDR
   - medium: ≥3 residues, or ≥2 with at least one in a CDR
   - low: any other emitted patch
   Without ``--cdr-json``, nothing is called high.

Run this only on an antibody model (Fv / Fab variable domains), not on the
antigen. A CA-only or backbone-only model is reported as ``coarse_model`` and
must not be used for a patch call.

Usage:
    uv run scripts/scan_surface_patches.py pred.model_idx_0.cif \\
        --chains A B --cdr-json cdrs.json [--json]

# /// script
# requires-python = ">=3.10"
# dependencies = ["biopython>=1.84", "numpy>=1.26"]
# ///
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict

from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB import MMCIFParser
from Bio.PDB.SASA import ShrakeRupley

import sys
from pathlib import Path

# uv does not put the script's directory on sys.path (uv >=0.13), so make the
# sibling _common module importable regardless of how this file is invoked.
sys.path.insert(0, str(Path(__file__).resolve().parent))

from _common import emit_json

_3to1 = {k.upper(): v for k, v in protein_letters_3to1.items()}

# Tien et al., PLOS ONE 2013, theoretical maximum ASA (Å²), ALLOWED region.
MAX_ASA = {
    "A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0,
    "Q": 225.0, "E": 223.0, "G": 104.0, "H": 224.0, "I": 197.0,
    "L": 201.0, "K": 236.0, "M": 224.0, "F": 240.0, "P": 159.0,
    "S": 155.0, "T": 172.0, "W": 285.0, "Y": 263.0, "V": 174.0,
}

# Kyte & Doolittle 1982 where the scale is positive; aromatic overrides for W/Y.
HYDROPHOBIC_WEIGHT = {
    "A": 1.8, "V": 4.2, "L": 3.8, "I": 4.5, "M": 1.9, "F": 2.8,
    "W": 1.0, "Y": 0.8,
}

POSITIVE = {"K": 1.0, "R": 1.0}
NEGATIVE = {"D": -1.0, "E": -1.0}
HIS_CHARGE = 0.5


def _chain_sequence(chain):
    seq, nums = [], []
    for residue in chain:
        if residue.id[0] != " ":
            continue
        seq.append(_3to1.get(residue.resname.upper(), "X"))
        nums.append(residue.id[1])
    return "".join(seq), nums


def _cdr_lookup(model, chain_ids, cdr_map):
    lookup: dict[tuple[str, int], str] = {}
    if not cdr_map:
        return lookup
    for cid in chain_ids:
        if cid not in model or cid not in cdr_map:
            continue
        seq, nums = _chain_sequence(model[cid])
        for cdr in cdr_map[cid]:
            cdr_seq = "".join(cdr.get("seq", "").split()).upper()
            if not cdr_seq:
                continue
            i = seq.find(cdr_seq)
            if i == -1:
                continue
            for resnum in nums[i:i + len(cdr_seq)]:
                lookup[(cid, resnum)] = cdr["name"]
    return lookup


def _sidechain_centroid(residue):
    atoms = [a for a in residue.get_atoms() if a.name not in ("N", "CA", "C", "O", "OXT")]
    if not atoms:
        atoms = [residue["CA"]] if "CA" in residue else list(residue.get_atoms())
    if not atoms:
        return None
    acc = [0.0, 0.0, 0.0]
    for atom in atoms:
        acc[0] += float(atom.coord[0])
        acc[1] += float(atom.coord[1])
        acc[2] += float(atom.coord[2])
    n = float(len(atoms))
    return (acc[0] / n, acc[1] / n, acc[2] / n)


def _dist(a, b) -> float:
    return ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2) ** 0.5


def _components(members, cutoff: float):
    parent = list(range(len(members)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i, j):
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[rj] = ri

    for i in range(len(members)):
        for j in range(i + 1, len(members)):
            if _dist(members[i]["centroid"], members[j]["centroid"]) <= cutoff:
                union(i, j)
    groups: dict[int, list] = defaultdict(list)
    for i, member in enumerate(members):
        groups[find(i)].append(member)
    return list(groups.values())


def _severity(n_residues: int, n_cdr: int, cdr_tagged: bool) -> str:
    if not cdr_tagged:
        return "medium" if n_residues >= 3 else "low"
    if n_residues >= 3 and n_cdr >= 1:
        return "high"
    if n_residues >= 3 or n_cdr >= 1:
        return "medium"
    return "low"


def scan(cif_path, chain_ids, exposure, cluster_cutoff, min_size, cdr_map, include_his):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("pred", cif_path)
    model = structure[0]
    available = [c.id for c in model]
    missing = [c for c in chain_ids if c not in model]
    if missing:
        raise SystemExit(
            f"[!] chain(s) {', '.join(missing)} not in {cif_path}; "
            f"available: {', '.join(available) or '(none)'}"
        )

    # Solvent accessibility must be computed on the antibody chains ALONE. If a
    # bound antigen (or any other chain) is left in the structure it buries the
    # paratope, those CDR residues drop below the exposure cutoff, and the
    # CDR-local patches this scan exists to find disappear. Detach everything
    # not selected before computing SASA.
    chains_in_file = list(available)
    for extra in [c.id for c in model if c.id not in set(chain_ids)]:
        model.detach_child(extra)

    ShrakeRupley(probe_radius=1.4).compute(model, level="A")
    cdr_lookup = _cdr_lookup(model, chain_ids, cdr_map)
    charge = dict(POSITIVE)
    charge.update(NEGATIVE)
    if include_his:
        charge["H"] = HIS_CHARGE

    members = {"positive_charge": [], "negative_charge": [], "hydrophobic": []}
    n_residues = 0
    n_no_sidechain = 0

    for cid in chain_ids:
        for residue in model[cid]:
            if residue.id[0] != " ":
                continue
            aa = _3to1.get(residue.resname.upper())
            if aa not in MAX_ASA:
                continue
            n_residues += 1
            sidechain = [a for a in residue.get_atoms() if a.name not in ("N", "CA", "C", "O", "OXT")]
            if not sidechain:
                n_no_sidechain += 1
            sasa = sum(float(getattr(atom, "sasa", 0.0) or 0.0) for atom in residue.get_atoms())
            rel = sasa / MAX_ASA[aa]
            if rel < exposure:
                continue
            centroid = _sidechain_centroid(residue)
            if centroid is None:
                continue
            region = cdr_lookup.get((cid, residue.id[1]), "FR" if cdr_map else "unassigned")
            base = {
                "chain": cid,
                "resname": residue.resname.upper(),
                "resnum": int(residue.id[1]),
                "region": region,
                "rel_sasa": round(rel, 3),
                "centroid": centroid,
            }
            if aa in charge and charge[aa] > 0:
                members["positive_charge"].append({**base, "contribution": round(charge[aa] * rel, 3)})
            elif aa in charge and charge[aa] < 0:
                members["negative_charge"].append({**base, "contribution": round(abs(charge[aa]) * rel, 3)})
            if aa in HYDROPHOBIC_WEIGHT:
                members["hydrophobic"].append({
                    **base,
                    "contribution": round(HYDROPHOBIC_WEIGHT[aa] * rel, 3),
                })

    patches = []
    cdr_tagged = cdr_map is not None
    for kind, group in members.items():
        for component in _components(group, cluster_cutoff):
            if len(component) < min_size:
                continue
            n_cdr = sum(1 for r in component if r["region"] not in ("FR", "unassigned"))
            residues = []
            for r in sorted(component, key=lambda x: (x["chain"], x["resnum"])):
                residues.append({k: r[k] for k in ("chain", "resname", "resnum", "region", "rel_sasa", "contribution")})
            patches.append({
                "kind": kind,
                "severity": _severity(len(component), n_cdr, cdr_tagged),
                "n_residues": len(component),
                "score": round(sum(r["contribution"] for r in component), 3),
                "cdr_residues": n_cdr,
                "residues": residues,
            })

    order = {"high": 0, "medium": 1, "low": 2}
    kind_order = {"positive_charge": 0, "negative_charge": 1, "hydrophobic": 2}
    patches.sort(key=lambda p: (order[p["severity"]], kind_order[p["kind"]], -p["score"]))

    coarse = n_residues > 0 and (n_no_sidechain / n_residues) > 0.5
    warnings = []
    if coarse:
        warnings.append(
            "More than half of the selected residues have no side-chain atoms. "
            "Patch calls on a CA/backbone-only model are not usable."
        )
    if not cdr_tagged:
        warnings.append(
            "No --cdr-json supplied, so no patch is flagged high. "
            "Pass find_cdrs.py output to distinguish CDR-local patches from framework patches."
        )
    stripped = [c for c in chains_in_file if c not in set(chain_ids)]
    if stripped:
        warnings.append(
            f"Chain(s) {', '.join(stripped)} were present in the file and were "
            "excluded before the SASA calculation, so exposure is that of the "
            "isolated antibody chains (apo), not the bound complex."
        )
    return {
        "cif": cif_path,
        "method": "surface-patch-v1",
        "not_tap": True,
        "not_sap": True,
        "exposure_cutoff": exposure,
        "cluster_cutoff_A": cluster_cutoff,
        "min_size": min_size,
        "include_his": include_his,
        "chains": chain_ids,
        "chains_in_file": chains_in_file,
        "chains_excluded_from_sasa": stripped,
        "coarse_model": coarse,
        "warnings": warnings,
        "patches": patches,
    }


def report(data) -> None:
    print(f"CIF: {data['cif']}")
    print(f"Method: {data['method']}  (not TAP PPC/PNC/PSH, not SAP)")
    print(f"Chains: {', '.join(data['chains'])}   "
          f"exposure ≥ {data['exposure_cutoff']}   "
          f"cluster ≤ {data['cluster_cutoff_A']} Å")
    for warning in data["warnings"]:
        print(f"[!] {warning}")
    if not data["patches"]:
        print("No surface patches at or above the minimum size.")
        return
    print(f"\n{'Sev':<8} {'Kind':<18} {'n':>3} {'CDR':>4} {'Score':>8}  Residues")
    for patch in data["patches"]:
        residues = ", ".join(
            f"{r['chain']}:{r['resname']}{r['resnum']}"
            + (f"[{r['region']}]" if r["region"] not in ("FR", "unassigned") else "")
            for r in patch["residues"]
        )
        print(f"{patch['severity']:<8} {patch['kind']:<18} {patch['n_residues']:>3} "
              f"{patch['cdr_residues']:>4} {patch['score']:>8.3f}  {residues}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Charged and hydrophobic surface patches from a predicted structure")
    parser.add_argument("cif", help="Path to a predicted mmCIF (Fv/Fab variable domains, not the antigen)")
    parser.add_argument("--chains", nargs="+", required=True, metavar="CHAIN",
                        help="Antibody chain IDs to scan. Do not include the antigen.")
    parser.add_argument("--cdr-json", help="JSON of {chain: [{name, seq}]} from find_cdrs.py")
    parser.add_argument("--exposure", type=float, default=0.20,
                        help="Minimum relative SASA to count as surface (default 0.20)")
    parser.add_argument("--cluster-cutoff", type=float, default=7.0,
                        help="Side-chain centroid cutoff in Å (default 7.0)")
    parser.add_argument("--min-size", type=int, default=2,
                        help="Minimum residues in a reported patch (default 2)")
    parser.add_argument("--include-his", action="store_true",
                        help="Count His as +0.5 in positive patches (default: exclude)")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    args = parser.parse_args()

    cdr_map = None
    if args.cdr_json:
        with open(args.cdr_json) as fh:
            cdr_map = json.load(fh)

    data = scan(
        args.cif, args.chains, args.exposure, args.cluster_cutoff,
        args.min_size, cdr_map, args.include_his,
    )
    if args.json:
        emit_json(data)
        return
    report(data)


if __name__ == "__main__":
    main()
