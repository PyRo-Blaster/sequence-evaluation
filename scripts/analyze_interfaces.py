"""Analyze antibody-antigen interface contacts from a Chai-1 / AlphaFold mmCIF file.

Uses a KD-tree (Bio.PD.NeighborSearch) for fast contact detection, reports
per-residue-pair contacts with minimum distance, lists paratope/epitope
residues, and can map paratope residues back to CDRs (output of find_cdrs.py)
to flag framework-mediated binding.

Usage:
    uv run scripts/analyze_interfaces.py <model.cif> \
        --antibody A B --antigen C [--cutoff 5.0] [--cdr-json cdrs.json] [--json]

The optional --cdr-json file maps antibody chain IDs to CDR lists, e.g.:
    {"A": [{"name": "CDR-H1", "seq": "GFTFSSYA"}, ...]}

# /// script
# requires-python = ">=3.10"
# dependencies = ["biopython>=1.84", "numpy>=1.26"]
# ///
"""

import argparse
import json
from collections import defaultdict
from dataclasses import dataclass

from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB import MMCIFParser, NeighborSearch

from _common import emit_json

_3to1 = {k.upper(): v for k, v in protein_letters_3to1.items()}


@dataclass
class ResiduePair:
    ab_chain: str
    ab_resname: str
    ab_resnum: int
    ag_chain: str
    ag_resname: str
    ag_resnum: int
    n_contacts: int
    min_distance: float


def _chain_sequence(chain):
    """Return (one-letter sequence, [residue numbers]) for standard residues."""
    seq, nums = [], []
    for residue in chain:
        if residue.id[0] != " ":
            continue
        seq.append(_3to1.get(residue.resname.upper(), "X"))
        nums.append(residue.id[1])
    return "".join(seq), nums


def _cdr_ranges(structure_model, ab_chains, cdr_map):
    """Map (chain_id, resnum) -> CDR name by locating CDR sequences in each chain."""
    lookup: dict[tuple[str, int], str] = {}
    for cid in ab_chains:
        if cid not in structure_model or cid not in cdr_map:
            continue
        seq, nums = _chain_sequence(structure_model[cid])
        for cdr in cdr_map[cid]:
            cdr_seq = cdr.get("seq", "")
            if not cdr_seq:
                continue
            i = seq.find(cdr_seq)
            if i == -1:
                continue
            for resnum in nums[i:i + len(cdr_seq)]:
                lookup[(cid, resnum)] = cdr["name"]
    return lookup


def analyze_interface(cif_path, ab_chains, ag_chains, cutoff=5.0, cdr_map=None):
    parser = MMCIFParser(QUIET=True)
    model = parser.get_structure("pred", cif_path)[0]

    ab_set, ag_set = set(ab_chains), set(ag_chains)
    ab_atoms, ag_atoms = [], []
    for chain in model:
        if chain.id in ab_set:
            target = ab_atoms
        elif chain.id in ag_set:
            target = ag_atoms
        else:
            continue
        for residue in chain:
            if residue.id[0] != " ":
                continue
            for atom in residue:
                target.append(atom)

    # KD-tree over antigen atoms; query each antibody atom's neighborhood.
    pairs: dict[tuple, list[float]] = defaultdict(list)
    total_contacts = 0
    if ab_atoms and ag_atoms:
        ns = NeighborSearch(ag_atoms)
        for ab_atom in ab_atoms:
            for ag_atom in ns.search(ab_atom.coord, cutoff):
                dist = float(ab_atom - ag_atom)  # Bio.PDB Atom.__sub__ = distance
                ab_res = ab_atom.get_parent()
                ag_res = ag_atom.get_parent()
                key = (
                    ab_res.get_parent().id, ab_res.resname, ab_res.id[1],
                    ag_res.get_parent().id, ag_res.resname, ag_res.id[1],
                )
                pairs[key].append(dist)
                total_contacts += 1

    cdr_lookup = _cdr_ranges(model, ab_chains, cdr_map) if cdr_map else {}

    residue_pairs = []
    seen_ab, seen_ag = {}, {}
    for key, dists in pairs.items():
        abc, abrn, abnum, agc, agrn, agnum = key
        residue_pairs.append(ResiduePair(abc, abrn, abnum, agc, agrn, agnum,
                                         len(dists), min(dists)))
        seen_ab[(abc, abnum, abrn)] = cdr_lookup.get((abc, abnum), "FR")
        seen_ag[(agc, agnum, agrn)] = None
    residue_pairs.sort(key=lambda p: p.min_distance)

    paratope = [(c, n, r, seen_ab[(c, n, r)]) for (c, n, r) in
                sorted(seen_ab, key=lambda x: (x[0], x[1]))]
    epitope = sorted(seen_ag, key=lambda x: (x[0], x[1]))
    return total_contacts, residue_pairs, paratope, epitope


def main() -> None:
    parser = argparse.ArgumentParser(description="Antibody-antigen interface analysis from mmCIF")
    parser.add_argument("cif", help="Path to mmCIF structure file")
    parser.add_argument("--antibody", nargs="+", required=True, metavar="CHAIN")
    parser.add_argument("--antigen", nargs="+", required=True, metavar="CHAIN")
    parser.add_argument("--cutoff", type=float, default=5.0, help="Distance cutoff in A (default 5.0)")
    parser.add_argument("--cdr-json", help="JSON of {chain: [{name, seq}]} to map paratope to CDRs")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    args = parser.parse_args()

    cdr_map = None
    if args.cdr_json:
        with open(args.cdr_json) as fh:
            cdr_map = json.load(fh)

    total, residue_pairs, paratope, epitope = analyze_interface(
        args.cif, args.antibody, args.antigen, args.cutoff, cdr_map
    )

    if args.json:
        emit_json({
            "cif": args.cif, "cutoff": args.cutoff,
            "antibody_chains": args.antibody, "antigen_chains": args.antigen,
            "total_atom_contacts": total,
            "residue_pairs": [p.__dict__ for p in residue_pairs],
            "paratope": [{"chain": c, "resnum": n, "resname": r, "region": reg}
                         for (c, n, r, reg) in paratope],
            "epitope": [{"chain": c, "resnum": n, "resname": r} for (c, n, r) in epitope],
        })
        return

    print(f"CIF file:        {args.cif}")
    print(f"Antibody chains: {', '.join(args.antibody)}   Antigen chains: {', '.join(args.antigen)}")
    print(f"Cutoff:          {args.cutoff} A")
    print(f"\nTotal atom-atom contacts: {total}")
    print(f"Residue-residue contact pairs: {len(residue_pairs)}")
    if residue_pairs:
        print(f"\n  {'min d':>6}  {'#':>3}  Antibody        Antigen")
        print(f"  {'-----':>6}  {'---':>3}  --------------  --------------")
        for p in residue_pairs:
            print(f"  {p.min_distance:>6.2f}  {p.n_contacts:>3}  "
                  f"{p.ab_chain}:{p.ab_resname}{p.ab_resnum:<8}  "
                  f"{p.ag_chain}:{p.ag_resname}{p.ag_resnum}")

    region_note = any(reg not in (None, "FR") for *_, reg in paratope)
    para_str = ", ".join(
        f"{c}:{r}{n}" + (f"[{reg}]" if reg and reg != "FR" else "")
        for (c, n, r, reg) in paratope
    )
    print(f"\nParatope residues (n={len(paratope)}): {para_str or 'None'}")
    epi_str = ", ".join(f"{c}:{r}{n}" for (c, n, r) in epitope)
    print(f"Epitope residues  (n={len(epitope)}):  {epi_str or 'None'}")
    if cdr_map and not region_note:
        print("[note] no paratope residue mapped to a CDR — check chain IDs / CDR seqs, "
              "or this may be a framework-mediated binding mode.")
    elif cdr_map:
        fr_contacts = [f"{c}:{r}{n}" for (c, n, r, reg) in paratope if reg == "FR"]
        if fr_contacts:
            print(f"[note] framework (non-CDR) paratope contacts: {', '.join(fr_contacts)}")


if __name__ == "__main__":
    main()
