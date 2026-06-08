"""Analyze antibody-antigen interface contacts from a Chai-1 / AlphaFold mmCIF file.

Reports atom-atom contacts within a distance cutoff, and lists paratope
(antibody) and epitope (antigen) residues.

Usage:
    uv run scripts/analyze_interfaces.py <model.cif> \
        --antibody A B \
        --antigen C \
        [--cutoff 5.0]

# /// script
# requires-python = ">=3.10"
# dependencies = ["biopython>=1.84", "numpy>=1.26"]
# ///
"""

import argparse
from dataclasses import dataclass

import numpy as np
from Bio.PDB import MMCIFParser


@dataclass
class Contact:
    ab_chain: str
    ab_resname: str
    ab_resnum: int
    ag_chain: str
    ag_resname: str
    ag_resnum: int
    distance: float


def analyze_interface(
    cif_path: str,
    ab_chains: list[str],
    ag_chains: list[str],
    cutoff: float = 5.0,
) -> tuple[list[Contact], list[tuple[str, int, str]], list[tuple[str, int, str]]]:
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("pred", cif_path)
    model = structure[0]

    ab_atoms: list[tuple[str, object, object]] = []
    ag_atoms: list[tuple[str, object, object]] = []

    for chain in model:
        bucket = ab_atoms if chain.id in ab_chains else (ag_atoms if chain.id in ag_chains else None)
        if bucket is None:
            continue
        for residue in chain:
            if residue.id[0] != " ":  # skip HETATMs and water
                continue
            for atom in residue:
                bucket.append((chain.id, residue, atom))

    contacts: list[Contact] = []
    seen_ab: set[tuple[str, int, str]] = set()
    seen_ag: set[tuple[str, int, str]] = set()

    for ab_cid, ab_res, ab_atom in ab_atoms:
        for ag_cid, ag_res, ag_atom in ag_atoms:
            diff = ab_atom.coord - ag_atom.coord
            dist = float(np.sqrt(np.dot(diff, diff)))
            if dist <= cutoff:
                contacts.append(Contact(
                    ab_chain=ab_cid,
                    ab_resname=ab_res.resname,
                    ab_resnum=ab_res.id[1],
                    ag_chain=ag_cid,
                    ag_resname=ag_res.resname,
                    ag_resnum=ag_res.id[1],
                    distance=dist,
                ))
                seen_ab.add((ab_cid, ab_res.id[1], ab_res.resname))
                seen_ag.add((ag_cid, ag_res.id[1], ag_res.resname))

    paratope = sorted(seen_ab, key=lambda x: (x[0], x[1]))
    epitope = sorted(seen_ag, key=lambda x: (x[0], x[1]))
    return contacts, paratope, epitope


def fmt_residues(residues: list[tuple[str, int, str]]) -> str:
    return ", ".join(f"{chain}:{resname}{num}" for chain, num, resname in residues)


def main() -> None:
    parser = argparse.ArgumentParser(description="Antibody-antigen interface analysis from mmCIF")
    parser.add_argument("cif", help="Path to mmCIF structure file (Chai-1 / AlphaFold output)")
    parser.add_argument("--antibody", nargs="+", required=True, metavar="CHAIN",
                        help="Chain IDs for the antibody (e.g. A B)")
    parser.add_argument("--antigen", nargs="+", required=True, metavar="CHAIN",
                        help="Chain IDs for the antigen (e.g. C)")
    parser.add_argument("--cutoff", type=float, default=5.0,
                        help="Distance cutoff in Ångströms (default: 5.0)")
    args = parser.parse_args()

    contacts, paratope, epitope = analyze_interface(
        args.cif, args.antibody, args.antigen, args.cutoff
    )

    print(f"CIF file:       {args.cif}")
    print(f"Antibody chains: {', '.join(args.antibody)}")
    print(f"Antigen chains:  {', '.join(args.antigen)}")
    print(f"Cutoff:         {args.cutoff} Å")
    print()
    print(f"Total atom-atom contacts: {len(contacts)}")
    print(f"Paratope residues (n={len(paratope)}): {fmt_residues(paratope) or 'None'}")
    print(f"Epitope  residues (n={len(epitope)}):  {fmt_residues(epitope) or 'None'}")


if __name__ == "__main__":
    main()
