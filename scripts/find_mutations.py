"""Map mutations in an IgG1 heavy chain constant region vs WT IGHG1*01 (EU numbering).

Usage:
    uv run scripts/find_mutations.py <full_heavy_chain_sequence>

The script auto-detects the start of the constant region by searching for the
conserved ASTKGPSVF motif that opens CH1.

# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""

import argparse
import sys

# WT human IgG1 constant region (IGHG1*01): CH1 + Hinge + CH2 + CH3
WT_IgG1_CONST = (
    "ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSS"
    "VVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTL"
    "MISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEY"
    "KCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQP"
    "ENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVLHEALHNHYTQKSLSLSPGK"
)

# EU numbering: CH1 118-215, Hinge 216-230, CH2 231-340, CH3 341-447
EU_NUMBERS: list[int] = (
    list(range(118, 216))   # CH1: 98 residues
    + list(range(216, 231)) # Hinge: 15 residues
    + list(range(231, 341)) # CH2: 110 residues
    + list(range(341, 448)) # CH3: 107 residues
)

DOMAIN_MAP = {
    **{n: "CH1"   for n in range(118, 216)},
    **{n: "Hinge" for n in range(216, 231)},
    **{n: "CH2"   for n in range(231, 341)},
    **{n: "CH3"   for n in range(341, 448)},
}


def find_mutations(chain_seq: str, label: str) -> None:
    anchor = "ASTKGPSVF"
    pos = chain_seq.find(anchor)
    if pos == -1:
        print(f"[!] Could not find constant region anchor (ASTKGPSVF) in {label}.", file=sys.stderr)
        sys.exit(1)

    const_seq = chain_seq[pos:]
    max_len = min(len(const_seq), len(WT_IgG1_CONST))

    mutations: list[tuple[int, str, str, str]] = []
    for idx in range(max_len):
        wt_aa = WT_IgG1_CONST[idx]
        mut_aa = const_seq[idx]
        if wt_aa != mut_aa:
            eu = EU_NUMBERS[idx]
            domain = DOMAIN_MAP.get(eu, "?")
            mutations.append((eu, wt_aa, mut_aa, domain))

    print(f"\n{label} — {len(mutations)} mutation(s) vs WT IGHG1*01:")
    if not mutations:
        print("  (none — wild-type constant region)")
        return

    print(f"  {'EU #':>6}  {'Domain':<7}  Mutation")
    print(f"  {'------':>6}  {'-------':<7}  --------")
    for eu, wt, mut, domain in mutations:
        print(f"  {eu:>6}  {domain:<7}  {wt}{eu}{mut}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Map EU-numbered mutations in an IgG1 heavy chain constant region"
    )
    parser.add_argument("sequence", help="Full heavy chain amino acid sequence (with or without signal peptide)")
    parser.add_argument("--label", default="Heavy Chain", help="Label for display (default: 'Heavy Chain')")
    args = parser.parse_args()

    find_mutations(args.sequence.strip().upper(), args.label)


if __name__ == "__main__":
    main()
