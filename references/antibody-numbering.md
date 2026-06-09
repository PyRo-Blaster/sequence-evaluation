# Antibody Numbering Schemes Reference

## CDR Definitions: Kabat vs IMGT vs Chothia

| Scheme | CDR1 (VH) | CDR2 (VH) | CDR3 (VH) | CDR1 (VL) | CDR2 (VL) | CDR3 (VL) |
|--------|-----------|-----------|-----------|-----------|-----------|-----------|
| Kabat  | 31–35     | 50–65     | 95–102    | 24–34     | 50–56     | 89–97     |
| Chothia| 26–32     | 52–56     | 95–102    | 26–32     | 50–52     | 91–96     |
| IMGT   | 27–38     | 56–65     | 105–117   | 27–38     | 56–65     | 105–117   |

*Residue numbers follow Kabat unless stated otherwise.*

The `find_cdrs.py` script uses conserved anchor residues (not fixed positions) and produces Kabat-compatible boundaries:

- **VH CDR1**: residues between C22...C + 4 residues and the Trp that opens FR2
- **VH CDR2**: between FR2-end and the conserved `RFTIS` motif
- **VH CDR3**: between the second conserved Cys (in YYC) and the WGQG motif

## EU Numbering for IgG1 Constant Region

EU numbering aligns the Fc constant region against a reference IgG1 sequence, independent of variable domain insertions/deletions. The `find_mutations.py` script uses this scheme.

| Domain | EU Range  | Residue Count | Notes                              |
|--------|-----------|---------------|------------------------------------|
| CH1    | 118–215   | 98            | Pairs with CL via disulfide C220   |
| Hinge  | 216–230   | 15            | Inter-chain disulfides at C226/C229|
| CH2    | 231–340   | 110           | N-glycosylation at N297            |
| CH3    | 341–447   | 107           | FcRn binding site                  |

## Common Therapeutic Fc Mutations (EU Numbering)

### Effector Silencing
| Mutation | Effect |
|----------|--------|
| L234A    | Abolishes FcγR binding (ADCC/ADCP) |
| L235A    | Part of LALA silencing pair         |
| P329G    | Completes LALA-PG triple mutant; abolishes C1q/CDC |
| L234F + L235E + P331S | FLES variant — alternative silencing |

### Half-Life Extension (FcRn binding)
| Mutation | Effect |
|----------|--------|
| M428L/N434S (LS) | +3–4× serum half-life via pH 6 FcRn affinity |
| N434S    | Partial LS; used alone in some therapeutics |
| T250Q/M428L (QM) | Alternative half-life extension pair |
| YTE (M252Y/S254T/T256E) | Alternative; used in RSV mAbs |

### Allotypic Variants (non-immunogenic, common in therapeutics)
| Mutation | Allotype |
|----------|---------|
| D356E + L358M | EEM |
| G357A         | EMV |

### Knobs-into-Holes (Heterodimerization)
| Chain  | Mutations              | Effect               |
|--------|------------------------|----------------------|
| Knob   | T366W                  | Bulky "knob" in CH3  |
| Hole   | T366S + L368A + Y407V  | Complementary "hole" |

### ADCP/CDC-Retaining Variants
| Name    | Mutations   |
|---------|-------------|
| GASDALIE| G236A/S239D/A330L/I332E |
| SELF    | S267E/L328F |

## VHH (Nanobody) Framework Residues

VHH domains differ from VH in Framework 2 (Kabat positions 44–47):
- VH (conventional): G44, L45, W47 (hydrophobic core)
- VHH (camelid): E44, R45, F47 (hydrophilic, enables VHH solubility without VL)

The `find_cdrs.py` script handles VHH identically to VH — the CDR boundaries are the same.

## Isotype Differences (IgG1 / IgG2 / IgG4)

`find_mutations.py` ships WT references for IgG1, IgG2 and IgG4 and auto-detects
the isotype by alignment. EU numbering is inherited from IgG1 via alignment, so
CH1/CH2/CH3 and the conserved CPxCP hinge core are reliable; the non-conserved
N-terminal hinge residues of IgG2/IgG4 are approximate.

| Isotype | Core hinge | Lower hinge (EU 233–238) | Notes |
|---------|-----------|--------------------------|-------|
| IgG1 | `...CPPCP` (C226, C229) | `ELLGGP` | Most common; strong effector function |
| IgG2 | `ERKCCVECPPCP` (4 hinge Cys) | `PVAGP` | Reduced effector; rigid hinge |
| IgG4 | `ESKYGPPCPSCP` | `EFLGGP` (F234) | Anti-inflammatory; **S228P** prevents Fab-arm exchange |

The famous IgG4 **S228P** stabilizing mutation sits in the CPSCP→CPPCP hinge core
and is numbered correctly (EU 228) because that motif aligns cleanly to IgG1.

> The built-in IgG2/IgG4 references are reconstructed from canonical isotype
> differences. For definitive work supply an authoritative sequence (IMGT
> IGHG2*01 / IGHG4*01, or UniProt P01859 / P01861) via `--reference-fasta`.

## Chain Types and Constant Domains

| Chain | Constant Region | Isotype marker |
|-------|----------------|----------------|
| IgG1 HC | CH1–Hinge–CH2–CH3 | ELLGG in CH2 hinge-proximal |
| κ LC | Cκ | FNRGEC at C-terminus |
| λ LC (λ2) | Cλ2 | SSELTQ in framework |
| VHH only | — | No constant domain in nanobody alone |
