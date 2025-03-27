# Understanding Immune Receptor Gene Notation in TRUST4 Output

## Overview

In TRUST4 output files, immune receptor genes are represented using a standardized nomenclature. Let's break down this example:

```
V_gene         D_gene    J_gene      C_gene
IGKV1-8*01     *         IGKJ4*01    IGKC
```

## Gene Segment Nomenclature

### V_gene: IGKV1-8*01

- **IG**: Immunoglobulin (antibody)
- **K**: Kappa light chain
- **V**: Variable region gene
- **1-8**: Family 1, member 8
- ***01**: Allele variant 01

### D_gene: *

- The asterisk indicates no D gene was identified
- This is normal for kappa light chains, which don't use D genes
- D (Diversity) genes only occur in heavy chains (IGH) and some T-cell receptors (TRB, TRD)

### J_gene: IGKJ4*01

- **IG**: Immunoglobulin
- **K**: Kappa light chain
- **J**: Joining region gene
- **4**: The 4th J gene in the kappa locus
- ***01**: Allele variant 01

### C_gene: IGKC

- **IG**: Immunoglobulin
- **K**: Kappa light chain
- **C**: Constant region gene
- Kappa chains have only one constant region gene (no number needed)

## Biological Context

This notation represents a complete kappa light chain from an antibody, showing which specific gene segments were recombined during B-cell development to create this particular immune receptor. The V and J regions (particularly their CDR regions) contribute to antigen-binding specificity, while the C region provides the structural framework.

## Other Common Notations

- **Heavy chains**: IGHV, IGHD, IGHJ, IGHC (uses all gene segments)
- **Lambda light chains**: IGLV, *, IGLJ, IGLC
- **T-cell receptor alpha**: TRAV, *, TRAJ, TRAC
- **T-cell receptor beta**: TRBV, TRBD, TRBJ, TRBC

The "*" always indicates absence of a gene segment, which is structurally normal for certain chain types.