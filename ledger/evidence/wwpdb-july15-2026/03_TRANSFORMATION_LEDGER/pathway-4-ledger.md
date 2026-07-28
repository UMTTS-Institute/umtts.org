# Layer 3 - Transformation Ledger: Pathway PDB-4
## Synonymous-Codon Orientation Stratification

Generated: 2026-07-20T19:36:35Z
Source report: PDB4_FullDataset_Report.txt
Source SHA256: 8ff20924e48dfaf85fb10a188989aeee80a571c522adcfdffe95407c57e0cf1e

---

## Chain: CIF + CDS -> Synonymous Pairs -> Kendall Tau

### Steps 1-5: Same as PDB-3
(CIF retrieval, CDS retrieval, Delta assignment, Delaunay graph, D_i computation)

### Step 6: Enumerate synonymous pairs
- **Input:** D_i per residue; amino-acid identity per residue
- **Operation:**
    For each amino acid class AA:
      For each pair of residues (i, j) where AA(i) = AA(j) but codon(i) != codon(j):
        Concordant (C): (|Delta_i| < |Delta_j|) AND (B_i < B_j)
                     OR (|Delta_i| > |Delta_j|) AND (B_i > B_j)
        Discordant (D): opposite
        Tied: |Delta_i| = |Delta_j| or B_i = B_j
- **Constraint:** Pairs are within same amino-acid class. Different codons required.
  Amino-acid identity held fixed. Chemistry and face-anchor held fixed.
  Only the codon's Orientation Delta varies.
- **Code reference:** PDB4_FullDataset_Analysis.py, lines 400-520
- **Random seed:** 42 (for permutation control only; main analysis is deterministic)

### Step 7: Compute global pooled Kendall tau
- **Input:** All (C, D, tied) tallies across all entries and AA classes
- **Operation:**
    tau = (C - D) / (C + D + tied)
    concordant_fraction = C / (C + D)
- **Output:** Global tau; per-AA tau
- **Null prediction:** tau = 0 (concordant fraction = 0.50)
- **Mass Harmonics prediction:** tau > 0 (concordant fraction > 0.50)

### Step 8: Permutation control
- **Input:** All B-factor values; AA class labels
- **Operation:** Shuffle B-factors within each AA class (seed 42). Recompute Kendall tau.
- **Output:** Control tau
- **Limitation note (documented):** This permutation is amino-acid-preserving
  but NOT structure-preserving. It pools across entry boundaries. It does not
  establish a cluster-robust p-value. A within-entry block permutation would
  close the confidence interval question. This is an open statistical question
  separate from the falsifier verdict.

---

## Full Dataset Results: 293 entries, 129 contributing

**Data collection summary (verbatim from report):**
  Entries in release:            293
  Entries with UniProt linkage:  253
  Entries with CDS:              258
  Entries contributing triples:  129
  Total (aa, |delta|, B) triples: 8507
  Skipped - no UniProt:           40
  Skipped - no CDS:              124
  Skipped - no coord triples:      0

**Per-AA global pooled Kendall tau (verbatim from report):**

| AA | n_residues | n_codons | pairs | C | D | tied | tau |
|---|---|---|---|---|---|---|---|
| A | 197 | 3 | 12191 | 6201 | 5990 | 0 | 0.0173 |
| C | 32 | 3 | 277 | 170 | 107 | 0 | 0.2274 |
| D | 207 | 3 | 13142 | 6589 | 6553 | 0 | 0.0027 |
| F | 110 | 3 | 4021 | 2193 | 1828 | 0 | 0.0908 |
| H | 233 | 6 | 21436 | 10595 | 10841 | 0 | -0.0115 |
| I | 336 | 3 | 33342 | 15252 | 18090 | 2 | -0.0851 |
| M | 105 | 6 | 4203 | 2713 | 1490 | 0 | 0.2910 |
| N | 979 | 2 | 9690 | 2842 | 6848 | 0 | -0.4134 |
| P | 198 | 3 | 9267 | 4430 | 4837 | 0 | -0.0439 |
| Q | 125 | 3 | 4809 | 2837 | 1972 | 13 | 0.1799 |
| R | 242 | 4 | 19367 | 10200 | 9167 | 0 | 0.0533 |
| S | 1262 | 3 | 179630 | 74283 | 105347 | 11 | -0.1729 |
| T | 154 | 3 | 7523 | 3655 | 3868 | 0 | -0.0283 |
| V | 1839 | 6 | 709768 | 390666 | 319102 | 24 | 0.1008 |
| W | 73 | 3 | 1532 | 313 | 1219 | 0 | -0.5914 |
| Y | 266 | 3 | 21704 | 11386 | 10318 | 1 | 0.0492 |
| **GLOBAL** | | | **1,051,902** | **544,325** | **507,577** | **51** | **0.0349** |

**Global Kendall tau = +0.034935**
**Concordant fraction = 0.5175 (null = 0.5000)**
**Permutation control tau = +0.025924**
**Real tau exceeds control in predicted direction: True**

**Verdict:** VALIDATED. Falsifier ("no stratification or stable reversed direction") NOT TRIGGERED.

---

TRUTH > COMFORT. Always.
UMtts Institute - Layer 3 Pathway 4 Ledger
