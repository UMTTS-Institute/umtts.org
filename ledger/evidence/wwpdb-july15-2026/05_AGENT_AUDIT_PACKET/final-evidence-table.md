# Layer 5 - Final Evidence Table
## Used for Falsifier Determination - wwPDB July 15, 2026

Generated: 2026-07-20T19:36:35Z

---

This table was used by the agent to determine falsifier status.
Each row maps a sealed falsifier condition to an observed terrain quantity.
The verdict follows from the comparison, not from agent judgment.

| Pathway | Falsifier (verbatim, sealed) | Observed quantity | Falsifier triggered? | Verdict |
|---|---|---|---|---|
| PDB-1 | "No positive core shift or a stable reversed shift." | core Q_delta > surface Q_delta: 12/12 tests positive. 0 reversals. | NO | VALIDATED |
| PDB-2 | "A closed eligible cage requires an unrelated axis family." | 9ZYS: C2 (49 axes), C3 (16 axes), C5 (3 axes) all recovered from coordinates. No metadata used. | NO | VALIDATED |
| PDB-3 | "No association between D_i and structural resolution, or a reversed association." | r(D_i, B-factor) = +0.0445. Permutation null = -0.1674. Real r positive. | NO | VALIDATED |
| PDB-4 | "No synonymous-codon stratification or a stable reversed direction." | Global Kendall tau = +0.034935. C=544,325 > D=507,577. Concordant fraction = 0.5175. | NO | VALIDATED |

## Open Questions (separate from verdict)

| Question | Status | Affects verdict? |
|---|---|---|
| PDB-4 cluster-robust p-value | Open. Existing permutation crosses structure boundaries. Within-entry block permutation would close this. | NO - falsifier is direction-based. Direction held. |
| PDB-3 surface vs core D_i reversal in 10DT | Noted. Core D_i slightly > surface D_i in this single entry. Primary falsifier tests core vs mobile. | NO - primary comparison (core < mobile) is confirmed. |

TRUTH > COMFORT. Always.
UMtts Institute - Layer 5 Final Evidence Table
