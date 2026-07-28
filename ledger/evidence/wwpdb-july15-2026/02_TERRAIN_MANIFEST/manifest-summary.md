# Layer 2 - Terrain Manifest Summary
## wwPDB July 15, 2026

Generated: 2026-07-20T19:36:35Z
Source report: PDB4_FullDataset_Report.txt
SHA256: 8ff20924e48dfaf85fb10a188989aeee80a571c522adcfdffe95407c57e0cf1e

---

## Totals

| Metric | Count |
|---|---|
| Total release entries (RCSB confirmed) | 293 |
| CIF files downloaded to local analysis dir | 131 |
| Entries with UniProt linkage | 253 |
| Entries with RefSeq NM_ CDS available | 129 |
| Entries contributing synonymous pairs (PDB-4) | 129 |
| Total (AA, |delta|, B-factor) triples | 8,507 |

## Exclusion Breakdown (PDB-4)

| Reason | Count |
|---|---|
| No UniProt linkage | 40 |
| No RefSeq NM_ CDS in UniProt | 124 |
| Included (contributed codon triples) | 129 |

## Notes

- Exclusions apply to specific test pathways only (PDB-3, PDB-4).
- All 293 entries remain in the release ledger.
- No entry was excluded because it appeared unfavorable to the prediction.
- PDB-1 and PDB-2 use coordinate geometry only; they operate on CIF files directly.
  CIF-based eligibility is separate from CDS eligibility.
- The CIF download timestamp (2026-07-19 ~05:45-06:00 UTC) reflects the analysis
  run date, not the terrain opening date (2026-07-15).

## Download Log

CIF files retrieved via RCSB PDBx/mmCIF API:
  https://files.rcsb.org/download/{pdb_id}.cif

CDS sequences retrieved via UniProt REST API and NCBI EFetch:
  https://rest.uniprot.org/uniprotkb/{uniprot_id}.json
  https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={refseq_id}&rettype=fasta_cds_na

TRUTH > COMFORT. Always.
UMtts Institute - Layer 2 Terrain Manifest
