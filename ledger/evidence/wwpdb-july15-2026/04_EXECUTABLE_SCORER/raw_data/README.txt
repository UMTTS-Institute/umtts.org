# Layer 4 - Raw Data Directory
## wwPDB July 15, 2026 Terrain Contact

This directory contains the raw input files used by the scorer.
These are the actual files fetched during the 2026-07-19 analysis run.

cds_fasta/     - RefSeq coding sequences (NM_* accessions)
                 Source: NCBI EFetch (https://eutils.ncbi.nlm.nih.gov/)
uniprot_json/  - UniProt entry records in JSON format
                 Source: UniProt REST API (https://rest.uniprot.org/)

All CIF coordinate files (the terrain itself) are in:
  02_TERRAIN_MANIFEST/cif_files/
  Source: RCSB PDB (https://files.rcsb.org/download/)

File integrity: verify SHA256 against 02_TERRAIN_MANIFEST/manifest.csv
Download timestamps: 2026-07-19 ~05:45-06:00 UTC

TRUTH > COMFORT. Always.
UMtts Institute
