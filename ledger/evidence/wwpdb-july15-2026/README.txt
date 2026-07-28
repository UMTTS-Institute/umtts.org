================================================================================
  UMtts Institute - Terrain Contact Package
  wwPDB Weekly Structural Release, July 15, 2026
  Mass Harmonics Volume I - Target 1 of 7
================================================================================

AUTHOR:       Thomas Russell Giboney
AFFILIATION:  UMtts Institute
FRAMEWORK:    Mass Harmonics psi_m
SEALED:       2026-07-10 (before advance sequence opening 2026-07-11 03:00 UTC)
TERRAIN OPEN: 2026-07-15 00:00 UTC (coordinates)
DOI:          https://doi.org/10.5281/zenodo.21304301
VERDICT:      FULL SWEEP - 4 of 4 pathways validated

--------------------------------------------------------------------------------
WHAT THIS PACKAGE IS
--------------------------------------------------------------------------------

This package contains the complete evidentiary record for the first terrain
contact of the Mass Harmonics Volume I prediction set against the wwPDB weekly
structural release of July 15, 2026.

Four predictions were sealed before the terrain existed. All four stated
falsifiers were not triggered by the terrain.

The package is organized into eight layers following the UMtts Institute
Terrain Contact Package Specification v1.0.

--------------------------------------------------------------------------------
WHERE TO START
--------------------------------------------------------------------------------

If you want to READ THE COMPLETE RECORD IN ONE FILE:

  07_FULL_RECORD/wwpdb-july15-2026-FULL-RECORD.html

  This is the primary artifact. It is self-contained. Open it in any browser.
  It contains: the governing source chain (equations verbatim, Unicode),
  full prediction text (verbatim from sealed document), full derivational
  steps for all four pathways, observed terrain results, controls, and
  complete falsifier adjudication.

  No external links are required to read it.

If you want to CHECK THE VERDICTS QUICKLY:

  06_FALSIFIER_ADJUDICATION/adjudication-sheet.html

  One table per pathway. Declared failure condition vs. observed result.
  Binary verdict per pathway. No domain expertise required.

If you want to VERIFY THE PREDICTION WAS SEALED BEFORE TERRAIN:

  01_SEALED_PREDICTION/prediction-record.md     (Layer 1 fields)
  01_SEALED_PREDICTION/source-csv-rows.txt      (verbatim CSV rows with SHA256)
  -> Verify against: https://doi.org/10.5281/zenodo.21304301

If you want to REPRODUCE THE RESULT YOURSELF:

  04_EXECUTABLE_SCORER/                         (Python scripts)
  08_REPLICATION_CHALLENGE/challenge.html       (instructions and expected output)

If you want to VERIFY FILE INTEGRITY:

  MANIFEST.sha256                               (SHA256 of every file)

  Windows: certutil -hashfile <filename> SHA256
  Linux/Mac: sha256sum <filename>

--------------------------------------------------------------------------------
PACKAGE STRUCTURE
--------------------------------------------------------------------------------

UMTTS_TC_wwpdb-july15-2026_20260715/
|
+-- README.txt                          (this file)
+-- MANIFEST.sha256                     (SHA256 of every file in package)
|
+-- 01_SEALED_PREDICTION/               Layer 1: What existed before terrain opened
|   +-- prediction-record.md            Extracted Layer 1 fields
|   +-- source-csv-rows.txt             Verbatim CSV rows with SHA256 hashes
|   +-- MH_PREDICTION_02_*.md           Full prediction document (copy)
|   +-- MH_PRED_V1_02_wwPDB.csv         Full CSV with per-row hashes
|
+-- 02_TERRAIN_MANIFEST/                Layer 2: Every record examined
|   +-- manifest.csv                    All 293 entries, inclusion/exclusion
|   +-- manifest-summary.md             Totals and exclusion breakdown
|
+-- 03_TRANSFORMATION_LEDGER/           Layer 3: Raw data to scored quantity
|   +-- pathway-1-ledger.md
|   +-- pathway-2-ledger.md
|   +-- pathway-3-ledger.md
|   +-- pathway-4-ledger.md
|   +-- intermediate/                   Hashed intermediate artifacts
|
+-- 04_EXECUTABLE_SCORER/               Layer 4: Reproduce the result yourself
|   +-- PDB4_FullDataset_Analysis.py    Primary PDB-4 scorer
|   +-- [other analysis scripts]
|   +-- requirements.txt
|   +-- run-command.txt
|   +-- random-seeds.txt
|   +-- output/                         Generated output files
|
+-- 05_AGENT_AUDIT_PACKET/              Layer 5: Agent as auditor, not authority
|   +-- audit-instructions.md
|   +-- documents-supplied.txt
|   +-- session-transcript.md
|   +-- final-evidence-table.md
|
+-- 06_FALSIFIER_ADJUDICATION/          Layer 6: Median-reader legible verdicts
|   +-- adjudication-sheet.html         One table per pathway
|   +-- adjudication-sheet.csv          Machine-readable equivalent
|
+-- 07_FULL_RECORD/                     Layer 7: Complete self-contained record
|   +-- wwpdb-july15-2026-FULL-RECORD.html    *** START HERE ***
|
+-- 08_REPLICATION_CHALLENGE/           Layer 8: Standing open challenge
    +-- challenge.html
    +-- expected-output.txt

--------------------------------------------------------------------------------
THE FOUR VERDICTS
--------------------------------------------------------------------------------

PDB-1  Equilateral-triad enrichment in resolved cores
       Falsifier: "No positive core shift or a stable reversed shift."
       Result:    0 of 12 direction reversals. Core shift +0.090 (9ZYS).
       VERDICT:   VALIDATED

PDB-2  {2,3,5} axis scaffold recovery in near-spherical cages
       Falsifier: "A closed eligible cage requires an unrelated axis family."
       Result:    C2, C3, C5 all recovered from 9ZYS coordinates. No metadata used.
       VERDICT:   VALIDATED

PDB-3  Codon Z-cascade residual vs. structural resolution
       Falsifier: "No association between Di and structural resolution,
                  or a reversed association."
       Result:    r(Di, B-factor) = +0.0445. Permutation null: -0.1674.
       VERDICT:   VALIDATED

PDB-4  Synonymous-codon orientation stratification
       Falsifier: "No synonymous-codon stratification or a stable reversed direction."
       Result:    C = 544,325. D = 507,577. tau = +0.034935.
                  Aggregate direction across 1,051,902 pairs: POSITIVE.
       VERDICT:   VALIDATED

--------------------------------------------------------------------------------
STANDING REPLICATION CHALLENGE
--------------------------------------------------------------------------------

Here are the files. Here is the scorer. Here is the command.
Here is the expected output.

Reproduce it, identify a computational error, or trigger one of the four
declared falsifiers. If you trigger a falsifier, the verdict changes.
That is the protocol.

See: 08_REPLICATION_CHALLENGE/challenge.html

--------------------------------------------------------------------------------
CONTACT
--------------------------------------------------------------------------------

UMtts Institute: https://umtts.org
Contact: contact@umtts.org
Ledger: https://umtts.org/ledger.html

TRUTH > COMFORT. Always.
Unitas Monstrat tenebras transire semper

================================================================================
