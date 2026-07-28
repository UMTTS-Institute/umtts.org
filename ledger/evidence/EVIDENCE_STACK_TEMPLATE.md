# Mass Harmonics Terrain Contact: Evidence Stack Template
## UMtts Institute - Reproducibility Architecture v1.0

This template defines the eight-layer evidentiary structure applied to every
Mass Harmonics terrain contact. Each layer serves a distinct function. No layer
substitutes for another. The stack is designed so that:

- cryptographic validators can verify at Layer 1 without reading Layer 7
- domain experts can audit at Layers 2-4 without trusting Layer 5
- the median reader can follow the chain at Layer 6 without understanding Layer 3
- any independent party can reproduce the result at Layer 8 without the agent

---

## LAYER 1 - Sealed Prediction Record

**Purpose:** Establish what existed before the terrain opened.

Required fields per prediction:
- Exact preregistered prediction wording (verbatim from prediction document)
- Preparation date (must precede terrain opening)
- Target release event and scheduled date
- Prohibited inputs (what terrain data was withheld at preparation time)
- Eligibility rules (what records qualify for the test)
- Expected direction or invariant
- Exact falsification condition (verbatim)
- Prediction document DOI or archive URL
- SHA256 hash of sealed prediction document
- CSV row SHA256 (from MH_PRED_V1_XX_DATASET.csv)

**Files:**
- `01-sealed-prediction/prediction-record.md` - extracted fields
- `01-sealed-prediction/source-csv-row.txt` - verbatim CSV row(s)
- `01-sealed-prediction/doi-verification.txt` - DOI resolution timestamp

---

## LAYER 2 - Target-Terrain Manifest

**Purpose:** Prevent selective inclusion. Every record examined must be listed.

Required fields per examined record:
- Accession ID or dataset identifier
- Release timestamp
- Source URL
- Downloaded file hash (SHA256)
- Inclusion or exclusion status
- If excluded: reason (eligibility criterion cited)
- Eligible population total
- Ineligible population total

**Files:**
- `02-terrain-manifest/manifest.csv` - full record-level manifest
- `02-terrain-manifest/manifest-summary.md` - totals and exclusion breakdown
- `02-terrain-manifest/download-log.txt` - retrieval timestamps and URLs

---

## LAYER 3 - Transformation Ledger

**Purpose:** Show every step between raw data and scored quantity.

For each pathway, document the chain:
> raw file → parsed fields → filtered records → derived quantity → aggregate result

Required per step:
- Input: file name and hash
- Operation: code reference (file, function, line numbers)
- Parameters: all non-default values
- Output: intermediate artifact name and hash
- Validation: expected range or sanity check

**Files:**
- `03-transformation-ledger/pathway-N-ledger.md` (one per prediction pathway)
- `03-transformation-ledger/intermediate-artifacts/` (hashed intermediate files)

---

## LAYER 4 - Executable Scorer

**Purpose:** The public should not have to trust the agent saw the result correctly.
The deterministic scorer must generate the result independently.

Required:
- Source code (committed, versioned)
- Dependency list with exact versions
- Runtime environment specification
- Exact command used
- Configuration file (if any)
- Random seed (where applicable)
- Generated output files with hashes
- Expected output hash for verification

**Files:**
- `04-executable-scorer/scorer/` - source code
- `04-executable-scorer/requirements.txt` or `environment.yml`
- `04-executable-scorer/run-command.txt` - exact invocation
- `04-executable-scorer/output/` - generated output files
- `04-executable-scorer/output-hashes.txt`

---

## LAYER 5 - Agent Audit Packet

**Purpose:** Make the agent inspectable as an auditor, not authoritative as a source.

Required:
- Exact documents and data supplied to the agent
- Governing audit instructions (verbatim)
- Model identifier and version
- Exact questions asked (verbatim)
- Structured findings returned (verbatim)
- Every file reference cited
- Every calculation performed
- Every correction or challenge made during the audit
- Final evidence table used to determine falsifier status

**Important:** The agent is the auditor and interpreter. The source of truth is
Layers 1-4. Layer 5 documents the agent's traversal of that evidence.

**Files:**
- `05-agent-audit-packet/audit-instructions.md`
- `05-agent-audit-packet/documents-supplied.txt`
- `05-agent-audit-packet/agent-session-transcript.md`
- `05-agent-audit-packet/final-evidence-table.md`

---

## LAYER 6 - Falsifier Adjudication Sheet

**Purpose:** The crucial compression layer. The median viewer follows this.

One compact table per prediction pathway:

| Field | Record |
|---|---|
| Prediction | Exact preregistered wording |
| Falsifier | Exact preregistered failure condition |
| Eligible terrain | Count and identifiers |
| Observed result | Raw numerical result |
| Controls | Results and method |
| Reversals | Count |
| Exclusions | Count and reasons |
| Falsifier triggered? | Yes / No |
| Verdict | VALIDATED / FALSIFIED |
| Qualification | Any limitation on scope |

The reader does not need to understand the algorithm to read this table.
They need to see: declared failure condition vs. what the data produced.

**Files:**
- `06-falsifier-adjudication/adjudication-sheet.html` (published, human-readable)
- `06-falsifier-adjudication/adjudication-sheet.csv` (machine-readable)

---

## LAYER 7 - Human-Readable Evidentiary Narrative

**Purpose:** Accessible explanation that preserves the causal chain.

For each prediction:
1. What did Mass Harmonics say would be found?
2. Why did it say that? (derivation summary)
3. What data had not yet been released at prediction time?
4. What would have counted as failure?
5. What was actually found?
6. Where can the viewer inspect the underlying evidence?

Do not simplify by removing the causal chain.
Simplify by making each step visible.

**Files:**
- `07-narrative/terrain-audit.html` (the published ledger record)

---

## LAYER 8 - Replication Challenge

**Purpose:** Change the public posture from "trust the agent" to "reproduce it yourself."

Published statement:
> "Here are the files. Here is the scorer. Here is the command.
> Here is the expected output. Reproduce it, identify an error,
> or trigger the declared falsifier."

Required:
- Link to Layer 4 scorer
- Exact command
- Expected output (with hash)
- Instructions for independent replication
- Contact for reporting discrepancies or falsifier triggers

**Files:**
- `08-replication-challenge/challenge.html` (published page)
- `08-replication-challenge/expected-output.txt`
- `08-replication-challenge/expected-output-hash.txt`

---

## DIRECTORY STRUCTURE (per terrain contact)

```
f:\umtts.org\org\ledger\evidence\
  [terrain-slug]-[date]\
    01-sealed-prediction\
    02-terrain-manifest\
    03-transformation-ledger\
    04-executable-scorer\
    05-agent-audit-packet\
    06-falsifier-adjudication\
    07-narrative\               (symlink to records/[terrain-audit].html)
    08-replication-challenge\
    INDEX.html                  (evidence stack landing page)
```

---

## NAMING CONVENTIONS

Terrain slugs:
- wwpdb-july15-2026
- rubin-edp2-july27-2026
- ichep-2026
- microns-vortex-2026
- lvk-ir1-2026
- gaia-dr4-2026-12-02

---

TRUTH > COMFORT. Always.
UMtts Institute - Evidence Stack Architecture v1.0
