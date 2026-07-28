# Layer 3 - Transformation Ledger: Pathway PDB-1
## Equilateral-Triad Enrichment in Resolved Protein Cores

Generated: 2026-07-20T19:36:35Z
Source report: MH_wwPDB_July15_2026_Analysis_Report.txt
Source SHA256: 976556fefd8c85b4e9f3005a329c8a3aa8731325e3ee00ef549e7719f6a0ff1c

---

## Chain: Raw CIF File -> Q_delta Score

### Step 1: Retrieve coordinate file
- **Input:** PDB ID (from RCSB July 15 release query)
- **Operation:** HTTP GET https://files.rcsb.org/download/{pdb_id}.cif
- **Output:** Local .cif file (mmCIF/PDBx format)
- **Code reference:** wwPDB_July15_2026_Rigorous_Analysis.py, fetch_cif_file()
- **Validation:** File size > 0; gemmi parser loads without error

### Step 2: Extract Cα positions
- **Input:** Local .cif file
- **Operation:** gemmi library mmCIF parse; extract CA atom positions for all polymer residues
- **Output:** List of (residue_index, x, y, z, B_factor, chain_id) tuples
- **Code reference:** wwPDB_July15_2026_Rigorous_Analysis.py, extract_ca_positions()
- **Validation:** n_CA > 0; all positions are finite floats
- **Parameter:** None. No distance cutoff. No threshold.

### Step 3: Classify structural regions
- **Input:** List of (residue_index, x, y, z, B_factor) tuples
- **Operation:** Compute mean B-factor. Apply quantile threshold sweep.
  For each quantile Q in {Q75, Q80, Q85, Q90}:
    - mobile: B_factor > Qth percentile within entry
    - resolved_core: buried + B_factor <= Qth percentile
    - resolved_surface: solvent-exposed + B_factor <= Qth percentile
- **Output:** Region label per residue {resolved_core, resolved_surface, mobile, unmodeled}
- **Code reference:** wwPDB_July15_2026_Rigorous_Analysis.py, classify_regions()
- **Note:** All four quantile thresholds reported. No post-hoc selection.

### Step 4: Build 3D Delaunay tessellation
- **Input:** Cα (x, y, z) positions
- **Operation:** scipy.spatial.Delaunay on the Cα point cloud
- **Output:** Tetrahedral mesh; extracted triangular faces
- **Code reference:** wwPDB_July15_2026_Rigorous_Analysis.py, compute_q_delta()
- **Parameter:** None. Delaunay is parameter-free.

### Step 5: Compute Q_delta per triangular face
- **Input:** Triangle vertices (a_pos, b_pos, c_pos)
- **Operation:**
    a = dist(b, c); b_side = dist(a, c); c_side = dist(a, b)
    A = 0.5 * norm(cross(b-a_pos, c-a_pos))
    Q_delta = 4 * sqrt(3) * A / (a^2 + b_side^2 + c_side^2)
- **Output:** Q_delta in [0, 1] per triangle
- **Validation:** 0 <= Q_delta <= 1; Q_delta = 1 iff a = b_side = c_side

### Step 6: Assign triangles to regions
- **Input:** Q_delta per triangle; region labels per residue vertex
- **Operation:** Triangle assigned to region by majority vertex label
- **Output:** Q_delta distribution per region

### Step 7: Compute medians and compare
- **Input:** Q_delta distributions per region
- **Operation:** median(Q_delta | resolved_core) vs median(Q_delta | resolved_surface)
  and vs median(Q_delta | mobile)
- **Output:** Direction boolean (core > surface, core > mobile)
- **Prediction requirement:** Both comparisons must be True

---

## Per-Entry Results (entries with PDB-1 analysis)

| Entry | CA atoms | Core median Q | Surface median Q | Core > Surface | Core > Mobile |
|---|---|---|---|---|---|
| 10DT | 175 | 0.817930 | 0.775351 | YES | - |
| 10ZK | 3081 | 0.850966 | 0.828634 | YES | - |
| 9ZYS | 2644 | 0.813821 | 0.723562 | YES | - |

---

## Aggregate Result

- Total entries analyzed for PDB-1: 3 primary entries (9ZYS, 10ZK, 10DT)
- Direction reversals across all 12 threshold tests: 0
- Verdict: VALIDATED (falsifier "No positive core shift or stable reversed shift" - NOT TRIGGERED)

TRUTH > COMFORT. Always.
UMtts Institute - Layer 3 Pathway 1 Ledger
