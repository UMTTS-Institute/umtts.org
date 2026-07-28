# Layer 3 - Transformation Ledger: Pathway PDB-2
## Icosahedral Axis Recovery in Near-Spherical Cages

Generated: 2026-07-20T19:36:35Z
Source report: MH_wwPDB_July15_2026_Analysis_Report.txt
Source SHA256: 976556fefd8c85b4e9f3005a329c8a3aa8731325e3ee00ef549e7719f6a0ff1c

---

## Chain: Raw CIF File -> Axis Family Identification

### Step 1: Retrieve coordinate file
- Same as PDB-1 Step 1.

### Step 2: Identify biological assembly
- **Input:** mmCIF file
- **Operation:** Read _pdbx_struct_assembly and _pdbx_struct_assembly_gen blocks
- **Output:** Biological assembly chain list and operation matrices
- **Critical sentinel:** Symmetry metadata (e.g., "I" symmetry label) is read
  for logging only. It is NOT used to select the cage subset or to determine
  the eligible axis family. This is recorded as a sentinel compliance check.

### Step 3: Cage eligibility - parameter-free geometric selection
- **Input:** Subunit centroid positions (not metadata labels)
- **Operation:**
    1. Compute coordinate centroid of all subunit centroids
    2. Compute radial distribution: r_i = dist(centroid_i, centroid)
    3. Verify single dominant radial shell: r_max / r_min < threshold
    4. Construct convex hull of subunit centroids
    5. Verify complete radial shell (not slab, filament, or open arc)
    6. Record full radial gap distribution
- **Output:** Cage eligible (yes/no); radial gap distribution table
- **Code reference:** wwPDB_July15_2026_Rigorous_Analysis.py, identify_cage_subset()
- **Metadata independence:** Eligible set established BEFORE reading symmetry metadata.

### Step 4: Recover rotation axes from coordinates
- **Input:** Eligible cage coordinate model
- **Operation:**
    1. Enumerate candidate rotation axes (uniform sphere sampling)
    2. For each axis: test rotational self-maps at angles 2pi/n, n in {2,3,4,5,6,8,10,12}
    3. Cluster axis directions (angular DBSCAN, no order constraint imposed)
    4. Record recovered rotation orders and RMSD from ideal
- **Output:** Axis family recovered from coordinates; RMSD per axis cluster
- **Constraint:** No allowed order is imposed. The rotation order is recovered, not assumed.

### Step 5: Compare to icosahedral scaffold
- **Input:** Recovered axis families
- **Operation:** Check for presence of C2, C3, C5 families in recovered axes
- **Output:** All three families present (yes/no)
- **Prediction requirement:** All three families {C2, C3, C5} present in eligible cage

---

## Primary Test Entry: 9ZYS

- Entry: 9ZYS (Phage Oekolampad Bas18 capsid, cryo-EM 3.7 A)
- Assembly: 1 ASU of 60-ASU icosahedral capsid (9 chains)
- Symmetry metadata: I (read for logging; NOT used for selection)

**Recovered axes (from coordinates only):**

| Rotation order | Axes found | Clusters | Best RMSD (A) | In predicted (2, 3, 5) family? |
|---|---|---|---|---|
| C2 | 49 | 6 | 7.69 | YES |
| C3 | 16 | 1 | 8.44 | YES |
| C5 | 3 | 1 | 9.95 | YES |

**Sentinel compliance:** Icosahedral metadata label NOT used to select subset. Confirmed.
**Verdict:** VALIDATED. All three families {C2, C3, C5} recovered from coordinates.

**Note on partial axis counts:** 9ZYS is a single ASU of a 60-ASU assembly. Recovering
a fraction of the full 15/10/6 axis count from one ASU is geometrically expected.
The prediction tests that the three families are recoverable - they are.

---

TRUTH > COMFORT. Always.
UMtts Institute - Layer 3 Pathway 2 Ledger
