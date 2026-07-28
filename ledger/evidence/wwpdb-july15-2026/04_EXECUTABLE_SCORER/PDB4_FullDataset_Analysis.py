"""
PDB-4 Full-Dataset Analysis - wwPDB Release 2026-07-15
=======================================================
Queries ALL 293 released entries, retrieves UniProt entity cross-references,
downloads coordinate files + CDS for every entry with a UniProt linkage,
then runs a GLOBAL POOLED PDB-4 Kendall tau concordance test.

PDB-4 Prediction (MH_Origin.md / MH_PREDICTION_02):
  Within any amino acid identity class, synonymous codons with smaller
  |Orientation Delta| (from the source-fixed MH Z-cascade ledger) correspond
  to lower B-factor (greater structural order).

Global pooled test:
  - Accumulate all (aa, |ori_delta|, B_factor) triples across all entries
  - For each AA, enumerate all pairs (i,j) with different synonymous codons
  - Count concordant pairs: |delta_i| < |delta_j| AND B_i < B_j
  - Count discordant pairs: |delta_i| < |delta_j| AND B_i > B_j
  - Kendall tau = (concordant - discordant) / total_pairs
  - Prediction: tau > 0

All Orientation Delta values verbatim from MH_Origin.md Appendix C.
No free parameters. No fitted coefficients.
"""

import os, sys, json, math, time, collections, itertools, statistics, urllib.request, urllib.parse

PHI = 1.618033988749895

# Source-fixed Orientation Delta ledger - MH_Origin.md Appendix C
# key: RNA codon, value: (z_final, aa_symbol, orientation_delta)
CODON_LEDGER = {
    'AAA':(1.021173,   'K', 0.000000000), 'AAC':(1.018997,   'W', 0.014435533),
    'AAG':(1.021248,   'Y',-0.000494282), 'AAU':(1.022292,   'N',-0.007422841),
    'ACA':(1.011980,   'W', 0.007418690), 'ACC':(1.009823,   'A', 0.021724520),
    'ACG':(1.012054,   'H', 0.006928849), 'ACU':(1.013089,   'M', 0.000062543),
    'AGA':(1.021488,   'Y',-0.000254021), 'AGC':(1.019311,   'H', 0.014185953),
    'AGG':(1.021563,   'I',-0.000748455), 'AGU':(1.022608,   'R',-0.007679146),
    'AUA':(1.025901,   'N',-0.003814737), 'AUC':(1.023714,   'M', 0.010687490),
    'AUG':(1.025975,   'R',-0.004311303), 'AUU':(1.027025,   'Q',-0.011271872),
    'CAA':(0.982706773,'W',-0.021854223), 'CAC':(0.980612085,'A',-0.007485997),
    'CAG':(0.982778497,'H',-0.022346201), 'CAU':(0.983783875,'M',-0.029242455),
    'CCA':(0.973859559,'A',-0.014238523), 'CCC':(0.971783729,'G', 0.000000000),
    'CCG':(0.973930637,'C',-0.014726060), 'CCU':(0.974926963,'S',-0.021560060),
    'CGA':(0.983009708,'H',-0.022114990), 'CGC':(0.980914374,'C',-0.007742322),
    'CGG':(0.983081453,'P',-0.022607120), 'CGU':(0.984087141,'V',-0.029505505),
    'CUA':(0.987256066,'M',-0.025770264), 'CUC':(0.985151681,'S',-0.011335343),
    'CUG':(0.987328121,'V',-0.026264525), 'CUU':(0.988338154,'T',-0.033192790),
    'GAA':(1.022490,   'Y', 0.000748303), 'GAC':(1.020311,   'H', 0.015186141),
    'GAG':(1.022565,   'I', 0.000253942), 'GAU':(1.023611,   'R',-0.006675723),
    'GCA':(1.013285,   'H', 0.008160248), 'GCC':(1.011125,   'C', 0.022468382),
    'GCG':(1.013359,   'P', 0.007670328), 'GCU':(1.014396,   'V', 0.000802916),
    'GGA':(1.022806,   'I', 0.000494513), 'GGC':(1.020625,   'P', 0.014936792),
    'GGG':(1.022880,   'L', 0.000000000), 'GGU':(1.023927,   'D',-0.006931796),
    'GUA':(1.027224,   'R',-0.003062965), 'GUC':(1.025034,   'V', 0.011441567),
    'GUG':(1.027299,   'D',-0.003559610), 'GUU':(1.028350,   'F',-0.010521285),
    'UAA':(1.040953,   'N', 0.011237578), 'UAC':(1.038734,   'M', 0.025707720),
    'UAG':(1.041029,   'R', 0.010742110), 'UAU':(1.042094,   'Q', 0.003796940),
    'UCA':(1.031581,   'M', 0.018554964), 'UCC':(1.029382,   'S', 0.032895403),
    'UCG':(1.031657,   'V', 0.018063938), 'UCU':(1.032712,   'T', 0.011181021),
    'UGA':(1.041274,   'R', 0.010987026), 'UGC':(1.039054,   'V', 0.025461609),
    'UGG':(1.041350,   'D', 0.010491406), 'UGU':(1.042415,   'F', 0.003544104),
    'UUA':(1.045772,   'Q', 0.007474932), 'UUC':(1.043543,   'T', 0.022011769),
    'UUG':(1.045848,   'F', 0.006977181), 'UUU':(1.046918,   'E', 0.000000000),
}
STOP_CODONS = {'UAA','UAG','UGA'}

# Unique sorted synonymous codons per AA (for pair enumeration)
AA_CODONS = collections.defaultdict(list)
for cod,(z,aa,delta) in CODON_LEDGER.items():
    if cod not in STOP_CODONS:
        AA_CODONS[aa].append((abs(delta), cod))
for aa in AA_CODONS:
    AA_CODONS[aa].sort()  # ascending |delta|

THREE_TO_ONE = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
    'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
    'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
    'MSE':'M','SEP':'S','TPO':'T','PTR':'Y','HYP':'P','CSO':'C','HSD':'H',
}

OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "..", "wwPDB_July15_analysis")
os.makedirs(OUTDIR, exist_ok=True)

LOG = []
def log(m=""):
    print(m, flush=True)
    LOG.append(str(m))

def http_get(url, timeout=45):
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "MH-pdb4/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.read()
    except Exception as e:
        return None

def download(url, dest, label=""):
    if os.path.exists(dest) and os.path.getsize(dest) > 100:
        return dest
    data = http_get(url)
    if not data:
        return None
    with open(dest, 'wb') as f:
        f.write(data)
    return dest

# ===========================================================================
# STEP 1: Get all 293 entry IDs from July 15, 2026 release
# ===========================================================================

RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_ENTITY = "https://data.rcsb.org/rest/v1/core/polymer_entity/{pdbid}/{entity_id}"
RCSB_ENTRY  = "https://data.rcsb.org/rest/v1/core/entry/{pdbid}"
RCSB_CIF    = "https://files.rcsb.org/download/{pdbid}.cif"
NCBI_ESEARCH= ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
               "db=nuccore&term={q}&retmax=3&retmode=json")
NCBI_EFETCH = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
               "db=nuccore&id={acc}&rettype=fasta_cds_na&retmode=text")
UNIPROT_JSON= "https://rest.uniprot.org/uniprotkb/{uid}.json"
RCSB_GRAPHQL= "https://data.rcsb.org/graphql"

def get_all_release_ids():
    """Retrieve all PDB IDs released on 2026-07-15."""
    log("Querying RCSB for all July 15, 2026 release IDs...")
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_accession_info.initial_release_date",
                "operator": "equals",
                "value": "2026-07-15T00:00:00Z"
            }
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": 500},
            "results_verbosity": "compact"
        }
    }
    # Use POST
    import urllib.request as ur
    try:
        req = ur.Request(RCSB_SEARCH,
                         data=json.dumps(query).encode(),
                         headers={"Content-Type":"application/json","User-Agent":"MH/1.0"})
        with ur.urlopen(req, timeout=30) as r:
            data = r.read()
        result = json.loads(data.decode('utf-8'))
        raw    = result.get('result_set', [])
        # result_set is a flat list of PDB ID strings e.g. ['10DT', '10EC', ...]
        ids    = [h if isinstance(h, str) else h.get('identifier', h.get('id',''))
                  for h in raw]
        ids    = [x for x in ids if x]
        total  = result.get('total_count', len(ids))
        log(f"  Total entries: {total}, retrieved: {len(ids)}")
        return ids
    except Exception as e:
        log(f"  POST search failed: {e}")


    # Fallback: use the RCSB search with GET and different result format
    try:
        qstr = urllib.parse.quote(json.dumps(query))
        data = http_get(f"{RCSB_SEARCH}?json={qstr}")
        if data:
            result = json.loads(data.decode('utf-8'))
            # Some API versions return identifiers directly as list
            raw = result.get('result_set', result.get('identifiers', []))
            ids = []
            for h in raw:
                if isinstance(h, dict):
                    v = h.get('identifier') or h.get('entry_id') or ''
                    ids.append(v)
                elif isinstance(h, str):
                    ids.append(h)
            ids = [x for x in ids if x]
            if ids:
                log(f"  GET fallback: retrieved {len(ids)} IDs")
                return ids
    except Exception as e2:
        log(f"  GET fallback also failed: {e2}")

    # Last resort: use the RCSB monthly holdings endpoint
    try:
        data = http_get("https://data.rcsb.org/rest/v1/holdings/released/entries?date=2026-07-15")
        if data:
            j = json.loads(data.decode('utf-8'))
            ids = j if isinstance(j, list) else j.get('ids', j.get('entry_ids', []))
            log(f"  Holdings endpoint: {len(ids)} IDs")
            return ids
    except Exception as e3:
        log(f"  Holdings fallback failed: {e3}")
    return []


def get_entry_entities(pdbid):
    """Get polymer entity IDs for a PDB entry."""
    data = http_get(RCSB_ENTRY.format(pdbid=pdbid))
    if not data:
        return []
    try:
        j = json.loads(data.decode('utf-8'))
        return j.get('rcsb_entry_container_identifiers', {}).get('polymer_entity_ids', [])
    except:
        return []


def get_uniprot_from_entity(pdbid, entity_id):
    """Get UniProt accession for a PDB polymer entity."""
    data = http_get(RCSB_ENTITY.format(pdbid=pdbid, entity_id=entity_id))
    if not data:
        return None
    try:
        j = json.loads(data.decode('utf-8'))
        # rcsb_polymer_entity_align -> uniprot_ids
        aligns = j.get('rcsb_polymer_entity_align', [])
        for align in aligns:
            uid = align.get('reference_database_accession')
            db  = align.get('reference_database_name', '')
            if db == 'UniProt' and uid:
                return uid
        # fallback: entity_src_nat
        for src in j.get('entity_src_nat', []):
            pass
        return None
    except:
        return None


def get_cds_from_uniprot(uniprot_id, pdbid):
    """UniProt -> RefSeq NM_ -> NCBI CDS FASTA."""
    dest = os.path.join(OUTDIR, f"up_{uniprot_id}.json")
    path = download(UNIPROT_JSON.format(uid=uniprot_id), dest)
    if not path:
        return []
    try:
        with open(path) as f:
            j = json.load(f)
    except:
        return []
    # Find RefSeq NM_
    nm = None
    for xref in j.get('uniProtKBCrossReferences', []):
        if xref.get('database') == 'RefSeq':
            for prop in xref.get('properties', []):
                if prop.get('key') == 'NucleotideSequenceId':
                    v = prop.get('value','')
                    if v.startswith('NM_'):
                        nm = v.split('.')[0]; break
        if nm: break
    if not nm:
        return []
    return fetch_ncbi_cds(nm, pdbid, f"UP:{uniprot_id}")


def fetch_ncbi_cds(accession, pdbid, label=""):
    """Fetch CDS nucleotide FASTA from NCBI Efetch."""
    dest = os.path.join(OUTDIR, f"cds_{accession}.fasta")
    url  = NCBI_EFETCH.format(acc=accession)
    path = download(url, dest)
    if not path:
        return []
    return parse_cds_fasta(path)


def search_ncbi_cds_gene(gene, organism, pdbid):
    """Search NCBI for gene CDS - bacteria do NOT have mRNA filter."""
    q   = urllib.parse.quote(f"{gene}[gene] AND {organism}[organism]")
    url = NCBI_ESEARCH.format(q=q)
    data = http_get(url)
    if not data: return []
    try:
        j   = json.loads(data.decode('utf-8'))
        ids = j.get('esearchresult',{}).get('idlist',[])
    except:
        return []
    if not ids: return []
    return fetch_ncbi_cds(ids[0], pdbid, f"NCBI:{gene}/{organism}")


def parse_cds_fasta(fasta_path):
    """Parse CDS FASTA -> list of (rna_codon, 1-indexed position)."""
    seq_lines = []
    started = False
    with open(fasta_path, errors='replace') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if started and seq_lines:
                    break  # first CDS only
                started = True
                continue
            if started:
                seq_lines.append(line.upper())
    seq = ''.join(seq_lines).replace('T','U')
    out = []
    for i in range(0, len(seq)-2, 3):
        c = seq[i:i+3]
        if len(c)==3:
            out.append((c, i//3+1))
    return out


def codons_to_posmap(codon_list):
    """position -> (rna_codon, orientation_delta)"""
    m = {}
    for c,pos in codon_list:
        rna = c.upper().replace('T','U')
        if rna in STOP_CODONS: continue
        e = CODON_LEDGER.get(rna)
        if e:
            m[pos] = (rna, e[2])
    return m


# ===========================================================================
# STEP 2: Parse mmCIF CA + B-factor (native, no external dep for speed)
# ===========================================================================

def parse_ca_bfactor(cif_path):
    """
    Extract (chain_id, seq_id, res_name, ca_xyz, bfactor) for all ATOM CA records.
    Returns list of dicts. Uses gemmi if available, else native parser.
    """
    try:
        import gemmi
        st = gemmi.read_structure(cif_path)
        rows = []
        for chain in st[0]:
            for res in chain:
                if res.entity_type != gemmi.EntityType.Polymer:
                    continue
                ca = res.find_atom('CA', '\0')
                if ca:
                    rows.append({
                        'chain': chain.name,
                        'seq_id': res.seqid.num,
                        'res_name': res.name,
                        'aa': gemmi.find_tabulated_residue(res.name).one_letter_code
                              if gemmi.find_tabulated_residue(res.name) else THREE_TO_ONE.get(res.name,'X'),
                        'xyz': (ca.pos.x, ca.pos.y, ca.pos.z),
                        'bfactor': ca.b_iso,
                    })
        return rows
    except ImportError:
        pass

    # native fallback
    rows = []
    header = {}
    in_loop = False
    i = 0
    with open(cif_path, encoding='utf-8', errors='replace') as f:
        lines = f.readlines()
    while i < len(lines):
        line = lines[i].strip()
        if line == 'loop_':
            j = i+1; cols={}; idx=0
            while j < len(lines):
                l = lines[j].strip()
                if l.startswith('_atom_site.'):
                    cols[l]=idx; idx+=1; j+=1
                else: break
            if '_atom_site.id' in cols or '_atom_site.type_symbol' in cols:
                header=cols; in_loop=True; i=j; continue
            in_loop=False
        if in_loop and line and not line.startswith(('_','#','loop_','data_','save_')):
            t = line.split()
            def g(k):
                ix=header.get(k); return t[ix] if ix is not None and ix<len(t) else ''
            if g('_atom_site.group_PDB')=='ATOM' and g('_atom_site.label_atom_id')=='CA':
                try:
                    rows.append({
                        'chain': g('_atom_site.auth_asym_id') or g('_atom_site.label_asym_id'),
                        'seq_id': int(g('_atom_site.auth_seq_id') or g('_atom_site.label_seq_id')),
                        'res_name': g('_atom_site.label_comp_id'),
                        'aa': THREE_TO_ONE.get(g('_atom_site.label_comp_id'),'X'),
                        'xyz': (float(g('_atom_site.Cartn_x')),
                                float(g('_atom_site.Cartn_y')),
                                float(g('_atom_site.Cartn_z'))),
                        'bfactor': float(g('_atom_site.B_iso_or_equiv')),
                    })
                except: pass
        i+=1
    return rows


# ===========================================================================
# STEP 3: Collect (aa, |ori_delta|, B_factor) triples from one entry
# ===========================================================================

def collect_triples(cif_path, codon_map, pdbid):
    """
    Returns list of (aa_sym, abs_ori_delta, bfactor, codon_rna) for all
    residues that have both a CA coordinate and a codon annotation.
    """
    ca_rows = parse_ca_bfactor(cif_path)
    if not ca_rows:
        return []
    triples = []
    for row in ca_rows:
        pos = row['seq_id']
        if pos not in codon_map:
            continue
        codon_rna, ori_delta = codon_map[pos]
        entry = CODON_LEDGER.get(codon_rna)
        if not entry:
            continue
        aa = entry[1]  # MH source AA symbol
        # Cross-check structure AA matches MH assignment (codon -> AA must agree with structure)
        struct_aa = row['aa']
        # Accept if MH AA matches structure AA (1-letter), or structure reports X/unknown
        if struct_aa != 'X' and struct_aa != aa:
            continue  # skip mismatch - frame error or selenomethionine etc
        triples.append((aa, abs(ori_delta), row['bfactor'], codon_rna))
    return triples


# ===========================================================================
# STEP 4: Global pooled Kendall tau PDB-4 test
# ===========================================================================

def kendall_tau_pdb4(all_triples):
    """
    Global pooled Kendall tau test for PDB-4.

    Within each amino acid identity:
      For every ordered pair (i, j) where codon_i != codon_j:
        concordant if |delta_i| < |delta_j| AND B_i < B_j  (prediction direction)
        discordant if |delta_i| < |delta_j| AND B_i > B_j
        tied if B_i == B_j (excluded)

    tau = (C - D) / (C + D)
    Prediction: tau > 0
    """
    # Group by AA
    aa_groups = collections.defaultdict(list)
    for aa, abs_delta, bfac, codon in all_triples:
        aa_groups[aa].append((abs_delta, bfac, codon))

    total_C = 0
    total_D = 0
    total_tied = 0
    per_aa = {}

    log(f"\n  Per-AA global pooled Kendall tau:")
    log(f"  {'AA':<4} {'n_residues':>10} {'n_codons':>9} {'pairs':>8} {'C':>8} {'D':>8} {'tied':>6} {'tau':>8}")
    log(f"  {'-'*70}")

    for aa in sorted(aa_groups.keys()):
        items = aa_groups[aa]
        if len(items) < 2:
            continue
        # Check if this AA actually has multiple synonymous codons observed
        codons_obs = set(it[2] for it in items)
        if len(codons_obs) < 2:
            continue

        C = D = tied = 0
        # All ordered pairs (i,j) where codon_i != codon_j
        for i in range(len(items)):
            for j in range(i+1, len(items)):
                d_i, b_i, c_i = items[i]
                d_j, b_j, c_j = items[j]
                if c_i == c_j:
                    continue  # same codon, not a synonymous contrast
                # Orient so d_i <= d_j (WLOG by symmetry)
                if d_i > d_j:
                    d_i, d_j = d_j, d_i
                    b_i, b_j = b_j, b_i
                # Prediction: d_i < d_j -> B_i < B_j (concordant)
                if d_i == d_j:
                    tied += 1
                elif b_i < b_j:
                    C += 1
                elif b_i > b_j:
                    D += 1
                else:
                    tied += 1

        if (C + D) == 0:
            continue
        tau_aa = (C - D) / (C + D)
        n_codons = len(codons_obs)
        log(f"  {aa:<4} {len(items):>10} {n_codons:>9} {C+D:>8} {C:>8} {D:>8} {tied:>6} {tau_aa:>8.4f}")
        total_C += C; total_D += D; total_tied += tied
        per_aa[aa] = {'n': len(items), 'n_codons': n_codons,
                      'C': C, 'D': D, 'tied': tied, 'tau': tau_aa}

    log(f"  {'-'*70}")
    total_pairs = total_C + total_D
    if total_pairs == 0:
        log("  No pairs available."); return None

    tau_global = (total_C - total_D) / total_pairs
    log(f"  {'GLOBAL':<4} {'':>10} {'':>9} {total_pairs:>8} {total_C:>8} {total_D:>8} {total_tied:>6} {tau_global:>8.4f}")
    log(f"\n  Global Kendall tau = {tau_global:.6f}")
    log(f"  Total concordant (C) = {total_C}")
    log(f"  Total discordant (D) = {total_D}")
    log(f"  Total tied           = {total_tied}")
    log(f"  C / (C+D)            = {total_C/total_pairs:.4f}  (0.50 = null; >0.50 = predicted)")
    log(f"  Prediction: tau > 0, i.e. concordant fraction > 0.50")
    log(f"  Observed:   tau = {tau_global:+.6f}  -> concordant fraction = {total_C/total_pairs:.4f}")

    if tau_global > 0:
        log(f"\n  -> PDB-4 VALIDATED: Global Kendall tau > 0  "
            f"(C={total_C} > D={total_D}, fraction={total_C/total_pairs:.4f})")
    elif tau_global < 0:
        log(f"\n  -> PDB-4 FALSIFIED: Global Kendall tau < 0  "
            f"(C={total_C} < D={total_D})")
    else:
        log(f"\n  -> PDB-4 EXACTLY NULL: tau = 0")

    # Permutation control: shuffle B-factors within each AA group
    import random; random.seed(42)
    ctrl_C = ctrl_D = 0
    for aa, items in aa_groups.items():
        bfacs = [it[1] for it in items]
        random.shuffle(bfacs)
        shuffled = [(items[k][0], bfacs[k], items[k][2]) for k in range(len(items))]
        codons_obs = set(it[2] for it in shuffled)
        if len(codons_obs) < 2: continue
        for i in range(len(shuffled)):
            for j in range(i+1, len(shuffled)):
                d_i,b_i,c_i = shuffled[i]; d_j,b_j,c_j = shuffled[j]
                if c_i==c_j: continue
                if d_i>d_j: d_i,d_j=d_j,d_i; b_i,b_j=b_j,b_i
                if d_i==d_j: continue
                if b_i<b_j: ctrl_C+=1
                elif b_i>b_j: ctrl_D+=1
    ctrl_total = ctrl_C+ctrl_D
    ctrl_tau = (ctrl_C-ctrl_D)/ctrl_total if ctrl_total else 0
    log(f"\n  Permutation control (shuffled B-factors within AA): tau = {ctrl_tau:+.6f}")
    log(f"  Real tau = {tau_global:+.6f}  vs  shuffled = {ctrl_tau:+.6f}")
    log(f"  Real tau exceeds control in predicted direction: {tau_global > ctrl_tau}")

    return {
        'tau_global': tau_global, 'C': total_C, 'D': total_D,
        'tied': total_tied, 'n_pairs': total_pairs,
        'concordant_fraction': total_C/total_pairs,
        'permutation_tau': ctrl_tau,
        'per_aa': per_aa
    }


# ===========================================================================
# STEP 5: Main - process full 293-entry release
# ===========================================================================

def main():
    log("="*70)
    log("PDB-4 FULL DATASET ANALYSIS - wwPDB July 15, 2026 (293 entries)")
    log("="*70)
    log(f"phi = {PHI}")
    log(f"Ledger: MH_Origin.md Appendix C (source-fixed, 64 codons)")
    log(f"Test: Global pooled Kendall tau across all synonymous pairs")
    log(f"Output dir: {OUTDIR}")
    log()

    # Step 1: Get all IDs
    all_ids = get_all_release_ids()
    if not all_ids:
        log("FATAL: Could not retrieve release IDs."); return
    log(f"Processing {len(all_ids)} entries...\n")

    all_triples = []
    processed   = 0
    skipped_no_uniprot = 0
    skipped_no_cds     = 0
    skipped_no_triples = 0
    entries_with_data  = 0

    for idx, pdbid in enumerate(all_ids):
        time.sleep(0.05)  # polite rate limiting
        log(f"[{idx+1:3d}/{len(all_ids)}] {pdbid}", )

        # Get entity list
        entity_ids = get_entry_entities(pdbid)
        if not entity_ids:
            log(f"  -> no entities")
            skipped_no_uniprot += 1
            continue

        # Find UniProt for first protein entity with a mapping
        uniprot_id = None
        for eid in entity_ids:
            uid = get_uniprot_from_entity(pdbid, eid)
            if uid:
                uniprot_id = uid
                break

        if not uniprot_id:
            log(f"  -> no UniProt linkage")
            skipped_no_uniprot += 1
            continue

        log(f"  UniProt: {uniprot_id}")

        # Get CDS
        codon_list = get_cds_from_uniprot(uniprot_id, pdbid)
        if not codon_list:
            log(f"  -> no CDS (no RefSeq NM_ in UniProt)")
            skipped_no_cds += 1
            continue

        codon_map = codons_to_posmap(codon_list)
        if not codon_map:
            log(f"  -> empty codon map")
            skipped_no_cds += 1
            continue

        # Download CIF
        cif_url  = RCSB_CIF.format(pdbid=pdbid.upper())
        cif_dest = os.path.join(OUTDIR, f"{pdbid.upper()}.cif")
        cif_path = download(cif_url, cif_dest)
        if not cif_path:
            log(f"  -> CIF download failed")
            skipped_no_triples += 1
            continue

        # Collect triples
        triples = collect_triples(cif_path, codon_map, pdbid)
        if not triples:
            log(f"  -> no codon-annotated residues")
            skipped_no_triples += 1
            continue

        n_aa = len(set(t[0] for t in triples))
        log(f"  -> {len(triples)} residues, {n_aa} unique AAs with codon annotation")
        all_triples.extend(triples)
        entries_with_data += 1
        processed += 1

    log(f"\n{'='*70}")
    log(f"Data collection summary:")
    log(f"  Entries in release:             {len(all_ids)}")
    log(f"  Entries with UniProt linkage:   {len(all_ids) - skipped_no_uniprot}")
    log(f"  Entries with CDS:               {processed + entries_with_data - processed + processed}")
    log(f"  Entries contributing triples:   {entries_with_data}")
    log(f"  Total (aa, |delta|, B) triples: {len(all_triples)}")
    log(f"  Skipped - no UniProt:           {skipped_no_uniprot}")
    log(f"  Skipped - no CDS:               {skipped_no_cds}")
    log(f"  Skipped - no coord triples:     {skipped_no_triples}")

    if not all_triples:
        log("\nFATAL: No triples accumulated. Cannot run PDB-4."); return

    # Run global pooled Kendall tau
    log(f"\n{'='*70}")
    log(f"PDB-4: GLOBAL POOLED KENDALL TAU - Full July 15, 2026 Dataset")
    log(f"{'='*70}")
    result = kendall_tau_pdb4(all_triples)

    # Save report
    rpt = os.path.join(OUTDIR, "PDB4_FullDataset_Report.txt")
    with open(rpt, 'w', encoding='utf-8') as f:
        f.write('\n'.join(LOG))
    log(f"\nReport saved: {rpt}")

    if result:
        log(f"\nFINAL PDB-4 VERDICT:")
        log(f"  tau = {result['tau_global']:+.6f}")
        log(f"  C/(C+D) = {result['concordant_fraction']:.4f}")
        if result['tau_global'] > 0:
            log(f"  PDB-4 VALIDATED across full July 15, 2026 wwPDB release dataset.")
        else:
            log(f"  PDB-4 NOT validated.")
        log(f"\n  TRUTH > COMFORT. Always.")


if __name__ == '__main__':
    main()
