"""
PDB-4 Sign-Corrected Test - Prediction 3 from PDB4_WN_Dissent_Derivation.md
=============================================================================
Re-runs the full-dataset Kendall tau using SIGNED Orientation Delta as the
predictor instead of |Delta|.

Standard PDB-4:  |Delta_i| < |Delta_j| -> B_i < B_j  (unsigned, tau > 0 predicted)
Corrected PDB-4: Delta_i  < Delta_j   -> B_i < B_j  (signed,   tau > 0 predicted)

The derivation (PDB4_WN_Dissent_Derivation.md) shows:
  - W (T2, face 5): CAA is sub-basin (Z<1, Delta=-0.022), should be MOST ordered
    (lowest B). With |Delta|, CAA gets ranked as LEAST ordered. Signed Delta fixes this.
  - N (G2, face 15): Only codons are AAU (Delta=-0.0074) and AUA (Delta=-0.0038).
    Both are below-anchor. The more-negative (AAU) should be more constrained.
    With |Delta|, AAU gets ranked as LESS ordered. Signed Delta fixes this.

Prediction 3: With signed Delta, tau for W and N should flip from negative to positive.
              Global tau should increase above the unsigned result (+0.035).

All data files are cached from the prior run - no network calls needed for most entries.
"""

import os, sys, json, time, collections, statistics
import urllib.request, urllib.parse

PHI = 1.618033988749895
BASIN_FLOOR_Z = 1.0  # Z = vx/(c/phi^9); sub-basin = Z < 1

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
THREE_TO_ONE = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
    'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
    'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
    'MSE':'M','SEP':'S','TPO':'T','PTR':'Y','HYP':'P','CSO':'C','HSD':'H',
}

OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "..", "wwPDB_July15_analysis")

LOG = []
def log(m=""):
    print(m, flush=True)
    LOG.append(str(m))

def http_get(url, timeout=45):
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "MH/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.read()
    except:
        return None

# ── Sub-basin flag ──────────────────────────────────────────────────────────
# A codon is sub-basin if its Z < 1.000 (vx < c/phi^9)
SUB_BASIN_CODONS = {cod for cod,(z,aa,d) in CODON_LEDGER.items() if z < 1.0}
log(f"Sub-basin codons ({len(SUB_BASIN_CODONS)}): {sorted(SUB_BASIN_CODONS)}")

# ── Reuse all cached data (no downloads if CIF+CDS already present) ─────────
RCSB_SEARCH  = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_ENTITY  = "https://data.rcsb.org/rest/v1/core/polymer_entity/{pdbid}/{entity_id}"
RCSB_ENTRY   = "https://data.rcsb.org/rest/v1/core/entry/{pdbid}"
RCSB_CIF     = "https://files.rcsb.org/download/{pdbid}.cif"
NCBI_EFETCH  = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
                "db=nuccore&id={acc}&rettype=fasta_cds_na&retmode=text")
UNIPROT_JSON = "https://rest.uniprot.org/uniprotkb/{uid}.json"

def download(url, dest):
    if os.path.exists(dest) and os.path.getsize(dest) > 100:
        return dest
    data = http_get(url)
    if not data: return None
    with open(dest,'wb') as f: f.write(data)
    return dest

def get_all_release_ids():
    query = {
        "query":{"type":"terminal","service":"text","parameters":{
            "attribute":"rcsb_accession_info.initial_release_date",
            "operator":"equals","value":"2026-07-15T00:00:00Z"}},
        "return_type":"entry",
        "request_options":{"paginate":{"start":0,"rows":500},"results_verbosity":"compact"}
    }
    import urllib.request as ur
    req = ur.Request(RCSB_SEARCH, data=json.dumps(query).encode(),
                     headers={"Content-Type":"application/json","User-Agent":"MH/1.0"})
    with ur.urlopen(req, timeout=30) as r:
        result = json.loads(r.read().decode())
    raw = result.get('result_set', [])
    return [h if isinstance(h,str) else h.get('identifier','') for h in raw]

def get_entry_entities(pdbid):
    data = http_get(RCSB_ENTRY.format(pdbid=pdbid))
    if not data: return []
    try:
        return json.loads(data).get('rcsb_entry_container_identifiers',{}).get('polymer_entity_ids',[])
    except: return []

def get_uniprot(pdbid, entity_ids):
    for eid in entity_ids:
        data = http_get(RCSB_ENTITY.format(pdbid=pdbid, entity_id=eid))
        if not data: continue
        try:
            for align in json.loads(data).get('rcsb_polymer_entity_align',[]):
                if align.get('reference_database_name')=='UniProt':
                    uid = align.get('reference_database_accession')
                    if uid: return uid
        except: continue
    return None

def get_cds(uniprot_id, pdbid):
    dest = os.path.join(OUTDIR, f"up_{uniprot_id}.json")
    path = download(UNIPROT_JSON.format(uid=uniprot_id), dest)
    if not path: return []
    try:
        with open(path) as f: j = json.load(f)
    except: return []
    nm = None
    for xref in j.get('uniProtKBCrossReferences',[]):
        if xref.get('database')=='RefSeq':
            for prop in xref.get('properties',[]):
                if prop.get('key')=='NucleotideSequenceId':
                    v = prop.get('value','')
                    if v.startswith('NM_'): nm = v.split('.')[0]; break
        if nm: break
    if not nm: return []
    dest2 = os.path.join(OUTDIR, f"cds_{nm}.fasta")
    path2 = download(NCBI_EFETCH.format(acc=nm), dest2)
    if not path2: return []
    return parse_cds(path2)

def parse_cds(fasta_path):
    seq_lines=[]; started=False
    with open(fasta_path, errors='replace') as f:
        for line in f:
            line=line.strip()
            if line.startswith('>'):
                if started and seq_lines: break
                started=True; continue
            if started: seq_lines.append(line.upper())
    seq=''.join(seq_lines).replace('T','U')
    return [(seq[i:i+3], i//3+1) for i in range(0,len(seq)-2,3) if len(seq[i:i+3])==3]

def codons_to_posmap(codon_list):
    m={}
    for c,pos in codon_list:
        rna=c.upper().replace('T','U')
        if rna in STOP_CODONS: continue
        e=CODON_LEDGER.get(rna)
        if e: m[pos]=(rna, e[2])  # (codon_rna, signed_delta)
    return m

def parse_ca_bfactor(cif_path):
    try:
        import gemmi
        st = gemmi.read_structure(cif_path)
        rows=[]
        for chain in st[0]:
            for res in chain:
                if res.entity_type != gemmi.EntityType.Polymer: continue
                ca = res.find_atom('CA','\0')
                if ca:
                    tr = gemmi.find_tabulated_residue(res.name)
                    aa = tr.one_letter_code if tr else THREE_TO_ONE.get(res.name,'X')
                    rows.append({'seq_id':res.seqid.num,'aa':aa,'bfactor':ca.b_iso})
        return rows
    except: pass
    # native fallback
    rows=[]; header={}; in_loop=False; i=0
    with open(cif_path,encoding='utf-8',errors='replace') as f: lines=f.readlines()
    while i<len(lines):
        line=lines[i].strip()
        if line=='loop_':
            j=i+1; cols={}; idx=0
            while j<len(lines):
                l=lines[j].strip()
                if l.startswith('_atom_site.'): cols[l]=idx; idx+=1; j+=1
                else: break
            if '_atom_site.id' in cols or '_atom_site.type_symbol' in cols:
                header=cols; in_loop=True; i=j; continue
            in_loop=False
        if in_loop and line and not line.startswith(('_','#','loop_','data_','save_')):
            t=line.split()
            def g(k): ix=header.get(k); return t[ix] if ix is not None and ix<len(t) else ''
            if g('_atom_site.group_PDB')=='ATOM' and g('_atom_site.label_atom_id')=='CA':
                try:
                    rows.append({'seq_id':int(g('_atom_site.auth_seq_id') or g('_atom_site.label_seq_id')),
                                 'aa':THREE_TO_ONE.get(g('_atom_site.label_comp_id'),'X'),
                                 'bfactor':float(g('_atom_site.B_iso_or_equiv'))})
                except: pass
        i+=1
    return rows

def collect_triples(cif_path, codon_map):
    """Returns (aa, signed_delta, bfactor, codon_rna) for each annotated residue."""
    rows = parse_ca_bfactor(cif_path)
    triples = []
    for row in rows:
        pos = row['seq_id']
        if pos not in codon_map: continue
        codon_rna, ori_delta = codon_map[pos]
        entry = CODON_LEDGER.get(codon_rna)
        if not entry: continue
        mh_aa = entry[1]
        if row['aa'] != 'X' and row['aa'] != mh_aa: continue
        triples.append((mh_aa, ori_delta, row['bfactor'], codon_rna))
    return triples

# ── Dual Kendall tau ─────────────────────────────────────────────────────────

def kendall_tau(all_triples, use_signed):
    """
    use_signed=False: concordant if |delta_i| < |delta_j| AND B_i < B_j  (original)
    use_signed=True : concordant if  delta_i  < delta_j  AND B_i < B_j  (corrected)
    """
    label = "SIGNED Delta" if use_signed else "UNSIGNED |Delta|"
    log(f"\n{'='*70}")
    log(f"Kendall tau - {label}")
    log(f"{'='*70}")

    aa_groups = collections.defaultdict(list)
    for aa, signed_d, bfac, codon in all_triples:
        predictor = signed_d if use_signed else abs(signed_d)
        aa_groups[aa].append((predictor, bfac, codon))

    total_C=0; total_D=0; total_tied=0
    per_aa={}

    log(f"  {'AA':<4} {'n':>8} {'n_cod':>6} {'pairs':>9} {'C':>9} {'D':>9} {'tau':>9}")
    log(f"  {'-'*65}")

    for aa in sorted(aa_groups.keys()):
        items = aa_groups[aa]
        codons_obs = set(it[2] for it in items)
        if len(codons_obs) < 2: continue
        C=D=tied=0
        for i in range(len(items)):
            for j in range(i+1,len(items)):
                d_i,b_i,c_i = items[i]; d_j,b_j,c_j = items[j]
                if c_i==c_j: continue
                if d_i > d_j: d_i,d_j=d_j,d_i; b_i,b_j=b_j,b_i
                if d_i==d_j: tied+=1
                elif b_i<b_j: C+=1
                elif b_i>b_j: D+=1
                else: tied+=1
        if (C+D)==0: continue
        tau_aa=(C-D)/(C+D)
        log(f"  {aa:<4} {len(items):>8} {len(codons_obs):>6} {C+D:>9} {C:>9} {D:>9} {tau_aa:>9.4f}")
        total_C+=C; total_D+=D; total_tied+=tied
        per_aa[aa]={'C':C,'D':D,'tau':tau_aa}

    log(f"  {'-'*65}")
    tp=total_C+total_D
    if tp==0: log("  No pairs."); return None
    tau_g=(total_C-total_D)/tp
    log(f"  {'GLOBAL':<4} {'':>8} {'':>6} {tp:>9} {total_C:>9} {total_D:>9} {tau_g:>9.6f}")
    log(f"\n  tau = {tau_g:+.6f}   C/(C+D) = {total_C/tp:.4f}  (null=0.5000)")
    log(f"  C={total_C}  D={total_D}  margin={total_C-total_D:+d}")

    # Permutation control
    import random; random.seed(42)
    ctrl_C=ctrl_D=0
    for aa, items in aa_groups.items():
        bfacs=[it[1] for it in items]; random.shuffle(bfacs)
        sh=[(items[k][0],bfacs[k],items[k][2]) for k in range(len(items))]
        codons_obs=set(it[2] for it in sh)
        if len(codons_obs)<2: continue
        for i in range(len(sh)):
            for j in range(i+1,len(sh)):
                d_i,b_i,c_i=sh[i]; d_j,b_j,c_j=sh[j]
                if c_i==c_j: continue
                if d_i>d_j: d_i,d_j=d_j,d_i; b_i,b_j=b_j,b_i
                if d_i==d_j: continue
                if b_i<b_j: ctrl_C+=1
                elif b_i>b_j: ctrl_D+=1
    ctrl_t=(ctrl_C-ctrl_D)/(ctrl_C+ctrl_D) if (ctrl_C+ctrl_D) else 0
    log(f"  Permutation control: tau = {ctrl_t:+.6f}")
    log(f"  Real tau {tau_g:+.6f} vs shuffled {ctrl_t:+.6f} -> real exceeds control: {tau_g > ctrl_t}")

    if tau_g > 0:
        log(f"\n  -> {label}: VALIDATED  (tau > 0)")
    else:
        log(f"\n  -> {label}: NOT VALIDATED  (tau <= 0)")

    return {'tau':tau_g,'C':total_C,'D':total_D,'frac':total_C/tp,'ctrl_tau':ctrl_t,'per_aa':per_aa}

# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    log("="*70)
    log("PDB-4 SIGN-CORRECTED TEST - wwPDB July 15, 2026")
    log("Prediction 3 from PDB4_WN_Dissent_Derivation.md")
    log("="*70)
    log(f"Basin floor Z = {BASIN_FLOOR_Z}")
    log(f"Sub-basin codons (Z < 1.0): {sorted(SUB_BASIN_CODONS)}")
    log()

    all_ids = get_all_release_ids()
    log(f"Release entries: {len(all_ids)}")

    all_triples=[]
    for idx, pdbid in enumerate(all_ids):
        time.sleep(0.02)
        entity_ids = get_entry_entities(pdbid)
        uniprot_id = get_uniprot(pdbid, entity_ids) if entity_ids else None
        if not uniprot_id: continue
        codon_list = get_cds(uniprot_id, pdbid)
        if not codon_list: continue
        codon_map = codons_to_posmap(codon_list)
        if not codon_map: continue
        cif_dest = os.path.join(OUTDIR, f"{pdbid.upper()}.cif")
        cif_path = download(RCSB_CIF.format(pdbid=pdbid.upper()), cif_dest)
        if not cif_path: continue
        triples = collect_triples(cif_path, codon_map)
        if triples:
            all_triples.extend(triples)
            log(f"[{idx+1:3d}] {pdbid}: +{len(triples)} triples  (total={len(all_triples)})")

    log(f"\nTotal triples: {len(all_triples)}")

    # Run both versions for direct comparison
    r_unsigned = kendall_tau(all_triples, use_signed=False)
    r_signed   = kendall_tau(all_triples, use_signed=True)

    log(f"\n{'='*70}")
    log(f"COMPARISON SUMMARY")
    log(f"{'='*70}")
    log(f"  Metric              Unsigned |Delta|    Signed Delta")
    log(f"  tau (global)        {r_unsigned['tau']:+.6f}         {r_signed['tau']:+.6f}")
    log(f"  C/(C+D)             {r_unsigned['frac']:.4f}              {r_signed['frac']:.4f}")
    log(f"  permutation ctrl    {r_unsigned['ctrl_tau']:+.6f}         {r_signed['ctrl_tau']:+.6f}")
    log()

    # W and N specifically
    for aa in ['W','N']:
        ru = r_unsigned['per_aa'].get(aa,{})
        rs = r_signed['per_aa'].get(aa,{})
        log(f"  {aa}: unsigned tau={ru.get('tau','N/A'):>+.4f}  signed tau={rs.get('tau','N/A'):>+.4f}  "
            f"flip={'YES' if ru and rs and rs.get('tau',0)>0 and ru.get('tau',0)<0 else 'NO'}")

    log(f"\n  Prediction 3: tau(W,signed) > 0 and tau(N,signed) > 0")
    w_flip = r_signed['per_aa'].get('W',{}).get('tau',0) > 0
    n_flip = r_signed['per_aa'].get('N',{}).get('tau',0) > 0
    log(f"  W flipped to positive: {w_flip}")
    log(f"  N flipped to positive: {n_flip}")
    log(f"  Both flipped: {w_flip and n_flip}")
    if w_flip and n_flip:
        log(f"\n  -> Prediction 3: CONFIRMED")
    elif w_flip or n_flip:
        log(f"\n  -> Prediction 3: PARTIAL")
    else:
        log(f"\n  -> Prediction 3: NOT CONFIRMED")

    log(f"\n  TRUTH > COMFORT. Always.")

    rpt=os.path.join(OUTDIR,"PDB4_SignCorrected_Report.txt")
    with open(rpt,'w',encoding='utf-8') as f:
        f.write('\n'.join(LOG))
    log(f"\nReport: {rpt}")

if __name__=='__main__':
    main()
