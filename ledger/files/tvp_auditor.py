#!/usr/bin/env python3
"""
TVP v1.9 Universal Forensic Auditor
Mass Harmonics psim Framework | UMTTS Institute
TRUTH > COMFORT. Always.

Takes ANY 2 of 3 focal variables {R_km, f_Hz, vx_kms} per object.
Derives the third. Runs all 9 calculations. Outputs 30-column ledger.

Engine only. No embedded datasets. Generic CSV pipeline.
Domain data (elements, nuclides, exoplanets, pulsars, anything)
lives in external CSVs and is fed in through --input.

Single-object usage (supply any 2 of 3 focal variables; engine derives the third):
    # Map mode -- R + f both independently measured; engine solves local substrate vx
    python tvp_auditor.py --name "Proton-map"  --f_Hz 2.27e23 --R_km 8.41e-19
    python tvp_auditor.py --name "Crab Pulsar" --R_km 10 --f_Hz 29.947

    # Predict mode -- R known, vx_kms explicitly supplied as bare/vac reference
    python tvp_auditor.py --name "Earth PCF"   --R_km 6371 --vx_kms 299792.458

    # Detect mode -- f known, vx_kms explicitly supplied as bare/vac reference
    python tvp_auditor.py --name "Schumann"    --f_Hz 7.83 --vx_kms 299792.458

Batch CSV usage (ordinary terrain records; no special dataset handling):
    python tvp_auditor.py --input objects.csv
    python tvp_auditor.py --input elements.csv
    python tvp_auditor.py --input nuclides.csv
    python tvp_auditor.py --input objects.csv --output results.csv

Input CSV columns (any 2 of 3 focal variables, leave missing column blank):
    name, R_km, f_Hz, vx_kms, note, fault_hint

Unit-helper CSV columns also supported (preserved verbatim):
    R_pm, R_fm, R_m, f_MHz, f_GHz, f_eV, f_MeV

Example element row:
    name,R_pm,f_eV,note,fault_hint
    Hydrogen (H Z=1),31,13.5984,"geo=measured Z=1 IE=13.5984eV","C1"

Example nuclide row:
    name,R_fm,f_MeV,note,fault_hint
    Uranium-238 (Z=92 A=238),5.8571,1801.697,"src=Angeli2013 BE=1801.697MeV","C1"

Unit-helper CLI flags (single-object mode):
    --R_pm   (picometers  -> km)
    --R_fm   (femtometers -> km)
    --R_m    (meters      -> km)
    --f_MHz  (MHz -> Hz)
    --f_GHz  (GHz -> Hz)
    --f_eV   (eV  -> Hz via h)
    --f_MeV  (MeV -> Hz via h)
"""

import argparse
import csv
import io
import math
import os
import sys
from datetime import datetime

# ── Substrate constants ───────────────────────────────────────────────────────
VX_VAC_REF = 299_792.458             # km/s | bare substrate ref at Z=1 vacuum limit | MFE-derived
KPSIM      = VX_VAC_REF / (2*math.pi) # 47713.4516 Hz·km | Giboney Gradient vacuum coupling constant
PHI        = (1 + math.sqrt(5)) / 2
ALPHA      = 7.2973525693e-3
EV_HZ      = 2.417989242e14           # 1 eV in Hz  (h = 4.135667696e-15 eV·s, exact)
MEV_HZ     = EV_HZ * 1e6

# ── Unit conversions ──────────────────────────────────────────────────────────
PM_KM = 1e-15
FM_KM = 1e-18
M_KM  = 1e-3

# ── Ratio factorization ───────────────────────────────────────────────────────
RATIO_SET = {
    # Fine-structure eigenvalues
    'alpha':            ALPHA,
    'alpha^2':          ALPHA**2,
    'alpha^3':          ALPHA**3,           # HeNe laser cavity signal
    '1/alpha':          1/ALPHA,
    '2/alpha':          2/ALPHA,
    '(1/alpha)^5':      (1/ALPHA)**5,       # heavy-element NMR signature
    'pi/(2alpha)':      math.pi/(2*ALPHA),
    # Circular constants
    'pi':               math.pi,
    'pi/2':             math.pi/2,
    'pi^2':             math.pi**2,
    '2pi':              2*math.pi,
    # Icosahedral eigenvalues
    'phi':              PHI,
    'phi^2':            PHI**2,
    'phi^3':            PHI**3,
    'phi^6':            PHI**6,
    'phi^9':            PHI**9,
    # Polyphonic composites
    '6pi^5':            6*math.pi**5,
    '(6pi^5)^3':        (6*math.pi**5)**3,            # max nuclear coherence (Fe)
    '(6pi^5)^3*phi^3':  (6*math.pi**5)**3 * PHI**3,   # Fe cross-product signature
    # Integer scaffolding
    '1': 1.0, '2': 2.0, '4': 4.0,
}
def factorize(ratio):
    if ratio is None or ratio <= 0:
        return "no_ratio" 
    best_delta, best_name = 1e9, None
    for name, val in RATIO_SET.items():
        if val <= 0:
            continue
        delta = abs(ratio - val) / val 
        if delta < best_delta:
            best_delta, best_name = delta, name 
    if best_name is None:
        return "no_basis_available" 
    delta_pct = best_delta * 100
    if best_delta < 0.0005:
        return f"basis_lock {best_name} (delta={delta_pct:.2f}%)"  
    if best_delta < 0.5:
        return f"basis_near {best_name} (delta={delta_pct:.2f}%)"
    return f"basis_scan nearest={best_name} (delta={delta_pct:.1f}%)"

def p3gg_active():
    """
    All five P3GG harmonic voices are simultaneously active in every substrate interaction.
    The substrate does not switch. There are zero regimes. No voice switches off.
    beta_n = phi^(3*(n-1)) — icosahedral eigenvalues, all present always.
    """
    b = [PHI**(3*(n-1)) for n in range(1, 6)]
    return (f'All P3GG voices active: '
            f'b1={b[0]:.3f} b2={b[1]:.3f} b3={b[2]:.3f} b4={b[3]:.3f} b5={b[4]:.3f}')

def fmt(v, p=6):
    if v is None:
        return 'null'
    if isinstance(v, float):
        return f'{v:.{p}g}'
    return str(v)

# ── Core engine ───────────────────────────────────────────────────────────────
def nine(R_km, f_Hz, vx_kms):
    """
    All 9 calculations. All three focal variables must be provided (one may be derived).
    TVP v1.9: 3 topologies × 3 modes = 9.
    Torus: T2-Torus first (R_minor excavation), T1-Torus back-populates.

    Guards: R, f, vx must be nonzero. Zero or negative inputs raise ValueError;
    fault codes annotate, they do not null cells, but a zero-divide is not a fault —
    it is missing data and must surface to the auditor.
    """
    vx = vx_kms
    R  = R_km
    f  = f_Hz

    if R is None or R <= 0:
        raise ValueError(f'R_km must be > 0; got {R}')
    if f is None or f <= 0:
        raise ValueError(f'f_Hz must be > 0; got {f}')
    if vx is None or vx <= 0:
        raise ValueError(f'vx_kms must be > 0; got {vx}')

    # ── T1: Predict f from R and vx ──────────────────────────────────────────
    f1s = vx / (2*math.pi*R)      # sphere
    f1l = vx / (2*R)              # slab
    Rm  = vx / (2*math.pi*f)      # T2-Torus first: R_minor excavated
    f1t = vx / (2*math.pi*Rm)     # T1-Torus back-populates (= f, always)

    # ── T2: Recover geometry from f and vx ───────────────────────────────────
    R2s = vx / (2*math.pi*f)      # sphere
    t2l = vx / (2*f)              # slab (recovers thickness t)
    # Rm already computed above (torus R_minor)

    # ── T3: Map vx from R and f ───────────────────────────────────────────────
    V3s = 2*math.pi*R*f           # sphere
    V3l = 2*R*f                   # slab
    V3t = 2*math.pi*Rm*f          # torus

    # ── Closure % ────────────────────────────────────────────────────────────
    def cl(obs, pred):
        return (obs/pred*100) if pred else 0.0

    return dict(
        # T1 predictions
        f1s=f1s, f1l=f1l, f1t=f1t,
        # T2 recovered geometry
        R2s=R2s, t2l=t2l, Rm=Rm,
        # T3 mapped vx
        V3s=V3s, V3l=V3l, V3t=V3t,
        # T1 closures: how well does f_input match each topology prediction?
        cl1s=cl(f, f1s), cl1l=cl(f, f1l), cl1t=cl(f, f1t),
        # T2 closures: how well does R_input match recovered geometry?
        # cl2t = cl(R, Rm): R/Rm relational ratio. Always populated.
        # For sphere entries, R is R_major vs excavated R_minor — a relational ratio.
        # For torus entries where R_km IS R_minor, it is a direct closure check.
        cl2s=cl(R, R2s), cl2l=cl(R, t2l), cl2t=cl(R, Rm),
        # T3 closures: how well does vx_input match mapped vx?
        cl3s=cl(vx, V3s), cl3l=cl(vx, V3l), cl3t=cl(vx, V3t),
    )

# ── Mode detection ────────────────────────────────────────────────────────────
MODE_PREDICT = 'Predict: R+vx->f'
MODE_DETECT  = 'Detect:  f+vx->R'
MODE_MAP     = 'Map:     f+R->vx'
MODE_VERIFY  = 'Verify: all 3 supplied'

def resolve(R_km, f_Hz, vx_kms):
    """
    Given exactly 2 of 3, derive the third via substrate law.
    Returns (R_km, f_Hz, vx_kms, mode, derived_var).
    Raises ValueError if fewer than 2 are provided.
    If all 3 are provided, returns Verify mode (full verification run, no derivation).
    """
    have = sum(x is not None for x in [R_km, f_Hz, vx_kms])
    if have == 3:
        return R_km, f_Hz, vx_kms, MODE_VERIFY, None
    if have < 2:
        raise ValueError(f'TVP v1.9 requires at least 2 of 3 focal variables. Got {have}.')

    if R_km is None:
        # Detect mode: f + vx -> R
        R_km = vx_kms / (2*math.pi*f_Hz)
        return R_km, f_Hz, vx_kms, MODE_DETECT, 'R_km'
    if f_Hz is None:
        # Predict mode: R + vx -> f (sphere topology gives the primary prediction)
        f_Hz = vx_kms / (2*math.pi*R_km)
        return R_km, f_Hz, vx_kms, MODE_PREDICT, 'f_Hz'
    if vx_kms is None:
        # Map mode: f + R -> vx
        vx_kms = 2*math.pi*R_km*f_Hz
        return R_km, f_Hz, vx_kms, MODE_MAP, 'vx_kms'

    raise ValueError('Unreachable resolve() state.')

# ── Ledger header ─────────────────────────────────────────────────────────────
HEADER = [
    'Object', 'Mode', 'Derived_Var',
    'R_km', 'f_Hz', 'vx_kms',
    'f_pred_Sphere_T1_Hz', 'f_pred_Slab_T1_Hz', 'f_pred_Torus_T1_Hz',
    'R_rec_Sphere_T2_km',  't_rec_Slab_T2_km',  'Rminor_Torus_T2_km',
    'Vx_map_Sphere_T3_kms','Vx_map_Slab_T3_kms','Vx_map_Torus_T3_kms',
    'Cl_T1_Sphere_%', 'Cl_T1_Slab_%', 'Cl_T1_Torus_%',
    'Cl_T2_Sphere_%', 'Cl_T2_Slab_%', 'R_over_Rminor_T2_Torus_%',
    'Cl_T3_Sphere_%', 'Cl_T3_Slab_%', 'Cl_T3_Torus_%',
    'T1_Cl_Summary', 'P3GG_Active',
    'Ratio_f_sphere_over_f_obs', 'Ratio_factorized',
    'Fault_Code',
    'Note',
]

def run_object(name, R_km, f_Hz, vx_kms, note='', fault_hint=''):
    """
    Run TVP v1.9 nine-calculation audit on one object. Returns ledger row.

    fault_hint : auditor-declared fault code string. Annotates the output row;
                 never nulls any cell. All nine cells always populated,
                 including cl2t (R/Rm relational ratio).
                 'C1' — frequency is a transition energy, not a boundary closure frequency.
                 'G1' — derived theoretical scale used as observed dimension.
                 '00' or blank — no declared fault.
                 Composite codes via '+': e.g. 'C1+G1'.
    """
    R, f, vx, mode, derived = resolve(R_km, f_Hz, vx_kms)
    c = dict(nine(R, f, vx))   # mutable copy

    # ── Fault code ────────────────────────────────────────────────────────────
    fault_code = str(fault_hint).strip().upper() if fault_hint else '00'

    # ── T1 closure summary: all three topologies reported ───────────────────────
    # All nine run unconditionally. The substrate has no regimes and no exclusions.
    # Analyst reads sphere, slab, and torus T1 closures directly from this field.
    # cl1t ≡ 100% in Map mode (algebraic identity — that IS the forensic data).
    cl1s = c['cl1s'] if c['cl1s'] is not None else 0.0
    cl1l = c['cl1l'] if c['cl1l'] is not None else 0.0
    cl1t = c['cl1t'] if c['cl1t'] is not None else 0.0
    t1_summary = f'Sph:{cl1s:.2f}% | Slab:{cl1l:.2f}% | Tor:{cl1t:.2f}%'

    # ── Ratio ─────────────────────────────────────────────────────────────────
    # Map mode: vx was derived from R and f, so f1s = vx/(2πR) = f always.
    # f1s/f = 1 is a trivial identity — forensically useless.
    # Use VX_VAC_REF/(2πRf) instead: the vacuum-reference sphere eigenfrequency
    # over f_obs. This is the P3GG fingerprint for every Map mode entry.
    # Non-Map modes: vx is independently supplied, so f1s/f_obs is meaningful.
    if mode == MODE_MAP:
        ratio = VX_VAC_REF / (2*math.pi*R*f) if (R and f) else None
    else:
        ratio = c['f1s']/f if (f and c['f1s'] is not None) else None
    ratio_fac = factorize(ratio) if ratio else 'n/a'
    h = p3gg_active()

    return [
        name, mode, derived or 'none',
        fmt(R), fmt(f), fmt(vx),
        fmt(c['f1s']), fmt(c['f1l']), fmt(c['f1t']),
        fmt(c['R2s']), fmt(c['t2l']), fmt(c['Rm']),
        fmt(c['V3s']), fmt(c['V3l']), fmt(c['V3t']),
        fmt(c['cl1s']), fmt(c['cl1l']), fmt(c['cl1t']),
        fmt(c['cl2s']), fmt(c['cl2l']), fmt(c['cl2t']),
        fmt(c['cl3s']), fmt(c['cl3l']), fmt(c['cl3t']),
        t1_summary, h,
        fmt(ratio) if ratio else 'n/a', ratio_fac,
        fault_code,
        note,
    ]

# ── CLI ───────────────────────────────────────────────────────────────────────
def build_parser():
    p = argparse.ArgumentParser(
        description='TVP v1.9 Universal Forensic Auditor | Mass Harmonics | UMTTS')
    p.add_argument('--input',      default='', help='Batch CSV input path.')
    p.add_argument('--output',     default='', help='Output CSV path.')
    p.add_argument('--stdout',     action='store_true', help='Echo CSV to stdout.')
    p.add_argument('--mode',       default='auto',
                   choices=['auto','map','predict','detect'],
                   help=('auto: use whatever 2 of 3 the CSV/CLI supplies (default). '
                         'map: force R+f->vx (ignores vx column; vx derived). '
                         'predict: force R+VX_VAC_REF->f (ignores f column). '
                         'detect: force f+VX_VAC_REF->R (ignores R column).'))
    p.add_argument('--name',       default='Object')
    p.add_argument('--note',       default='')
    p.add_argument('--fault_hint', default='')
    p.add_argument('--R_km',       type=float, default=None)
    p.add_argument('--f_Hz',       type=float, default=None)
    p.add_argument('--vx_kms',     type=float, default=None)
    p.add_argument('--R_pm',       type=float, default=None)
    p.add_argument('--R_fm',       type=float, default=None)
    p.add_argument('--R_m',        type=float, default=None)
    p.add_argument('--f_MHz',      type=float, default=None)
    p.add_argument('--f_GHz',      type=float, default=None)
    p.add_argument('--f_eV',       type=float, default=None)
    p.add_argument('--f_MeV',      type=float, default=None)
    return p

def apply_unit_helpers(args):
    R  = args.R_km
    f  = args.f_Hz
    vx = args.vx_kms
    if args.R_pm  is not None: R = args.R_pm  * PM_KM
    if args.R_fm  is not None: R = args.R_fm  * FM_KM
    if args.R_m   is not None: R = args.R_m   * M_KM
    if args.f_MHz is not None: f = args.f_MHz * 1e6
    if args.f_GHz is not None: f = args.f_GHz * 1e9
    if args.f_eV  is not None: f = args.f_eV  * EV_HZ
    if args.f_MeV is not None: f = args.f_MeV * MEV_HZ
    return R, f, vx

def read_input_csv(path, mode='auto'):
    """
    Read batch CSV and apply mode override.
    mode='auto'    : use whatever 2-of-3 the CSV supplies (normal).
    mode='map'     : force vx=None so engine derives vx from R+f.
    mode='predict' : drop f, inject vx=VX_VAC_REF → engine derives f from R+vx.
    mode='detect'  : drop R, inject vx=VX_VAC_REF → engine derives R from f+vx.
    """
    rows = []
    with open(path, newline='', encoding='utf-8') as fp:
        for row in csv.DictReader(fp):
            def g(k):
                v = row.get(k, '')
                v = v.strip() if isinstance(v, str) else v
                return float(v) if v not in (None, '') else None

            R  = g('R_km')
            f  = g('f_Hz')
            vx = g('vx_kms')

            if g('R_pm')  is not None: R = g('R_pm')  * PM_KM
            if g('R_fm')  is not None: R = g('R_fm')  * FM_KM
            if g('R_m')   is not None: R = g('R_m')   * M_KM
            if g('f_MHz') is not None: f = g('f_MHz') * 1e6
            if g('f_GHz') is not None: f = g('f_GHz') * 1e9
            if g('f_eV')  is not None: f = g('f_eV')  * EV_HZ
            if g('f_MeV') is not None: f = g('f_MeV') * MEV_HZ

            # Apply mode override AFTER unit conversion
            if mode == 'map':
                vx = None                   # force engine to derive vx from R+f
            elif mode == 'predict':
                f  = None                   # drop f; engine derives f from R+VX_VAC_REF
                vx = VX_VAC_REF
            elif mode == 'detect':
                R  = None                   # drop R; engine derives R from f+VX_VAC_REF
                vx = VX_VAC_REF

            fault_hint = (row.get('fault_hint') or row.get('fault_code') or '').strip()
            rows.append((row.get('name', 'Object'), R, f, vx, row.get('note', ''), fault_hint))
    return rows

def write_output(rows, out_path, to_stdout=False):
    out_dir = os.path.dirname(out_path) if os.path.dirname(out_path) else '.'
    os.makedirs(out_dir, exist_ok=True)
    with open(out_path, 'w', newline='', encoding='utf-8') as fp:
        w = csv.writer(fp)
        w.writerow(HEADER)
        w.writerows(rows)
    print(f'TVP v1.9 -> {out_path}  ({len(rows)} records)')
    if to_stdout:
        buf = io.StringIO()
        w = csv.writer(buf)
        w.writerow(HEADER)
        w.writerows(rows)
        print(buf.getvalue())

def auto_out(tag):
    ts = datetime.now().strftime('%Y%m%d_%H%M%S')
    d  = './data'
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, f'tvp_{tag}_{ts}.csv')

def main():
    parser = build_parser()
    args   = parser.parse_args()

    if args.input:
        objects = read_input_csv(args.input, mode=args.mode)
        stem = os.path.splitext(os.path.basename(args.input))[0].replace('_input', '')
        tag  = f'{stem}_{args.mode}' if args.mode != 'auto' else stem
        ledger  = []
        for (name, R, f, vx, note, fh) in objects:
            try:
                ledger.append(run_object(name, R, f, vx, note, fh))
            except Exception as e:
                ledger.append([name, 'ERROR', str(e)] + [''] * (len(HEADER) - 3))
        write_output(ledger, args.output or auto_out(tag), args.stdout)
        return

    R, f, vx = apply_unit_helpers(args)
    if sum(x is not None for x in [R, f, vx]) < 2:
        parser.print_help()
        print('\nSupply at least 2 of 3 focal variables (R, f, vx).')
        sys.exit(1)

    row = run_object(args.name, R, f, vx, args.note, args.fault_hint)
    for h, v in zip(HEADER, row):
        print(f'  {h:<30s}: {v}')
    write_output([row], args.output or auto_out(args.name.replace(' ', '_')[:20]))

if __name__ == '__main__':
    main()
