#!/usr/bin/env python3
"""
asic_results.py - every headline number, computed from data/ and structures/.

Usage:
    python3 codes/asic_results.py            prints the result dictionary
    from asic_results import results         for the other scripts

Values are computed at call time; nothing is read from a derived table.

Exits non-zero on any fatal condition.
"""
import csv, itertools, math, os, re, sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import picode_paths as P

OXY = {'ASP': ('OD1', 'OD2'), 'GLU': ('OE1', 'OE2')}
POCKET_H = [(237, 351), (238, 347), (219, 409), (79, 418)]
POCKET_C = [(238, 350), (239, 346), (220, 408), (80, 417)]
CLUSTER_H = [237, 238, 347, 351]
CLUSTER_C = [238, 239, 346, 350]
STRUCTS = [('9E4A.pdb1', 8.5, 'human'), ('9E4F.pdb1', 5.7, 'human'),
           ('6VTL.pdb1', 8.0, 'chicken'), ('2QTS.pdb1', 6.5, 'chicken')]
NEUTRAL = {a: 0 for a in 'AVFMGLISTCPW'}
RNG_SEED = 20261007
MC_DRAWS = 20000         # permutation draws where exact enumeration is impossible


def die(msg):
    sys.exit('FATAL asic_results.py: %s' % msg)


def _rows(path, col):
    out = []
    for r in csv.DictReader(open(path)):
        if r['mutation'] == 'WT' or r.get(col) in (None, ''):
            continue
        out.append(r)
    if not out:
        die('no usable rows in %s for column %s' % (path, col))
    return out


def _spearman(x, y):
    def rank(v):
        o = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v); i = 0
        while i < len(o):
            j = i
            while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]:
                j += 1
            for k in range(i, j + 1):
                r[o[k]] = (i + j) / 2.0 + 1.0
            i = j + 1
        return r
    rx, ry = rank(x), rank(y)
    n = len(x); mx = sum(rx) / n; my = sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return num / (dx * dy) if dx and dy else float('nan')


def _perm_p(x, y, obs):
    n = len(x)
    if n > 9:
        die('exact permutation requested for n=%d' % n)
    hits = tot = 0
    for p in itertools.permutations(range(n)):
        tot += 1
        if abs(_spearman(x, [y[i] for i in p])) >= abs(obs) - 1e-12:
            hits += 1
    return hits / float(tot)


def _mc_p(x, y, obs, draws=MC_DRAWS):
    """Two-sided permutation p by Monte Carlo, for n too large to enumerate.

    Shuffling y against a fixed x is the exact null for 'no monotone relation',
    and the (hits+1)/(draws+1) form keeps p strictly positive so it can never be
    reported as zero.
    """
    if len(x) != len(y):
        die('_mc_p got unequal lengths %d and %d' % (len(x), len(y)))
    if len(x) < 6:
        die('_mc_p called with n=%d; use the exact test' % len(x))
    rng = np.random.default_rng(RNG_SEED)
    yv = np.asarray(y, dtype=float)
    hits = 0
    for _ in range(draws):
        if abs(_spearman(x, rng.permutation(yv).tolist())) >= abs(obs) - 1e-12:
            hits += 1
    return (hits + 1) / float(draws + 1)


def _partial(x, y, z):
    """Spearman partial rank correlation of x with y, controlling for z."""
    rxy, rxz, ryz = _spearman(x, y), _spearman(x, z), _spearman(y, z)
    den = math.sqrt(max(1e-12, (1 - rxz ** 2) * (1 - ryz ** 2)))
    return (rxy - rxz * ryz) / den


def _ols(y, X):
    XtX = X.T @ X
    if np.linalg.matrix_rank(XtX) < X.shape[1]:
        die('rank-deficient design matrix')
    Xi = np.linalg.inv(XtX)
    b = Xi @ (X.T @ y)
    r = y - X @ b
    S = (X * r[:, None]).T @ (X * r[:, None])
    cov = Xi @ S @ Xi
    return b, np.sqrt(np.diag(cov))


def _demean(y, X, g):
    y2, X2 = y.astype(float).copy(), X.astype(float).copy()
    keep = np.zeros(len(y), bool)
    for v in np.unique(g):
        m = g == v
        if m.sum() < 2:
            continue
        keep |= m
        y2[m] -= y2[m].mean()
        X2[m] -= X2[m].mean(axis=0)
    return y2[keep], X2[keep], keep


def _props():
    import openpyxl
    x = P.project('data', 'sources', '42003_2025_9090_MOESM4_ESM.xlsx')
    if not os.path.exists(x):
        die('missing %s' % x)
    wb = openpyxl.load_workbook(x, read_only=True, data_only=True)
    out, hdr = {}, None
    for sh in wb.sheetnames:
        if not re.match(r'^[A-Z]\d{1,3}_(act|des)$', sh):
            continue
        for i, row in enumerate(wb[sh].iter_rows(values_only=True)):
            if i == 0:
                hdr = [str(c) if c else '' for c in row]
                continue
            if not row or not isinstance(row[0], str) or len(row[0]) != 1:
                continue
            if not row[0].isupper():
                continue
            v = row[2:12]
            if all(isinstance(q, (int, float)) for q in v):
                out.setdefault(row[0], [float(q) for q in v])
    if not out or hdr is None or hdr[2] != 'size':
        die('could not harvest the side-chain property block')
    return out, hdr[2:12].index('size')


def _pdb(path):
    if not os.path.isfile(path):
        die('%s is not present. The deposited assemblies are not carried in '
            'this repository.\nRun: python3 codes/fetch_structures.py'
            % os.path.relpath(path, P.project()))
    at = {}
    for l in open(path, errors='replace'):
        if l[:4] != 'ATOM':
            continue
        at.setdefault((l[21], int(l[22:26])), [l[17:20].strip(), {}])[1][
            l[12:16].strip()] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    if not at:
        die('no ATOM records in %s' % path)
    return at


def _ox(at, k):
    e = at.get(k)
    return [e[1][o] for o in OXY.get(e[0], ()) if o in e[1]] if e else []


def _mind(A, B):
    return min(math.dist(a, b) for a in A for b in B) if A and B else float('nan')


def struct_descriptors(rows, struct='9E4F.pdb1'):
    """Per-position burial, distance to the pocket and largest measured shift.

    Computed on the low-pH human trimer 9E4F with the accessibility algorithm
    and maximum-ASA table the method folder uses. Exposed as a function so the
    values can be re-ranked without re-implementing the geometry.

    Returns rsa, dist, shift, lpo, nsub, mshift: relative side-chain solvent
    accessibility, distance to the pocket centre, the largest measured shift at
    the position, the burial-threshold module, the number of substitutions
    measured there, and the mean shift over them.
    """
    at = _pdb(P.project('structures', struct))
    pocket_pts = []
    for a_, b_ in POCKET_H:
        for rr in (a_, b_):
            o = _ox(at, ('A', rr))
            if o:
                pocket_pts.extend(o)
    if not pocket_pts:
        die('no pocket carboxyl oxygens found in %s chain A' % struct)
    centre = np.mean(np.array(pocket_pts), axis=0)

    lpo = P.load_lpo()
    from Bio.PDB import PDBParser
    from Bio.PDB.SASA import ShrakeRupley
    import io as _io
    txt = ''.join(l for l in open(P.project('structures', struct),
                                  errors='replace')
                  if l[:6] in ('ATOM  ', 'HETATM', 'TER   ') and l[17:20] != 'HOH')
    st = PDBParser(QUIET=True).get_structure('x', _io.StringIO(txt))
    mdl = st[0]
    ShrakeRupley().compute(mdl, level='R')
    sasa = {(c.id, r.id[1]): r.sasa for c in mdl for r in c if r.id[0] == ' '}

    dist, rsa, shift, nsub, mshift = [], [], [], [], []
    seenp = set()
    for r in rows:
        ri = int(r['resi'])
        if ri in seenp:
            continue
        seenp.add(ri)
        res = at.get(('A', ri))
        if not res:
            continue
        side = [v for k, v in res[1].items() if k not in ('N', 'CA', 'C', 'O')]
        if not side:
            continue
        d = float(np.linalg.norm(np.mean(np.array(side), axis=0) - centre))
        mx = lpo.MAXASA.get(r['wt_aa'])
        if mx is None or ('A', ri) not in sasa:
            continue
        here = [abs(float(x['d_act'])) for x in rows if int(x['resi']) == ri]
        dist.append(d)
        rsa.append(100.0 * sasa[('A', ri)] / mx)
        shift.append(max(here))
        # The maximum rises with the number of tries, so the count and the
        # mean are returned too and that confound can be tested.
        nsub.append(len(here))
        mshift.append(float(np.mean(here)))
    if len(dist) < 10:
        die('only %d positions could be placed in the structure' % len(dist))
    return rsa, dist, shift, lpo, nsub, mshift


def results():
    R = {}
    L = P.load_ladder()
    props, si = _props()
    core = set(L.CORE)
    impl = set(L.CORE) | set(L.ARG)
    allg = impl | set(L.EXT)

    # ---- cohorts
    pri = P.project('data', 'asic1a_ph50.csv')
    sec = P.project('data', 'asic1a_ph50_canessa.csv')
    for p in (pri, sec):
        if not os.path.exists(p):
            die('missing %s' % p)
    wt = next(r for r in csv.DictReader(open(pri)) if r['mutation'] == 'WT')
    rows = _rows(pri, 'd_act')
    R['n_subs'] = len(rows)
    R['n_pos'] = len({r['resi'] for r in rows})
    R['wt_act'] = float(wt['ph50_act']); R['wt_act_n'] = int(wt['n_act'])
    R['wt_ssd'] = float(wt['ph50_ssd']); R['wt_ssd_n'] = int(wt['n_ssd'])
    R['n_on_edge'] = sum(1 for r in rows if (r['wt_aa'], r['mut_aa']) in allg)

    crows = _rows(sec, 'd_act')
    cwt = next(r for r in csv.DictReader(open(sec)) if r['mutation'] == 'WT')
    R['c_n_subs'] = len(crows)
    R['c_n_pos'] = len({r['resi'] for r in crows})
    R['c_wt_act'] = float(cwt['ph50_act'])

    # ---- edge coverage
    per = defaultdict(list)
    for r in rows:
        per[(r['wt_aa'], r['mut_aa'])].append(abs(float(r['d_act'])))
    cper = defaultdict(list)
    for r in crows:
        cper[(r['wt_aa'], r['mut_aa'])].append(abs(float(r['d_act'])))
    R['core_covered'] = sum(1 for p in L.CORE if p in per)
    R['core_missing'] = ['%s>%s' % p for p in L.CORE if p not in per]
    R['core_covered_both'] = sum(1 for p in L.CORE if p in per or p in cper)
    R['filled_by_second'] = ['%s>%s' % p for p in L.CORE
                             if p not in per and p in cper]

    # ---- edge costs
    ec = P.project('results', 'edge_costs.csv')
    if not os.path.exists(ec):
        die('missing %s; run codes/extract_edge_costs.py first' % ec)
    cliff, hums = {}, {}
    for r in csv.DictReader(open(ec)):
        cliff[(r['source'], r['target'])] = float(r['cliff_delta'])
        hums[(r['source'], r['target'])] = float(r['humsavar_pct'])

    edges = [p for p in L.CORE if p in per]
    shift = [float(np.mean(per[p])) for p in edges]
    # Substitutions behind each of the eight edge means.
    _ns = sorted(len(per[p]) for p in edges)
    R['edge_n_min'] = _ns[0]
    R['edge_n_max'] = _ns[-1]
    R['edge_singletons'] = sum(1 for n in _ns if n == 1)
    R['edge_singleton_names'] = ['%s>%s' % p for p in edges if len(per[p]) == 1]
    cl = [cliff[p] for p in edges]
    hu = [hums[p] for p in edges]
    dq = [abs(L.CHARGE[p[1]] - L.CHARGE[p[0]]) for p in edges]
    R['rho_cliff'] = _spearman(shift, cl)
    R['p_cliff'] = _perm_p(shift, cl, R['rho_cliff'])
    R['rho_hums'] = _spearman(shift, hu)
    R['p_hums'] = _perm_p(shift, hu, R['rho_hums'])
    R['rho_charge'] = _spearman(shift, dq)
    R['p_charge'] = _perm_p(shift, dq, R['rho_charge'])
    R['n_charge_levels'] = len(set(dq))

    # ---- fixed-effects model
    # Indicator: the twelve core edges. The wider sets follow as a sensitivity.
    def _fe_fit(eset, colname):
        rr = _rows(pri, colname)
        y = np.array([abs(float(r[colname])) for r in rr])
        g = np.array([int(r['resi']) for r in rr])
        X = np.column_stack([
            [abs(L.CHARGE.get(r['mut_aa'], 0) - L.CHARGE.get(r['wt_aa'], 0)) for r in rr],
            [abs(props[r['mut_aa']][si] - props[r['wt_aa']][si]) for r in rr],
            [1.0 if (r['wt_aa'], r['mut_aa']) in eset else 0.0 for r in rr]])
        yc, Xc, keep = _demean(y, X, g)
        b, se = _ols(yc, Xc)
        rng = np.random.default_rng(RNG_SEED)
        lad0 = X[keep, 2].copy()
        pg = g[keep]
        hits = 0
        NP = 20000
        for _ in range(NP):
            sh = lad0.copy()
            for v in np.unique(pg):
                m = pg == v
                sh[m] = rng.permutation(sh[m])
            c = sh.astype(float).copy()
            for v in np.unique(pg):
                m = pg == v
                c[m] -= c[m].mean()
            if np.allclose(c, 0):
                continue
            Xp = Xc.copy(); Xp[:, 2] = c
            bp, _ = _ols(yc, Xp)
            if abs(bp[2]) >= abs(b[2]) - 1e-12:
                hits += 1
        return b, se, keep, g, (hits + 1) / float(NP + 1)

    for tag, colname in (('act', 'd_act'), ('ssd', 'd_ssd')):
        b, se, keep, g, pperm = _fe_fit(core, colname)
        R['fe_%s_n' % tag] = int(keep.sum())
        R['fe_%s_npos' % tag] = len(set(g[keep]))
        R['fe_%s_charge' % tag] = b[0]
        R['fe_%s_charge_se' % tag] = se[0]
        R['fe_%s_charge_t' % tag] = b[0] / se[0]
        R['fe_%s_volume' % tag] = b[1]
        R['fe_%s_volume_se' % tag] = se[1]
        R['fe_%s_ladder' % tag] = b[2]
        R['fe_%s_ladder_se' % tag] = se[2]
        R['fe_%s_ladder_p' % tag] = pperm

    # Sensitivity to the edge set, activation readout.
    for suf, eset in (('_impl', impl), ('_ext', allg)):
        b, se, _k, _g, pperm = _fe_fit(eset, 'd_act')
        R['fe_act_ladder%s' % suf] = b[2]
        R['fe_act_ladder_se%s' % suf] = se[2]
        R['fe_act_ladder_p%s' % suf] = pperm
    R['fe_n_edges_core'] = len(core)
    R['fe_n_edges_impl'] = len(impl)
    R['fe_n_edges_ext'] = len(allg)

    # ---- within-position charge-matched pairs
    # Ladder against non-ladder, same position, same charge change.
    CH = dict(L.CHARGE); CH.update(NEUTRAL)
    byp = defaultdict(list)
    for r in rows:
        byp[(r['wt_aa'], int(r['resi']))].append(r)

    def _pair_diffs(eset):
        out = []
        for (w, pz), v in byp.items():
            for r in v:
                if (w, r['mut_aa']) not in eset:
                    continue
                q = CH[r['mut_aa']] - CH[w]
                for s in v:
                    if s is r or (w, s['mut_aa']) in eset:
                        continue
                    if abs((CH[s['mut_aa']] - CH[w]) - q) > 1e-9:
                        continue
                    out.append(abs(float(r['d_act'])) - abs(float(s['d_act'])))
        return np.array(out)

    from scipy import stats

    def _pair_stats(d, suf):
        n = len(d)
        sm = int((d < 0).sum())
        k = min(sm, n - sm)
        R['pairs_n' + suf] = n
        R['pairs_ladder_smaller' + suf] = sm
        R['pairs_sign_p' + suf] = min(1.0, 2 * sum(math.comb(n, i)
                                      for i in range(k + 1)) / 2.0 ** n)
        R['pairs_wilcoxon_p' + suf] = float(stats.wilcoxon(d)[1])
        if suf:
            return
        R['pairs_mean_diff'] = float(d.mean())
        R['pairs_win_by'] = float(np.mean([-x for x in d if x < 0]))
        R['pairs_lose_by'] = float(np.mean([x for x in d if x > 0]))

    dcore = _pair_diffs(core)
    if not len(dcore):
        die('no charge-matched pairs built on the twelve core edges')
    _pair_stats(dcore, '')
    for suf, eset in (('_impl', impl), ('_ext', allg)):
        dd = _pair_diffs(eset)
        if not len(dd):
            die('no charge-matched pairs built on the %s edge set' % suf.strip('_'))
        _pair_stats(dd, suf)
    # What the proposed extensions add to the pair set.
    R['pairs_ext_added'] = R['pairs_n_ext'] - R['pairs_n']
    R['pairs_ext_added_wins'] = (R['pairs_ladder_smaller_ext']
                                 - R['pairs_ladder_smaller'])
    if R['pairs_ext_added'] < 0 or R['pairs_ext_added_wins'] < 0:
        die('widening the edge set removed pairs, which cannot happen; the edge '
            'sets are not nested as assumed')

    # ---- agreement between the two cost scales, twelve core edges.
    # The scales run opposite ways, so agreement is a negative correlation.
    ec = {}
    with open(P.project('results', 'edge_costs.csv')) as fh:
        for row in csv.DictReader(fh):
            ec[(row['source'], row['target'])] = row
    miss = [e for e in L.CORE if e not in ec]
    if miss:
        die('edge_costs.csv is missing %d core edges: %s'
            % (len(miss), ', '.join('%s>%s' % e for e in miss)))
    cl = [float(ec[e]['cliff_delta']) for e in L.CORE]
    pc = [float(ec[e]['humsavar_pct']) for e in L.CORE]
    R['rho_cost_agree'] = float(stats.spearmanr(cl, pc).statistic)
    R['n_cost_agree'] = len(cl)

    # ---- buried uncompensated carboxylates, monomer frame against assembly
    ic = P.project('results', 'interface_carboxylates.csv')
    if not os.path.exists(ic):
        die('missing %s; run codes/extract_interface_carboxylates.py first'
            % os.path.relpath(ic, P.project()))
    with open(ic) as fh:
        icr = list(csv.DictReader(fh))
    want = {f.split('.')[0] for f, _ph, _sp in STRUCTS}
    have = {r['structure'] for r in icr}
    if want - have:
        die('interface_carboxylates.csv is missing %s'
            % ', '.join(sorted(want - have)))
    add = [int(r['added_tier0']) for r in icr if r['structure'] in want]
    R['iface_added_min'] = min(add)
    R['iface_added_max'] = max(add)
    R['iface_n_struct'] = len(add)
    R['iface_added'] = ['%s %+d' % (r['structure'], int(r['added_tier0']))
                        for r in icr if r['structure'] in want]
    R['iface_mono_tier0'] = [int(r['monomer_tier0']) for r in icr
                             if r['structure'] in want]

    # ---- geometry
    geo, palm = {}, {}
    for f, ph, sp in STRUCTS:
        at = _pdb(P.project('structures', f))
        chs = sorted({c for c, _ in at})
        prs = POCKET_H if sp == 'human' else POCKET_C
        clu = CLUSTER_H if sp == 'human' else CLUSTER_C
        for i, (r1, r2) in enumerate(prs):
            ia, ie = [], []
            for ca in chs:
                for cb in chs:
                    dd = _mind(_ox(at, (ca, r1)), _ox(at, (cb, r2)))
                    if dd == dd:
                        (ia if ca == cb else ie).append(dd)
            if ia and ie:
                geo[(f, i)] = (min(ia), min(ie))
        for rr in (prs[2][0], prs[2][1]):
            s_, a_ = [], []
            for ca in chs:
                A = _ox(at, (ca, rr))
                for cb in chs:
                    for q in clu:
                        dd = _mind(A, _ox(at, (cb, q)))
                        if dd == dd:
                            (s_ if ca == cb else a_).append(dd)
            if s_ and a_:
                palm[(f, rr)] = (min(s_), min(a_))
    if not geo or not palm:
        die('geometry could not be computed')
    R['geo_all_intra'] = all(i < e for i, e in geo.values())
    R['geo_min_margin'] = min(e - i for i, e in geo.values())
    R['geo_n_comparisons'] = len(geo)
    R['palm_same_min'] = min(s for s, _ in palm.values())
    R['palm_same_max'] = max(s for s, _ in palm.values())
    R['palm_adj_min'] = min(a for _, a in palm.values())
    R['palm_adj_max'] = max(a for _, a in palm.values())
    R['collapse_lo'] = geo[('9E4F.pdb1', 0)][0]
    R['collapse_hi'] = geo[('9E4A.pdb1', 0)][0]

    # ---- signed direction of the measured shift
    sa = [float(r['d_act']) for r in rows]
    R['signed_down'] = sum(1 for v in sa if v < 0)
    R['signed_up'] = sum(1 for v in sa if v > 0)
    R['signed_mean'] = float(np.mean(sa))

    # ---- coupling between the two readouts
    both = [(float(r['d_act']), float(r['d_ssd'])) for r in rows
            if r['d_ssd'] not in ('', None)]
    R['coupling_n'] = len(both)
    R['coupling_rho'] = _spearman([x for x, _ in both], [y for _, y in both])
    R['coupling_p'] = _mc_p([x for x, _ in both], [y for _, y in both],
                            R['coupling_rho'])
    # Sign test: independent readouts would disagree in sign half the time.
    R['coupling_same_sign'] = sum(1 for a, b in both if a * b > 0)
    R['coupling_sign_n'] = sum(1 for a, b in both if a * b != 0)
    k = min(R['coupling_same_sign'], R['coupling_sign_n'] - R['coupling_same_sign'])
    R['coupling_sign_p'] = min(1.0, 2 * sum(math.comb(R['coupling_sign_n'], i)
                                            for i in range(k + 1))
                               / 2.0 ** R['coupling_sign_n'])

    # ---- per-position structural context, from the low-pH human trimer
    rsa, dist, shift, lpo, nsub, mshift = struct_descriptors(rows)
    R['struct_n'] = len(dist)
    R['rho_dist'] = _spearman(dist, shift)
    R['p_dist'] = _mc_p(dist, shift, R['rho_dist'])
    R['rho_rsa'] = _spearman(rsa, shift)
    R['p_rsa'] = _mc_p(rsa, shift, R['rho_rsa'])
    R['dist_min'] = float(min(dist)); R['dist_max'] = float(max(dist))
    R['rsa_min'] = float(min(rsa)); R['rsa_max'] = float(max(rsa))
    # Burial and pocket distance are correlated, so the partial rank
    # correlation separates them.
    R['rho_rsa_dist'] = _spearman(rsa, dist)
    R['rho_rsa_partial'] = _partial(rsa, shift, dist)

    # The maximum rises with the count, so both controls are reported.
    R['nsub_min'] = int(min(nsub))
    R['nsub_max'] = int(max(nsub))
    R['nsub_median'] = float(np.median(nsub))
    R['rho_n_shift'] = _spearman(nsub, shift)
    R['rho_n_rsa'] = _spearman(nsub, rsa)
    R['rho_n_mshift'] = _spearman(nsub, mshift)
    R['rho_rsa_partial_n'] = _partial(rsa, shift, nsub)
    R['rho_rsa_mean'] = _spearman(rsa, mshift)
    R['p_rsa_mean'] = _mc_p(rsa, mshift, R['rho_rsa_mean'])
    if R['rho_rsa_mean'] * R['rho_rsa'] <= 0:
        die('burial correlates with the per-position maximum at %+.3f and with '
            'the mean at %+.3f. The two estimators disagree in sign, so the '
            'result depends on which one is chosen and cannot be reported as a '
            'burial effect.' % (R['rho_rsa'], R['rho_rsa_mean']))

    # The resting structure of the same construct gives a different answer,
    # so the correlation belongs to one conformer.
    rsa_r, _d_r, shift_r = struct_descriptors(rows, '9E4A.pdb1')[:3]
    if len(rsa_r) != len(rsa):
        die('the two conformers place a different number of positions (%d on '
            '9E4F, %d on 9E4A), so the burial correlations are not comparable.'
            % (len(rsa), len(rsa_r)))
    R['rho_rsa_resting'] = _spearman(rsa_r, shift_r)
    R['p_rsa_resting'] = _mc_p(rsa_r, shift_r, R['rho_rsa_resting'])

    # Published so the manuscript can state the cutoffs.
    R['rsa_buried_t'] = float(lpo.BURIED_RSA)
    R['rsa_partial_t'] = float(lpo.PARTIAL_RSA)
    R['salt_bridge_t'] = float(lpo.SALT_BRIDGE)

    R['buried_n'] = sum(1 for v in rsa if v < lpo.BURIED_RSA)
    R['exposed_n'] = sum(1 for v in rsa if v >= lpo.PARTIAL_RSA)
    bu = [s for v, s in zip(rsa, shift) if v < lpo.BURIED_RSA]
    ex = [s for v, s in zip(rsa, shift) if v >= lpo.PARTIAL_RSA]
    if not bu or not ex:
        die('burial stratification is empty at the %g and %g percent relative '
            'accessibility thresholds' % (lpo.BURIED_RSA, lpo.PARTIAL_RSA))
    R['buried_mean'] = float(np.mean(bu))
    R['exposed_mean'] = float(np.mean(ex))
    R['burial_mwu_p'] = float(stats.mannwhitneyu(bu, ex,
                                                 alternative='two-sided')[1])

    # ---- predicted pKa: state dependence, and whether it ranks the measurement
    # Parsed from the PROPKA outputs with asic_second_shell.py's own parser.
    import importlib.util as _il
    _sp = _il.spec_from_file_location(
        'asic_second_shell', os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                          'asic_second_shell.py'))
    _ss = _il.module_from_spec(_sp)
    try:
        _sp.loader.exec_module(_ss)
    except Exception as e:
        die('cannot load asic_second_shell.py: %s' % e)
    pk_lo, _ = _ss.parse_pka(P.project('structures', '9E4F.pka'))
    pk_hi, _ = _ss.parse_pka(P.project('structures', '9E4A.pka'))
    CH = 'A'
    moves = []
    for resi in sorted(_ss.POCKET):
        a_, b_ = (resi, CH), (resi, CH)
        if a_ not in pk_lo or b_ not in pk_hi:
            die('pocket position %d absent from one of the PROPKA outputs' % resi)
        moves.append(abs(pk_lo[a_][1] - pk_hi[b_][1]))
    R['pka_n_pocket'] = len(moves)
    R['pka_max_move'] = float(max(moves))
    R['pka_moved_over_2'] = sum(1 for m in moves if m > 2.0)
    # Calibration: does |pKa - model compound| rank the measured sensitivity?
    byres = defaultdict(list)
    for r in rows:
        byres[int(r['resi'])].append(r)
    xs, ya = [], []
    for resi, rs in sorted(byres.items()):
        k = (resi, CH)
        if k not in pk_lo:
            continue
        rn, v, _b = pk_lo[k]
        if rn not in _ss.MODEL:
            continue                       # calibrate on carboxylates only
        xs.append(abs(v - _ss.MODEL[rn]))
        ya.append(max(abs(float(x['d_act'])) for x in rs))
    if len(xs) < 10:
        die('only %d carboxylate positions overlap the PROPKA output' % len(xs))
    R['cal_n'] = len(xs)
    R['cal_rho'] = _spearman(xs, ya)
    R['cal_p'] = _mc_p(xs, ya, R['cal_rho'])

    # ---- the interval the activation null actually excludes
    R['fe_act_ladder_ci'] = 1.96 * R['fe_act_ladder_se']
    R['fe_act_ci_ratio'] = R['fe_act_charge'] / R['fe_act_ladder_ci']

    # ---- domain stratification of the charge effect
    from collections import defaultdict as _dd
    byloc = _dd(list)
    for r in rows:
        byloc[r.get('location', '')].append(abs(float(r['d_act'])))
    if '' in byloc:
        die('%d substitutions carry no domain label' % len(byloc['']))
    order = sorted(byloc, key=lambda k: -float(np.mean(byloc[k])))
    R['dom_counts'] = {k: len(byloc[k]) for k in order}
    R['dom_means'] = {k: round(float(np.mean(byloc[k])), 3) for k in order}
    R['dom_rank'] = ['%s %.3f' % (k, float(np.mean(byloc[k]))) for k in order]
    R['dom_top'] = order[0]
    R['dom_bottom'] = order[-1]
    R['dom_top_mean'] = float(np.mean(byloc[order[0]]))
    R['dom_bottom_mean'] = float(np.mean(byloc[order[-1]]))
    R['dom_top_n'] = len(byloc[order[0]])
    R['dom_bottom_n'] = len(byloc[order[-1]])
    big = sorted(byloc, key=lambda k: -len(byloc[k]))[:2]
    R['dom_largest'] = ['%s %d' % (k, len(byloc[k])) for k in big]
    R['dom_ratio'] = R['dom_top_mean'] / R['dom_bottom_mean']
    R['dom_kw_p'] = float(stats.kruskal(*[byloc[k] for k in order])[1])
    R['dom_n'] = len(order)

    # ---- the asparagine to histidine edge, quoted on its own
    # Unmeasured in both cohorts, scored well by the predictor.
    NH, DH = ('N', 'H'), ('D', 'H')
    for e in (NH, DH):
        if e not in cliff:
            die('edge %s>%s is missing from edge_costs.csv' % e)
    if NH in per or NH in cper:
        die('N>H now carries a measured pH50, so the unmeasured-edge passage '
            'no longer holds and must be rewritten')
    R['nh_cliff'] = cliff[NH]
    R['nh_hums'] = hums[NH]
    R['dh_cliff'] = cliff[DH]
    R['dh_hums'] = hums[DH]
    R['nh_dq'] = abs(L.CHARGE['H'] - L.CHARGE['N'])
    R['nh_pred'] = R['nh_dq'] * R['fe_act_charge']
    R['nh_better_count'] = sum(1 for e in L.CORE
                               if e != NH and cliff[e] < cliff[NH])
    meas = {}
    for e in L.CORE:
        v = per.get(e) or cper.get(e)
        if v:
            meas[e] = float(np.mean(v))
    if len(meas) != R['core_covered_both']:
        die('%d core edges carry a mean shift but coverage says %d'
            % (len(meas), R['core_covered_both']))
    R['min_edge_shift'] = min(meas.values())
    R['min_edge_name'] = '%s>%s' % min(meas, key=meas.get)
    return R


def main():
    R = results()
    w = max(len(k) for k in R)
    for k in sorted(R):
        v = R[k]
        if isinstance(v, float):
            print('%-*s  %.4f' % (w, k, v))
        else:
            print('%-*s  %s' % (w, k, v))


if __name__ == '__main__':
    main()
