#!/usr/bin/env python3
"""
asic_second_shell.py - second-shell positions for the ASIC1a acidic pocket.

Usage:
    python3 asic_second_shell.py --pka <low_pH.pka> --pka-resting <resting.pka> \
        --ph50 asic1a_ph50.csv

The second shell is taken from the PROPKA decomposition: for each titratable
group, the residues it names as contributing to the predicted pKa through
side-chain hydrogen bonding, backbone hydrogen bonding and Coulombic
interaction. Candidates are then filtered on pI-code reachability using the
ladder definition.

Calibration compares the predicted pKa perturbation at every titratable
position against the measured activation and SSD shift at that position, over
the positions the measured set covers.

Exits non-zero on any fatal condition.
"""
import argparse, csv, itertools, math, os, re, sys
from collections import defaultdict

POCKET = {79: 'E', 219: 'E', 237: 'D', 238: 'E', 242: 'E',
          347: 'D', 351: 'D', 409: 'D', 418: 'E'}
TITRATABLE = {'ASP': 'D', 'GLU': 'E', 'HIS': 'H', 'LYS': 'K', 'ARG': 'R',
              'TYR': 'Y', 'CYS': 'C'}
MODEL = {'ASP': 3.80, 'GLU': 4.50}

DET = re.compile(
    r'^(?P<rn>[A-Z]{3})\s+(?P<ri>\d+)\s+(?P<ch>[A-Z])\s+'
    r'(?:(?P<pka>-?\d+\.\d+)\*?\s+(?P<bur>\d+)\s*%\s+\S+\s+\S+\s+\S+\s+\S+\s+)?'
    r'(?P<rest>.*)$')
PARTNER = re.compile(r'(-?\d+\.\d+)\s+([A-Z]{3})\s+(\d+)\s+([A-Z])')


def die(msg):
    sys.exit('FATAL asic_second_shell.py: %s' % msg)


def load_ladder(path=None):
    """The ladder edges and charges, by way of codes/picode_paths.py.

    An explicit --ladder path wins. Otherwise the pi-code method folder is used
    where it is present and the snapshot committed in data/ where it is not, so
    the definitions keep a single home and a clone of this repository without
    the method folder still resolves the same published constants.
    """
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import picode_paths as P
    return P.load_ladder(path=path)


def parse_pka(path):
    """Return (pka, determinants). pka[(resi,chain)] = (resname, pka, buried%).
    determinants[(resi,chain)] = list of (value, partner_resname, partner_resi,
    partner_chain)."""
    try:
        lines = open(path).read().splitlines()
    except OSError as e:
        die('cannot read %s: %s' % (path, e))
    try:
        stop = next(i for i, l in enumerate(lines)
                    if l.startswith('SUMMARY OF THIS PREDICTION'))
    except StopIteration:
        die('%s has no SUMMARY block; not a PROPKA .pka file?' % path)

    pka, det = {}, defaultdict(list)
    for l in lines[:stop]:
        m = DET.match(l.rstrip())
        if not m:
            continue
        key = (int(m.group('ri')), m.group('ch'))
        if m.group('pka') is not None:
            pka[key] = (m.group('rn'), float(m.group('pka')),
                        int(m.group('bur')))
        for val, prn, pri, pch in PARTNER.findall(m.group('rest') or ''):
            if prn == 'XXX':
                continue
            if (int(pri), pch) == key:
                continue                      # self term, not a second-shell partner
            det[key].append((float(val), prn, int(pri), pch))
    if not pka:
        die('parsed no pKa values from %s' % path)
    return pka, det


def spearman(x, y):
    n = len(x)
    if n < 4:
        die('calibration needs at least 4 paired positions, got %d' % n)

    def rank(v):
        o = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(o):
            j = i
            while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                r[o[k]] = avg
            i = j + 1
        return r
    rx, ry = rank(x), rank(y)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return num / (dx * dy) if dx and dy else float('nan')


def perm_p(x, y, obs, trials=200000):
    """Two-sided permutation p. Random for n>9, exact below."""
    n = len(x)
    if n <= 8:
        hits = tot = 0
        for p in itertools.permutations(range(n)):
            tot += 1
            if abs(spearman(x, [y[i] for i in p])) >= abs(obs) - 1e-12:
                hits += 1
        return hits / float(tot)
    import random
    rnd = random.Random(20261007)
    hits = 0
    idx = list(range(n))
    for _ in range(trials):
        rnd.shuffle(idx)
        if abs(spearman(x, [y[i] for i in idx])) >= abs(obs) - 1e-12:
            hits += 1
    return hits / float(trials)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pka', required=True, help='low-pH / collapsed state')
    ap.add_argument('--pka-resting', required=True)
    ap.add_argument('--ph50', required=True)
    ap.add_argument('--ladder', help='path to ladder_edge_cost.py; by default the method folder, then data/ladder_definition.json')
    ap.add_argument('--chain', default='A')
    a = ap.parse_args()

    L = load_ladder(a.ladder)
    edges = set(L.CORE) | set(L.ARG) | set(L.EXT)
    reach = defaultdict(list)
    for s, t in edges:
        reach[s].append(t)

    pk_lo, det_lo = parse_pka(a.pka)
    pk_hi, _ = parse_pka(a.pka_resting)

    rows = [r for r in csv.DictReader(open(a.ph50)) if r['mutation'] != 'WT']
    if not rows:
        die('no substitution rows in %s' % a.ph50)
    meas = defaultdict(list)
    for r in rows:
        meas[int(r['resi'])].append(r)

    # ---- 1. the pocket itself, in both states
    print('--- 1. acidic pocket carboxylates, predicted pKa, chain %s ---' % a.chain)
    print('%-7s %-9s %-9s %-9s %-7s %s'
          % ('residue', 'pKa low', 'pKa rest', 'delta', 'buried', 'shift vs model (low)'))
    missing = []
    for resi in sorted(POCKET):
        k = (resi, a.chain)
        if k not in pk_lo or k not in pk_hi:
            missing.append(resi)
            continue
        rn, lo, bur = pk_lo[k]
        _, hi, _ = pk_hi[k]
        mod = MODEL.get(rn)
        if mod is None:
            die('no model pKa for %s' % rn)
        print('%-7s %-9.2f %-9.2f %+-9.2f %-7s %+.2f'
              % ('%s%d' % (POCKET[resi], resi), lo, hi, lo - hi, '%d%%' % bur, lo - mod))
    if missing:
        die('pocket positions %s absent from the PROPKA output; wrong chain or a '
            'different construct' % missing)

    # ---- 2. calibration against the measured set
    print('\n--- 2. calibration: does the predicted pKa rank measured pH sensitivity? ---')
    for state, pk, tag in (('low-pH', pk_lo, a.pka.split('/')[-1]),
                           ('resting', pk_hi, a.pka_resting.split('/')[-1])):
        xs, ya, ys, used = [], [], [], []
        for resi, rs in sorted(meas.items()):
            k = (resi, a.chain)
            if k not in pk:
                continue
            rn, v, _ = pk[k]
            if rn not in MODEL:
                continue              # calibrate on carboxylates only
            xs.append(abs(v - MODEL[rn]))
            ya.append(max(abs(float(r['d_act'])) for r in rs))
            ys.append(max(abs(float(r['d_ssd'])) for r in rs))
            used.append(resi)
        if len(xs) < 4:
            die('only %d carboxylate positions overlap between %s and the measured '
                'set' % (len(xs), tag))
        ra, rs_ = spearman(xs, ya), spearman(xs, ys)
        print('  %s state (%s), %d carboxylate positions with both a predicted pKa '
              'and a measured shift' % (state, tag, len(xs)))
        print('    |pKa - model| vs max |delta activation pH50|  rho = %+.3f  p = %.3f'
              % (ra, perm_p(xs, ya, ra)))
        print('    |pKa - model| vs max |delta SSD pHD50|        rho = %+.3f  p = %.3f'
              % (rs_, perm_p(xs, ys, rs_)))

    # ---- 3. the second shell, from PROPKA's own decomposition
    print('\n--- 3. second shell of the pocket, from the PROPKA determinants ---')
    contrib = defaultdict(float)
    absc = defaultdict(float)
    acts = defaultdict(list)
    for resi in POCKET:
        for val, prn, pri, pch in det_lo.get((resi, a.chain), []):
            if pri in POCKET:
                continue                      # pocket carboxylate, gate 1 forbids it
            contrib[(pri, pch, prn)] += val
            absc[(pri, pch, prn)] += abs(val)
            acts[(pri, pch, prn)].append('%s%d%s%+.2f'
                                         % (POCKET[resi], resi, '' if pch == a.chain
                                            else '(x)', val))
    if not contrib:
        die('no second-shell determinants found; the determinant parse is wrong')

    print('%-9s %-5s %-8s %-8s %-28s %s'
          % ('residue', 'chain', 'sum', 'sum|.|', 'pI-code reachable to', 'acts on'))
    cands = []
    for (pri, pch, prn), v in sorted(absc.items(), key=lambda x: -x[1]):
        aa = TITRATABLE.get(prn) or {'GLN': 'Q', 'ASN': 'N', 'SER': 'S', 'THR': 'T',
                                      'TRP': 'W', 'GLY': 'G', 'ALA': 'A', 'PHE': 'F',
                                      'MET': 'M', 'LEU': 'L', 'ILE': 'I', 'VAL': 'V',
                                      'PRO': 'P'}.get(prn, '?')
        targets = sorted(reach.get(aa, []))
        if not targets:
            continue                          # not reachable by the ladder at all
        cands.append((v, pri, pch, prn, aa, targets, contrib[(pri, pch, prn)]))
        print('%-9s %-5s %+-8.2f %-8.2f %-28s %s'
              % ('%s%d' % (aa, pri), pch, contrib[(pri, pch, prn)], v,
                 ', '.join('%s>%s' % (aa, t) for t in targets),
                 ' '.join(sorted(acts[(pri, pch, prn)]))))

    print('\n  (x) marks a determinant acting across the subunit interface')
    print('  %d position-chain entries are pI-code reachable, over %d distinct '
          'sequence positions' % (len(cands), len({c[1] for c in cands})))

    print('\n--- 4. candidates ranked, with whatever measurement exists at the position ---')
    print('%-9s %-5s %-8s %-22s %s'
          % ('residue', 'chain', 'sum|.|', 'ladder edges', 'measured at this position'))
    for v, pri, pch, prn, aa, targets, sm in cands[:14]:
        m = meas.get(pri, [])
        mtxt = ('none' if not m else
                ', '.join('%s dAct %+.2f' % (r['mutation'], float(r['d_act']))
                          for r in sorted(m, key=lambda r: r['mut_aa'])[:4]))
        print('%-9s %-5s %-8.2f %-22s %s'
              % ('%s%d' % (aa, pri), pch, v,
                 ','.join('%s>%s' % (aa, t) for t in targets), mtxt))


if __name__ == '__main__':
    main()
