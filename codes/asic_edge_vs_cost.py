#!/usr/bin/env python3
"""
asic_edge_vs_cost.py - rank the measured pH50 shift against pI-code edge cost.

Usage:
    python3 asic_edge_vs_cost.py asic1a_ph50.csv \
        --ladder <ladder_edge_cost.py> --nontm <nontm_substitutions.csv.gz> \
        --humsavar-script <ladder_humsavar.py> --humsavar <humsavar.txt>

Both cost columns are recomputed on each run by executing the two scoring
scripts as subprocesses and parsing their stdout. Those scripts are part of the
unpublished pI-code method and are not in this repository, so this script does
not run from a clone; its output is committed as results/edge_costs.csv.

Reports three analyses: the raw rank correlation per edge, the same partialled
on charge-change magnitude, and a within-position paired test at equal charge
change. Significance is by exact permutation over all orderings of the eight
edges.

Exits non-zero on any fatal condition.
"""
import argparse, csv, itertools, math, os, re, subprocess, sys
from collections import defaultdict

# Charges for source-paper substitutions the ladder table omits, all neutral.
# For charge matching only; not added to the method's own CHARGE table.
NEUTRAL_EXTRA = {aa: 0 for aa in 'AVFMGLISTCPW'}


def die(msg):
    sys.exit('FATAL asic_edge_vs_cost.py: %s' % msg)


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


def run(script, data, tag):
    """Execute a canonical script and return its stdout. Non-zero exit is fatal."""
    try:
        p = subprocess.run([sys.executable, script, data], capture_output=True,
                           text=True, timeout=1800)
    except FileNotFoundError:
        die('%s not found' % script)
    except subprocess.TimeoutExpired:
        die('%s timed out on %s' % (script, data))
    if p.returncode != 0:
        die('%s exited %d on %s:\n%s' % (script, p.returncode, data,
                                         (p.stderr or p.stdout)[-800:]))
    if not p.stdout.strip():
        die('%s produced no output for %s' % (tag, data))
    return p.stdout


EDGE_CLIFF = re.compile(r'^([A-Z])>([A-Z])\s+\S+\s+\d+\s+\d+\s+[\d.]+\s+[\d.]+\s+'
                        r'([+-][\d.]+)')
EDGE_HUMS = re.compile(r'^([A-Z])>([A-Z])\s+\S+\s+\d+\s+\d+\s+\d+\s+([\d.]+)\s')


def parse_costs(text, rx, what):
    out = {}
    for line in text.splitlines():
        m = rx.match(line.strip())
        if m:
            out[(m.group(1), m.group(2))] = float(m.group(3))
    if not out:
        die('parsed no %s values; the output format of the canonical script has '
            'changed' % what)
    return out


def spearman(x, y):
    n = len(x)
    if n != len(y) or n < 3:
        die('spearman needs equal-length vectors of at least 3')

    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r

    rx, ry = rank(x), rank(y)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    if dx == 0 or dy == 0:
        return float('nan')
    return num / (dx * dy)


def perm_p(x, y, observed):
    """Exact two-sided permutation p for Spearman. n<=9 only."""
    n = len(x)
    if n > 9:
        die('exact permutation requested for n=%d' % n)
    hits = tot = 0
    for perm in itertools.permutations(range(n)):
        r = spearman(x, [y[i] for i in perm])
        tot += 1
        if abs(r) >= abs(observed) - 1e-12:
            hits += 1
    return hits / float(tot)


def partial(x, y, z):
    """Spearman of x,y controlling z, via the first-order partial formula."""
    rxy, rxz, ryz = spearman(x, y), spearman(x, z), spearman(y, z)
    den = math.sqrt(max(0.0, (1 - rxz ** 2) * (1 - ryz ** 2)))
    if den < 1e-12:
        return float('nan'), rxy, rxz, ryz
    return (rxy - rxz * ryz) / den, rxy, rxz, ryz


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('csv')
    ap.add_argument('--ladder', help='path to ladder_edge_cost.py; by default the method folder, then data/ladder_definition.json')
    ap.add_argument('--nontm')
    ap.add_argument('--humsavar-script')
    ap.add_argument('--humsavar')
    ap.add_argument('--readout', choices=['act', 'ssd'], default='act')
    a = ap.parse_args()

    missing = [n for n, v in (('--nontm', a.nontm),
                              ('--humsavar-script', a.humsavar_script),
                              ('--humsavar', a.humsavar)) if not v]
    if missing:
        die('%s not given. This script recomputes the edge costs by executing '
            'the pI-code scoring\nscripts, which are unpublished and not part of '
            'this repository, so it does not run\nfrom a clone. The per-edge '
            'output it produces is committed as results/edge_costs.csv,\nand '
            'that file is what every figure and table downstream of it reads.'
            % ', '.join(missing))

    L = load_ladder(a.ladder)
    charge = dict(NEUTRAL_EXTRA)
    charge.update(L.CHARGE)

    cliff = parse_costs(run(a.ladder, a.nontm, 'ladder_edge_cost.py'),
                        EDGE_CLIFF, "Cliff's delta")
    hums = parse_costs(run(a.humsavar_script, a.humsavar, 'ladder_humsavar.py'),
                       EDGE_HUMS, 'humsavar pathogenic fraction')
    print('recomputed %d edge Cliff deltas and %d humsavar fractions from source'
          % (len(cliff), len(hums)))

    col = 'd_act' if a.readout == 'act' else 'd_ssd'
    label = 'activation pH50' if a.readout == 'act' else 'SSD pHD50'
    rows = [r for r in csv.DictReader(open(a.csv)) if r['mutation'] != 'WT']
    if not rows:
        die('no substitution rows in %s' % a.csv)

    core = set(L.CORE)
    per_edge = defaultdict(list)
    for r in rows:
        p = (r['wt_aa'], r['mut_aa'])
        if p in core:
            per_edge[p].append(r)

    missing = [p for p in L.CORE if p not in per_edge]
    print('\nreadout: %s   core edges measured: %d of 12' % (label, len(per_edge)))
    if missing:
        print('not measured: %s' % ', '.join('%s>%s' % p for p in missing))

    edges = [p for p in L.CORE if p in per_edge]
    if len(edges) < 5:
        die('only %d core edges measured; a rank correlation on this is not worth '
            'reporting' % len(edges))

    shift, cl, hu, dq, ns = [], [], [], [], []
    print('\n%-7s %-4s %-10s %-10s %-9s %-8s %s'
          % ('edge', 'n', 'mean|d|', 'median|d|', 'cliff', 'hums%', '|dQ|'))
    for p in edges:
        v = [abs(float(r[col])) for r in per_edge[p]]
        v.sort()
        md = v[len(v) // 2] if len(v) % 2 else (v[len(v) // 2 - 1] + v[len(v) // 2]) / 2
        if p not in cliff:
            die('no Cliff delta recomputed for %s>%s' % p)
        if p not in hums:
            die('no humsavar fraction recomputed for %s>%s' % p)
        q = abs(charge[p[1]] - charge[p[0]])
        shift.append(sum(v) / len(v)); cl.append(cliff[p]); hu.append(hums[p])
        dq.append(q); ns.append(len(v))
        print('%-7s %-4d %-10.3f %-10.3f %+-9.3f %-8.1f %.1f'
              % (p[0] + '>' + p[1], len(v), sum(v) / len(v), md, cliff[p], hums[p], q))

    print('\n--- 1. raw rank correlation, mean |delta %s| per edge ---' % label)
    for nm, cost, sign in (("AlphaMissense Cliff delta", cl, 'negative'),
                           ('humsavar pathogenic fraction', hu, 'positive')):
        r = spearman(shift, cost)
        print('  vs %-30s rho = %+.3f   exact p = %.4f   (transfer predicts %s)'
              % (nm, r, perm_p(shift, cost, r), sign))

    print('\n--- 2. partialled on |delta charge| ---')
    print('  |delta charge| takes %d distinct values over %d edges: %s'
          % (len(set(dq)), len(dq), sorted(set(dq))))
    r_sq = spearman(shift, dq)
    print('  mean |delta %s| vs |delta charge|      rho = %+.3f   exact p = %.4f'
          % (label, r_sq, perm_p(shift, dq, r_sq)))
    for nm, cost in (("AlphaMissense Cliff delta", cl),
                     ('humsavar pathogenic fraction', hu)):
        pr, rxy, rxz, ryz = partial(shift, cost, dq)
        print('  %-30s  rho(raw) %+.3f  rho(cost,|dQ|) %+.3f  PARTIAL %+.3f'
              % (nm, rxy, ryz, pr))

    print('\n--- 3. within-position, charge-matched pairs ---')
    print('  one position at a time, a ladder edge against a non-ladder '
          'substitution of equal charge change')
    bypos = defaultdict(list)
    for r in rows:
        bypos[(r['wt_aa'], int(r['resi']))].append(r)
    allowed = set(L.CORE) | set(L.ARG) | set(L.EXT)
    pairs = []
    for (waa, resi), v in sorted(bypos.items(), key=lambda x: x[0][1]):
        for r in v:
            if (waa, r['mut_aa']) not in allowed:
                continue
            if r['mut_aa'] not in charge:
                continue
            ql = charge[r['mut_aa']] - charge[waa]
            for s in v:
                if s is r or (waa, s['mut_aa']) in allowed:
                    continue
                if s['mut_aa'] not in charge:
                    continue
                if abs((charge[s['mut_aa']] - charge[waa]) - ql) > 1e-9:
                    continue
                pairs.append((waa, resi, r['mutation'], abs(float(r[col])),
                              s['mutation'], abs(float(s[col])), ql))
    if not pairs:
        die('no charge-matched within-position pairs found; the matching logic '
            'is not doing what it claims')
    won = sum(1 for p in pairs if p[3] < p[5])
    tied = sum(1 for p in pairs if abs(p[3] - p[5]) < 1e-9)
    print('  %-9s %-9s %-9s %-9s %-7s %s'
          % ('position', 'ladder', '|d|', 'non-ladder', '|d|', 'dQ'))
    for waa, resi, lm, lv, sm, sv, q in pairs:
        print('  %s%-8d %-9s %-9.3f %-9s %-7.3f %+.1f  %s'
              % (waa, resi, lm, lv, sm, sv, q,
                 'ladder smaller' if lv < sv else 'ladder larger'))
    def sign_p(w, n):
        if n == 0:
            return float('nan')
        k = min(w, n - w)
        return min(1.0, 2 * sum(math.comb(n, i) for i in range(k + 1)) / 2.0 ** n)

    n = len(pairs) - tied
    print('\n  ladder substitution perturbs less in %d of %d charge-matched pairs'
          % (won, len(pairs)))
    print('  pair-level sign test, two-sided p = %.3f' % sign_p(won, n))

    # The sign test discards magnitude; the paired difference does not agree.
    diffs = [p[3] - p[5] for p in pairs]          # ladder minus non-ladder
    md = sum(diffs) / len(diffs)
    wins = [-d for d in diffs if d < 0]
    loss = [d for d in diffs if d > 0]
    print('  mean paired difference %+.4f pH (negative means the ladder perturbs '
          'less)' % md)
    if wins and loss:
        print('  when the ladder wins it wins by %.3f pH on average (n=%d); when it '
              'loses it loses by %.3f (n=%d)'
              % (sum(wins) / len(wins), len(wins), sum(loss) / len(loss), len(loss)))
        if sum(loss) / len(loss) > sum(wins) / len(wins):
            print('  THE SIGN TEST AND THE MAGNITUDE DISAGREE: the ladder wins more '
                  'often and by less than it loses. Quote the magnitude, and see '
                  'asic_fixed_effects.py, which is the properly powered form of this '
                  'comparison.')
    print('  NOTE: these pairs are not independent. One ladder mutant can appear in')
    print('  several pairs at its position, so the pair-level p overstates the')
    print('  evidence. The position-clustered test below is the one to quote.')

    # One vote per position.
    votes = defaultdict(lambda: [0, 0])
    for waa, resi, lm, lv, sm, sv, q in pairs:
        v = votes[(waa, resi)]
        if abs(lv - sv) < 1e-9:
            continue
        v[0 if lv < sv else 1] += 1
    pw = sum(1 for v in votes.values() if v[0] > v[1])
    pl = sum(1 for v in votes.values() if v[1] > v[0])
    print('\n  by position, one vote each: ladder smaller at %d positions, larger at '
          '%d, split at %d' % (pw, pl, len(votes) - pw - pl))
    print('  position-clustered sign test, two-sided p = %.3f' % sign_p(pw, pw + pl))
    against = sorted('%s%d' % k for k, v in votes.items() if v[1] > v[0])
    if against:
        print('  positions where the ladder substitution perturbs MORE: %s'
              % ', '.join(against))


if __name__ == '__main__':
    main()
