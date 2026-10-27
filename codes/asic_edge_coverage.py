#!/usr/bin/env python3
"""
asic_edge_coverage.py - which ladder edges the measured pH50 set covers.

Usage:
    python3 asic_edge_coverage.py asic1a_ph50.csv

Edge definitions are imported from the ladder definition, not restated. Exits
non-zero if the import fails or the constants are not the expected shape.
"""
import argparse, csv, os, sys
from collections import defaultdict


def die(msg):
    sys.exit('FATAL asic_edge_coverage.py: %s' % msg)


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('csv')
    ap.add_argument('--ladder', help='path to ladder_edge_cost.py; by default the method folder, then data/ladder_definition.json')
    a = ap.parse_args()
    L = load_ladder(a.ladder)

    rows = [r for r in csv.DictReader(open(a.csv)) if r['mutation'] != 'WT']
    if not rows:
        die('no substitution rows in %s' % a.csv)

    kinds = {}
    for p in L.CORE:
        kinds[p] = 'core'
    for p in L.ARG:
        kinds[p] = 'arginine'
    for p in L.EXT:
        kinds.setdefault(p, 'extension')

    hit = defaultdict(list)
    for r in rows:
        p = (r['wt_aa'], r['mut_aa'])
        if p in kinds:
            hit[p].append(r)

    print('measured substitutions: %d over %d positions' % (len(rows),
          len({r['resi'] for r in rows})))
    print('of these, on a pI-code edge: %d\n' % sum(len(v) for v in hit.values()))

    for tag, pairs in (('12 core edges', L.CORE), ('5 arginine edges', L.ARG),
                       ('proposed extensions', L.EXT)):
        print('--- %s ---' % tag)
        print('%-7s %-8s %-5s %-9s %-9s %s' % ('edge', 'dir', 'n', 'mean dAct',
                                               'mean dSSD', 'positions'))
        cov = 0
        for p in pairs:
            v = hit.get(p, [])
            if not v:
                print('%-7s %-8s %-5s %-9s %-9s %s'
                      % (p[0] + '>' + p[1],
                         'basic' if p in L.BASICDIR else
                         'acidic' if (p[1], p[0]) in L.BASICDIR else '-',
                         0, '-', '-', 'NOT MEASURED'))
                continue
            cov += 1
            da = sum(float(r['d_act']) for r in v) / len(v)
            ds = sum(float(r['d_ssd']) for r in v) / len(v)
            print('%-7s %-8s %-5d %+-9.3f %+-9.3f %s'
                  % (p[0] + '>' + p[1],
                     'basic' if p in L.BASICDIR else
                     'acidic' if (p[1], p[0]) in L.BASICDIR else '-',
                     len(v), da, ds,
                     ','.join(sorted(r['mutation'] for r in v))))
        print('covered: %d of %d\n' % (cov, len(pairs)))

    # Positions with more than one measured edge allow a within-position test.
    bypos = defaultdict(list)
    for p, v in hit.items():
        for r in v:
            bypos[(r['wt_aa'], int(r['resi']))].append((p, r))
    multi = {k: v for k, v in bypos.items() if len(v) > 1}
    print('positions carrying 2 or more measured ladder edges: %d' % len(multi))
    for (aa, resi), v in sorted(multi.items(), key=lambda x: x[0][1]):
        print('  %s%-4d  %s' % (aa, resi, '  '.join(
            '%s>%s %+.2f' % (p[0], p[1], float(r['d_act'])) for p, r in
            sorted(v, key=lambda x: x[0][1]))))


if __name__ == '__main__':
    main()
