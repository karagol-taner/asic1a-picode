#!/usr/bin/env python3
"""
extract_edge_costs.py - write results/edge_costs.csv.

Usage:
    python3 codes/extract_edge_costs.py

Executes targets/ladder_edge_cost.py over nontm_substitutions.csv.gz and
targets/ladder_humsavar.py over humsavar.txt and parses their stdout. Those
scripts are part of the unpublished pI-code method and are not in this
repository, so this script does not run from a clone; its output is committed
as results/edge_costs.csv, which is what everything downstream reads.

Exits non-zero if either script fails, if its output cannot be parsed, or if
any of the twelve core edges is missing.
"""
import csv, os, re, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import picode_paths as P

CLIFF = re.compile(r'^([A-Z])>([A-Z])\s+(\S+)\s+(\d+)\s+(\d+)\s+([\d.]+)\s+[\d.]+\s+'
                   r'([+-][\d.]+)')
HUMS = re.compile(r'^([A-Z])>([A-Z])\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([\d.]+)\s')


def die(msg):
    sys.exit('FATAL extract_edge_costs.py: %s' % msg)


def run(script, data):
    try:
        p = subprocess.run([sys.executable, script, data], capture_output=True,
                           text=True, timeout=3600)
    except subprocess.TimeoutExpired:
        die('%s timed out' % script)
    if p.returncode != 0:
        die('%s exited %d:\n%s' % (script, p.returncode, (p.stderr or p.stdout)[-900:]))
    return p.stdout


def main():
    L = P.load_ladder()
    out_c = run(P.ladder_edge_cost(), P.nontm_substitutions())
    out_h = run(P.ladder_humsavar(), P.humsavar())

    cliff, hums = {}, {}
    for line in out_c.splitlines():
        m = CLIFF.match(line.strip())
        if m:
            cliff[(m.group(1), m.group(2))] = dict(
                direction=m.group(3), n=int(m.group(4)), snv_pct=int(m.group(5)),
                mean=float(m.group(6)), cliff_delta=float(m.group(7)))
    for line in out_h.splitlines():
        m = HUMS.match(line.strip())
        if m:
            hums[(m.group(1), m.group(2))] = dict(
                lp=int(m.group(4)), lb=int(m.group(5)), us=int(m.group(6)),
                pct=float(m.group(7)))
    if not cliff:
        die('parsed no Cliff deltas; ladder_edge_cost.py output format has changed')
    if not hums:
        die('parsed no humsavar fractions; ladder_humsavar.py output format has changed')

    missing = [p for p in L.CORE if p not in cliff or p not in hums]
    if missing:
        die('no cost recomputed for core edges %s'
            % ', '.join('%s>%s' % p for p in missing))

    kinds = {}
    for p in L.CORE:
        kinds[p] = 'core'
    for p in L.ARG:
        kinds[p] = 'arginine'
    for p in L.EXT:
        kinds.setdefault(p, 'extension')

    path = P.project('results', 'edge_costs.csv')
    os.makedirs(os.path.dirname(path), exist_ok=True)
    cols = ['source', 'target', 'kind', 'direction', 'delta_charge',
            'alphamissense_n', 'alphamissense_mean', 'cliff_delta',
            'humsavar_lp_p', 'humsavar_lb_b', 'humsavar_pct']
    nrow = 0
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for p in list(L.CORE) + list(L.ARG) + list(L.EXT):
            if p not in cliff or p not in hums:
                continue
            c, h = cliff[p], hums[p]
            w.writerow(dict(source=p[0], target=p[1], kind=kinds[p],
                            direction=c['direction'],
                            delta_charge=round(L.CHARGE[p[1]] - L.CHARGE[p[0]], 3),
                            alphamissense_n=c['n'], alphamissense_mean=c['mean'],
                            cliff_delta=c['cliff_delta'],
                            humsavar_lp_p=h['lp'], humsavar_lb_b=h['lb'],
                            humsavar_pct=h['pct']))
            nrow += 1
    print('recomputed %d edge costs from source' % nrow)
    print('  AlphaMissense: %s' % P.nontm_substitutions())
    print('  humsavar:      %s' % P.humsavar())
    print('wrote %s' % path)


if __name__ == '__main__':
    main()
