#!/usr/bin/env python3
"""
extract_interface_carboxylates.py - write results/interface_carboxylates.csv.

Usage:
    python3 codes/extract_interface_carboxylates.py

Executes targets/interface_carboxylates.py over the four deposited assemblies
and parses its table. That script is part of the unpublished pI-code method and
is not in this repository, so this script does not run from a clone; its output
is committed as results/interface_carboxylates.csv, which is what everything
downstream reads.

Counts are buried, charge-uncompensated carboxylates at the thresholds of
data/burial_thresholds.json, scored once in the monomer frame and once in the
assembly frame. Exits non-zero if the script fails, if its output cannot be
parsed, or if any expected assembly is missing from the table.
"""
import csv, os, re, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import picode_paths as P

ENTRIES = ['9E4F', '9E4A', '6VTL', '2QTS']

ROW = re.compile(r'^(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+\|\s+'
                 r'(\d+)/(\d+)/(\d+)/(\d+)\s+\|\s+'
                 r'(\d+)/(\d+)/(\d+)/(\d+)\s+\|\s+(\d+)\s*(.*)$')


def die(msg):
    sys.exit('FATAL extract_interface_carboxylates.py: %s' % msg)


def main():
    script = P.interface_carboxylates()
    paths = []
    for e in ENTRIES:
        p = P.project('structures', '%s.pdb1' % e)
        if not os.path.isfile(p):
            die('%s is missing; fetch the assemblies first with\n'
                '  python3 codes/fetch_structures.py' % os.path.relpath(p, P.project()))
        paths.append(p)

    env = dict(os.environ)
    env['PICODE_LPO_DIR'] = os.path.dirname(P.lpo_design())
    try:
        r = subprocess.run([sys.executable, script] + paths, capture_output=True,
                           text=True, timeout=3600, env=env)
    except subprocess.TimeoutExpired:
        die('%s timed out' % script)
    if r.returncode != 0:
        die('%s exited %d:\n%s' % (script, r.returncode,
                                   (r.stderr or r.stdout)[-900:]))

    rows = {}
    for line in r.stdout.splitlines():
        m = ROW.match(line.strip())
        if not m:
            continue
        g = m.groups()
        name = g[0]
        mono = [int(x) for x in g[5:9]]
        asm = [int(x) for x in g[9:13]]
        rows[name] = dict(structure=name, chains=int(g[2]), residues=int(g[3]),
                          carboxylates=int(g[4]),
                          monomer_tier0=mono[0], assembly_tier0=asm[0],
                          added_tier0=asm[0] - mono[0],
                          interface_tier0=int(g[13]),
                          interface_positions=' '.join(g[14].split()))
    if not rows:
        die('parsed no rows; interface_carboxylates.py output format has changed')
    missing = [e for e in ENTRIES if e not in rows]
    if missing:
        die('no counts recomputed for %s' % ', '.join(missing))
    bad = [e for e in ENTRIES if rows[e]['added_tier0'] < 0]
    if bad:
        die('assembly tier 0 below monomer tier 0 for %s, which cannot happen; '
            'check the thresholds' % ', '.join(bad))

    path = P.project('results', 'interface_carboxylates.csv')
    os.makedirs(os.path.dirname(path), exist_ok=True)
    cols = ['structure', 'chains', 'residues', 'carboxylates', 'monomer_tier0',
            'assembly_tier0', 'added_tier0', 'interface_tier0',
            'interface_positions']
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for e in ENTRIES:
            w.writerow(rows[e])
    added = [rows[e]['added_tier0'] for e in ENTRIES]
    print('recomputed buried uncompensated carboxylates for %d assemblies'
          % len(ENTRIES))
    for e in ENTRIES:
        r_ = rows[e]
        print('  %-6s monomer %d, assembly %d, added %d'
              % (e, r_['monomer_tier0'], r_['assembly_tier0'], r_['added_tier0']))
    print('the assembly frame adds %d to %d' % (min(added), max(added)))
    print('wrote %s' % path)


if __name__ == '__main__':
    main()
