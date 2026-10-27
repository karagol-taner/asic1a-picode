#!/usr/bin/env python3
"""
fetch_structures.py - download the deposited assemblies into structures/.

    python3 codes/fetch_structures.py [--force]

The four assemblies are public RCSB depositions and are fetched rather than
carried here; the PROPKA output beside them is derived and is committed. Each
download is checked for chain count, ATOM count and residue span before it is
written.

Exits non-zero if any download or check fails.
"""
import argparse, gzip, io, os, sys, urllib.error, urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import picode_paths as P

URL = 'https://files.rcsb.org/download/%s.pdb1.gz'

# entry: (chains, ATOM records, first residue, last residue)
ENTRIES = {'9E4F': (3, 10842, 15, 465),
           '9E4A': (3, 10842, 15, 465),
           '6VTL': (3, 10158, 17, 462),
           '2QTS': (3, 10064, 40, 461)}


def die(msg):
    sys.exit('FATAL fetch_structures.py: %s' % msg)


def check(text, entry):
    want_ch, want_at, want_lo, want_hi = ENTRIES[entry]
    atoms = [l for l in text.splitlines() if l.startswith('ATOM')]
    if len(atoms) != want_at:
        return '%d ATOM records, expected %d' % (len(atoms), want_at)
    chains = sorted({l[21] for l in atoms})
    if len(chains) != want_ch:
        return 'chains %s, expected %d' % (chains, want_ch)
    res = [int(l[22:26]) for l in atoms]
    if min(res) != want_lo or max(res) != want_hi:
        return ('residues %d-%d, expected %d-%d'
                % (min(res), max(res), want_lo, want_hi))
    return None


def fetch(entry, dest, force):
    out = os.path.join(dest, entry + '.pdb1')
    if os.path.isfile(out) and not force:
        with open(out) as fh:
            bad = check(fh.read(), entry)
        if bad:
            die('%s is already present but fails its check: %s\n'
                'Delete it and rerun, or rerun with --force.' % (out, bad))
        print('  %-12s present' % (entry + '.pdb1'))
        return
    url = URL % entry
    try:
        with urllib.request.urlopen(url, timeout=120) as r:
            raw = r.read()
    except urllib.error.URLError as e:
        die('cannot download %s: %s' % (url, e))
    try:
        text = gzip.GzipFile(fileobj=io.BytesIO(raw)).read().decode('utf-8')
    except Exception as e:
        die('cannot decompress %s: %s' % (url, e))
    bad = check(text, entry)
    if bad:
        die('%s downloaded but fails its check: %s' % (entry, bad))
    with open(out, 'w') as fh:
        fh.write(text)
    print('  %-12s %6.1f KB' % (entry + '.pdb1', len(text) / 1024.0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--force', action='store_true',
                    help='download again even where the file is present')
    a = ap.parse_args()
    dest = P.project('structures')
    if not os.path.isdir(dest):
        os.makedirs(dest)
    print('structures/')
    for entry in sorted(ENTRIES):
        fetch(entry, dest, a.force)
    missing = [e for e in ENTRIES
               if not os.path.isfile(os.path.join(dest, e + '.pdb1'))]
    if missing:
        die('still missing: %s' % ', '.join(sorted(missing)))
    print('\n%d assemblies in %s' % (len(ENTRIES), dest))


if __name__ == '__main__':
    main()
