#!/usr/bin/env python3
"""
asic_numbering_map.py - residue mapping between structure numbering and human
ASIC1a, from the UniProt sequences.

Usage:
    python3 asic_numbering_map.py --chick Q1XA76.fasta --human P78348.fasta \
        [--rat P55926.fasta] [--pdb <assembly.pdb> --pdb-species chick|human]

The chicken to human offset is not constant across the extracellular domain, so
the mapping is derived from a global Needleman-Wunsch alignment (BLOSUM62) and
printed position by position, then checked residue by residue against a
deposited structure where one is given.

Exits non-zero if a position fails to map to the expected residue type, or if a
structure disagrees with the sequence.
"""
import argparse, sys

# Pocket pairs as the reviews state them, plus every measured position.
# Chicken numbering.
POCKET_CHICK = [238, 350, 239, 346, 220, 408, 80, 417]
# Human positions of interest: the paper's pH sensors plus the interface position.
HUMAN_OF_INTEREST = [73, 78, 211, 219, 237, 238, 242, 243, 347, 351, 375, 409,
                     413, 418, 433, 434]


def die(msg):
    sys.exit('FATAL asic_numbering_map.py: %s' % msg)


def read_fasta(p):
    try:
        lines = open(p).read().split('\n')
    except OSError as e:
        die('cannot read %s: %s' % (p, e))
    hdr = lines[0] if lines and lines[0].startswith('>') else ''
    seq = ''.join(l.strip() for l in lines[1:] if l and not l.startswith('>'))
    if not seq:
        die('no sequence in %s' % p)
    return hdr, seq


def align(a, b):
    from Bio import Align
    al = Align.PairwiseAligner()
    al.mode = 'global'
    try:
        from Bio.Align import substitution_matrices
        al.substitution_matrix = substitution_matrices.load('BLOSUM62')
    except Exception as e:
        die('cannot load BLOSUM62: %s' % e)
    al.open_gap_score, al.extend_gap_score = -11, -1
    res = al.align(a, b)
    if len(res) == 0:
        die('alignment produced no result')
    return res[0]


def index_map(aln):
    """Return {pos_in_a: pos_in_b} (1-based), gaps omitted."""
    m, ident = {}, 0
    ta, tb = aln[0], aln[1]
    ia = ib = 0
    for ca, cb in zip(ta, tb):
        if ca != '-':
            ia += 1
        if cb != '-':
            ib += 1
        if ca != '-' and cb != '-':
            m[ia] = ib
            if ca == cb:
                ident += 1
    return m, ident


def pdb_residues(path, chain=None):
    out = {}
    T3 = {'ASP': 'D', 'GLU': 'E', 'HIS': 'H', 'ASN': 'N', 'GLN': 'Q', 'LYS': 'K',
          'ARG': 'R', 'GLY': 'G', 'ALA': 'A', 'VAL': 'V', 'LEU': 'L', 'ILE': 'I',
          'PHE': 'F', 'TYR': 'Y', 'TRP': 'W', 'SER': 'S', 'THR': 'T', 'CYS': 'C',
          'MET': 'M', 'PRO': 'P'}
    for l in open(path, errors='replace'):
        if l[:6] != 'ATOM  ':
            continue
        ch = l[21]
        if chain and ch != chain:
            continue
        try:
            ri = int(l[22:26])
        except ValueError:
            continue
        out.setdefault(ri, T3.get(l[17:20].strip(), 'X'))
    if not out:
        die('no ATOM records parsed from %s' % path)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--chick', required=True)
    ap.add_argument('--human', required=True)
    ap.add_argument('--rat')
    ap.add_argument('--pdb')
    ap.add_argument('--pdb-species', choices=['chick', 'human'], default='chick')
    a = ap.parse_args()

    hc, sc = read_fasta(a.chick)
    hh, sh = read_fasta(a.human)
    print('chicken %s  length %d' % (hc.split()[0].lstrip('>'), len(sc)))
    print('human   %s  length %d' % (hh.split()[0].lstrip('>'), len(sh)))

    aln = align(sc, sh)
    c2h, ident = index_map(aln)
    print('chicken to human: %d aligned pairs, %.1f%% identity over the alignment'
          % (len(c2h), 100.0 * ident / max(1, len(c2h))))

    h2c = {v: k for k, v in c2h.items()}

    rat = None
    if a.rat:
        hr, sr = read_fasta(a.rat)
        alr = align(sh, sr)
        h2r, identr = index_map(alr)
        rat = h2r
        print('rat     %s  length %d' % (hr.split()[0].lstrip('>'), len(sr)))
        print('human to rat: %d aligned pairs, %.1f%% identity'
              % (len(h2r), 100.0 * identr / max(1, len(h2r))))

    print('\n--- acidic pocket, chicken structure numbering to human ASIC1a ---')
    print('%-14s %-14s %-8s %s' % ('chicken', 'human', 'offset', 'rat'))
    bad = []
    for cpos in POCKET_CHICK:
        if cpos not in c2h:
            die('chicken %d does not align to human; mapping incomplete' % cpos)
        hpos = c2h[cpos]
        ca, ha = sc[cpos - 1], sh[hpos - 1]
        if ca != ha:
            bad.append('chicken %s%d aligns to human %s%d' % (ca, cpos, ha, hpos))
        rr = ''
        if rat and hpos in rat:
            rpos = rat[hpos]
            rr = '%s%d' % (read_fasta(a.rat)[1][rpos - 1], rpos)
        print('%-14s %-14s %-8s %s' % ('%s%d' % (ca, cpos), '%s%d' % (ha, hpos),
                                       '%+d' % (hpos - cpos), rr))
    if bad:
        die('residue type changes across the alignment at: ' + '; '.join(bad))

    offsets = sorted({c2h[p] - p for p in POCKET_CHICK})
    print('\ndistinct chicken-to-human offsets over the pocket: %s'
          % ', '.join('%+d' % o for o in offsets))
    if len(offsets) > 1:
        print('A SINGLE GLOBAL OFFSET IS WRONG for this protein. The offset changes '
              'sign inside the extracellular domain.')

    print('\n--- human positions used in this work, back-mapped to chicken ---')
    print('%-10s %-10s %s' % ('human', 'chicken', 'offset'))
    for hpos in HUMAN_OF_INTEREST:
        if hpos not in h2c:
            print('%-10s %-10s %s' % ('%s%d' % (sh[hpos - 1], hpos), 'no align', '-'))
            continue
        cpos = h2c[hpos]
        print('%-10s %-10s %+d' % ('%s%d' % (sh[hpos - 1], hpos),
                                   '%s%d' % (sc[cpos - 1], cpos), hpos - cpos))

    if a.pdb:
        seq = sc if a.pdb_species == 'chick' else sh
        pres = pdb_residues(a.pdb, 'A')
        checked = mism = 0
        for ri, aa in sorted(pres.items()):
            if 1 <= ri <= len(seq):
                checked += 1
                if seq[ri - 1] != aa:
                    mism += 1
                    if mism <= 8:
                        print('  structure %s%d vs sequence %s%d'
                              % (aa, ri, seq[ri - 1], ri))
        print('\nstructure check %s (%s): %d residues compared, %d mismatches'
              % (a.pdb.split('/')[-1], a.pdb_species, checked, mism))
        if mism:
            die('%d residues in %s do not match the %s sequence, so the deposited '
                'numbering is not UniProt numbering' % (mism, a.pdb, a.pdb_species))


if __name__ == '__main__':
    main()
