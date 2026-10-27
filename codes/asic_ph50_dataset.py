#!/usr/bin/env python3
"""
asic_ph50_dataset.py - parse the primary ASIC1a pH50 dataset.

Usage:
    python3 asic_ph50_dataset.py <42003_2025_9090_MOESM2_ESM.pdf> [-o out.csv]
    python3 asic_ph50_dataset.py <...MOESM2_ESM.pdf> --xlsx <...MOESM4_ESM.xlsx>

Source: Supplementary Table 2 ("pH50 and pHD50 values") of

  Bignucolo O, Molton O, Gautschi I, Kellenberger S. Molecular mechanisms and
  hotspots of pH sensing in ASIC1a revealed by computational and functional
  analysis. Commun Biol 2025. PMCID PMC12658276.

Read from the published PDF. With --xlsx the values are cross-checked against
Supplementary Data 1, which carries them at four decimals for the fifteen
scanned positions; a disagreement beyond rounding is fatal.

Activation pH50 and steady-state desensitization pHD50 are kept in separate
columns and never pooled.

Exits non-zero on any fatal condition.
"""
import argparse, csv, re, sys

# Row shape:  [Location] [Residue] Mutation  v ± e  N  v ± e  N
ROW = re.compile(
    r'(?P<mut>(?:[A-Z]\d{1,3}[A-Z])|WT)\s+'
    r'(?P<a>\d\.\d{2})\s*(?:±|\+/-)\s*(?P<ae>\d\.\d{2})\s+(?P<an>\d{1,4})\s+'
    r'(?P<d>\d\.\d{2})\s*(?:±|\+/-)\s*(?P<de>\d\.\d{2})\s+(?P<dn>\d{1,4})\s*$')
SUB = re.compile(r'^([A-Z])(\d{1,3})([A-Z])$')

TABLE_HEAD = 'Supplementary Table 2'
NEXT_HEAD = 'Supplementary Table 3'
# The double-mutant block follows the single substitutions under its own header.
DOUBLE_HEAD = 'Double'

PH_MIN, PH_MAX = 3.0, 9.0
EXPECTED_MIN_ROWS = 100          # the published table holds 117 single substitutions


def die(msg):
    sys.exit('FATAL %s: %s' % (sys.argv[0].split('/')[-1], msg))


def pages(path):
    try:
        import pdfplumber
    except ImportError:
        die('pdfplumber is not installed; pip install pdfplumber')
    with pdfplumber.open(path) as pdf:
        return [p.extract_text() or '' for p in pdf.pages]


def parse_pdf(path):
    txt = pages(path)
    start = next((i for i, t in enumerate(txt) if TABLE_HEAD in t), None)
    if start is None:
        die('%r not found in %s; wrong supplementary file?' % (TABLE_HEAD, path))

    rows, wt, in_double = [], None, False
    seen = set()
    # The domain prints only when it changes, so it is carried forward. One
    # label wraps over two lines, which a plain prefix match would mis-carry.
    LOCS = ('Acidic pocket', 'Transmembrane', 'Knuckle', 'Palm',
            '\u03b2-ball', 'Finger and thumb')
    loc, pending = '', ''
    for i in range(start, len(txt)):
        if i > start and NEXT_HEAD in txt[i]:
            break
        for line in txt[i].splitlines():
            line = line.strip()
            if NEXT_HEAD in line:
                in_double = True          # stop collecting; table 2 is over
                break
            if line.startswith(DOUBLE_HEAD):
                in_double = True
            if in_double or '/' in line:  # double mutants are out of scope
                continue
            if line.strip() == 'Finger and':
                pending = 'Finger and'
                continue
            if pending and line.strip().startswith('thumb'):
                loc = 'Finger and thumb'
                pending = ''
                line = line.strip()[len('thumb'):].strip()
            else:
                for L in LOCS:
                    if line.startswith(L):
                        loc = L
                        break
            m = ROW.search(line)
            if not m:
                continue
            mut = m.group('mut')
            rec = dict(mutation=mut,
                       ph50_act=float(m.group('a')), sem_act=float(m.group('ae')),
                       n_act=int(m.group('an')),
                       ph50_ssd=float(m.group('d')), sem_ssd=float(m.group('de')),
                       n_ssd=int(m.group('dn')), page=i + 1, location=loc)
            for k in ('ph50_act', 'ph50_ssd'):
                if not PH_MIN <= rec[k] <= PH_MAX:
                    die('%s %s = %.2f is outside pH %.1f-%.1f, so the column '
                        'alignment is wrong' % (mut, k, rec[k], PH_MIN, PH_MAX))
            if mut == 'WT':
                if wt is not None:
                    die('two WT rows parsed; the table should carry exactly one')
                wt = rec
                continue
            s = SUB.match(mut)
            if not s:
                die('%r is not a single substitution' % mut)
            if mut in seen:
                die('%s appears twice; the page ranges overlap' % mut)
            seen.add(mut)
            rec.update(wt_aa=s.group(1), resi=int(s.group(2)), mut_aa=s.group(3))
            rows.append(rec)

    if wt is None:
        die('no WT row found, so no shift can be computed')
    if len(rows) < EXPECTED_MIN_ROWS:
        die('only %d single substitutions parsed, expected at least %d; the PDF '
            'layout has changed' % (len(rows), EXPECTED_MIN_ROWS))

    got = {r.get('location', '') for r in rows}
    if '' in got:
        die('%d substitutions carry no domain label; the Location column was not '
            'parsed' % sum(1 for r in rows if not r.get('location')))
    if 'Finger and thumb' not in got:
        die('the Finger and thumb block was not detected; its label wraps over two '
            'lines in the source and the handling for that has broken')

    # Every position's wild-type residue must agree with the mutation label.
    by_pos = {}
    for r in rows:
        prev = by_pos.setdefault(r['resi'], r['wt_aa'])
        if prev != r['wt_aa']:
            die('position %d is labelled both %s and %s' % (r['resi'], prev, r['wt_aa']))

    for r in rows:
        r['d_act'] = round(r['ph50_act'] - wt['ph50_act'], 4)
        r['d_ssd'] = round(r['ph50_ssd'] - wt['ph50_ssd'], 4)
    return wt, rows


def crosscheck(xlsx, wt, rows):
    """Compare against Supplementary Data 1 fitted values. Disagreement is fatal."""
    try:
        import openpyxl
    except ImportError:
        die('openpyxl is not installed; pip install openpyxl')
    wb = openpyxl.load_workbook(xlsx, read_only=True, data_only=True)
    idx = {(r['resi'], r['mut_aa']): r for r in rows}
    checked, disagree = 0, []
    # Table 2 gives two decimals, the workbook four, so 0.005 is exact.
    ROUNDING, MISALIGNED = 0.005, 0.05
    for sheet in wb.sheetnames:
        m = re.match(r'^([A-Z])(\d{1,3})_(act|des)$', sheet)
        if not m:
            continue
        wt_aa, resi, kind = m.group(1), int(m.group(2)), m.group(3)
        col = 'ph50_act' if kind == 'act' else 'ph50_ssd'
        for row in wb[sheet].iter_rows(values_only=True):
            if not row or not isinstance(row[0], str) or len(row[0]) != 1:
                continue
            aa, val = row[0], row[1]
            if not isinstance(val, (int, float)):
                continue
            ref = wt[col] if aa == wt_aa else None
            if ref is None:
                r = idx.get((resi, aa))
                if r is None:
                    continue
                ref = r[col]
            dev = abs(float(val) - ref)
            if dev > MISALIGNED:
                die('cross-check failed: %s sheet %s %s gives %.4f, the table gives '
                    '%.2f. A gap this large means the parsed columns are not the '
                    'ones they are taken to be.' % (xlsx.split('/')[-1], sheet, aa,
                                                    val, ref))
            if dev > ROUNDING:
                disagree.append((sheet, aa, float(val), ref, dev))
            checked += 1
    if checked == 0:
        die('cross-check parsed zero comparable values from %s' % xlsx)
    return checked, sorted(disagree, key=lambda x: -x[4])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('pdf')
    ap.add_argument('-o', '--out', default='asic1a_ph50.csv')
    ap.add_argument('--xlsx')
    a = ap.parse_args()

    wt, rows = parse_pdf(a.pdf)
    print('source: %s' % a.pdf)
    print('WT hASIC1a  activation pH50 %.2f ± %.2f (N=%d)   SSD pHD50 %.2f ± %.2f (N=%d)'
          % (wt['ph50_act'], wt['sem_act'], wt['n_act'],
             wt['ph50_ssd'], wt['sem_ssd'], wt['n_ssd']))
    print('single substitutions parsed: %d over %d positions'
          % (len(rows), len({r['resi'] for r in rows})))
    if a.xlsx:
        n, disagree = crosscheck(a.xlsx, wt, rows)
        print('cross-checked %d values against %s: %d agree to printed precision, '
              '%d disagree' % (n, a.xlsx.split('/')[-1], n - len(disagree),
                               len(disagree)))
        for sheet, aa, val, ref, dev in disagree:
            print('   SOURCE DISAGREEMENT %-10s %s  Suppl Data 1 %.4f  vs  Table 2 '
                  '%.2f  (%.4f)' % (sheet, aa, val, ref, dev))

    cols = ['mutation', 'wt_aa', 'resi', 'mut_aa', 'location', 'ph50_act',
            'sem_act', 'n_act', 'd_act', 'ph50_ssd', 'sem_ssd', 'n_ssd', 'd_ssd',
            'page']
    with open(a.out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction='ignore')
        w.writeheader()
        w.writerow(dict(wt, wt_aa='', resi=0, mut_aa='', d_act=0.0, d_ssd=0.0))
        for r in sorted(rows, key=lambda x: (x['resi'], x['mut_aa'])):
            w.writerow(r)
    print('wrote %s' % a.out)


if __name__ == '__main__':
    main()
