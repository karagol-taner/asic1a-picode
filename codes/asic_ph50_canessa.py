#!/usr/bin/env python3
"""
asic_ph50_canessa.py - parse the second ASIC1a pH50 cohort.

Usage:
    python3 asic_ph50_canessa.py PMC6497441.xml -o asic1a_ph50_canessa.csv

Source: Table 1 of

  Wu Y, Chen Z, Canessa CM. A valve-like mechanism controls desensitization of
  functional mammalian isoforms of acid-sensing ion channels. eLife 2019;8:e45851.
  PMCID PMC6497441.

Read from the Europe PMC full-text XML. This cohort is never pooled with the
primary one: different laboratory, different solutions and SD rather than SEM,
so only within-cohort shifts against its own wild type are comparable.

The columns headed N and n in Table 1 are Hill coefficients, not sample sizes.
Sample size is given only at table level, as five to eight oocytes, and is
stored as a range.

Exits non-zero on any fatal condition.
"""
import argparse, csv, re, sys
from xml.etree import ElementTree as ET

CENSORED = re.compile(r'^\s*[<>]')
VALUE = re.compile(r'^\s*(\d\.\d{1,2})\s*±\s*(\d\.\d{1,3})\s*$')
SUB = re.compile(r'^([A-Z])(\d{1,3})([A-Z])$')

N_MIN, N_MAX = 5, 8          # Methods: "five to eight independent measurements"
PH_MIN, PH_MAX = 3.0, 9.0


def die(msg):
    sys.exit('FATAL asic_ph50_canessa.py: %s' % msg)


def table1(path):
    try:
        raw = open(path, encoding='utf-8', errors='replace').read()
    except OSError as e:
        die('cannot read %s: %s' % (path, e))
    blocks = re.findall(r'<table-wrap\b.*?</table-wrap>', raw, re.S)
    if not blocks:
        die('no <table-wrap> in %s; not a Europe PMC full-text XML?' % path)
    for blk in blocks:
        lab = re.search(r'<label>(.*?)</label>', blk, re.S)
        lab = re.sub('<[^>]+>', '', lab.group(1)) if lab else ''
        if 'Table 1' not in lab:
            continue
        cap = re.search(r'<caption>(.*?)</caption>', blk, re.S)
        cap = re.sub(r'\s+', ' ', re.sub('<[^>]+>', ' ', cap.group(1))).strip() if cap else ''
        if 'Hill coefficient' not in cap:
            die('Table 1 caption does not mention Hill coefficients; the column '
                'meaning cannot be confirmed, so parsing is unsafe. Caption: %r' % cap[:200])
        try:
            root = ET.fromstring(re.sub(r'xmlns="[^"]+"', '', blk))
        except ET.ParseError as e:
            die('cannot parse Table 1 XML: %s' % e)
        return root, cap
    die('Table 1 not found among %d tables in %s' % (len(blocks), path))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('xml')
    ap.add_argument('-o', '--out', default='asic1a_ph50_canessa.csv')
    a = ap.parse_args()

    root, cap = table1(a.xml)
    rows, wt, seen = [], None, set()
    for tr in root.iter('tr'):
        cells = [''.join(c.itertext()).strip() for c in tr]
        cells = [c for c in cells if c != '']
        if len(cells) < 3:
            continue
        name = cells[0]
        if name in ('pH50a', 'pH50A'):
            continue
        rest = cells[1:]

        def take(i):
            """Return (value, hill, censored_text) for the i-th measurement block."""
            if i >= len(rest):
                return None, None, None
            v = rest[i]
            h = rest[i + 1] if i + 1 < len(rest) else ''
            if CENSORED.match(v):
                return None, None, v.strip()
            m = VALUE.match(v)
            if not m:
                return None, None, None
            return (float(m.group(1)), float(m.group(2))), h.strip(), None

        act, hill_a, cens_a = take(0)
        ssd, hill_s, cens_s = take(2)
        if act is None and cens_a is None:
            continue

        rec = dict(mutation=name,
                   ph50_act=act[0] if act else '', sd_act=act[1] if act else '',
                   hill_act=hill_a or '', act_censored=cens_a or '',
                   ph50_ssd=ssd[0] if ssd else '', sd_ssd=ssd[1] if ssd else '',
                   hill_ssd=hill_s or '', ssd_censored=cens_s or '',
                   n_min=N_MIN, n_max=N_MAX)
        for k in ('ph50_act', 'ph50_ssd'):
            if rec[k] != '' and not PH_MIN <= rec[k] <= PH_MAX:
                die('%s %s = %s is outside pH %.1f-%.1f' % (name, k, rec[k], PH_MIN, PH_MAX))
        if name.upper() == 'WT':
            if wt is not None:
                die('two WT rows parsed')
            wt = rec
            continue
        s = SUB.match(name)
        if not s:
            continue                      # section headers and non-substitution rows
        if name in seen:
            die('%s appears twice' % name)
        seen.add(name)
        rec.update(wt_aa=s.group(1), resi=int(s.group(2)), mut_aa=s.group(3))
        rows.append(rec)

    if wt is None:
        die('no WT row in Table 1, so no shift can be computed')
    if wt['ph50_act'] == '' or wt['ph50_ssd'] == '':
        die('the WT row is censored in one readout; shifts cannot be formed')
    if len(rows) < 30:
        die('only %d substitutions parsed, expected at least 30; the table layout '
            'has changed' % len(rows))

    for r in rows:
        r['d_act'] = round(r['ph50_act'] - wt['ph50_act'], 4) if r['ph50_act'] != '' else ''
        r['d_ssd'] = round(r['ph50_ssd'] - wt['ph50_ssd'], 4) if r['ph50_ssd'] != '' else ''

    pos = sorted({r['resi'] for r in rows})
    print('source: %s' % a.xml)
    print('caption confirms the N columns are Hill coefficients, not sample sizes')
    print('WT hASIC1a  activation %.2f ± %.2f (SD)   SSD %.2f ± %.2f (SD)   n = %d to %d'
          % (wt['ph50_act'], wt['sd_act'], wt['ph50_ssd'], wt['sd_ssd'], N_MIN, N_MAX))
    print('substitutions parsed: %d over %d positions %s'
          % (len(rows), len(pos), pos))
    ncens = sum(1 for r in rows if r['ssd_censored'])
    print('SSD censored (reported as an inequality, unusable): %d of %d' % (ncens, len(rows)))

    cols = ['mutation', 'wt_aa', 'resi', 'mut_aa', 'ph50_act', 'sd_act', 'hill_act',
            'act_censored', 'd_act', 'ph50_ssd', 'sd_ssd', 'hill_ssd', 'ssd_censored',
            'd_ssd', 'n_min', 'n_max']
    with open(a.out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction='ignore')
        w.writeheader()
        w.writerow(dict(wt, wt_aa='', resi=0, mut_aa='', d_act=0.0, d_ssd=0.0))
        for r in sorted(rows, key=lambda x: (x['resi'], x['mut_aa'])):
            w.writerow(r)
    print('wrote %s' % a.out)


if __name__ == '__main__':
    main()
