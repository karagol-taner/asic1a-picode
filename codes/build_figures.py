#!/usr/bin/env python3
"""
build_figures.py - regenerate the figures for the ASIC1a pI-code paper.

Usage:
    python3 codes/build_figures.py          (run from the project root)

Writes figures/Figure1 to Figure6 and FigureS1 to FigureS2, each as PNG at
400 dpi and as SVG. Every number plotted is read from data/ or structures/ at
run time.

Figure 1  the ladder as two codon triangles, and which edges have a measured pH50
Figure 2  measured shift against edge cost, and the charge confound behind it
Figure 3  the within-position charge-matched comparison
Figure 4  burial, distance and domain against the measured shift
Figure 5  the acidic pocket: pairs within a subunit, cavity across two
Figure 6  the second shell, and what sets the pocket pKa

Supplementary Figure S3 is written by build_md_figure.py.

Exits non-zero on any fatal condition.
"""
import csv, importlib.util, math, os, sys
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle
from matplotlib.transforms import Bbox

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import picode_paths as P

ROOT = P.ROOT
DATA = os.path.join(ROOT, 'data')
STRUCT = os.path.join(ROOT, 'structures')
FIGS = os.path.join(ROOT, 'figures')

INK = '#1a1a1a'
GREY = '#9a9a9a'
LAD = '#2166ac'
NON = '#b2182b'
ACC = '#4d4d4d'
FILL = '#d9d9d9'

plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 7.5,
    'axes.edgecolor': INK, 'axes.labelcolor': INK, 'text.color': INK,
    'xtick.color': INK, 'ytick.color': INK,
    'axes.linewidth': 0.7, 'xtick.major.width': 0.7, 'ytick.major.width': 0.7,
    'axes.spines.top': False, 'axes.spines.right': False,
    'figure.dpi': 120, 'savefig.bbox': 'tight',
    # Fixed salt, so SVG element ids are the same on every run.
    'svg.hashsalt': 'asic1a-picode',
})


def die(msg):
    sys.exit('FATAL build_figures.py: %s' % msg)


def need(p):
    if not os.path.exists(p):
        die('missing input %s; run the analysis scripts first' % p)
    return p


def load_ladder():
    return P.load_ladder()


def read_ph50(path, col='d_act'):
    rows = []
    for r in csv.DictReader(open(need(path))):
        if r['mutation'] == 'WT' or r.get(col) in (None, ''):
            continue
        rows.append(r)
    if not rows:
        die('no usable rows in %s' % path)
    return rows


def save(fig, name):
    os.makedirs(FIGS, exist_ok=True)
    place_panel_letters(fig)
    # No build date in the SVG, so a rebuild is byte-identical.
    for ext, kw in (('png', dict(dpi=400)),
                    ('svg', dict(metadata={'Date': None}))):
        fig.savefig(os.path.join(FIGS, '%s.%s' % (name, ext)), **kw)
    plt.close(fig)
    print('  wrote figures/%s.png and .svg' % name)


def place_labels(ax, xs, ys, labels, fontsize=6, dy=7):
    """Point labels nudged apart until no two overlap."""
    fig = ax.figure
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    ts = [ax.annotate(l, (x, y), textcoords='offset points', xytext=(0, dy),
                      ha='center', fontsize=fontsize, color=ACC)
          for x, y, l in zip(xs, ys, labels)]
    for _ in range(120):
        fig.canvas.draw()
        pts = [ax.transData.transform((x, y)) for x, y in zip(xs, ys)]
        # A 2 px cushion, so labels separate instead of settling flush
        # against one another.
        boxes = [Bbox.from_extents(b.x0 - 2, b.y0 - 2, b.x1 + 2, b.y1 + 2)
                 for b in (q.get_window_extent(rend) for q in ts)]
        moved = False
        for i in range(len(ts)):
            for j in range(i + 1, len(ts)):
                if not boxes[i].overlaps(boxes[j]):
                    continue
                ci = ((boxes[i].x0 + boxes[i].x1) / 2.0, (boxes[i].y0 + boxes[i].y1) / 2.0)
                cj = ((boxes[j].x0 + boxes[j].x1) / 2.0, (boxes[j].y0 + boxes[j].y1) / 2.0)
                vx, vy = ci[0] - cj[0], ci[1] - cj[1]
                if vx == 0 and vy == 0:
                    vx, vy = 1.0, 0.6
                m = (vx * vx + vy * vy) ** 0.5
                sx, sy = 1.4 * vx / m, 1.4 * vy / m
                ox, oy = ts[i].xyann; ts[i].xyann = (ox + sx, oy + sy)
                ox, oy = ts[j].xyann; ts[j].xyann = (ox - sx, oy - sy)
                moved = True
        # keep every label clear of every marker, not only its own
        for k, q in enumerate(ts):
            b = q.get_window_extent(rend)
            for px, py in pts:
                if b.x0 - 3 <= px <= b.x1 + 3 and b.y0 - 3 <= py <= b.y1 + 3:
                    ox, oy = q.xyann; q.xyann = (ox, oy + 1.5)
                    moved = True
                    break
        if not moved:
            break
    return ts


def panel(ax, letter, stat=None, center=False):
    # Queue the letter; save() places it once the layout pass has settled the
    # axes positions. Placing it here would freeze a position that
    # tight_layout then moves, which drops the letter onto panel content.
    fig = ax.figure
    if not hasattr(fig, 'panel_letters'):
        fig.panel_letters = []
    fig.panel_letters.append((ax, letter))
    if stat:
        ax.set_title(stat, loc='center' if center else 'right',
                     fontsize=8 if center else 7, color=INK, pad=6)


def place_panel_letters(fig):
    # Outside the panel's own left edge, so the letter clears the axis label
    # and the tick labels rather than starting inside the plot. The gridspec
    # cell is checked as well because an equal aspect ratio shrinks the axes
    # box and centres it inside its cell, which would carry the letter
    # inwards with it.
    letters = getattr(fig, 'panel_letters', ())
    if not letters:
        return
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    dy = 9.0 / 72.0 / fig.get_figheight()
    placed = []
    for ax, letter in letters:
        spec = getattr(ax, 'get_subplotspec', lambda: None)()
        box = ax.get_position()
        tight = inv.transform_bbox(ax.get_tightbbox(rend))
        left = min(tight.x0, box.x0)
        col = None
        if spec is not None:
            cell = spec.get_position(fig)
            left = min(left, cell.x0)
            col = round(cell.x0, 3)
        placed.append([letter, max(left, 0.006), box.y1 + dy, col])
    # One left edge per column, so letters stack flush whether or not a given
    # panel carries an axis label.
    for c in {p[3] for p in placed if p[3] is not None}:
        x = min(p[1] for p in placed if p[3] == c)
        for p in placed:
            if p[3] == c:
                p[1] = x
    for letter, x, y, _ in placed:
        fig.text(x, y, letter, fontsize=9, fontweight='bold', color=INK,
                 ha='left', va='bottom')


# structures
OXY = {'ASP': ('OD1', 'OD2'), 'GLU': ('OE1', 'OE2')}


def load_pdb(path):
    at = {}
    for l in open(need(path), errors='replace'):
        if l[:4] != 'ATOM':
            continue
        key = (l[21], int(l[22:26]))
        at.setdefault(key, [l[17:20].strip(), {}])[1][l[12:16].strip()] = (
            float(l[30:38]), float(l[38:46]), float(l[46:54]))
    if not at:
        die('no ATOM records in %s' % path)
    return at


def carboxyl(at, key):
    e = at.get(key)
    if not e:
        return []
    return [e[1][o] for o in OXY.get(e[0], ()) if o in e[1]]


def mind(A, B):
    return min(math.dist(a, b) for a in A for b in B) if A and B else float('nan')


def ca_trace(at, chain):
    pts = [(ri, v[1]['CA']) for (c, ri), v in sorted(at.items())
           if c == chain and 'CA' in v[1]]
    return np.array([p for _, p in pts]), [ri for ri, _ in pts]


def principal_plane(P):
    """Project onto the two principal axes, so the trimer is seen down its axis."""
    C = P - P.mean(axis=0)
    _, _, Vt = np.linalg.svd(C, full_matrices=False)
    return C @ Vt[:2].T, Vt


# Figure 1
def figure1(L):
    rows = read_ph50(os.path.join(DATA, 'asic1a_ph50.csv'))
    can = read_ph50(os.path.join(DATA, 'asic1a_ph50_canessa.csv'))
    seen = defaultdict(int)
    for r in rows:
        seen[(r['wt_aa'], r['mut_aa'])] += 1
    seen2 = defaultdict(int)
    for r in can:
        seen2[(r['wt_aa'], r['mut_aa'])] += 1

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.9))
    tri = [('a', 'NAY triangle, third base U or C', ['D', 'H', 'N']),
           ('b', 'NAR triangle, third base A or G', ['E', 'Q', 'K'])]
    for ax, (tag, title, res) in zip(axes, tri):
        panel(ax, tag, title, center=True)
        ang = {r: a for r, a in zip(res, [90, 210, 330])}
        pos = {r: (math.cos(math.radians(a)), math.sin(math.radians(a)))
               for r, a in ang.items()}
        for a in res:
            for b in res:
                if a == b:
                    continue
                n1, n2 = seen.get((a, b), 0), seen2.get((a, b), 0)
                tot = n1 + n2
                x1, y1 = pos[a]
                x2, y2 = pos[b]
                dx, dy = x2 - x1, y2 - y1
                nx, ny = -dy, dx
                nrm = math.hypot(nx, ny) or 1
                off = 0.085
                sx, sy = x1 + nx / nrm * off, y1 + ny / nrm * off
                ex, ey = x2 + nx / nrm * off, y2 + ny / nrm * off
                sh = 0.22
                sx, sy = sx + dx * sh, sy + dy * sh
                ex, ey = ex - dx * sh, ey - dy * sh
                col = LAD if tot else GREY
                ax.add_patch(FancyArrowPatch(
                    (sx, sy), (ex, ey), arrowstyle='-|>', mutation_scale=7,
                    lw=1.9 if tot else 0.8, color=col,
                    linestyle='-' if tot else (0, (2, 2)), zorder=2))
                mx, my = ((sx + ex) / 2 + nx / nrm * 0.17,
                          (sy + ey) / 2 + ny / nrm * 0.17)
                lbl = ('n=%d' % n1) if n1 else ('n=%d*' % n2 if n2 else 'none')
                ax.text(mx, my, lbl, ha='center', va='center', fontsize=6,
                        color=col if tot else GREY, zorder=3)
        for r in res:
            x, y = pos[r]
            ax.add_patch(Circle((x, y), 0.17, facecolor='white',
                                edgecolor=INK, lw=1.0, zorder=4))
            ax.text(x, y, r, ha='center', va='center', fontsize=11,
                    weight='bold', color=INK, zorder=5)
        ax.set_xlim(-1.62, 1.62); ax.set_ylim(-1.30, 1.42)
        ax.set_aspect('equal'); ax.axis('off')
    fig.subplots_adjust(bottom=0.06, top=0.90, wspace=0.02)
    save(fig, 'Figure1')


# Figure 2
def figure2(L, costs):
    rows = read_ph50(os.path.join(DATA, 'asic1a_ph50.csv'))
    per = defaultdict(list)
    for r in rows:
        p = (r['wt_aa'], r['mut_aa'])
        if p in set(L.CORE):
            per[p].append(abs(float(r['d_act'])))
    edges = [p for p in L.CORE if p in per]
    y = np.array([np.mean(per[p]) for p in edges])
    n = np.array([len(per[p]) for p in edges])
    cl = np.array([costs['cliff'][p] for p in edges])
    hu = np.array([costs['hums'][p] for p in edges])
    dq = np.array([abs(L.CHARGE[p[1]] - L.CHARGE[p[0]]) for p in edges])

    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.5))
    panels = [(cl, "AlphaMissense Cliff delta", 'a'),
              (hu, 'humsavar pathogenic fraction (%)', 'b'),
              (dq, '|charge change| (units)', 'c')]
    for ax, (x, xl, tag) in zip(axes, panels):
        ax.scatter(x, y, s=18 + 7 * n, facecolor=LAD, edgecolor=INK,
                   linewidth=0.5, alpha=0.9, zorder=3)
        if len(x) > 2:
            b = np.polyfit(x, y, 1)
            xs = np.linspace(x.min(), x.max(), 50)
            ax.plot(xs, np.polyval(b, xs), color=GREY, lw=0.9, ls='--', zorder=1)
        rho = spearman(list(x), list(y))
        ax.set_xlabel(xl, fontsize=7)
        panel(ax, tag, 'rho = %+.3f' % rho)
        ax.set_ylim(-0.05, max(y) * 1.32)
    axes[0].set_ylabel('mean |delta activation pH50|', fontsize=7)
    fig.tight_layout()
    for ax, (x, _xl, _tag) in zip(axes, panels):
        place_labels(ax, x, y, ['%s>%s' % p for p in edges])
    save(fig, 'Figure2')


def spearman(x, y):
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
    n = len(x)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return num / (dx * dy) if dx and dy else float('nan')


# Figure 3
def figure3(L):
    rows = read_ph50(os.path.join(DATA, 'asic1a_ph50.csv'))
    # The twelve core edges, matching the primary analysis.
    edges = set(L.CORE)
    CH = dict(L.CHARGE); CH.update({a: 0 for a in 'AVFMGLISTCPW'})
    byp = defaultdict(list)
    for r in rows:
        byp[(r['wt_aa'], int(r['resi']))].append(r)
    pairs = []
    for (w, p), v in byp.items():
        for r in v:
            if (w, r['mut_aa']) not in edges:
                continue
            ql = CH[r['mut_aa']] - CH[w]
            for s in v:
                if s is r or (w, s['mut_aa']) in edges:
                    continue
                if abs((CH[s['mut_aa']] - CH[w]) - ql) > 1e-9:
                    continue
                pairs.append((p, r['mutation'], abs(float(r['d_act'])),
                              s['mutation'], abs(float(s['d_act']))))
    if not pairs:
        die('no charge-matched pairs built for Figure 3')
    pairs.sort(key=lambda x: x[2] - x[4])
    d = np.array([x[2] - x[4] for x in pairs])

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2),
                             gridspec_kw={'width_ratios': [1.55, 1]})
    ax = axes[0]
    yy = np.arange(len(pairs))
    for i, (p, lm, lv, sm, sv) in enumerate(pairs):
        ax.plot([sv, lv], [i, i], color=GREY, lw=0.7, zorder=1)
    ax.scatter([x[4] for x in pairs], yy, s=13, facecolor=NON,
               edgecolor=INK, lw=0.4, label='non-ladder, same charge change', zorder=3)
    ax.scatter([x[2] for x in pairs], yy, s=13, facecolor=LAD,
               edgecolor=INK, lw=0.4, label='pI-code edge', zorder=3)
    ax.set_yticks(yy)
    ax.set_yticklabels(['%s / %s' % (x[1], x[3]) for x in pairs], fontsize=5.4)
    ax.set_xlabel('|delta activation pH50|', fontsize=7)
    ax.legend(frameon=False, fontsize=6.3, loc='upper center',
              bbox_to_anchor=(0.5, -0.12), ncol=2)
    panel(ax, 'a')
    ax.set_ylim(-1, len(pairs))

    ax = axes[1]
    ax.axvline(0, color=INK, lw=0.7, zorder=1)
    ax.hist(d, bins=np.arange(-0.32, 0.33, 0.04), color=FILL,
            edgecolor=INK, linewidth=0.6, zorder=2)
    ax.axvline(d.mean(), color=NON, lw=1.4, zorder=3)
    ax.annotate('mean %+.3f' % d.mean(), (d.mean(), ax.get_ylim()[1] * 0.92),
                textcoords='offset points', xytext=(5, 0), fontsize=6.5, color=NON)
    ax.set_xlabel('ladder minus non-ladder, |delta pH50|', fontsize=7)
    ax.set_ylabel('pairs', fontsize=7)
    panel(ax, 'b')
    fig.tight_layout()
    save(fig, 'Figure3')



# Figure 4, structure
def figure4():
    """Where a substitution sits predicts its effect better than what it is."""
    import picode_paths as _P
    lpo = _P.load_lpo()
    rows = read_ph50(os.path.join(DATA, 'asic1a_ph50.csv'))
    at = load_pdb(os.path.join(STRUCT, '9E4F.pdb1'))
    chains = sorted({c for c, _ in at})

    pts = []
    for a_, b_ in POCKET_H:
        for rr in (a_, b_):
            o = carboxyl(at, ('A', rr))
            if o:
                pts.extend(o)
    centre = np.mean(np.array(pts), axis=0)

    from Bio.PDB import PDBParser
    from Bio.PDB.SASA import ShrakeRupley
    import io as _io
    txt = ''.join(l for l in open(os.path.join(STRUCT, '9E4F.pdb1'),
                                  errors='replace')
                  if l[:6] in ('ATOM  ', 'HETATM', 'TER   ') and l[17:20] != 'HOH')
    mdl = PDBParser(QUIET=True).get_structure('x', _io.StringIO(txt))[0]
    ShrakeRupley().compute(mdl, level='R')
    sasa = {(c.id, r.id[1]): r.sasa for c in mdl for r in c if r.id[0] == ' '}

    per = defaultdict(list)
    loc = {}
    for r in rows:
        per[int(r['resi'])].append(abs(float(r['d_act'])))
        loc[int(r['resi'])] = r.get('location', '')
    D, S, B, Lb, P_ = [], [], [], [], []
    for ri, v in sorted(per.items()):
        res = at.get(('A', ri))
        if not res or ('A', ri) not in sasa:
            continue
        side = [x for k, x in res[1].items() if k not in ('N', 'CA', 'C', 'O')]
        if not side:
            continue
        aa = {'ASP': 'D', 'GLU': 'E', 'HIS': 'H', 'LYS': 'K'}.get(res[0])
        if aa is None or aa not in lpo.MAXASA:
            continue
        D.append(float(np.linalg.norm(np.mean(np.array(side), axis=0) - centre)))
        B.append(100.0 * sasa[('A', ri)] / lpo.MAXASA[aa])
        S.append(max(v)); Lb.append(loc[ri]); P_.append(ri)
    if len(D) < 10:
        die('too few positions placed in the structure for Figure 4')

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.2))

    # (a) measured shift painted onto the trimer
    ax = axes[0, 0]
    B2, ctr = c3_axis(at, chains)
    proj = lambda q: (np.asarray(q) - ctr) @ B2.T
    # Different weights per subunit; a raw CA trace down the three-fold axis
    # is otherwise unreadable.
    shades = {chains[0]: ('#7fa8cc', 0.95, 3), chains[1]: ('#dde8f2', 0.45, 1),
              chains[2]: ('#dde8f2', 0.45, 1)}
    for c in chains:
        Q = proj(ca_trace(at, c)[0])
        col, lw, z = shades[c]
        # Smooths the trace only; markers below use unsmoothed coordinates.
        if len(Q) >= 5:
            ker = np.ones(5) / 5.0
            Q = np.column_stack([np.convolve(Q[:, 0], ker, mode='valid'),
                                 np.convolve(Q[:, 1], ker, mode='valid')])
        ax.plot(Q[:, 0], Q[:, 1], color=col, lw=lw, zorder=z, alpha=0.9,
                solid_capstyle='round')
    smax = max(S)
    cm = plt.get_cmap('YlOrRd')
    for ri, sv in zip(P_, S):
        res = at.get(('A', ri))
        side = [x for k, x in res[1].items() if k not in ('N', 'CA', 'C', 'O')]
        q = proj(np.mean(np.array(side), axis=0))
        ax.scatter([q[0]], [q[1]], s=16 + 90 * (sv / smax), facecolor=cm(sv / smax),
                   edgecolor=INK, lw=0.4, zorder=4)
    for a_, b_ in POCKET_H:
        for rr in (a_, b_):
            o = carboxyl(at, ('A', rr))
            if o:
                q = proj(np.mean(o, axis=0))
                ax.scatter([q[0]], [q[1]], s=9, marker='x', color=LAD, lw=0.8,
                           zorder=5)
    ax.set_aspect('equal'); ax.axis('off')
    panel(ax, 'a')
    # (b) burial against shift
    ax = axes[0, 1]
    ax.scatter(B, S, s=20, facecolor=LAD, edgecolor=INK, lw=0.4, zorder=3)
    b = np.polyfit(B, S, 1)
    xs = np.linspace(min(B), max(B), 40)
    ax.plot(xs, np.polyval(b, xs), color=GREY, lw=0.9, ls='--', zorder=1)
    ax.set_xlabel('relative solvent accessibility (%)', fontsize=7)
    ax.set_ylabel('largest |delta activation pH50|', fontsize=7)
    panel(ax, 'b', 'rho = %+.3f' % spearman(B, S))

    # (c) distance to the pocket against shift
    ax = axes[1, 0]
    ax.scatter(D, S, s=20, facecolor=NON, edgecolor=INK, lw=0.4, zorder=3)
    b = np.polyfit(D, S, 1)
    xs = np.linspace(min(D), max(D), 40)
    ax.plot(xs, np.polyval(b, xs), color=GREY, lw=0.9, ls='--', zorder=1)
    ax.set_xlabel('distance to the acidic pocket centre (A)', fontsize=7)
    ax.set_ylabel('largest |delta activation pH50|', fontsize=7)
    panel(ax, 'c', 'rho = %+.3f' % spearman(D, S))

    # (d) by domain, as the source labels them
    ax = axes[1, 1]
    byd = defaultdict(list)
    for r in rows:
        byd[r.get('location', '')].append(abs(float(r['d_act'])))
    if '' in byd:
        die('%d substitutions carry no domain label; rebuild data/asic1a_ph50.csv'
            % len(byd['']))
    order = sorted(byd, key=lambda k: -np.mean(byd[k]))
    xs = np.arange(len(order))
    # Headroom for the n= labels.
    ytop = max(max(v) for v in byd.values()) * 1.16
    ax.set_ylim(0, ytop)
    ax.bar(xs, [np.mean(byd[k]) for k in order], 0.62, color=FILL,
           edgecolor=INK, lw=0.6, zorder=2)
    for i, k in enumerate(order):
        ax.scatter(np.full(len(byd[k]), i) + np.random.default_rng(7).normal(
            0, 0.07, len(byd[k])), byd[k], s=5, facecolor=ACC, alpha=0.55,
            edgecolor='none', zorder=3)
        ax.text(i, ytop * 0.955, 'n=%d' % len(byd[k]), ha='center', va='top',
                fontsize=6.0, color=INK)
    ax.set_xticks(xs)
    ax.set_xticklabels([k.replace('Finger and thumb', 'finger,\nthumb')
                        .replace('Acidic pocket', 'acidic\npocket')
                        .replace('Transmembrane', 'trans-\nmembrane')
                        for k in order], fontsize=6)
    ax.set_ylabel('mean |delta activation pH50|', fontsize=7)
    panel(ax, 'd')
    fig.tight_layout(h_pad=2.0, w_pad=1.8)
    save(fig, 'Figure4')


# Figure 5
POCKET_H = [(237, 351), (238, 347), (219, 409), (79, 418)]
POCKET_C = [(238, 350), (239, 346), (220, 408), (80, 417)]
CLUSTER_H = [237, 238, 347, 351]
CLUSTER_C = [238, 239, 346, 350]
PAIRNAME = ['D237-D351 / D238-D350', 'E238-D347 / E239-D346',
            'E219-D409 / E220-D408', 'E79-E418 / E80-E417']


def c3_axis(at, chains):
    """The three-fold axis of a C3 trimer: the normal to the plane of the three
    chain centroids. Projecting onto the perpendicular plane views it down the axis."""
    cen = []
    for c in chains:
        P, _ = ca_trace(at, c)
        cen.append(P.mean(axis=0))
    cen = np.array(cen)
    if len(cen) != 3:
        die('expected 3 chains to define the three-fold axis, found %d' % len(cen))
    n = np.cross(cen[1] - cen[0], cen[2] - cen[0])
    nn = np.linalg.norm(n)
    if nn < 1e-6:
        die('the three chain centroids are collinear; no three-fold axis')
    n = n / nn
    e1 = cen[0] - cen.mean(axis=0)
    e1 = e1 - n * (e1 @ n)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(n, e1)
    return np.vstack([e1, e2]), cen.mean(axis=0)


def figure5():
    files = [('9E4A.pdb1', 8.5, 'human'), ('9E4F.pdb1', 5.7, 'human'),
             ('6VTL.pdb1', 8.0, 'chicken'), ('2QTS.pdb1', 6.5, 'chicken')]
    loaded = {f: load_pdb(os.path.join(STRUCT, f)) for f, _, _ in files}

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.4))

    # (a) trimer down the true three-fold
    ax = axes[0, 0]
    at = loaded['9E4F.pdb1']
    chains = sorted({c for c, _ in at})
    B, ctr = c3_axis(at, chains)
    proj = lambda p: (np.asarray(p) - ctr) @ B.T
    for c, col in zip(chains, ['#6baed6', '#9ecae1', '#c6dbef']):
        P, _ = ca_trace(at, c)
        Q = proj(P)
        ax.plot(Q[:, 0], Q[:, 1], color=col, lw=0.5, zorder=1)
    for c in chains:
        for a_, b_ in POCKET_H:
            for rr in (a_, b_):
                o = carboxyl(at, (c, rr))
                if o:
                    q = proj(np.mean(o, axis=0))
                    ax.scatter([q[0]], [q[1]], s=13, facecolor=NON,
                               edgecolor=INK, lw=0.35, zorder=4)
    ax.set_aspect('equal'); ax.axis('off')
    panel(ax, 'a')
    # (b) every pair is intra-subunit
    ax = axes[0, 1]
    labs, intra, inter = [], [], []
    for f, ph, sp in files:
        a2 = loaded[f]
        chs = sorted({c for c, _ in a2})
        prs = POCKET_H if sp == 'human' else POCKET_C
        for i, (r1, r2) in enumerate(prs):
            ia, ie = [], []
            for ca in chs:
                for cb in chs:
                    d = mind(carboxyl(a2, (ca, r1)), carboxyl(a2, (cb, r2)))
                    if d == d:
                        (ia if ca == cb else ie).append(d)
            if ia and ie:
                labs.append('%s %s' % (f.split('.')[0],
                                       PAIRNAME[i].split(' / ')[0 if sp == 'human' else 1]))
                intra.append(min(ia)); inter.append(min(ie))
    y = np.arange(len(labs))
    for i in range(len(labs)):
        ax.plot([intra[i], inter[i]], [i, i], color=GREY, lw=0.6, zorder=1)
    ax.scatter(intra, y, s=14, facecolor=LAD, edgecolor=INK, lw=0.4,
               label='within one subunit', zorder=3)
    ax.scatter(inter, y, s=14, facecolor=NON, edgecolor=INK, lw=0.4,
               label='closest across subunits', zorder=3)
    ax.set_xscale('log')
    ax.set_xticks([2, 3, 5, 10, 20, 30, 50])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_xlim(2.0, 70)
    ax.set_yticks(y); ax.set_yticklabels(labs, fontsize=5.0)
    ax.set_xlabel('min carboxyl O to O (A, log scale)', fontsize=7)
    ax.legend(frameon=False, fontsize=6.2, loc='upper center',
              bbox_to_anchor=(0.5, -0.17), ncol=2)
    panel(ax, 'b')

    # (c) the palm pair belongs to the neighbouring subunit
    ax = axes[1, 0]
    rows, same, adj = [], [], []
    for f, ph, sp in files:
        a2 = loaded[f]
        chs = sorted({c for c, _ in a2})
        pal = (219, 409) if sp == 'human' else (220, 408)
        clu = CLUSTER_H if sp == 'human' else CLUSTER_C
        for rr in pal:
            s_, d_ = [], []
            for ca in chs:
                A = carboxyl(a2, (ca, rr))
                for cb in chs:
                    for q in clu:
                        dd = mind(A, carboxyl(a2, (cb, q)))
                        if dd == dd:
                            (s_ if ca == cb else d_).append(dd)
            if s_ and d_:
                rows.append('%s %s%d' % (f.split('.')[0],
                                         'E' if rr in (219, 220) else 'D', rr))
                same.append(min(s_)); adj.append(min(d_))
    y = np.arange(len(rows))
    ax.barh(y - 0.2, same, 0.4, color=GREY, edgecolor=INK, lw=0.5,
            label='to its own subunit cluster')
    ax.barh(y + 0.2, adj, 0.4, color=NON, edgecolor=INK, lw=0.5,
            label='to the adjacent subunit cluster')
    ax.set_yticks(y); ax.set_yticklabels(rows, fontsize=5.4)
    ax.set_xlabel('min carboxyl O to O (A)', fontsize=7)
    ax.legend(frameon=False, fontsize=6.2, loc='upper center',
              bbox_to_anchor=(0.5, -0.17), ncol=2)
    panel(ax, 'c')

    # (d) collapse on acidification
    ax = axes[1, 1]
    marks = {'human': 'o', 'chicken': 's'}
    for idx, colr, nm in ((0, NON, PAIRNAME[0]), (1, LAD, PAIRNAME[1])):
        for f, ph, sp in files:
            a2 = loaded[f]
            pr = (POCKET_H if sp == 'human' else POCKET_C)[idx]
            vals = [mind(carboxyl(a2, (c, pr[0])), carboxyl(a2, (c, pr[1])))
                    for c in sorted({c for c, _ in a2})]
            vals = [v for v in vals if v == v]
            if vals:
                ax.scatter([ph], [min(vals)], s=28, marker=marks[sp],
                           facecolor=colr, edgecolor=INK, lw=0.5, zorder=3)
    ax.scatter([], [], s=28, marker='o', facecolor=NON, edgecolor=INK,
               label=PAIRNAME[0])
    ax.scatter([], [], s=28, marker='o', facecolor=LAD, edgecolor=INK,
               label=PAIRNAME[1])
    ax.scatter([], [], s=28, marker='s', facecolor='white', edgecolor=INK,
               label='square, chicken')
    ax.set_xlabel('structure pH', fontsize=7)
    ax.set_ylabel('min carboxyl O to O (A)', fontsize=7)
    ax.legend(frameon=False, fontsize=6.2, loc='upper left')
    panel(ax, 'd')

    fig.tight_layout(h_pad=3.0, w_pad=1.8)
    save(fig, 'Figure5')


# Figure 6
def figure6():
    pka = os.path.join(STRUCT, '9E4F.pka')
    rest = os.path.join(STRUCT, '9E4A.pka')
    lo = parse_pka_summary(need(pka))
    hi = parse_pka_summary(need(rest))
    order = [79, 219, 237, 238, 242, 347, 351, 409, 418]
    names = {79: 'E79', 219: 'E219', 237: 'D237', 238: 'E238', 242: 'E242',
             347: 'D347', 351: 'D351', 409: 'D409', 418: 'E418'}
    MODEL = {'ASP': 3.80, 'GLU': 4.50}

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9),
                             gridspec_kw={'width_ratios': [1.25, 1]})
    ax = axes[0]
    xs = np.arange(len(order))
    L_ = [lo[(r, 'A')][1] for r in order]
    H_ = [hi[(r, 'A')][1] for r in order]
    M_ = [MODEL[lo[(r, 'A')][0]] for r in order]
    ax.bar(xs - 0.2, np.array(L_) - np.array(M_), 0.4, color=NON,
           edgecolor=INK, lw=0.5, label='pH 5.7 (9E4F)')
    ax.bar(xs + 0.2, np.array(H_) - np.array(M_), 0.4, color=LAD,
           edgecolor=INK, lw=0.5, label='pH 8.5 (9E4A)')
    ax.axhline(0, color=INK, lw=0.7)
    ax.set_xticks(xs); ax.set_xticklabels([names[r] for r in order], fontsize=6.5)
    ax.set_ylabel('predicted pKa minus model value', fontsize=7)
    ax.legend(frameon=False, fontsize=6.3)
    panel(ax, 'a')

    ax = axes[1]
    det = parse_pka_determinants(need(pka))
    tot = defaultdict(float)
    for r in order:
        for val, prn, pri, pch in det.get((r, 'A'), []):
            if pri in order:
                continue
            tot[(pri, prn)] += abs(val)
    top = sorted(tot.items(), key=lambda kv: -kv[1])[:8][::-1]
    lbl = ['%s%d' % ({'ARG': 'R', 'GLU': 'E', 'ASP': 'D', 'LYS': 'K',
                      'HIS': 'H', 'GLN': 'Q', 'ASN': 'N'}.get(k[1], k[1]), k[0])
           for k, _ in top]
    ax.barh(range(len(top)), [v for _, v in top], color=FILL,
            edgecolor=INK, lw=0.6)
    ax.set_yticks(range(len(top))); ax.set_yticklabels(lbl, fontsize=6.5)
    ax.set_xlabel('summed |contribution| to pocket pKa', fontsize=7)
    panel(ax, 'b')
    fig.tight_layout()
    save(fig, 'Figure6')


def parse_pka_summary(path):
    import re
    lines = open(path).read().splitlines()
    try:
        s = next(i for i, l in enumerate(lines)
                 if l.startswith('SUMMARY OF THIS PREDICTION'))
    except StopIteration:
        die('%s has no SUMMARY block' % path)
    out = {}
    for l in lines[s + 1:]:
        m = re.match(r'\s*([A-Z]{3})\s+(\d+)\s+([A-Z])\s+(-?\d+\.\d+)', l)
        if m:
            out[(int(m.group(2)), m.group(3))] = (m.group(1), float(m.group(4)))
    if not out:
        die('parsed no pKa values from %s' % path)
    return out


def parse_pka_determinants(path):
    import re
    DET = re.compile(r'^(?P<rn>[A-Z]{3})\s+(?P<ri>\d+)\s+(?P<ch>[A-Z])\s+(?P<rest>.*)$')
    PAR = re.compile(r'(-?\d+\.\d+)\s+([A-Z]{3})\s+(\d+)\s+([A-Z])')
    lines = open(path).read().splitlines()
    stop = next((i for i, l in enumerate(lines)
                 if l.startswith('SUMMARY OF THIS PREDICTION')), len(lines))
    det = defaultdict(list)
    for l in lines[:stop]:
        m = DET.match(l.rstrip())
        if not m:
            continue
        key = (int(m.group('ri')), m.group('ch'))
        for val, prn, pri, pch in PAR.findall(m.group('rest')):
            if prn == 'XXX' or (int(pri), pch) == key:
                continue
            det[key].append((float(val), prn, int(pri), pch))
    if not det:
        die('parsed no determinants from %s' % path)
    return det


# supplementary figures
def uniprot_transmem(name):
    """TRANSMEM spans and sequence length from a UniProt flat file.

    Read from data/sources so the shading in Figure S2 traces to a file in the
    repository. A missing or featureless entry is fatal: an unshaded figure
    would silently assert the protein has no transmembrane segment.
    """
    p = need(os.path.join(DATA, 'sources', name))
    spans, ln = [], None
    for l in open(p, errors='replace'):
        if l.startswith('FT   TRANSMEM'):
            f = l.split()
            if len(f) < 3 or '..' not in f[2]:
                die('cannot read a TRANSMEM range from %r in %s' % (l.strip(), p))
            a, b = f[2].split('..')
            spans.append((int(a.lstrip('<')), int(b.lstrip('>'))))
        elif l.startswith('SQ   SEQUENCE'):
            ln = int(l.split()[2])
    if not spans:
        die('%s carries no TRANSMEM feature; ASIC1a has two' % p)
    if ln is None:
        die('%s carries no SQ line, so the sequence length is unknown' % p)
    return spans, ln


def figureS1():
    """The two readouts are coupled, so they are not two independent tests."""
    rows = read_ph50(os.path.join(DATA, 'asic1a_ph50.csv'))
    xy = [(float(r['d_act']), float(r['d_ssd'])) for r in rows
          if r.get('d_ssd') not in (None, '')]
    if len(xy) < 20:
        die('only %d substitutions carry both readouts' % len(xy))
    X = np.array([a for a, _ in xy]); Y = np.array([b for _, b in xy])
    same = int(np.sum(X * Y > 0)); nz = int(np.sum(X * Y != 0))

    fig, ax = plt.subplots(figsize=(4.4, 4.0))
    lim = max(abs(X).max(), abs(Y).max()) * 1.08
    ax.axhline(0, color=GREY, lw=0.6, zorder=1)
    ax.axvline(0, color=GREY, lw=0.6, zorder=1)
    # The quadrants carry the claim, so they are shaded rather than described.
    for sx, sy in ((1, 1), (-1, -1)):
        ax.add_patch(plt.Rectangle((0 if sx > 0 else -lim, 0 if sy > 0 else -lim),
                                   lim, lim, facecolor=LAD, alpha=0.055,
                                   edgecolor='none', zorder=0))
    ax.scatter(X, Y, s=16, facecolor=LAD, edgecolor=INK, lw=0.3, alpha=0.8,
               zorder=3)
    b, a = np.polyfit(X, Y, 1)
    xs = np.linspace(-lim, lim, 2)
    ax.plot(xs, b * xs + a, '--', color=GREY, lw=0.9, zorder=2)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_aspect('equal')
    ax.set_xlabel('shift in activation pH50')
    ax.set_ylabel('shift in steady-state desensitization pHD50')
    ax.set_title('rho = %+.3f   %d of %d move the same way'
                 % (spearman(list(X), list(Y)), same, nz), fontsize=7.5,
                 loc='left')
    fig.tight_layout()
    save(fig, 'FigureS1')


def figureS2():
    """Where along the sequence the scan looked, and what it found."""
    rows = read_ph50(os.path.join(DATA, 'asic1a_ph50.csv'))
    byloc = defaultdict(list)
    for r in rows:
        byloc[r['location']].append((int(r['resi']), abs(float(r['d_act']))))
    if '' in byloc:
        die('%d substitutions carry no domain label; rebuild data/asic1a_ph50.csv'
            % len(byloc['']))
    order = sorted(byloc, key=lambda k: -np.mean([v for _, v in byloc[k]]))
    cm = plt.get_cmap('tab10')

    fig, ax = plt.subplots(figsize=(7.2, 3.4))
    # Transmembrane spans read from the UniProt flat file, not written in here.
    tm, seqlen = uniprot_transmem('up_P78348.txt')
    top = max(v for vs in byloc.values() for _, v in vs) * 1.16
    for lo, hi in tm:
        ax.axvspan(lo, hi, color=FILL, alpha=0.8, zorder=0)
        ax.text((lo + hi) / 2.0, top * 0.985, 'TM', ha='center', va='top',
                fontsize=6, color=ACC)
    for i, k in enumerate(order):
        xs = [p for p, _ in byloc[k]]
        ys = [v for _, v in byloc[k]]
        ax.scatter(xs, ys, s=20, facecolor=cm(i), edgecolor=INK, lw=0.3,
                   alpha=0.85, label='%s (n=%d)' % (k, len(xs)), zorder=3)
    npos = len({int(r['resi']) for r in rows})
    for p in sorted({int(r['resi']) for r in rows}):
        ax.plot([p, p], [0, -0.05 * top], color=ACC, lw=0.5, zorder=2)
    ax.set_xlim(0, seqlen + 12); ax.set_ylim(-0.10 * top, top)
    ax.set_xlabel('position in human ASIC1a (UniProt P78348, %d residues); ticks '
                  'below the axis mark the %d scanned positions' % (seqlen, npos))
    ax.set_ylabel('|shift in activation pH50|')
    # Legend under the axis; every band of the plot area carries data.
    ax.legend(fontsize=6, ncol=6, frameon=False, loc='upper center',
              bbox_to_anchor=(0.5, -0.20), handletextpad=0.3, columnspacing=1.2)
    fig.tight_layout()
    save(fig, 'FigureS2')


def main():
    L = load_ladder()
    costs = read_costs()
    print('building figures into %s' % FIGS)
    figure1(L)
    figure2(L, costs)
    figure3(L)
    figure4()
    figure5()
    figure6()
    figureS1()
    figureS2()
    print('done')


def read_costs():
    """Edge costs, from the results file written by the analysis step."""
    p = os.path.join(ROOT, 'results', 'edge_costs.csv')
    need(p)
    cliff, hums = {}, {}
    for r in csv.DictReader(open(p)):
        k = (r['source'], r['target'])
        cliff[k] = float(r['cliff_delta'])
        hums[k] = float(r['humsavar_pct'])
    if not cliff:
        die('no edge costs in %s' % p)
    return dict(cliff=cliff, hums=hums)


if __name__ == '__main__':
    main()
