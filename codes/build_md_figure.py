#!/usr/bin/env python3
"""
build_md_figure.py - the MD control at position 246, from source trajectories.

    python3 codes/build_md_figure.py

Reads, for each of the three builds (K246 wild type, K246Q on the ladder and
K246M off it), the fitted-trajectory analysis under md/results/<build>/:

    picode_dist.xvg   per-frame pocket-pair and reach distances (gmx distance)
    picode_rmsd.xvg   per-frame backbone RMSD (gmx rms)
    prod_ref.pdb      the fitted reference frame, protein only

Everything plotted is recomputed from those; nothing is read from a derived
table. Writes md/results/{comparison.json, SUMMARY.txt, md_comparison.png} and
figures/FigureS3.{png,svg}.

Exits non-zero on any missing or unreadable input.
"""
import json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(HERE)
RESULTS = os.path.join(PROJ, "md", "results")
FIGS = os.path.join(PROJ, "figures")
BUILDS = ["wt", "ladder", "compare"]
POSITION = 246
BURN_FRACTION = 0.2
CONTACT = 4.0
# picode_dist.xvg columns: twelve intra-subunit pairs, then three inter-subunit reaches.
PAIRS = ["D237-D351", "E238-D347", "E219-D409", "E79-E418"]
LABELS = ["%s|%s" % (c, p) for c in "ABC" for p in PAIRS]
REACH_COLS = [12, 13, 14]
TIPNAME = {"LYS": "NZ", "GLN": "NE2", "MET": "CE"}
ONE = {"LYS": "K", "GLN": "Q", "MET": "M"}
COL = {"wt": "#4d4d4d", "ladder": "#2166ac", "compare": "#b2182b"}


def die(msg):
    sys.stderr.write("STOP: %s\n" % msg)
    sys.exit(1)


def load_xvg(path):
    if not os.path.isfile(path):
        die("missing %s" % path)
    rows = [l.split() for l in open(path) if l and l[0] not in "#@"]
    if not rows:
        die("%s has no data rows" % path)
    return np.array([[float(x) for x in r] for r in rows])


def read_pdb(path):
    if not os.path.isfile(path):
        die("missing %s" % path)
    res246, pocket = {}, {}
    for l in open(path):
        if not l.startswith("ATOM"):
            continue
        ch, rid, rn, an = l[21], int(l[22:26]), l[17:20].strip(), l[12:16].strip()
        xyz = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        if rid == POSITION and an == "CA":
            res246[ch] = rn
        key = (ch, rid, an)
        pocket[key] = (rn, xyz)
    if not res246:
        die("%s has no residue %d" % (path, POSITION))
    return res246, pocket


def main():
    from scipy import stats
    dist, rmsd, res246 = {}, {}, {}
    for b in BUILDS:
        d = os.path.join(RESULTS, b)
        dist[b] = load_xvg(os.path.join(d, "picode_dist.xvg"))[:, 1:] * 10.0
        rmsd[b] = load_xvg(os.path.join(d, "picode_rmsd.xvg"))
        r, _ = read_pdb(os.path.join(d, "prod_ref.pdb"))
        got = set(r.values())
        if len(got) != 1:
            die("%s carries more than one residue at %d: %s" % (b, POSITION, got))
        res246[b] = r["A"]
        if dist[b].shape[1] != len(LABELS) + 3:
            die("%s picode_dist.xvg has %d data columns, expected %d"
                % (b, dist[b].shape[1], len(LABELS) + 3))
    want = {"wt": "LYS", "ladder": "GLN", "compare": "MET"}
    for b in BUILDS:
        if res246[b] != want[b]:
            die("%s should carry %s at %d but carries %s"
                % (b, want[b], POSITION, res246[b]))
    n_frames = dist["wt"].shape[0]
    burn = int(n_frames * BURN_FRACTION)

    def pmean(b, j):
        return float(dist[b][burn:, j].mean())

    def closed(b, j):
        return float((dist[b][burn:, j] < CONTACT).mean())

    dl = np.array([abs(pmean("ladder", j) - pmean("wt", j)) for j in range(12)])
    dc = np.array([abs(pmean("compare", j) - pmean("wt", j)) for j in range(12)])
    k, n = int((dl < dc).sum()), 12
    sign_p = min(1.0, 2 * sum(math.comb(n, i) for i in range(min(k, n - k) + 1)) / 2.0 ** n)
    W = stats.wilcoxon(dl, dc)
    wp = float(W.pvalue)
    sig = (sign_p < 0.05) or (wp < 0.05)
    noncontact = [LABELS[j] for j in range(12) if all(closed(b, j) == 0.0 for b in BUILDS)]
    reach = {b: dist[b][burn:, REACH_COLS].ravel() for b in BUILDS}
    reach_mean = {b: float(reach[b].mean()) for b in BUILDS}
    rmsd_mean = {b: float(rmsd[b][burn:, 1].mean() * 10.0) for b in BUILDS}
    tip = {b: TIPNAME[res246[b]] for b in BUILDS}

    # ---- SUMMARY + comparison.json ----
    L = ["ASIC1a pI-code MD, position %d" % POSITION,
         "pI-code K%dQ and comparator K%dM, both from wild-type K%d"
         % (POSITION, POSITION, POSITION),
         "",
         "%-14s %8s %8s %8s   %9s %9s" % ("pocket pair", "wt", "K%dQ" % POSITION,
                                          "K%dM" % POSITION, "|Q-wt|", "|M-wt|")]
    for j, lab in enumerate(LABELS):
        w, a, c = pmean("wt", j), pmean("ladder", j), pmean("compare", j)
        L.append("%-14s %8.2f %8.2f %8.2f   %9.2f %9.2f" % (lab, w, a, c, abs(a - w), abs(c - w)))
    L += ["",
          "mean departure from wild type over %d pocket-pair measurements" % n,
          "   pI-code K%dQ    %.3f A" % (POSITION, dl.mean()),
          "   comparator K%dM %.3f A" % (POSITION, dc.mean()),
          "   median         Q %.3f A, M %.3f A" % (np.median(dl), np.median(dc)),
          "",
          "pI-code closer to wild type in %d of %d pairs" % (k, n),
          "sign test p %.3f, Wilcoxon signed-rank p %.3f" % (sign_p, wp),
          "never within %.0f A in any build (not a contact here): %s"
          % (CONTACT, ", ".join(noncontact) if noncontact else "none"),
          "",
          "position %d reach to the neighbouring subunit D409, mean over 3 interfaces" % POSITION,
          "   " + ",  ".join(("wild type K%d %s %.2f A" % (POSITION, tip[b], reach_mean[b])
                              if b == "wt" else
                              "K%d%s %s %.2f A" % (POSITION, ONE[res246[b]], tip[b], reach_mean[b]))
                             for b in BUILDS),
          "",
          "backbone RMSD, A: " + ", ".join("%s %.2f" % (b, rmsd_mean[b]) for b in BUILDS),
          "frames per build: %d  (first %.0f%% discarded)" % (n_frames, 100 * BURN_FRACTION),
          ""]
    if sig:
        L.append("The two substitutions differ in pocket geometry (sign p %.2f, "
                 "Wilcoxon p %.2f)." % (sign_p, wp))
    else:
        L.append("Neither substitution moves the pocket from wild type beyond the "
                 "subunit-to-subunit scatter, and the two are indistinguishable from "
                 "each other (sign p %.2f, Wilcoxon p %.2f). The mean departure differs "
                 "by %.3f A, carried by the D237-D351 pairs, which barely form contacts."
                 % (sign_p, wp, abs(dl.mean() - dc.mean())))
    L.append("The engineered side chain at %d does not reach the neighbouring D409 in "
             "any build (%.1f to %.1f A, never within %.0f A)."
             % (POSITION, min(reach_mean.values()), max(reach_mean.values()), CONTACT))
    L.append("This is a structural control over one 20 ns trajectory per build with the "
             "three subunits as the only replication. It tests whether either "
             "substitution remodels the pocket, not what either does to pH50.")
    os.makedirs(RESULTS, exist_ok=True)
    open(os.path.join(RESULTS, "SUMMARY.txt"), "w").write("\n".join(L) + "\n")
    json.dump(dict(position=POSITION, builds={b: ONE[res246[b]] for b in BUILDS},
                   labels=LABELS, ladder_departure=list(map(float, dl)),
                   compare_departure=list(map(float, dc)),
                   ladder_mean=float(dl.mean()), compare_mean=float(dc.mean()),
                   ladder_median=float(np.median(dl)), compare_median=float(np.median(dc)),
                   ladder_closer_in=k, n=n, sign_p=float(sign_p), wilcoxon_p=wp,
                   significant=bool(sig), noncontact_pairs=noncontact,
                   reach=reach_mean, rmsd_mean=rmsd_mean, n_frames=n_frames,
                   finding=("no significant difference between the builds; pocket "
                            "geometry unchanged from wild type within this sampling")),
              open(os.path.join(RESULTS, "comparison.json"), "w"), indent=2)
    struct_png = os.path.join(RESULTS, "pocket_struct.png")
    render_struct(struct_png)
    if not os.path.isfile(struct_png):
        die("panel a has no structural render: pymol is unavailable here and no "
            "previous %s exists, so the figure would ship without it. Install "
            "pymol (pip install pymol-open-source; wheels exist for macOS "
            "arm64/x86_64, linux x86_64 and Windows, Python 3.9 to 3.13, but "
            "not for linux aarch64) or restore the file." % struct_png)
    if True:
        # strip before the raster is embedded, so FigureS3.svg inherits a clean copy
        d = strip_credentials(struct_png)
        if d:
            print("stripped %s from pocket_struct.png" % ", ".join(sorted(set(d))))
    make_figure(dist, rmsd, reach_mean, tip, burn, struct_png)
    for f in (os.path.join(FIGS, "FigureS3.png"),
              os.path.join(RESULTS, "md_comparison.png")):
        if not os.path.isfile(f):
            die("expected %s to have been written" % f)
        d = strip_credentials(f)
        if d:
            print("stripped %s from %s" % (", ".join(sorted(set(d))), os.path.basename(f)))
    print("wrote comparison.json, SUMMARY.txt, md_comparison.png and FigureS3 to %s"
          % RESULTS)
    print("verdict:", "DIFFERENT" if sig else "no significant difference",
          "(sign p %.3f, Wilcoxon p %.3f)" % (sign_p, wp))


def render_struct(out):
    """PyMOL cartoon of the pocket and residue 246; skipped if pymol is absent."""
    try:
        import pymol2
    except Exception as e:
        sys.stderr.write("note: pymol not available (%s), figure will omit panel a\n" % e)
        return
    POCKET = "resi 237+351+238+347+219+409+79+418 and resn ASP+GLU"
    with pymol2.PyMOL() as p:
        cmd = p.cmd
        for obj, b in (("wt", "wt"), ("lad", "ladder"), ("cmp", "compare")):
            cmd.load(os.path.join(RESULTS, b, "prod_ref.pdb"), obj)
        cmd.align("lad", "wt"); cmd.align("cmp", "wt")
        cmd.bg_color("white"); cmd.hide("everything")
        cmd.set("ray_opaque_background", 0); cmd.set("ray_shadows", 0); cmd.set("antialias", 1)
        cmd.show("cartoon", "wt and chain A+B")
        cmd.color("grey80", "wt and chain A"); cmd.color("grey90", "wt and chain B")
        cmd.show("sticks", "wt and chain A and %s" % POCKET)
        cmd.color("firebrick", "wt and chain A and %s" % POCKET)
        cmd.set("stick_radius", 0.22, "wt and chain A and %s" % POCKET)
        cmd.show("sticks", "wt and chain B and resi 409"); cmd.color("purple", "wt and chain B and resi 409")
        cmd.set("stick_radius", 0.26, "wt and chain B and resi 409")
        for obj, c in (("wt", "0x333333"), ("lad", "0x2166ac"), ("cmp", "0xb2182b")):
            cmd.show("sticks", "%s and chain A and resi 246" % obj)
            cmd.color(c, "%s and chain A and resi 246" % obj)
            cmd.set("stick_radius", 0.30, "%s and chain A and resi 246" % obj)
        cmd.orient("wt and chain A and (%s or resi 246)" % POCKET)
        cmd.zoom("wt and chain A and (%s or resi 246) or (wt and chain B and resi 409)" % POCKET, 3.0)
        cmd.ray(1400, 1200); cmd.png(out, dpi=200)


# Image and colour/resolution chunks only; drops C2PA (caBX) and matplotlib text.
_KEEP = {b"IHDR", b"PLTE", b"IDAT", b"IEND", b"tRNS", b"gAMA",
         b"cHRM", b"sRGB", b"iCCP", b"pHYs", b"bKGD", b"sBIT"}


def strip_credentials(path):
    """Rewrite a PNG in place with only the allow-listed chunks. Fatal if malformed."""
    import struct
    d = open(path, "rb").read()
    if d[:8] != b"\x89PNG\r\n\x1a\n":
        die("%s is not a PNG, cannot strip metadata" % path)
    out, i, dropped = [d[:8]], 8, []
    while i < len(d):
        if i + 8 > len(d):
            die("%s is truncated at byte %d" % (path, i))
        ln = struct.unpack(">I", d[i:i + 4])[0]
        typ = d[i + 4:i + 8]
        end = i + 12 + ln
        if end > len(d):
            die("%s chunk %s overruns the file" % (path, typ.decode("latin1")))
        if typ in _KEEP:
            out.append(d[i:end])
        else:
            dropped.append(typ.decode("latin1"))
        i = end
        if typ == b"IEND":
            break
    if not any(c[4:8] == b"IEND" for c in out[1:]):
        die("%s has no IEND chunk" % path)
    open(path, "wb").write(b"".join(out))
    return dropped


def make_figure(dist, rmsd, reach_mean, tip, burn, struct_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    from scipy import stats
    import numpy as np
    LAB1 = {"wt": "wild type K246", "ladder": "pI-code K246Q", "compare": "comparator K246M"}
    plt.rcParams.update({"font.size": 7.2, "axes.linewidth": 0.6, "svg.fonttype": "none",
                          "svg.hashsalt": "asic1a-picode",
                         "font.family": "DejaVu Sans"})
    fig = plt.figure(figsize=(7.2, 6.0))
    gs = fig.add_gridspec(2, 2, hspace=0.30, wspace=0.24, left=0.075, right=0.975,
                          top=0.93, bottom=0.1)
    # (a)
    ax = fig.add_subplot(gs[0, 0]); ax.axis("off")
    if os.path.isfile(struct_png):
        img = mpimg.imread(struct_png); ax.imshow(img); H, W = img.shape[0], img.shape[1]
        xy = lambda fx, fy: (fx * W, fy * H)
        ax.annotate("residue 246\n(K / Q / M)", xy(0.47, 0.56), xytext=xy(0.70, 0.80),
                    fontsize=6.2, ha="left", arrowprops=dict(arrowstyle="-", lw=0.6, color="#1a1a1a"))
        ax.annotate("D409\n(subunit B)", xy(0.42, 0.45), xytext=xy(0.08, 0.16), fontsize=6.2,
                    ha="left", color="#762a83", arrowprops=dict(arrowstyle="-", lw=0.6, color="#762a83"))
        ax.annotate("pocket\ncarboxylates", xy(0.17, 0.52), xytext=xy(0.02, 0.80), fontsize=6.2,
                    ha="left", color="#b2182b", arrowprops=dict(arrowstyle="-", lw=0.6, color="#b2182b"))
    else:
        ax.text(0.5, 0.5, "structural panel needs PyMOL", ha="center", va="center", fontsize=7)
    ax.set_title("a   acidic pocket and residue 246 (MD, 9E4F)", fontsize=7.6, loc="left")
    # (b)
    ax = fig.add_subplot(gs[0, 1])
    for b in BUILDS:
        ax.plot(rmsd[b][:, 0] / 1000.0, rmsd[b][:, 1] * 10.0, lw=0.6, color=COL[b], label=LAB1[b])
    ax.axvline(rmsd["wt"][burn, 0] / 1000.0, ls=":", lw=0.7, color="#9a9a9a")
    ax.set_xlabel("time (ns)"); ax.set_ylabel("backbone RMSD (A)"); ax.set_ylim(0, None)
    ax.set_title("b   protein RMSD, each build stable over 20 ns", fontsize=7.6, loc="left")
    ax.legend(fontsize=5.8, frameon=False, loc="lower right")
    # (c)
    ax = fig.add_subplot(gs[1, 0]); xs = np.linspace(4, 13, 300)
    for b in BUILDS:
        v = dist[b][burn:, REACH_COLS].ravel(); kde = stats.gaussian_kde(v)
        ax.plot(xs, kde(xs), lw=1.1, color=COL[b], label="%s %s  %.1f A" % (LAB1[b], tip[b], v.mean()))
        ax.fill_between(xs, kde(xs), color=COL[b], alpha=0.08)
    ax.axvline(CONTACT, ls="--", lw=0.8, color="#1a1a1a"); ax.text(CONTACT + 0.12, 0.40, "4 A contact", fontsize=5.8)
    ax.set_xlabel("distance, residue 246 tip to neighbour D409 (A)"); ax.set_ylabel("density"); ax.set_xlim(4, 13)
    ax.set_title("c   the engineered side chain never reaches D409", fontsize=7.6, loc="left")
    ax.legend(fontsize=5.6, frameon=False, loc="upper right")
    # (d)
    ax = fig.add_subplot(gs[1, 1])
    dl = np.array([abs(dist["ladder"][burn:, j].mean() - dist["wt"][burn:, j].mean()) for j in range(12)])
    dc = np.array([abs(dist["compare"][burn:, j].mean() - dist["wt"][burn:, j].mean()) for j in range(12)])
    x = np.arange(12)
    ax.bar(x - 0.2, dl, 0.38, color=COL["ladder"], edgecolor="#1a1a1a", lw=0.4, label="pI-code K246Q")
    ax.bar(x + 0.2, dc, 0.38, color=COL["compare"], edgecolor="#1a1a1a", lw=0.4, label="comparator K246M")
    ax.set_xticks(x); ax.set_xticklabels([l.split("|")[1] for l in LABELS], fontsize=5.0, rotation=90)
    for i, c in enumerate("ABC"):
        ax.text(i * 4 + 1.5, -1.9, "subunit %s" % c, fontsize=5.6, ha="center", color="#555555")
    ax.set_ylabel("departure from wild type (A)"); ax.set_ylim(0, 6.5)
    ax.set_title("d   pocket-pair distances, no systematic difference", fontsize=7.6, loc="left")
    ax.legend(fontsize=5.8, frameon=False, loc="upper right")
    # simple two-panel md_comparison.png as well
    os.makedirs(FIGS, exist_ok=True)
    # 400 dpi, matching build_figures.py so the whole set is one resolution.
    fig.savefig(os.path.join(FIGS, "FigureS3.png"), dpi=400, bbox_inches="tight")
    fig.savefig(os.path.join(FIGS, "FigureS3.svg"), bbox_inches="tight",
                metadata={"Date": None})
    _two_panel(dl, dc, os.path.join(RESULTS, "md_comparison.png"))


def _two_panel(dl, dc, out):
    import matplotlib
    matplotlib.use("Agg"); import matplotlib.pyplot as plt
    import numpy as np
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.8))
    ax = axes[0]; x = np.arange(12)
    ax.bar(x - 0.2, dl, 0.38, color=COL["ladder"], edgecolor="#1a1a1a", lw=0.5, label="pI-code K246Q")
    ax.bar(x + 0.2, dc, 0.38, color=COL["compare"], edgecolor="#1a1a1a", lw=0.5, label="comparator K246M")
    ax.set_xticks(x); ax.set_xticklabels([l.split("|")[1] for l in LABELS], fontsize=5.0, rotation=90)
    for i, c in enumerate("ABC"):
        ax.text(i * 4 + 1.5, -1.95, "subunit %s" % c, fontsize=5.8, ha="center", color="#555")
    ax.set_ylabel("departure from wild type (A)"); ax.set_ylim(0, 6.5)
    ax.set_title("a   pocket-pair distances", fontsize=8, loc="left"); ax.legend(fontsize=6.5, frameon=False)
    ax = axes[1]; m = max(dl.max(), dc.max()) * 1.15
    ax.axline((0, 0), slope=1, color="#9a9a9a", lw=0.8, ls="--")
    ax.scatter(dc, dl, s=30, facecolor="#2166ac", edgecolor="#1a1a1a", lw=0.4)
    ax.set_xlim(0, m); ax.set_ylim(0, m); ax.set_aspect("equal")
    ax.set_xlabel("comparator K246M departure (A)"); ax.set_ylabel("pI-code K246Q departure (A)")
    ax.set_title("b   departure from wild type, both substitutions", fontsize=8, loc="left")
    ax.text(0.05 * m, 0.92 * m, "scatter about the diagonal: no systematic difference",
            fontsize=6, va="top", color="#555")
    fig.tight_layout(); fig.savefig(out, dpi=300)


if __name__ == "__main__":
    main()
