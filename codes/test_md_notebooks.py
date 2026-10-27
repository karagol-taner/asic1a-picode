#!/usr/bin/env python3
"""
test_md_notebooks.py - regression tests for the MD notebooks and templates.

    python3 codes/test_md_notebooks.py [--templates <folder>]

Each test pulls the function out of the notebook or template that will be run
rather than a copy kept here, and the file parsers run against the CHARMM-GUI
downloads in md/. The templates are a separate notebook library and are not
part of this repository, so this runs only where MD_INFO or --templates points
at one.

Exits non-zero on any failure.
"""
import ast, glob, json, os, re, shutil, subprocess, sys, tempfile, textwrap

HERE = os.path.dirname(os.path.abspath(__file__))
PRJ = os.path.dirname(HERE)
def find_templates():
    """The MD_info notebooks folder, wherever this is being run from.

    The templates are a separate notebook library and are not part of this
    repository, so this check runs only where that library is available. MD_INFO
    in the environment wins, then the argument --templates.
    """
    cands = []
    if os.environ.get("MD_INFO"):
        cands.append(os.environ["MD_INFO"])
    for i, a in enumerate(sys.argv):
        if a == "--templates" and i + 1 < len(sys.argv):
            cands.append(sys.argv[i + 1])
    for c in cands:
        if os.path.isfile(os.path.join(c, "MD_CHARMMGUI_GROMACS_template_2026.ipynb")):
            return c
    raise SystemExit("FATAL: could not find the MD notebook templates, which are "
                     "not part of this repository.\nLooked in:"
                     "\n  " + ("\n  ".join(cands) or "(nothing configured)") +
                     "\nSet MD_INFO or pass --templates <folder>.")


TPL = find_templates()
SIM = os.path.join(TPL, "MD_CHARMMGUI_GROMACS_template_2026.ipynb")
ANA = os.path.join(TPL, "MD_analysis_pH_pair_template_2026.ipynb")
RUNNB = os.path.join(PRJ, "md", "ASIC1a_picode_MD_2_wt.ipynb")
FAIL = []


def note(ok, what, detail=""):
    print("   %-4s %s%s" % ("ok" if ok else "FAIL", what, "  " + detail if detail else ""))
    if not ok:
        FAIL.append(what)


def demagic(s):
    s = re.sub(r"^(\s*)!.*$", r"\1pass", s, flags=re.M)
    return re.sub(r"^%.*$", "pass", s, flags=re.M)


def pull(path, names, env=None):
    """exec only the named top-level defs from a notebook: no side effects."""
    g = {"os": os, "re": re, "json": json, "subprocess": subprocess, "glob": glob,
         "shutil": shutil, "PREFIX": "/usr/local/gromacs"}
    g.update(env or {})
    for c in json.load(open(path))["cells"]:
        if c["cell_type"] != "code":
            continue
        for n in ast.parse(demagic("".join(c["source"]))).body:
            if isinstance(n, ast.FunctionDef) and n.name in names:
                exec(compile(ast.Module(body=[n], type_ignores=[]), path, "exec"), g)
    missing = [n for n in names if n not in g]
    if missing:
        raise SystemExit("FATAL: %s has no %s" % (os.path.basename(path), missing))
    return g


class Res(object):
    def __init__(self, out, rc):
        self.stdout, self.returncode = out, rc


WHOLE = ("Checking file prod_01.xtc\nReading frame       0 time    0.000\n"
         "# Atoms  127552\nPrecision 0.001 (nm)\nReading frame     500 time 5000.000\n"
         "Last frame        500 time 5000.000\n\nItem        #frames Timestep (ps)\n"
         "Step             501    10\nTime             501    10\nCoords           501    10\n"
         "Velocities         0\nBox              501    10\n")
TRUNC = ("Checking file prod_03.xtc\nReading frame       0 time    0.000\n"
         "# Atoms  127552\nReading frame     160 time 1600.000\n"
         "Last frame        160 time 1600.000\n\nItem        #frames Timestep (ps)\n"
         "Coords           161    10\n")
NOTAB = ("Checking file prod_04.xtc\nReading frame       0 time    0.000\n"
         "Last frame        499 time 4990.000\n")
BROKEN = "Fatal error:\nCould not open prod_02.xtc\n"


def test_frame_parsing():
    print("frame counting, from the run notebook's own xtc_frames and xtc_last_time")
    g = pull(RUNNB, ["xtc_frames", "xtc_last_time"])
    real = subprocess.run
    g["os"] = type("O", (), {"path": type("P", (), {"isfile": staticmethod(lambda p: True),
                                                    "basename": staticmethod(os.path.basename),
                                                    "join": staticmethod(os.path.join)})})
    try:
        for name, out, rc, want_n, want_t in [("whole chunk", WHOLE, 0, 501, 5000.0),
                                              ("truncated chunk", TRUNC, 0, 161, 1600.0),
                                              ("no summary table", NOTAB, 0, 500, 4990.0),
                                              ("gmx could not read it", BROKEN, 1, None, None)]:
            subprocess.run = lambda *a, **k: Res(out, rc)
            got_n, got_t = g["xtc_frames"]("x.xtc"), g["xtc_last_time"]("x.xtc")
            note(got_n == want_n, "%-22s frames" % name, "%s, want %s" % (got_n, want_n))
            note(got_t == want_t, "%-22s last time" % name, "%s, want %s" % (got_t, want_t))
        subprocess.run = lambda *a, **k: Res("Checking file x\n", 0)
        try:
            g["xtc_frames"]("x.xtc")
            note(False, "unparseable output stops the run", "it returned instead")
        except SystemExit:
            note(True, "unparseable output stops the run")
    finally:
        subprocess.run = real


def test_analysis_frames_in():
    print("frame counting, from the analysis template's own frames_in")
    s = "".join(json.load(open(ANA))["cells"][7]["source"])
    blk = textwrap.dedent(s[s.index("    def frames_in"):s.index('    mode = "reused"')])
    state = {"rc": 0, "text": ""}
    g = {"gmx": lambda *a, **k: (state["rc"], state["text"]), "re": re, "os": os}
    exec(blk, g)
    for name, out, rc, want in [("whole", WHOLE, 0, 501), ("truncated", TRUNC, 0, 161),
                                ("no summary table", NOTAB, 0, 500),
                                ("gmx failed", BROKEN, 1, None)]:
        state["rc"], state["text"] = rc, out
        got = g["frames_in"]("fixed.xtc")
        note(got == want, "%-22s" % name, "%s, want %s" % (got, want))
    state["rc"], state["text"] = 0, "Checking file x\n"
    try:
        g["frames_in"]("fixed.xtc")
        note(False, "unparseable output stops the run", "it returned instead")
    except SystemExit:
        note(True, "unparseable output stops the run")


def test_resume_logic():
    print("which chunk the template's production cell starts from")
    src = "".join(json.load(open(SIM))["cells"][15]["source"])
    block = src[src.index("RESUME_FROM = 1"):src.index("t_all = time.time()")]
    sound = {}
    tmp = tempfile.mkdtemp()

    def frames(p):
        return sound.get(p) if os.path.isfile(p) else None

    def scen(name, files, snd, want, want_left=None):
        d = os.path.join(tmp, re.sub(r"\W+", "_", name))
        os.makedirs(d)
        cwd = os.getcwd()
        os.chdir(d)
        try:
            for f in files:
                open(f, "w").write("x")
            sound.clear()
            sound.update(snd)
            g = {"os": os, "re": re, "glob": glob, "xtc_frames": frames, "CHUNKS": 4,
                 "PER_CHUNK": 500, "EXT": "xtc", "print": lambda *a, **k: None}
            exec(block, g)
            left = sorted(os.listdir(d))
            ok = g["RESUME_FROM"] == want and (want_left is None or left == want_left)
            note(ok, "%-46s" % name, "starts at chunk %d" % g["RESUME_FROM"])
        finally:
            os.chdir(cwd)

    three = ["prefix"] and ["step5_%d.%s" % (k, e) for k in (1, 2, 3)
                            for e in ("gro", "xtc", "cpt")]
    scen("fresh run, nothing on disk", [], {}, 1, [])
    scen("chunks 1-3 whole, 4 never run", three,
         {"step5_1.xtc": 500, "step5_2.xtc": 500, "step5_3.xtc": 500}, 4)
    scen("chunk 2 truncated, so 2 and 3 both re-run", three,
         {"step5_1.xtc": 500, "step5_2.xtc": 101, "step5_3.xtc": 500}, 2,
         ["step5_1.cpt", "step5_1.gro", "step5_1.xtc"])
    scen("chunk 3 died mid-chunk, its cpt must survive", three[:6] +
         ["step5_3.cpt", "step5_3.xtc"],
         {"step5_1.xtc": 500, "step5_2.xtc": 500, "step5_3.xtc": 220}, 3,
         ["step5_1.cpt", "step5_1.gro", "step5_1.xtc", "step5_2.cpt", "step5_2.gro",
          "step5_2.xtc", "step5_3.cpt", "step5_3.xtc"])
    scen("chunk 1 unreadable, so all four re-run",
         ["step5_%d.%s" % (k, e) for k in (1, 2, 3, 4) for e in ("gro", "xtc", "cpt")],
         {}, 1, [])
    scen("499 of 500 frames still counts as whole", three,
         {"step5_1.xtc": 499, "step5_2.xtc": 499, "step5_3.xtc": 499}, 4)
    shutil.rmtree(tmp, ignore_errors=True)


def test_real_files():
    print("the templates' file parsers, against the CHARMM-GUI downloads in md/")
    runs = sorted(glob.glob(os.path.join(PRJ, "md", "*-charmm-gui-*", "gromacs")))
    if not runs:
        note(False, "found no CHARMM-GUI run folders under md/")
        return
    g = pull(SIM, ["read_pdb_atoms", "read_gro", "mdp_values", "frames_per_chunk"])
    counts = []
    for run in runs:
        b = os.path.basename(os.path.dirname(run)).split("-")[0]
        pdb = g["read_pdb_atoms"](os.path.join(run, "step3_input.pdb"))
        gro = g["read_gro"](os.path.join(run, "step3_input.gro"))
        note(len(pdb) == len(gro), "%-8s pdb and gro agree" % b,
             "%d atoms" % len(pdb))
        segs = sorted({s for _, s, _, _, _ in pdb})
        prot = [s for s in segs if s.startswith("PRO")]
        note(len(prot) == 3, "%-8s three protein segments" % b, ",".join(prot))
        dt, nsteps, nx, nf = g["mdp_values"](os.path.join(run, "step5_production.mdp"))
        note(dt > 0 and nsteps > 0 and nx > 0, "%-8s mdp reads" % b,
             "dt=%g nsteps=%d every=%d" % (dt, nsteps, nx))
        counts.append(len(pdb))
    note(max(counts) - min(counts) < 500, "atom counts comparable across builds",
         "spread %d" % (max(counts) - min(counts)))


def test_chunk_arithmetic():
    print("chunk length, write interval and frames per chunk, end to end")
    g = pull(RUNNB, ["mdp_get", "mdp_set", "want_frames"])
    runs = sorted(glob.glob(os.path.join(PRJ, "md", "*-charmm-gui-*", "gromacs")))
    if not runs:
        note(False, "found no CHARMM-GUI run folders under md/")
        return
    tmp = tempfile.mkdtemp()
    cp = os.path.join(tmp, "step5_production.mdp")
    shutil.copy2(os.path.join(runs[0], "step5_production.mdp"), cp)
    dt = float(g["mdp_get"](cp, "dt"))
    ns_build, target, ns_chunk = 20, 2000, 5
    want = int(round(ns_build * 1000.0 / dt / target))
    want = max(100, want - want % 100)
    g["mdp_set"](cp, "nstxout-compressed", want)
    per_chunk = g["want_frames"](cp, ns_chunk)
    nsteps = int(round(ns_chunk * 1000.0 / dt))
    note(nsteps // want == per_chunk, "mdrun -nsteps agrees with the frame target",
         "%d steps / %d = %d, want_frames says %d" % (nsteps, want, nsteps // want, per_chunk))
    note(abs(ns_build * 1000.0 / dt / want - target) < 1, "the build holds TARGET_FRAMES",
         "%d over %d ns" % (ns_build * 1000.0 / dt / want, ns_build))
    shutil.rmtree(tmp, ignore_errors=True)


def test_group_selections():
    print("every call that selects a CHARMM-GUI group passes an index file")
    import ast as _ast
    paths = sorted(glob.glob(os.path.join(PRJ, "md", "*.ipynb"))) + [SIM, ANA]
    bad = 0
    checked = 0
    for p in paths:
        for c in json.load(open(p))["cells"]:
            if c["cell_type"] != "code":
                continue
            src = demagic("".join(c["source"]))
            try:
                t = _ast.parse(src)
            except SyntaxError:
                continue
            for n in _ast.walk(t):
                if not isinstance(n, _ast.Call):
                    continue
                f = n.func
                name = getattr(f, "id", None) or getattr(f, "attr", None)
                if name not in ("gmx", "sh"):
                    continue
                kw = {k.arg: k.value for k in n.keywords if k.arg}
                si = kw.get("stdin_text")
                if not isinstance(si, _ast.Constant) or not isinstance(si.value, str):
                    continue
                groups = [g for g in si.value.split("\n") if g.strip()]
                # SOLU and SOLV live only in index.ndx; System is a default group
                if not any(g.upper() in ("SOLU", "SOLV") for g in groups):
                    continue
                checked += 1
                flat = []
                for a in n.args:
                    if isinstance(a, _ast.Constant):
                        flat.append(a.value)
                    elif isinstance(a, _ast.Starred):
                        fn = getattr(a.value.func, "id", "") if isinstance(a.value, _ast.Call) else ""
                        if fn == "ndx_args":
                            flat.append("-n")
                if "-n" not in flat:
                    bad += 1
                    note(False, "%s line %d selects %s with no index file"
                         % (os.path.basename(p), n.lineno, groups[0]))
    note(bad == 0, "%d SOLU/SOLV selections, all with -n" % checked)


def test_names():
    print("every name resolves on a Run All")
    paths = sorted(glob.glob(os.path.join(PRJ, "md", "*.ipynb"))) + [SIM, ANA]
    r = subprocess.run([sys.executable, os.path.join(HERE, "check_notebooks.py")] + paths,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    note(r.returncode == 0, "%d notebooks and templates" % len(paths),
         "" if r.returncode == 0 else r.stdout[-600:])


def main():
    for t in (test_frame_parsing, test_analysis_frames_in, test_resume_logic,
              test_real_files, test_chunk_arithmetic, test_group_selections,
              test_names):
        t()
        print("")
    if FAIL:
        raise SystemExit("%d test(s) failed:\n  %s" % (len(FAIL), "\n  ".join(FAIL)))
    print("all tests pass")


if __name__ == "__main__":
    main()
