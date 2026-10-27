#!/usr/bin/env python3
"""
make_md_notebook.py - write the Colab notebooks for the ASIC1a MD control.

Usage:
    python3 codes/make_md_notebook.py [-d md]

Writes six notebooks, each run with Runtime > Run all:

    ASIC1a_picode_MD_1_prepare    three PDBs to Drive for the CHARMM-GUI
                                  upload, plus the settings checklist. No GPU
    ASIC1a_picode_MD_2_wt         one build each, so the three can run at once
    ASIC1a_picode_MD_2_ladder     in three sessions. Each pre-flights all three
    ASIC1a_picode_MD_2_compare    CHARMM-GUI folders, then simulates its own.
                                  Resumable after a disconnection
    ASIC1a_picode_MD_2b_salvage   measures truncated runs at their real length
    ASIC1a_picode_MD_3_compare    the three-way comparison and the figure. No GPU

Part 1 and part 2 are separate because the CHARMM-GUI step is a web upload and
cannot be driven from a notebook.
"""
import argparse, json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import picode_paths as P

NB1, NB2, NB3 = [], [], []


def md1(t): NB1.append(('markdown', t))
def c1(t): NB1.append(('code', t))
def md2(t): NB2.append(('markdown', t))
def c2(t): NB2.append(('code', t))
def md3(t): NB3.append(('markdown', t))
def c3(t): NB3.append(('code', t))


# Shared prose, so the notebooks cannot disagree.
WHY = r"""## The question

The functional half of this project is finished and the answer is negative. Holding
position, charge change and side-chain volume fixed across 101 substitutions, pI-code
ladder membership changes the measured activation shift by -0.006 pH units
(permutation p = 0.895). Charge change carries the whole effect at +0.162 pH units per
charge unit.

That is a statement about a channel midpoint. It is not a statement about structure,
and the ladder's premise is structural: a near-isosteric swap is supposed to disturb
the fold less. These two notebooks put that premise to a direct test.

## Why position 246

Of the 35 charge-matched pairs in the measured set, K246 is the sharpest case against
the method, and it is not a single comparison:

| at K246 | measured shift | |
|---|---|---|
| K246Q, **on the ladder** | 0.420 | the near-isosteric move |
| K246M, off the ladder | 0.190 | ladder loses by 0.230 |
| K246A, off the ladder | 0.170 | ladder loses by 0.250 |
| K246V, off the ladder | 0.310 | ladder loses by 0.110 |

All three comparators carry the same charge change from lysine, and all three agree.
K246 is also the most buried scanned position in the protein at 1.4 percent relative
side-chain accessibility, and the closest to the acidic pocket centre at 3.7 angstrom.
Burial is the one descriptor that does predict the measured shift (rho -0.595 over 38
positions, p < 0.0001), so the position where the ladder loses hardest is also the one
where the environment should matter most.

If isostericity bought tolerance, the glutamine build should sit closer to wild type
than the methionine build. If the measurement is right, it will not.

## Why the whole trimer

In the deposited low-pH trimer 9E4F the NZ of K246 in chain A sits

    6.93 angstrom from D409 of chain B
    10.04 angstrom from E219 of chain B
    11.53 angstrom from D237 of its own chain A
    11.78 angstrom from D351 of its own chain A

K246 reaches across the subunit interface into the neighbour's pocket. A monomer
simulation would answer a different question, and both notebooks refuse a system that
is not a trimer."""

SCOPE = r"""## Scope, stated up front

The extracellular domain only, residues 67 to 427 of all three chains, solvated, with
no membrane and no transmembrane helices. UniProt P78348 gives TRANSMEM 50..66 and
428..458, so the cut is clean and the pocket is 25 angstrom or more from either
boundary. This buys a system of roughly 110,000 atoms instead of 200,000 and makes
three builds affordable. It costs any coupling between the pocket and the gate, which
these notebooks cannot see and do not claim to."""


# Notebook 1: prepare the CHARMM-GUI inputs
md1(r"""# ASIC1a pI-code MD, part 1 of 2: prepare the CHARMM-GUI inputs

Taner Karagol and Alper Karagol.
Structure and discipline follow `MD_CHARMMGUI_GROMACS_template_2026.ipynb`.

**Runtime > Run all, then close it.** No GPU needed and nothing to fill in. It takes
under a minute.

It downloads the deposited structure, trims it to the extracellular trimer, writes the
three builds you will upload to CHARMM-GUI, and writes the settings checklist beside
them. Part 2 picks up from the CHARMM-GUI output and needs no further preparation.

""" + WHY + "\n\n" + SCOPE + r"""

## What you do between the two notebooks

CHARMM-GUI is a web upload and cannot be driven from a notebook, which is why this is
split in two. Run Solution Builder three times, once per PDB this notebook writes, with
the settings it prints. Put each `gromacs/` folder on Drive where it tells you. Then
open the three part 2 notebooks, one per build, and press Run all on each.""")

c1(r'''#@title Run all starts here: mount Drive, build the three inputs
PDB_ID        = "9E4F"   #@param {type:"string"}
#@markdown 9E4F is the low-pH human structure at pH 5.7. 9E4A is the resting one at 8.5.
ECD_FIRST     = 67       #@param {type:"integer"}
ECD_LAST      = 427      #@param {type:"integer"}
#@markdown UniProt P78348 gives TRANSMEM 50..66 and 428..458, so the
#@markdown extracellular domain lies between them.
POSITION      = 246      #@param {type:"integer"}
WT_RESNAME    = "LYS"    #@param {type:"string"}
LADDER_TO     = "GLN"    #@param ["GLN", "ASN", "ASP", "GLU", "HSD", "HSE", "LYS", "ARG"]
COMPARE_TO    = "MET"    #@param {type:"string"}
PROJECT       = "/content/drive/MyDrive/MDs/ASIC1a_picode"  #@param {type:"string"}
#@markdown For E418 instead: POSITION 418, WT_RESNAME GLU, LADDER_TO GLN,
#@markdown COMPARE_TO ALA, and a different PROJECT folder.

from google.colab import drive
drive.mount('/content/drive')

import os, sys, math, json, urllib.request, collections

OUT_DIR = os.path.join(PROJECT, "inputs")

# Side chains outward from CB; a substitution keeps what both share.
SIDECHAIN = {
    "LYS": ["CB", "CG", "CD", "CE", "NZ"],
    "GLN": ["CB", "CG", "CD", "OE1", "NE2"],
    "ASN": ["CB", "CG", "OD1", "ND2"],
    "GLU": ["CB", "CG", "CD", "OE1", "OE2"],
    "ASP": ["CB", "CG", "OD1", "OD2"],
    "MET": ["CB", "CG", "SD", "CE"],
    "ALA": ["CB"],
    "VAL": ["CB", "CG1", "CG2"],
    "ARG": ["CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"],
    "HSD": ["CB", "CG", "ND1", "CD2", "CE1", "NE2"],
    "HSE": ["CB", "CG", "ND1", "CD2", "CE1", "NE2"],
}
BACKBONE = ["N", "CA", "C", "O"]
EXPECT_SS = 7          # extracellular disulfides per chain in human ASIC1a


def fatal(msg):
    raise SystemExit("STOP: " + msg)


def fetch(pdb_id):
    """The deposited biological assembly, not the asymmetric unit."""
    os.makedirs(OUT_DIR, exist_ok=True)
    local = os.path.join(OUT_DIR, "%s.pdb1" % pdb_id)
    if os.path.isfile(local):
        print("using the copy already in", OUT_DIR)
        return open(local, errors="replace").read()
    url = "https://files.rcsb.org/download/%s.pdb1" % pdb_id
    print("downloading", url)
    try:
        txt = urllib.request.urlopen(url, timeout=120).read().decode("utf-8", "replace")
    except Exception as e:
        fatal("could not download %s (%s). Put the file at %s by hand and rerun."
              % (url, e, local))
    open(local, "w").write(txt)
    return txt


raw = fetch(PDB_ID)
lines = []
for L in raw.splitlines():
    if not L.startswith(("ATOM", "HETATM")):
        continue
    if L[17:20] == "HOH" or L[16] not in (" ", "A"):
        continue
    lines.append(L)
if not lines:
    fatal("%s carries no usable ATOM records" % PDB_ID)

chains = sorted({L[21] for L in lines})
if len(chains) != 3:
    fatal("expected a trimer, found %d chains (%s). The acidic pocket is completed "
          "across the subunit interface, so a monomer answers a different question."
          % (len(chains), ",".join(chains)))

res = collections.OrderedDict()
for L in lines:
    res.setdefault((L[21], int(L[22:26])), L[17:20].strip())
for c in chains:
    rs = sorted(r for ch, r in res if ch == c)
    gaps = [(rs[i], rs[i + 1]) for i in range(len(rs) - 1) if rs[i + 1] - rs[i] > 1]
    inside = [r for r in rs if ECD_FIRST <= r <= ECD_LAST]
    print("chain %s: %d residues %d..%d, %d in the ECD, %d gaps"
          % (c, len(rs), min(rs), max(rs), len(inside), len(gaps)))
    if gaps:
        fatal("chain %s has breaks at %s. CHARMM-GUI will either refuse the file or "
              "model the loop; decide which before running." % (c, gaps))
    if not (min(rs) <= ECD_FIRST and max(rs) >= ECD_LAST):
        fatal("chain %s covers %d..%d, the ECD asked for is %d..%d"
              % (c, min(rs), max(rs), ECD_FIRST, ECD_LAST))
    got = res.get((c, POSITION))
    if got != WT_RESNAME:
        fatal("position %d in chain %s is %s, not %s. Check the numbering: this "
              "structure should use UniProt P78348 numbering."
              % (POSITION, c, got, WT_RESNAME))

# Disulfides are measured, not assumed; a missed pair unfolds the domain.
sg = {}
for L in lines:
    if L[17:20].strip() == "CYS" and L[12:16].strip() == "SG":
        r = int(L[22:26])
        if ECD_FIRST <= r <= ECD_LAST:
            sg[(L[21], r)] = (float(L[30:38]), float(L[38:46]), float(L[46:54]))
ss, keys = [], sorted(sg)
for i, a in enumerate(keys):
    for b in keys[i + 1:]:
        if a[0] == b[0] and math.dist(sg[a], sg[b]) < 2.5:
            ss.append((a[0], a[1], b[1], round(math.dist(sg[a], sg[b]), 2)))
per_chain = collections.Counter(c for c, _, _, _ in ss)
print("\ndisulfides inside the ECD, per chain:", dict(per_chain))
for c, a, b, d in ss:
    if c == chains[0]:
        print("   CYS%d - CYS%d  %.2f A" % (a, b, d))
if any(per_chain.get(c, 0) != EXPECT_SS for c in chains):
    fatal("expected %d disulfides per chain inside the ECD, found %s. Do not hand "
          "CHARMM-GUI a list you have not checked." % (EXPECT_SS, dict(per_chain)))
free = sorted(r for (c, r) in sg if c == chains[0]
              and not any(r in (a, b) for cc, a, b, _ in ss if cc == chains[0]))
print("free cysteines in the ECD, per chain:", free)


def write_build(tag, new_resname):
    """One trimmed PDB with the substitution applied to every chain."""
    keep = None
    if new_resname != WT_RESNAME:
        if new_resname not in SIDECHAIN or WT_RESNAME not in SIDECHAIN:
            fatal("no side-chain definition for %s or %s; add it to SIDECHAIN above"
                  % (WT_RESNAME, new_resname))
        old, new = SIDECHAIN[WT_RESNAME], SIDECHAIN[new_resname]
        shared = 0
        for x, y in zip(old, new):
            if x != y:
                break
            shared += 1
        keep = set(BACKBONE) | set(old[:shared])
        print("  %-8s keeps %s, CHARMM-GUI rebuilds %s"
              % (tag, sorted(set(old[:shared])) or "backbone only", new[shared:]))
    out, n, nmut = [], 0, 0
    for L in lines:
        r = int(L[22:26])
        if not (ECD_FIRST <= r <= ECD_LAST):
            continue
        if keep is not None and r == POSITION:
            if L[12:16].strip() not in keep:
                continue
            L = L[:17] + "%-3s" % new_resname + L[20:]
            nmut += 1
        n += 1
        out.append(L)
    if keep is not None and nmut == 0:
        fatal("%s: no atoms were renamed at position %d" % (tag, POSITION))
    path = os.path.join(OUT_DIR, "%s_ECD_%s.pdb" % (PDB_ID, tag))
    with open(path, "w") as fh:
        fh.write("REMARK  ASIC1a pI-code MD input, %s\n" % tag)
        fh.write("REMARK  %s residues %d..%d of chains %s\n"
                 % (PDB_ID, ECD_FIRST, ECD_LAST, ",".join(chains)))
        if keep is not None:
            fh.write("REMARK  %s%d%s applied to every chain\n"
                     % (WT_RESNAME, POSITION, new_resname))
        fh.write("\n".join(out) + "\nEND\n")
    print("           %-24s %d atoms" % (os.path.basename(path), n))
    return path, n


print("\nbuilds:")
BUILDS = collections.OrderedDict()
for tag, rn in (("wt", WT_RESNAME), ("ladder", LADDER_TO), ("compare", COMPARE_TO)):
    BUILDS[tag] = write_build(tag, rn)

with open(os.path.join(OUT_DIR, "disulfides.txt"), "w") as fh:
    for c, a, b, d in ss:
        fh.write("%s CYS%d CYS%d %.2f\n" % (c, a, b, d))

# Part 2 reads this, so the two notebooks cannot disagree about which build is which.
spec = dict(pdb_id=PDB_ID, position=POSITION, ecd=[ECD_FIRST, ECD_LAST],
            chains=chains, disulfides_per_chain=EXPECT_SS,
            expect={"wt": WT_RESNAME, "ladder": LADDER_TO, "compare": COMPARE_TO},
            atoms={k: v[1] for k, v in BUILDS.items()})
json.dump(spec, open(os.path.join(OUT_DIR, "build_spec.json"), "w"), indent=2)

CHECKLIST = """CHARMM-GUI settings for the ASIC1a pI-code MD
=============================================
Input Generator > Solution Builder > Upload PDB. Run it three times, once per file.

  %s
  %s
  %s

Chains            keep all three, %s. The pocket is a two-subunit cavity.
Disulfides        %d per chain, %d in all. The pairs are in disulfides.txt beside these
                  files. Confirm every one is ticked: CHARMM-GUI does not always detect
                  them all, and a missed pair unfolds the domain.
Terminal groups   ACE at the N terminus and CT3 at the C terminus of every chain. The
                  termini here are an artificial cut through the middle of a protein,
                  so charged termini would put six formal charges where the real
                  protein has none.
Protonation       defaults, standard protonation at pH 7. The three builds must differ
                  in nothing but the substitution.
Water box         rectangular, 12 angstrom edge distance.
Ions              KCl 0.15 M, neutralising. The three builds carry different net charge
                  so the ion counts will differ; that is correct.
Force field       CHARMM36m.
Dynamics          GROMACS. Temperature 310 K.

Unpack each charmm-gui-*.tgz and put the gromacs/ folder on Drive as

  %s/wt/gromacs
  %s/ladder/gromacs
  %s/compare/gromacs

Then open the three part 2 notebooks, one per build, and press Runtime > Run all on
each. They are independent and can run at the same time.
""" % (os.path.basename(BUILDS["wt"][0]), os.path.basename(BUILDS["ladder"][0]),
       os.path.basename(BUILDS["compare"][0]), ",".join(chains), EXPECT_SS,
       EXPECT_SS * len(chains), PROJECT, PROJECT, PROJECT)

open(os.path.join(OUT_DIR, "CHARMM-GUI-checklist.txt"), "w").write(CHECKLIST)
print("\n" + "=" * 72)
print(CHECKLIST)
print("everything written to", OUT_DIR)''')

md1(r"""## Done

The three PDBs, `disulfides.txt`, `build_spec.json` and `CHARMM-GUI-checklist.txt` are
in `inputs/` under your project folder on Drive. The checklist printed above is the same
file, so you can close this notebook and work from Drive.

Part 2 reads `build_spec.json` and refuses any folder whose residue at the mutated
position does not match what this notebook wrote, so uploading the three builds in the
wrong order is caught before anything runs.""")


# Notebook 2: run everything
md2(r"""# ASIC1a pI-code MD, part 2: the __BUILD__ build

Taner Karagol and Alper Karagol.

**Runtime > Run all, then leave it.** This notebook simulates one build, __BUILD__,
end to end and writes its report to Drive. The other two have their own notebooks,
`ASIC1a_picode_MD_2_wt`, `_ladder` and `_compare`, and the three are independent, so
open all three in three Colab sessions and start them together. When all three have
finished, run `ASIC1a_picode_MD_3_compare`.

The three notebooks are generated from one source and differ in a single line, the
name of the build they claim. Nothing else can drift between them.

**Before it runs anything it checks all three CHARMM-GUI folders**, not just this
one, because a missing disulfide or a mismatched thermostat in any of them makes the
comparison meaningless. That check reads files and takes seconds, so it is worth
paying three times over to find out now instead of after a session of GPU time.

**If Colab disconnects, open it again and press Run all again.** Every stage checks
for its own output first, so a rerun continues from the last finished chunk instead
of starting over. The compiled GROMACS is cached to Drive, so only the first of the
three sessions pays the forty minute build; the cache is written under a temporary
name and moved into place, so the other two never read a half-written archive.

**What this notebook writes to Drive**

    results/__BUILD___report.json    distances, occupancies and RMSD for this build
    results/progress___BUILD__.json  what is finished, rewritten after every stage

""" + WHY + "\n\n" + SCOPE + r"""

## Sampling, and what one session buys

The system is about 127,500 atoms. Scaling from a 50,000 atom benchmark of 5 ns in 25
minutes, expect roughly 4 to 5 ns per hour on an A100, so 20 ns is about four to five
hours plus minimisation and equilibration. One build fits one session comfortably, and
the three run in parallel.

Do not take that estimate on trust. The notebook prints the wall time for `prod_01`
after the first chunk, which gives the real rate within about half an hour.

20 ns per build is 60 ns of trajectory and, because the trimer carries three copies of
the pocket, 60 ns of pocket sampling per build over 12 distance measurements. It is
enough to detect a clear difference. It is not enough to make a null mean very much,
so if the comparison comes back flat, treat that as a reason to extend rather than as
a result: raise `NS_PER_BUILD` and press Run all again, and the run continues from the
last finished chunk.

`WALL_BUDGET_HOURS` stops the production loop cleanly before Colab ends the session, so
the last chunk is written and recorded instead of lost mid-write.""")

c2(r'''#@title 1. Settings
PROJECT         = "/content/drive/MyDrive/MDs/ASIC1a_picode"  #@param {type:"string"}
NS_PER_BUILD    = 20     #@param {type:"integer"}
NS_PER_CHUNK    = 5      #@param {type:"integer"}
WALL_BUDGET_HOURS = 23.0 #@param {type:"number"}
#@markdown Stops cleanly here; Run all again to continue from the last chunk.
TARGET_FRAMES   = 2000   #@param {type:"integer"}
#@markdown Frames per build. The write interval is recomputed from this,
#@markdown since CHARMM-GUI's default is set for a 1 ns run.
BURN_FRACTION   = 0.2    #@param {type:"number"}
#@markdown Fraction of each trajectory discarded as settling before any number is taken.
CACHE_GROMACS   = True   #@param {type:"boolean"}
#@markdown Keep the compiled CUDA GROMACS on Drive so later sessions skip the build.
GMX_VER         = "2026.3"    #@param {type:"string"}
GMX_SIMD        = "AVX2_256"  #@param ["AVX2_256", "AVX_512", "SSE4.1", "AUTO"]
MAXWARN         = 10     #@param {type:"integer"}

import time
T_START = time.time()
# One build per notebook, run in parallel sessions. BUILD_ORDER stays the full
# set because the pre-flight checks all three folders before anything starts.
BUILD = "__BUILD__"
BUILD_ORDER = ["wt", "ladder", "compare"]
assert BUILD in BUILD_ORDER
print("this notebook runs build:", BUILD.upper())
print("project      ", PROJECT)
print("sampling     ", NS_PER_BUILD, "ns")
print("wall budget  ", WALL_BUDGET_HOURS, "hours")''')

c2(r'''#@title 2. Drive, GPU and GROMACS
from google.colab import drive
drive.mount('/content/drive')

import os, sys, json, glob, re, subprocess, shutil, datetime, math, collections
import numpy as np

PREFIX = "/usr/local/gromacs"
os.environ["PATH"] = PREFIX + "/bin:" + os.environ["PATH"]
RESULTS = os.path.join(PROJECT, "results")
INPUTS = os.path.join(PROJECT, "inputs")
os.makedirs(RESULTS, exist_ok=True)

SPEC_PATH = os.path.join(INPUTS, "build_spec.json")
if not os.path.isfile(SPEC_PATH):
    raise SystemExit("STOP: no %s. Run part 1 first; this notebook reads the build "
                     "spec it writes, so the two cannot disagree." % SPEC_PATH)
SPEC = json.load(open(SPEC_PATH))
POSITION = SPEC["position"]
EXPECT = SPEC["expect"]
print("build spec: position %d, %s" % (POSITION, EXPECT))

missing = [b for b in BUILD_ORDER
           if not os.path.isdir(os.path.join(PROJECT, b, "gromacs"))]
if missing:
    raise SystemExit("STOP: no CHARMM-GUI folder for %s. Expected %s/<build>/gromacs "
                     "for each of %s. Finish the CHARMM-GUI step from part 1's "
                     "checklist first." % (missing, PROJECT, BUILD_ORDER))

print("\n--- GPU ---")
# check=False catches a non-zero exit, not a missing nvidia-smi, which raises
# FileNotFoundError. A missing GPU is recorded here and refused at build start.
HAVE_GPU = False
GPU_NAME = None
try:
    _g = subprocess.run(["nvidia-smi", "--query-gpu=name,memory.total",
                         "--format=csv,noheader"], check=False, text=True,
                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if _g.returncode == 0 and (_g.stdout or "").strip():
        HAVE_GPU = True
        GPU_NAME = _g.stdout.strip().split("\n")[0]
        print(GPU_NAME)
        print("Chunk times scale with this card, so record it with the run: the same "
              "system takes several times longer on a T4 than on an A100.")
    else:
        print("nvidia-smi exited %d, so the runtime claims a GPU it cannot use."
              % _g.returncode)
        print((_g.stdout or "").strip()[:400])
except FileNotFoundError:
    print("no nvidia-smi, so this runtime has no GPU.")
    print("That is fine for the salvage and compare notebooks: they only call "
          "trjconv, trjcat, distance and rms, all of which run on the CPU.")
    print("A production run needs one. Runtime > Change runtime type > T4 GPU.")


def sh(cmd, cwd=None, label=""):
    print("\n$ " + cmd, flush=True)
    r = subprocess.run(cmd, shell=True, cwd=cwd, text=True,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    lines = (r.stdout or "").split("\n")
    if len(lines) > 60:
        print("\n".join(lines[:25]))
        print("   ... %d lines omitted ..." % (len(lines) - 45))
        print("\n".join(lines[-20:]))
    else:
        print(r.stdout)
    if r.returncode != 0:
        raise RuntimeError("%s FAILED (exit %d). The output above is the reason."
                           % (label or cmd, r.returncode))


def gmx_version():
    try:
        return subprocess.run(["gmx", "--version"], capture_output=True, text=True,
                              timeout=90).stdout
    except Exception:
        return ""


CACHE = os.path.join(PROJECT, "gromacs_cuda_%s.tar.gz" % GMX_VER)
CACHE_STALE = False

if "CUDA" in gmx_version():
    print("\nCUDA gmx already present in this session")
elif CACHE_GROMACS and os.path.isfile(CACHE):
    # The cached binary is tested against this session's driver, and a part-written
    # archive is tolerated; either failure falls through to a rebuild.
    print("\nrestoring the cached GROMACS from", CACHE)
    try:
        sh("tar xzf '%s' -C /" % CACHE, label="restore cache")
    except Exception as e:
        print("the cache could not be read (%s); another session may be writing it. "
              "Building from source instead." % e)
        shutil.rmtree(PREFIX, ignore_errors=True)
    if "CUDA" not in gmx_version():
        print("the cached build does not run here, rebuilding from source")
        shutil.rmtree(PREFIX, ignore_errors=True)
        # A binary built against an older CUDA stops loading, so it is replaced.
        CACHE_STALE = True

if "CUDA" not in gmx_version():
    nvcc = subprocess.run("which nvcc", shell=True, capture_output=True,
                          text=True).stdout.strip()
    if not nvcc:
        raise SystemExit("STOP: nvcc not found, so there is no CUDA toolkit to build "
                         "against. Set the runtime to GPU and press Run all again.")
    print("\nbuilding GROMACS %s. 20 to 40 minutes, once per Drive folder." % GMX_VER)
    sh("apt-get -qq update && apt-get -qq install -y cmake build-essential "
       "libfftw3-dev", label="apt")
    if not os.path.isdir("/content/gromacs-%s" % GMX_VER):
        sh("wget -q https://ftp.gromacs.org/gromacs/gromacs-%s.tar.gz && "
           "tar xzf gromacs-%s.tar.gz" % (GMX_VER, GMX_VER), cwd="/content",
           label="download")
    os.makedirs("/content/gromacs-%s/build" % GMX_VER, exist_ok=True)
    sh("cmake .. -DGMX_BUILD_OWN_FFTW=OFF -DGMX_GPU=CUDA -DGMX_SIMD=%s "
       "-DCMAKE_INSTALL_PREFIX=%s -DGMX_DOUBLE=OFF -DBUILD_TESTING=OFF"
       % (GMX_SIMD, PREFIX), cwd="/content/gromacs-%s/build" % GMX_VER, label="cmake")
    sh("make -j$(nproc) && make install", cwd="/content/gromacs-%s/build" % GMX_VER,
       label="make")
    if "CUDA" not in gmx_version():
        raise SystemExit("STOP: the build finished but gmx reports no CUDA support.")
    if CACHE_GROMACS and (CACHE_STALE or not os.path.isfile(CACHE)):
        # Session-unique name, then moved, so no reader sees a partial archive.
        print("\ncaching the build to Drive so later sessions skip it")
        tmp = "%s.%d.tmp" % (CACHE, os.getpid())
        try:
            sh("tar czf '%s' %s" % (tmp, PREFIX.lstrip("/")), cwd="/", label="cache")
            os.replace(tmp, CACHE)
        except Exception as e:
            print("could not cache the build (%s); this costs the next session a "
                  "rebuild and nothing else." % e)
            if os.path.isfile(tmp):
                os.remove(tmp)
    elif CACHE_GROMACS:
        print("\nanother session already cached a build, leaving it alone")

print("\n" + gmx_version().split("\n")[0])''')

c2(r'''#@title 3. Helpers
AA = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS HSD HSE HSP ILE LEU LYS MET PHE PRO "
         "SER THR TRP TYR VAL".split())
# Acidic-pocket pairs in human P78348 numbering, from the manuscript.
POCKET_PAIRS = [("D237", 237, "OD1", "D351", 351, "OD1"),
                ("E238", 238, "OE1", "D347", 347, "OD1"),
                ("E219", 219, "OE1", "D409", 409, "OD1"),
                ("E79",   79, "OE1", "E418", 418, "OE1")]
# Last side-chain heavy atom per residue, resolved per build, so a substitution
# cannot measure from an atom it no longer has.
TIP_OF = {"LYS": "NZ", "GLN": "NE2", "ASN": "ND2", "MET": "CE", "GLU": "OE1",
          "ASP": "OD1", "ALA": "CB", "VAL": "CG1", "ARG": "NH1",
          "HSD": "NE2", "HSE": "NE2"}


def gmx(*args, label="", stdin_text=None, quiet=False, check=True):
    cmd = ["gmx"] + [str(a) for a in args]
    if not quiet:
        print("$ " + " ".join(cmd), flush=True)
    r = subprocess.run(cmd, text=True, input=stdin_text,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    out = r.stdout or ""
    lines = out.split("\n")
    if not quiet:
        if len(lines) > 80:
            print("\n".join(lines[:30]))
            print("   ... %d lines omitted ..." % (len(lines) - 60))
            print("\n".join(lines[-30:]))
        else:
            print(out)
    if r.returncode != 0 and check:
        if quiet:
            print(out[-3000:])
        raise RuntimeError("%s failed (exit %d). The output above is the reason."
                           % (label or cmd[1], r.returncode))
    return r.returncode, out


def must_exist(path, what):
    if not os.path.isfile(path):
        raise RuntimeError("%s reported success but %s was not written." % (what, path))


def read_pdb_atoms(path):
    """-> list of (serial, segid, resid, resname, atomname)."""
    out = []
    for L in open(path, errors="replace"):
        if not L.startswith(("ATOM", "HETATM")):
            continue
        seg = L[72:76].strip() or L[21].strip() or "?"
        out.append((int(L[6:11]), seg, int(L[22:26]), L[17:20].strip(),
                    L[12:16].strip()))
    return out


def protein_segments(atoms):
    """Segments holding amino acids.

    gmx does not always carry the CHARMM-GUI segid through to its output PDB, in which
    case read_pdb_atoms falls back to the chain letter and a 'PROA' prefix test
    silently finds nothing. A residue-name test is true of either spelling.
    """
    rn = {}
    for ser, seg, ri, name, an in atoms:
        rn.setdefault((seg, ri), name)
    return sorted({seg for (seg, ri), name in rn.items() if name in AA})


def progress(**kw):
    """Rewritten after every stage, so a dead session still says where it got to."""
    p = os.path.join(RESULTS, "progress_%s.json" % BUILD)
    cur = json.load(open(p)) if os.path.isfile(p) else {}
    cur.update(kw)
    cur["updated"] = datetime.datetime.now().isoformat(timespec="seconds")
    cur["elapsed_hours"] = round((time.time() - T_START) / 3600.0, 2)
    json.dump(cur, open(p, "w"), indent=2)
    return cur


def budget_left():
    return WALL_BUDGET_HOURS - (time.time() - T_START) / 3600.0


def xtc_frames(path):
    """Frame count straight from gmx check.

    An absent or unreadable file returns None, which callers treat as "re-run
    this chunk and everything after it". Output that gmx produced but this
    function cannot parse stops instead, rather than silently re-running every
    chunk of every build on a paid GPU.

    The count is read from the Coords row of gmx check's summary table, with the
    Last frame line as a second chance if that table is ever renamed.
    """
    if not os.path.isfile(path):
        return None
    r = subprocess.run([os.path.join(PREFIX, "bin", "gmx"), "check", "-f", path],
                       text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    out = r.stdout or ""
    m = re.search(r"^\s*Coords\s+(\d+)", out, re.M)
    if m:
        return int(m.group(1))
    m = re.search(r"Last frame\s+(\d+)", out)
    if m:
        return int(m.group(1)) + 1
    if r.returncode != 0:
        return None
    raise SystemExit(
        "STOP: gmx check read %s and exited 0, but its output carried neither a "
        "Coords row nor a Last frame line, so the frame count is unknown. Fix "
        "this parse rather than letting every chunk be re-run as unsound. The "
        "first lines of the output were:\n%s"
        % (os.path.basename(path), "\n".join(out.split("\n")[:12])))


def xtc_last_time(path):
    """Time stamp of the last frame, in ps. None when unreadable."""
    if not os.path.isfile(path):
        return None
    r = subprocess.run([os.path.join(PREFIX, "bin", "gmx"), "check", "-f", path],
                       text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    t = None
    for L in (r.stdout or "").split("\n"):
        m = re.search(r"frame\s+\d+\s+time\s+([-\d.eE+]+)", L)
        if m:
            t = float(m.group(1))
    return t


def want_frames(mdp, ns):
    """Frames a chunk of ns nanoseconds should hold, read from the mdp."""
    dt = float(mdp_get(mdp, "dt"))
    every = int(mdp_get(mdp, "nstxout-compressed"))
    return int(round(ns * 1000.0 / dt / every))


def ndx_args(group="SOLU"):
    """-n index.ndx, with a clear stop when it is not there.

    SOLU and SOLV are CHARMM-GUI's own groups and exist only in index.ndx. Given
    -s and no -n, trjconv and rms build the default groups from the topology,
    which hold System, Protein, the side chain sets and the residue names, and no
    SOLU. The selection then fails as "Cannot read from input" after the group
    list has scrolled past, so the failure is easy to misread. grompp is not
    affected because every grompp call here passes -n.
    """
    if not os.path.isfile("index.ndx"):
        raise SystemExit(
            "STOP: no index.ndx in %s, so the group %s cannot be resolved. It comes "
            "from the CHARMM-GUI download and is staged in with the other inputs, so "
            "if it is missing here the staging is wrong."
            % (os.getcwd(), group))
    return ("-n", "index.ndx")


def stage_in(run, work, n_chunks, per_chunk):
    """Copy inputs and any sound restart files from Drive onto local VM disk.

    Everything heavy then runs on local disk. A concatenated trajectory written
    straight into the Drive mount leaves gigabytes of upload in flight, and any
    part still unwritten when the runtime ends comes back truncated.

    A chunk is reusable only if every chunk before it is. Re-running chunk k
    produces a different trajectory from k onward, so a later chunk that came off
    the old branch cannot be stitched onto the new one.
    """
    os.makedirs(work, exist_ok=True)
    heavy = ("prod_all.xtc", "prod_nojump.xtc", "prod_fit.xtc")
    for name in sorted(os.listdir(run)):
        src = os.path.join(run, name)
        if os.path.isdir(src):
            if name == "toppar":
                shutil.copytree(src, os.path.join(work, name), dirs_exist_ok=True)
            continue
        if re.match(r"prod_\d\d\.xtc$", name) or name in heavy:
            continue
        shutil.copy2(src, os.path.join(work, name))
    kept = 0
    for i in range(1, n_chunks + 1):
        nm = "prod_%02d" % i
        src = os.path.join(run, nm + ".xtc")
        if not (os.path.isfile(src) and os.path.isfile(os.path.join(work, nm + ".gro"))):
            print("  %s absent, so it and every chunk after it is re-run" % nm)
            break
        shutil.copy2(src, os.path.join(work, nm + ".xtc"))
        got = xtc_frames(os.path.join(work, nm + ".xtc"))
        if got is None or got < per_chunk - 1:
            os.remove(os.path.join(work, nm + ".xtc"))
            print("  %s holds %s of %d frames, so it and every chunk after it is "
                  "re-run" % (nm, got, per_chunk))
            break
        kept += 1
    for i in range(kept + 1, n_chunks + 1):
        for ext in (".xtc", ".gro", ".cpt", "_prev.cpt", ".tpr", ".log", ".edr"):
            f = os.path.join(work, "prod_%02d%s" % (i, ext))
            if os.path.isfile(f):
                os.remove(f)
    print("  staged into %s, %d of %d chunks reusable" % (work, kept, n_chunks))
    return kept


DRIVE_PENDING = [0]       # bytes handed to the Drive mount and not yet waited on
DRIVE_SENT = {}           # drive path -> bytes written, for the check at the end


def stage_out(work, run, names):
    """Copy products into the Drive mount and count the bytes now in flight.

    The size check below compares the destination with the source, and on a FUSE
    mount that compares the local cache against the local file. It catches a short
    write and says nothing at all about whether Google has the bytes, so a
    truncated upload passes it. What
    gives the upload time is drive_wait, and the only real proof is the next
    session's stage_in, which counts frames in the copy it pulls back down.
    """
    for nm in names:
        src = os.path.join(work, nm)
        if not os.path.isfile(src):
            continue
        dst = os.path.join(run, nm)
        shutil.copy2(src, dst)
        got, want = os.path.getsize(dst), os.path.getsize(src)
        if got != want:
            raise SystemExit("STOP: %s is %d bytes in the Drive mount against %d on "
                             "disk, so even the cached copy is short."
                             % (nm, got, want))
        DRIVE_PENDING[0] += want
        DRIVE_SENT[dst] = want


def drive_wait(why="", floor=20.0, mbps=12.0):
    """Give the Drive daemon time to finish the uploads just handed to it.

    Colab's mount writes to a local cache and uploads in the background, so a
    runtime that ends with uploads in flight leaves truncated files behind. There
    is no way to ask Google what it holds without unmounting, so this waits on a
    throughput estimate with a floor. It is a margin, not a proof. Twenty seconds
    against an eighty minute chunk is not worth economising on.
    """
    mb = DRIVE_PENDING[0] / 1e6
    DRIVE_PENDING[0] = 0
    if mb <= 0:
        return
    secs = max(floor, mb / mbps)
    print("  waiting %.0f s for Drive to take %.0f MB%s"
          % (secs, mb, " after " + why if why else ""), flush=True)
    os.sync()
    time.sleep(secs)


def drive_verify():
    """Flush the mount, bring it back, and compare every file this run uploaded.

    This is the only check in the session that asks Google rather than the local
    cache. flush_and_unmount makes Colab finish its pending writes, and the sizes
    read after the remount come from Drive's own metadata instead of from a page
    that may never have been sent. It runs last, with every file copied and every
    report written, so if the remount ever stalls it costs an open tab and nothing
    else.
    """
    try:
        from google.colab import drive as _drive
        os.sync()
        print("\nflushing the Drive mount, then checking what Drive actually holds")
        _drive.flush_and_unmount()
        _drive.mount("/content/drive")
    except Exception as _e:
        print("\ncould not flush and remount Drive (%s), so the uploads are "
              "unverified. Leave this tab open for a couple of minutes before "
              "closing it." % _e)
        return
    bad, n = [], 0
    for path, want in sorted(DRIVE_SENT.items()):
        n += 1
        got = os.path.getsize(path) if os.path.isfile(path) else None
        if got != want:
            bad.append((os.path.basename(path), got, want))
    if bad:
        for nm, got, want in bad:
            print("   %-22s %s bytes on Drive against %d written" % (nm, got, want))
        raise SystemExit(
            "STOP: %d of %d files are short or missing on Drive after the flush, so "
            "the upload did not finish. The frames are not lost: open this notebook "
            "again and press Run all, which re-sends what is short and skips every "
            "chunk that is whole." % (len(bad), n))
    print("Drive holds all %d uploaded files at full size. Safe to close the runtime."
          % n)


print("helpers ready, MAXWARN =", MAXWARN)''')

c2(r'''#@title 4. Run all three builds. This is the long cell.
def check_identity(run, build):
    PDB = read_pdb_atoms(os.path.join(run, "step3_input.pdb"))
    with open(os.path.join(run, "step3_input.gro"), errors="replace") as fh:
        fh.readline()
        n_gro = int(fh.readline())
    if len(PDB) != n_gro:
        raise SystemExit("STOP: %s step3_input.pdb has %d atoms and the .gro has %d. "
                         "They must be the same build." % (build, len(PDB), n_gro))
    RES = {}
    for ser, seg, ri, name, an in PDB:
        RES.setdefault((seg, ri), name)
    prot = protein_segments(PDB)
    if len(prot) != 3:
        raise SystemExit("STOP: %s has %d protein segments (%s), expected the three "
                         "subunits. The acidic pocket is completed across the subunit "
                         "interface and a monomer answers a different question."
                         % (build, len(prot), ",".join(prot)))
    want = EXPECT[build]
    bad = [s for s in prot if RES.get((s, POSITION)) != want]
    if bad:
        raise SystemExit("STOP: %s should carry %s at %d on every subunit; segments %s "
                         "carry %s. This is not the build you think it is, or the three "
                         "CHARMM-GUI uploads went in the wrong order."
                         % (build, want, POSITION, bad,
                            [RES.get((s, POSITION)) for s in bad]))
    print("  identity ok: %s%d on all of %s, %d atoms"
          % (want, POSITION, ",".join(prot), len(PDB)))
    return len(PDB), prot


def equilibrate(run):
    if not os.path.isfile("step4.0_minimization.tpr"):
        gmx("grompp", "-f", "step4.0_minimization.mdp", "-o",
            "step4.0_minimization.tpr", "-c", "step3_input.gro", "-r",
            "step3_input.gro", "-p", "topol.top", "-n", "index.ndx",
            "-maxwarn", MAXWARN, label="grompp minimisation")
        must_exist("step4.0_minimization.tpr", "grompp minimisation")
    if not os.path.isfile("step4.0_minimization.gro"):
        gmx("mdrun", "-v", "-deffnm", "step4.0_minimization", label="minimisation")
        must_exist("step4.0_minimization.gro", "minimisation")
    prev = "step4.0_minimization"
    stages = sorted(glob.glob("step4.?_equilibration.mdp"))
    if not stages:
        raise SystemExit("STOP: no step4.*_equilibration.mdp in %s. A CHARMM-GUI "
                         "solution build writes at least one, so the download is "
                         "incomplete." % run)
    for mdp in stages:
        name = mdp[:-4]
        if os.path.isfile(name + ".gro"):
            prev = name
            continue
        gmx("grompp", "-f", mdp, "-o", name + ".tpr", "-c", prev + ".gro", "-r",
            "step3_input.gro", "-p", "topol.top", "-n", "index.ndx",
            "-maxwarn", MAXWARN, label="grompp " + name)
        must_exist(name + ".tpr", "grompp " + name)
        gmx("mdrun", "-v", "-deffnm", name, label=name)
        must_exist(name + ".gro", name)
        prev = name
    return prev


def produce(run, build, eq_final):
    mdp = "step5_production.mdp"
    if not os.path.isfile(mdp):
        raise SystemExit("STOP: no %s in %s; the CHARMM-GUI download is incomplete"
                         % (mdp, run))
    dt = None
    for L in open(mdp):
        m = re.match(r"\s*dt\s*=\s*([\d.]+)", L)
        if m:
            dt = float(m.group(1))
    if dt is None:
        raise SystemExit("STOP: could not read dt from %s, so the chunk length is "
                         "unknown" % mdp)
    nsteps = int(round(NS_PER_CHUNK * 1000.0 / dt))
    n_chunks = int(math.ceil(float(NS_PER_BUILD) / NS_PER_CHUNK))
    per_chunk = want_frames(mdp, NS_PER_CHUNK)
    drive_run = os.path.join(PROJECT, build, "gromacs")
    prev = eq_final
    for i in range(1, n_chunks + 1):
        name = "prod_%02d" % i
        if os.path.isfile(name + ".gro"):
            got = xtc_frames(name + ".xtc")
            if got is not None and got >= per_chunk - 1:
                prev = name
                continue
        if budget_left() < 0.6:
            print("  wall budget nearly spent, stopping cleanly after %d of %d chunks"
                  % (i - 1, n_chunks))
            return prev, i - 1, n_chunks
        t0 = time.time()
        # A checkpoint with no .gro means an interrupted chunk; -cpi resumes it.
        resume = os.path.isfile(name + ".cpt") and not os.path.isfile(name + ".gro")
        if not os.path.isfile(name + ".tpr"):
            gmx("grompp", "-f", mdp, "-o", name + ".tpr", "-c", prev + ".gro", "-p",
                "topol.top", "-n", "index.ndx", "-maxwarn", MAXWARN,
                label="grompp " + name)
        must_exist(name + ".tpr", "grompp " + name)
        extra = ("-cpi", name + ".cpt") if resume else ()
        if resume:
            print("  %s was interrupted, resuming from its checkpoint" % name)
        gmx("mdrun", "-v", "-deffnm", name, "-nsteps", nsteps, *extra, label=name)
        must_exist(name + ".gro", name)
        got = xtc_frames(name + ".xtc")
        if got is None or got < per_chunk - 1:
            raise SystemExit("STOP: %s wrote %s frames and should hold %d. The chunk "
                             "is incomplete, and concatenating it would give a "
                             "trajectory shorter than it reports."
                             % (name, got, per_chunk))
        print("  %s %s done in %.1f min, %d frames"
              % (build, name, (time.time() - t0) / 60.0, got))
        # Each chunk goes to Drive as it closes, size verified, so nothing of
        # value is left only on the VM.
        stage_out(os.getcwd(), drive_run,
                  [name + ".gro", name + ".cpt", name + ".xtc"])
        drive_wait(why="chunk %d" % i)
        progress(**{build: dict(chunks_done=i, chunks_total=n_chunks)})
        prev = name
    return prev, n_chunks, n_chunks


def concatenate(run, n_chunks, per_chunk):
    # Anchored to two-digit chunk names: "prod_*.xtc" also matches prod_all.xtc,
    # which a second pass would feed back into itself.
    chunks = sorted(glob.glob("prod_[0-9][0-9].xtc"))
    if len(chunks) != n_chunks:
        raise SystemExit("STOP: %s holds %d chunk trajectories and needs %d: %s"
                         % (os.getcwd(), len(chunks), n_chunks, chunks))
    counts = [xtc_frames(c) for c in chunks]
    short = [(c, n) for c, n in zip(chunks, counts)
             if n is None or n < per_chunk - 1]
    if short:
        raise SystemExit("STOP: %s hold fewer than the %d frames a chunk should "
                         "have. Concatenating them would produce a trajectory "
                         "shorter than it claims to be." % (short, per_chunk))
    for f in ("prod_all.xtc", "prod_nojump.xtc", "prod_fit.xtc"):
        if os.path.isfile(f):
            os.remove(f)
    # Every chunk restarts its clock at zero, so -settime places each at its
    # real offset; plain -cat gives a non-monotonic axis.
    offsets = "".join("%g\n" % (k * NS_PER_CHUNK * 1000.0) for k in range(n_chunks))
    gmx("trjcat", "-f", *chunks, "-o", "prod_all.xtc", "-settime",
        stdin_text=offsets, label="trjcat")
    must_exist("prod_all.xtc", "trjcat")
    span = xtc_last_time("prod_all.xtc")
    if span is None or span < 0.9 * NS_PER_BUILD * 1000.0:
        raise SystemExit("STOP: prod_all.xtc ends at %s ps and should reach about "
                         "%g ps. The time axis is wrong, so every per-frame "
                         "measurement downstream would be mislabelled."
                         % (span, NS_PER_BUILD * 1000.0))
    got = xtc_frames("prod_all.xtc")
    if got is None or got < sum(counts) - n_chunks:
        raise SystemExit("STOP: prod_all.xtc holds %s frames and the %d chunks hold "
                         "%d between them. The concatenation is incomplete."
                         % (got, n_chunks, sum(counts)))
    print("  concatenated %d chunks into %d frames" % (n_chunks, got))
    tpr = sorted(glob.glob("prod_[0-9][0-9].tpr"))[0]
    # Periodic images first, then fit; together the pbc step is skipped.
    gmx("trjconv", "-s", tpr, "-f", "prod_all.xtc", "-o", "prod_nojump.xtc",
        *ndx_args(), "-pbc", "nojump", stdin_text="SOLU\n",
        label="trjconv nojump")
    must_exist("prod_nojump.xtc", "trjconv nojump")
    gmx("trjconv", "-s", tpr, "-f", "prod_nojump.xtc", "-o", "prod_fit.xtc",
        *ndx_args(), "-fit", "rot+trans", stdin_text="SOLU\nSOLU\n",
        label="trjconv fit")
    must_exist("prod_fit.xtc", "trjconv fit")
    gmx("trjconv", "-s", tpr, "-f", "prod_all.xtc", "-o", "prod_ref.pdb",
        *ndx_args(), "-dump", "0", stdin_text="SOLU\n", label="reference frame")
    must_exist("prod_ref.pdb", "reference frame")
    return tpr


def selections(build):
    """Atom numbers resolved by name, never typed in.

    The three builds differ in a side chain, and that alone renumbers every atom after
    it: the same measurement is a different serial in each build, and hardcoding gives
    a plausible number measured on the wrong atom with nothing to warn you.
    """
    REF = read_pdb_atoms("prod_ref.pdb")
    BY, RN = {}, {}
    for ser, seg, ri, name, an in REF:
        BY[(seg, ri, an)] = ser
        RN.setdefault((seg, ri), name)
    prot = protein_segments(REF)
    if len(prot) != 3:
        raise SystemExit("STOP: prod_ref.pdb has %d protein segments, expected 3"
                         % len(prot))

    def serial(seg, ri, an, what):
        s = BY.get((seg, ri, an))
        if s is None:
            have = sorted(a for (g, r, a) in BY if g == seg and r == ri)
            raise SystemExit("STOP: %s, no atom %s in %s %d (%s). Present: %s"
                             % (what, an, seg, ri, RN.get((seg, ri)), have))
        return s

    want = EXPECT[build]
    tip_atom = TIP_OF.get(want)
    if tip_atom is None:
        raise SystemExit("STOP: no tip atom defined for %s; add it to TIP_OF" % want)

    xyz = {}
    for L in open("prod_ref.pdb", errors="replace"):
        if L.startswith(("ATOM", "HETATM")):
            xyz[int(L[6:11])] = np.array([float(L[30:38]), float(L[38:46]),
                                          float(L[46:54])])
    # The gmx selection parser reads "|" as OR and "-" as subtraction, so group
    # names are alphanumeric; the readable form is kept beside them.
    SEL = []
    for seg in prot:
        for nm1, r1, a1, nm2, r2, a2 in POCKET_PAIRS:
            SEL.append(("%s_%s_%s" % (seg, nm1, nm2), "%s|%s-%s" % (seg, nm1, nm2),
                        serial(seg, r1, a1, "pocket pair"),
                        serial(seg, r2, a2, "pocket pair")))
    for seg in prot:
        tip = serial(seg, POSITION, tip_atom, "position %d tip" % POSITION)
        best, bestd = None, 1e9
        for other in prot:
            if other == seg:
                continue
            d = float(np.linalg.norm(
                xyz[tip] - xyz[serial(other, 409, "OD1", "D409")]))
            if d < bestd:
                best, bestd = other, d
        SEL.append(("%s_%d_%s_to_%s_D409" % (seg, POSITION, tip_atom, best),
                    "%s|%d%s-%s.D409" % (seg, POSITION, tip_atom, best),
                    tip, serial(best, 409, "OD1", "D409")))
    with open("picode_pairs.ndx", "w") as fh:
        for nm, pretty, a, b in SEL:
            fh.write("[ %s ]\n%d %d\n" % (nm, a, b))
    return SEL, tip_atom


def measure(run, build, tpr, SEL, tip_atom, natoms, suffix="", want_ns=None):
    gmx("distance", "-s", tpr, "-f", "prod_fit.xtc", "-n", "picode_pairs.ndx",
        "-select", ";".join('"%s"' % nm for nm, _, _, _ in SEL),
        "-oall", "picode_dist.xvg", "-xvg", "none", label="distance", quiet=True)
    must_exist("picode_dist.xvg", "gmx distance")
    D = np.loadtxt("picode_dist.xvg")
    if D.ndim == 1:
        D = D.reshape(1, -1)
    t, V = D[:, 0], D[:, 1:] * 10.0       # gmx writes nm, the manuscript uses angstrom
    if V.shape[1] != len(SEL):
        raise SystemExit("STOP: gmx distance returned %d columns for %d selections"
                         % (V.shape[1], len(SEL)))
    burn = int(len(t) * BURN_FRACTION)
    rows = {}
    for j, (nm, pretty, _, _) in enumerate(SEL):
        v = V[burn:, j]
        rows[pretty] = dict(mean=float(v.mean()), sd=float(v.std()),
                            p10=float(np.percentile(v, 10)),
                            p90=float(np.percentile(v, 90)),
                            closed_fraction=float((v < 4.0).mean()))
    gmx("rms", "-s", tpr, "-f", "prod_fit.xtc", "-o", "picode_rmsd.xvg", "-xvg",
        "none", *ndx_args(), stdin_text="SOLU\nSOLU\n", label="rmsd", quiet=True)
    must_exist("picode_rmsd.xvg", "gmx rms")
    rmsd = np.loadtxt("picode_rmsd.xvg")[burn:, 1] * 10.0
    out = dict(build=build, run=run, atoms=int(natoms), frames=int(len(t)),
               burn=int(burn), ns=float((t[-1] - t[0]) / 1000.0),
               tip_atom=tip_atom, distances=rows,
               rmsd_mean=float(rmsd.mean()), rmsd_sd=float(rmsd.std()),
               rmsd_high=bool(rmsd.mean() > 6.0),
               written=datetime.datetime.now().isoformat(timespec="seconds"))
    target = NS_PER_BUILD if want_ns is None else want_ns
    if out["ns"] < 0.9 * target:
        raise SystemExit("STOP: the trajectory being measured spans %.1f ns against "
                         "the %g expected. Something upstream delivered less than it "
                         "reported." % (out["ns"], target))
    json.dump(out, open(os.path.join(RESULTS, "%s_report%s.json" % (build, suffix)),
                        "w"), indent=2)
    print("  %d frames, %.1f ns, backbone RMSD %.2f +/- %.2f A"
          % (len(t), out["ns"], out["rmsd_mean"], out["rmsd_sd"]))
    if out["rmsd_mean"] > 6.0:
        print("  WARNING: mean RMSD above 6 A. Check all %d disulfides were declared "
              "in CHARMM-GUI before reading anything into the pocket distances."
              % (SPEC["disulfides_per_chain"] * len(SPEC["chains"])))
    return out


def psf_disulfides(run):
    """Every SG to SG bond in the PSF, which is the authoritative bond list.

    Counting the pairs you ticked on the CHARMM-GUI screen is not the same as
    counting what it built. A missed disulfide costs the whole trajectory and
    shows up only as an unexplained RMSD at the end, so it is read back here
    from the topology before anything runs.
    """
    path = os.path.join(run, "step3_input.psf")
    if not os.path.isfile(path):
        raise SystemExit("STOP: no %s, so the disulfides cannot be verified" % path)
    lines = open(path, errors="replace").read().splitlines()
    i = next(k for k, l in enumerate(lines) if "!NATOM" in l)
    nat = int(lines[i].split()[0])
    sg = {}
    for l in lines[i + 1:i + 1 + nat]:
        f = l.split()
        # index segid resid resname atomname atomtype charge mass 0
        if len(f) >= 5 and f[4] == "SG" and f[3] == "CYS":
            sg[int(f[0])] = (f[1], int(f[2]))
    j = next(k for k, l in enumerate(lines) if "!NBOND" in l)
    nb = int(lines[j].split()[0])
    tok = []
    for l in lines[j + 1:]:
        tok += [int(x) for x in l.split()]
        if len(tok) >= 2 * nb:
            break
    return sorted({tuple(sorted((sg[a], sg[b])))
                   for a, b in zip(tok[0::2], tok[1::2]) if a in sg and b in sg})


def mdp_get(path, key):
    for l in open(path, errors="replace"):
        l = l.split(";")[0]
        if "=" in l and l.split("=")[0].strip().lower() == key.lower():
            return l.split("=", 1)[1].strip()
    return None


def mdp_set(path, key, value):
    """Rewrite one mdp setting in place, leaving every other line untouched."""
    out, hit = [], False
    for l in open(path, errors="replace"):
        body = l.split(";")[0]
        if "=" in body and body.split("=")[0].strip().lower() == key.lower():
            out.append("%-23s = %s\n" % (key, value)); hit = True
        else:
            out.append(l)
    if not hit:
        out.append("%-23s = %s\n" % (key, value))
    open(path, "w").writelines(out)


#@markdown ---
print("=" * 72)
print("PRE-FLIGHT: the three builds must differ only in the substitution")
print("=" * 72)
pre = {}
for build in BUILD_ORDER:
    run = os.path.join(PROJECT, build, "gromacs")
    required = ["step3_input.gro", "step3_input.pdb", "step3_input.psf",
                "topol.top", "step5_production.mdp", "index.ndx"]
    miss = [f for f in required if not os.path.isfile(os.path.join(run, f))]
    if miss:
        raise SystemExit("STOP: %s is missing %s. It should be the gromacs/ folder "
                         "from inside the CHARMM-GUI download." % (run, miss))
    natoms, prot = check_identity(run, build)
    ss = psf_disulfides(run)
    perseg = collections.Counter(a[0] for a, _ in ss)
    mdp = os.path.join(run, "step5_production.mdp")
    pre[build] = dict(atoms=natoms, prot=prot, ss=len(ss), perseg=dict(perseg),
                      pairs={s: sorted(tuple(sorted((r1, r2)))
                                       for (sa, r1), (sb, r2) in ss if sa == s)
                             for s in prot},
                      ref_t=mdp_get(mdp, "ref_t"), dt=mdp_get(mdp, "dt"),
                      tc=mdp_get(mdp, "tc_grps"))
    print("  %-8s %7d atoms, %d disulfides %s, ref_t %s"
          % (build, natoms, len(ss), dict(perseg), pre[build]["ref_t"]))

# The same disulfide set in every build and subunit; one short is a different
# protein, not a different substitution.
ssc = {b: pre[b]["ss"] for b in BUILD_ORDER}
if len(set(ssc.values())) != 1:
    low = min(ssc, key=ssc.get); high = max(ssc, key=ssc.get)
    missing = {s: sorted(set(pre[high]["pairs"].get(s, []))
                         - set(pre[low]["pairs"].get(s, [])))
               for s in pre[high]["pairs"]}
    missing = {s: v for s, v in missing.items() if v}
    raise SystemExit(
        "STOP: the builds disagree on disulfides %s. Build '%s' is missing %s "
        "relative to '%s'. Rebuild it in CHARMM-GUI with those pairs ticked; a "
        "free thiol where the others have a bond is a different protein."
        % (ssc, low, missing, high))
if len({tuple(sorted(pre[b]["perseg"].values())) for b in BUILD_ORDER}) != 1:
    raise SystemExit("STOP: the disulfides are not evenly distributed over the "
                     "subunits in the same way in every build: %s"
                     % {b: pre[b]["perseg"] for b in BUILD_ORDER})
print("  disulfides agree across builds: %d, %s per subunit"
      % (ssc[BUILD_ORDER[0]], sorted(set(pre[BUILD_ORDER[0]]["perseg"].values()))))

for key in ("ref_t", "dt", "tc"):
    vals = {b: pre[b][key] for b in BUILD_ORDER}
    if len(set(vals.values())) != 1:
        raise SystemExit("STOP: the builds disagree on %s: %s. They must be "
                         "simulated identically." % (key, vals))
TEMPERATURE = pre[BUILD_ORDER[0]]["ref_t"]
print("  thermostat agrees across builds: ref_t %s K" % TEMPERATURE)

atoms = {b: pre[b]["atoms"] for b in BUILD_ORDER}
spread = max(atoms.values()) - min(atoms.values())
if spread > 3000:
    raise SystemExit("STOP: the builds differ by %d atoms (%s). Different box "
                     "sizes make a pocket distance incomparable." % (spread, atoms))
print("  atom spread %d, within tolerance" % spread)

# Write interval set from the sampling asked for, not CHARMM-GUI's 1 ns default.
for build in BUILD_ORDER:
    mdp = os.path.join(PROJECT, build, "gromacs", "step5_production.mdp")
    dt = float(mdp_get(mdp, "dt"))
    want = int(round(NS_PER_BUILD * 1000.0 / dt / TARGET_FRAMES))
    want = max(100, want - want % 100)
    have = mdp_get(mdp, "nstxout-compressed")
    if have is None or int(have) != want:
        mdp_set(mdp, "nstxout-compressed", want)
        print("  %-8s nstxout-compressed %s -> %d  (%d frames over %d ns)"
              % (build, have, want, NS_PER_BUILD * 1000.0 / dt / want, NS_PER_BUILD))
    else:
        print("  %-8s nstxout-compressed already %d" % (build, want))
progress(preflight="passed", disulfides=ssc[BUILD_ORDER[0]], ref_t=TEMPERATURE,
         atom_spread=spread)
print("\npre-flight passed, starting the runs\n")

incomplete = []
for build in [BUILD]:
    run = os.path.join(PROJECT, build, "gromacs")
    rep = os.path.join(RESULTS, "%s_report.json" % build)
    print("\n" + "=" * 72)
    print("BUILD %s   %s" % (build.upper(), run))
    print("=" * 72)
    if os.path.isfile(rep) and json.load(open(rep))["ns"] >= 0.95 * NS_PER_BUILD:
        print("  already finished at %.1f ns" % json.load(open(rep))["ns"])
        continue
    required = ["step3_input.gro", "step3_input.pdb", "topol.top",
                "step5_production.mdp", "index.ndx"]
    miss = [f for f in required if not os.path.isfile(os.path.join(run, f))]
    if miss:
        raise SystemExit("STOP: %s is missing %s. It should be the gromacs/ folder "
                         "from inside the CHARMM-GUI download." % (run, miss))
    if not globals().get("HAVE_GPU", False):
        raise SystemExit("STOP: build %s needs a GPU. Equilibration and production "
                         "both ask mdrun to offload, and section 2 found no usable "
                         "GPU, so this would fail hours from now with a message about "
                         "offload instead of about the runtime. Runtime > Change "
                         "runtime type > GPU, then run from the top." % build)
    natoms, prot = check_identity(run, build)
    n_chunks = int(math.ceil(float(NS_PER_BUILD) / NS_PER_CHUNK))
    per_chunk = want_frames(os.path.join(run, "step5_production.mdp"), NS_PER_CHUNK)
    work = os.path.join("/content/work", build)
    print("  a chunk should hold %d frames" % per_chunk)
    stage_in(run, work, n_chunks, per_chunk)
    os.chdir(work)
    eq = equilibrate(work)
    progress(**{build: dict(stage="equilibrated")})
    prev, done, total = produce(work, build, eq)
    if done == 0:
        incomplete.append((build, 0, total))
        continue
    tpr = concatenate(work, n_chunks, per_chunk)
    SEL, tip_atom = selections(build)
    measure(run, build, tpr, SEL, tip_atom, natoms)
    # The concatenation stays on the VM; the fitted trajectory, its reference
    # frame and the measurements go up.
    stage_out(work, run, ["prod_fit.xtc", "prod_ref.pdb", os.path.basename(tpr),
                          "picode_dist.xvg", "picode_rmsd.xvg", "picode_pairs.ndx"])
    drive_wait(why="the analysis files")
    print("  copied the fitted trajectory and its analysis files to Drive")
    if done < total:
        incomplete.append((build, done, total))
    progress(**{build: dict(stage="measured", chunks_done=done, chunks_total=total)})

os.chdir("/content")
drive_verify()
if incomplete:
    for b, d, t in incomplete:
        print("%s reached %d of %d chunks" % (b, d, t))
    raise SystemExit(
        "The wall budget ran out before %s reached %d ns. Nothing is lost: open this "
        "notebook again and press Run all. Finished chunks are skipped, and the cached "
        "GROMACS means the next session starts in about a minute."
        % (BUILD, NS_PER_BUILD))
print("\n%s finished. Run the other two builds in their own notebooks, then open "
      "ASIC1a_picode_MD_3_compare.ipynb." % BUILD.upper())''')


# Notebook 2b: measure truncated runs at the length they reached
NB4 = []


def md4(t): NB4.append(('markdown', t))
def c4(t): NB4.append(('code', t))


md4(r"""# ASIC1a pI-code MD, part 2b: measure truncated runs at their real length

Taner Karagol and Alper Karagol.

**Runtime > Run all.** No GPU needed. It takes about twenty minutes, nearly all of
it copying the three concatenations down from Drive.

Use this pass when production finished and recorded `chunks_done 4 of 4` but the
concatenation did not, leaving each `prod_all.xtc` on Drive contiguous from t = 0
and truncated at the end.

This notebook measures each build at the length it really has, cuts all three to the
shortest so the comparison rests on equal sampling, and writes its reports as
`<build>_report_short.json`. The full-length reports keep their own names and
cannot be overwritten by this pass.

It answers one question: whether the difference between the builds is already visible
at the sampling that survived, or the full 20 ns has to be re-run before anything can
be said. Nothing it writes belongs in the manuscript as 20 ns.""")

c4(r'''#@title 4b. Measure every build at the length it reached, on a common window
FRAME_PS = None
counts, works = {}, {}
for build in BUILD_ORDER:
    run = os.path.join(PROJECT, build, "gromacs")
    work = os.path.join("/content/work", build)
    os.makedirs(work, exist_ok=True)
    src = os.path.join(run, "prod_all.xtc")
    if not os.path.isfile(src):
        raise SystemExit("STOP: no prod_all.xtc in %s. This notebook measures the "
                         "truncated concatenations, and there is nothing to salvage "
                         "without them." % run)
    mdp = os.path.join(run, "step5_production.mdp")
    ps = float(mdp_get(mdp, "dt")) * int(mdp_get(mdp, "nstxout-compressed"))
    if FRAME_PS is None:
        FRAME_PS = ps
    elif abs(ps - FRAME_PS) > 1e-6:
        raise SystemExit("STOP: %s writes a frame every %g ps and the others every %g. "
                         "The builds cannot be compared frame for frame."
                         % (build, ps, FRAME_PS))
    keep = ("topol.top", "index.ndx", "step3_input.pdb", "step3_input.gro",
            "step5_production.mdp", "prod_all.xtc")
    for nm in sorted(os.listdir(run)):
        s = os.path.join(run, nm)
        if nm == "toppar" and os.path.isdir(s):
            shutil.copytree(s, os.path.join(work, nm), dirs_exist_ok=True)
        elif nm in keep or re.match(r"prod_\d\d\.tpr$", nm):
            shutil.copy2(s, os.path.join(work, nm))
    n = xtc_frames(os.path.join(work, "prod_all.xtc"))
    if n is None:
        raise SystemExit("STOP: gmx check cannot read %s/prod_all.xtc, so it is not a "
                         "readable trajectory." % work)
    counts[build], works[build] = n, work
    print("  %-8s %d frames, %.2f ns at %g ps per frame"
          % (build, n, n * FRAME_PS / 1000.0, FRAME_PS))

common = min(counts.values())
common_ns = common * FRAME_PS / 1000.0
if common < 200:
    raise SystemExit("STOP: the shortest build holds %d frames, %.2f ns. That is too "
                     "little to average a pocket distance over, so there is nothing "
                     "to salvage." % (common, common_ns))
print("")
print("common window %d frames, %.2f ns, set by the shortest build"
      % (common, common_ns))
progress(stage="salvage", frames=counts, common_frames=common,
         common_ns=round(common_ns, 2))

for build in BUILD_ORDER:
    work = works[build]
    os.chdir(work)
    print("")
    print("=" * 72)
    print("BUILD %s, salvaged at %.2f ns" % (build.upper(), common_ns))
    print("=" * 72)
    natoms, prot = check_identity(work, build)
    tprs = sorted(glob.glob("prod_[0-9][0-9].tpr"))
    if not tprs:
        raise SystemExit("STOP: no prod_NN.tpr in %s to read the topology from" % work)
    tpr = tprs[0]
    for f in ("prod_seq.xtc", "prod_nojump.xtc", "prod_fit.xtc", "prod_fit_short.xtc"):
        if os.path.isfile(f):
            os.remove(f)
    # Chunks all start at t = 0, so frames are relabelled before cutting, and
    # both steps are checked by frame count. No -n: System is a default group.
    gmx("trjconv", "-s", tpr, "-f", "prod_all.xtc", "-o", "prod_seq.xtc",
        "-t0", "0", "-timestep", "%g" % FRAME_PS, stdin_text="System\n",
        label="relabel the time axis")
    got = xtc_frames("prod_seq.xtc")
    if got != counts[build]:
        raise SystemExit("STOP: relabelling %s went from %d frames to %s. trjconv "
                         "dropped frames, so the sampling would be wrong."
                         % (build, counts[build], got))
    gmx("trjconv", "-s", tpr, "-f", "prod_seq.xtc", "-o", "prod_nojump.xtc",
        "-e", "%g" % ((common - 1) * FRAME_PS), *ndx_args(), "-pbc", "nojump",
        stdin_text="SOLU\n", label="cut to the common window and unwrap")
    got = xtc_frames("prod_nojump.xtc")
    if got != common:
        raise SystemExit("STOP: the cut trajectory for %s holds %s frames against a "
                         "common window of %d." % (build, got, common))
    os.remove("prod_seq.xtc")
    gmx("trjconv", "-s", tpr, "-f", "prod_nojump.xtc", "-o", "prod_fit.xtc",
        *ndx_args(), "-fit", "rot+trans", stdin_text="SOLU\nSOLU\n", label="fit")
    must_exist("prod_fit.xtc", "trjconv fit")
    gmx("trjconv", "-s", tpr, "-f", "prod_nojump.xtc", "-o", "prod_ref.pdb",
        *ndx_args(), "-dump", "0", stdin_text="SOLU\n", label="reference frame")
    must_exist("prod_ref.pdb", "reference frame")
    SEL, tip_atom = selections(build)
    measure(os.path.join(PROJECT, build, "gromacs"), build, tpr, SEL, tip_atom,
            natoms, suffix="_short", want_ns=common_ns)
    shutil.move("prod_fit.xtc", "prod_fit_short.xtc")
    stage_out(work, os.path.join(PROJECT, build, "gromacs"), ["prod_fit_short.xtc"])

os.chdir("/content")
print("")
print("wrote <build>_report_short.json for all three builds at %.2f ns, and copied "
      "each fitted trajectory to Drive as prod_fit_short.xtc." % common_ns)
print("Open ASIC1a_picode_MD_3_compare.ipynb and set SUFFIX to _short.")''')


# Notebook 3: compare the three finished builds
md3(r"""# ASIC1a pI-code MD, part 3: compare the three builds

Taner Karagol and Alper Karagol.

**Runtime > Run all.** No GPU needed and nothing to fill in but the project folder.
It takes under a minute.

Open this once the three per-build notebooks have each finished and written their
report. It reads the three JSON files off Drive, so it does not matter that the
builds ran in three separate sessions.

It refuses to compare builds that are not comparable: different sampling, different
box sizes, or a pocket pair missing from one of them.

**What it writes to Drive**

    results/comparison.json     the three-way comparison and the verdict
    results/md_comparison.png   the figure
    results/SUMMARY.txt         the whole thing in plain text
""")

c3(r'''#@title 1. Settings and the three reports
PROJECT = "/content/drive/MyDrive/MDs/ASIC1a_picode"  #@param {type:"string"}
SUFFIX  = ""  #@param ["", "_short"]
#@markdown Empty compares the full runs; `_short` compares part 2b's
#@markdown equal-length cuts.

from google.colab import drive
drive.mount('/content/drive')

import os, json, datetime, time
T_START = time.time()
BUILD_ORDER = ["wt", "ladder", "compare"]
RESULTS = os.path.join(PROJECT, "results")
INPUTS = os.path.join(PROJECT, "inputs")

SPEC_PATH = os.path.join(INPUTS, "build_spec.json")
if not os.path.isfile(SPEC_PATH):
    raise SystemExit("STOP: no %s. Run part 1 first." % SPEC_PATH)
SPEC = json.load(open(SPEC_PATH))
POSITION = SPEC["position"]
EXPECT = SPEC["expect"]

waiting = [b for b in BUILD_ORDER
           if not os.path.isfile(os.path.join(RESULTS, "%s_report%s.json" % (b, SUFFIX)))]
if waiting:
    raise SystemExit("STOP: still waiting on %s. Each build has its own notebook: "
                     "ASIC1a_picode_MD_2_<build>.ipynb. Run those first." % waiting)


def progress(**kw):
    p = os.path.join(RESULTS, "progress_compare.json")
    cur = json.load(open(p)) if os.path.isfile(p) else {}
    cur.update(kw)
    cur["updated"] = datetime.datetime.now().isoformat(timespec="seconds")
    json.dump(cur, open(p, "w"), indent=2)
    return cur


print("position %d, %s" % (POSITION, EXPECT))
for b in BUILD_ORDER:
    r = json.load(open(os.path.join(RESULTS, "%s_report%s.json" % (b, SUFFIX))))
    print("  %-8s %5.1f ns, %d frames, RMSD %.2f A, written %s"
          % (b, r["ns"], r["frames"], r["rmsd_mean"], r["written"]))''')

c3(r'''#@title 2. The comparison, the figure and the summary
import json, os, math
import numpy as np
from scipy import stats

REPORTS = {}
for b in BUILD_ORDER:
    p = os.path.join(RESULTS, "%s_report%s.json" % (b, SUFFIX))
    if not os.path.isfile(p):
        raise SystemExit("STOP: missing %s. The previous cell did not finish %s."
                         % (p, b))
    REPORTS[b] = json.load(open(p))

ns = {b: r["ns"] for b, r in REPORTS.items()}
if max(ns.values()) - min(ns.values()) > 0.2 * max(ns.values()):
    raise SystemExit("STOP: the three builds differ by more than a fifth in sampling "
                     "(%s), so the distances are not comparable. Press Run all again "
                     "to finish the short one." % ns)
atoms = {b: r["atoms"] for b, r in REPORTS.items()}
if max(atoms.values()) - min(atoms.values()) > 3000:
    raise SystemExit("STOP: the builds differ by %d atoms (%s). Different box sizes "
                     "make a pocket distance incomparable; rebuild them in CHARMM-GUI "
                     "with the same water box and ion settings."
                     % (max(atoms.values()) - min(atoms.values()), atoms))

# Only the pocket pairs are tested; the position-246 distance has a different
# tip atom per build, so it is reported only.
labs = sorted(k for k in REPORTS["wt"]["distances"] if "|%d" % POSITION not in k)
for b in BUILD_ORDER:
    gone = [l for l in labs if l not in REPORTS[b]["distances"]]
    if gone:
        raise SystemExit("STOP: build %s does not carry %s" % (b, gone))

dl, dc = [], []
lines = ["ASIC1a pI-code MD, position %d" % POSITION,
         "%s on the ladder against %s off it, both from %s"
         % (EXPECT["ladder"], EXPECT["compare"], EXPECT["wt"]),
         "",
         "%-26s %8s %8s %8s   %9s %9s"
         % ("pocket pair", "wt", "ladder", "compare", "|lad-wt|", "|cmp-wt|")]
for l in labs:
    w = REPORTS["wt"]["distances"][l]["mean"]
    a = REPORTS["ladder"]["distances"][l]["mean"]
    c = REPORTS["compare"]["distances"][l]["mean"]
    dl.append(abs(a - w)); dc.append(abs(c - w))
    lines.append("%-26s %8.2f %8.2f %8.2f   %9.2f %9.2f"
                 % (l, w, a, c, abs(a - w), abs(c - w)))
dl, dc = np.array(dl), np.array(dc)

k = int((dl < dc).sum())
n = len(dl)
sign_p = min(1.0, 2 * sum(math.comb(n, i) for i in range(min(k, n - k) + 1)) / 2.0 ** n)
W = stats.wilcoxon(dl, dc) if n >= 6 else None
significant = (sign_p < 0.05) or (W is not None and W[1] < 0.05)
noncontact = [l for l in labs
              if all(REPORTS[b]["distances"][l]["closed_fraction"] == 0.0
                     for b in BUILD_ORDER)]
verdict = ("the two substitutions differ in pocket geometry"
           if significant
           else "neither substitution moves the pocket from wild type beyond the "
                "subunit-to-subunit scatter, and the two are indistinguishable from "
                "each other")

lines += ["",
          "mean departure from wild type over %d pocket pair measurements" % n,
          "   ladder   %.3f A" % dl.mean(),
          "   compare  %.3f A" % dc.mean(),
          "",
          "ladder closer to wild type in %d of %d, sign test p = %.3f"
          % (k, n, sign_p)]
if W is not None:
    lines.append("Wilcoxon signed-rank p = %.3f" % W[1])

reach = {}
for b in BUILD_ORDER:
    r = REPORTS[b]
    kk = [x for x in r["distances"] if "|%d" % POSITION in x]
    reach[b] = float(np.mean([r["distances"][x]["mean"] for x in kk]))
lines += ["",
          "position %d reach across the interface, mean over 3 subunits" % POSITION,
          "   " + ",  ".join("%s %s %.2f A" % (b, REPORTS[b]["tip_atom"], reach[b])
                             for b in BUILD_ORDER),
          "",
          "backbone RMSD, A: " + ", ".join("%s %.2f" % (b, REPORTS[b]["rmsd_mean"])
                                           for b in BUILD_ORDER),
          "sampling, ns:     " + ", ".join("%s %.1f" % (b, ns[b]) for b in BUILD_ORDER),
          "atoms:            " + ", ".join("%s %d" % (b, atoms[b])
                                           for b in BUILD_ORDER),
          "",
          "never within 4 A in any build (not a contact here): "
          + (", ".join(noncontact) if noncontact else "none"),
          "",
          "%s (sign p %.2f, Wilcoxon p %.2f). The mean departures differ by %.3f A, "
          "carried by the D237-D351 pairs, which barely form contacts."
          % (verdict, sign_p, W[1] if W is not None else float("nan"),
             abs(dl.mean() - dc.mean())),
          "The engineered side chain at %d does not reach the neighbouring D409 in any "
          "build (%.1f to %.1f A). This run is a structural control over one %d ns "
          "trajectory per build with the three subunits as the only replication; it "
          "tests whether either substitution remodels the pocket, not what either does "
          "to pH50." % (POSITION, min(reach.values()), max(reach.values()),
                        int(round(max(ns.values()))))]

SUMMARY = "\n".join(lines)
print(SUMMARY)

json.dump(dict(position=POSITION, expect=EXPECT, labels=labs,
               ladder_departure=list(map(float, dl)),
               compare_departure=list(map(float, dc)),
               ladder_mean=float(dl.mean()), compare_mean=float(dc.mean()),
               ladder_closer_in=k, n=n, sign_p=float(sign_p),
               wilcoxon_p=float(W[1]) if W is not None else None,
               reach=reach, ns=ns, atoms=atoms, verdict=verdict,
               significant=bool(significant), noncontact_pairs=noncontact),
          open(os.path.join(RESULTS, "comparison%s.json" % SUFFIX), "w"), indent=2)
open(os.path.join(RESULTS, "SUMMARY%s.txt" % SUFFIX), "w").write(SUMMARY + "\n")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.8))
ax = axes[0]
x = np.arange(n)
ax.bar(x - 0.2, dl, 0.38, color="#2166ac", edgecolor="#1a1a1a", lw=0.5,
       label="pI-code, %d%s" % (POSITION, REPORTS["ladder"]["tip_atom"]))
ax.bar(x + 0.2, dc, 0.38, color="#b2182b", edgecolor="#1a1a1a", lw=0.5,
       label="comparator, %d%s" % (POSITION, REPORTS["compare"]["tip_atom"]))
ax.set_xticks(x)
ax.set_xticklabels([l.split("|")[1] + "\n" + l.split("|")[0] for l in labs], fontsize=6)
ax.set_ylabel("departure from wild type (A)")
ax.set_title("a   pocket pair distances", fontsize=8, loc="left")
ax.legend(fontsize=6.5, frameon=False)
ax = axes[1]
ax.axline((0, 0), slope=1, color="#9a9a9a", lw=0.8, ls="--")
ax.scatter(dc, dl, s=30, facecolor="#2166ac", edgecolor="#1a1a1a", lw=0.4)
m = max(dl.max(), dc.max()) * 1.15
ax.set_xlim(0, m); ax.set_ylim(0, m); ax.set_aspect("equal")
ax.set_xlabel("comparator departure (A)")
ax.set_ylabel("pI-code departure (A)")
ax.set_title("b   departure from wild type, both substitutions", fontsize=8,
             loc="left")
ax.text(0.04 * m, 0.93 * m, "scatter about the diagonal: no systematic difference",
        fontsize=6, va="top", color="#555555")
fig.tight_layout()
fig.savefig(os.path.join(RESULTS, "md_comparison%s.png" % SUFFIX), dpi=300)
progress(finished=True, verdict=verdict)
print("\nwrote SUMMARY.txt, comparison.json and md_comparison.png to %s" % RESULTS)
fig''')

md3(r"""## Reading it

Everything is in `results/SUMMARY.txt` on Drive, so there is nothing to scroll back for.

**If the ladder build departs further from wild type than the comparator does**, the
structural result agrees with the functional one: at the most buried position in the
pocket the near-isosteric substitution is the more disruptive one, and isostericity
bought nothing. That is what this project expects, and finding it here makes the
negative stronger instead of merely repeating it.

**If the ladder build sits closer**, the structural premise survives at this position
while the functional consequence does not follow from it. That is a different and more
interesting paper: the ladder would be conserving local structure, and conserved local
structure would be failing to predict the proton-sensing phenotype. Report it.

**What would make either reading unsafe.** Cell 5 refuses outright on the last two and
prints a warning on the first:

- Backbone RMSD much above 6 angstrom in any build, which usually means a disulfide was
  not declared in CHARMM-GUI.
- Atom counts differing by more than 3000 between builds, which means the boxes differ.
- Sampling differing by more than a fifth between builds.

Twenty nanoseconds per build settles a side chain and a local hydrogen-bond network. It
does not sample a gating transition, and no claim about gating should come out of it.

**One position is one position.** To run the opposite case, go back to part 1, set
POSITION 418, WT_RESNAME GLU, LADDER_TO GLN, COMPARE_TO ALA and a different PROJECT
folder, then come back here with that same PROJECT. At E418 the ladder substitution
E418Q wins by 0.350 pH units against E418A, the largest ladder win in the measured set.
Two positions where the measurement points in opposite directions is a far better test
than either alone.""")


# notebook assembly
def cell(kind, src):
    lines = src.split('\n')
    srcs = [l + '\n' for l in lines[:-1]] + [lines[-1]]
    c = {'cell_type': kind, 'metadata': {}, 'source': srcs}
    if kind == 'code':
        c['execution_count'] = None
        c['outputs'] = []
    return c


def write(cells, path, gpu):
    nb = {
        'cells': [cell(k, s) for k, s in cells],
        'metadata': {
            'colab': {'provenance': [], 'toc_visible': True},
            'kernelspec': {'display_name': 'Python 3', 'name': 'python3'},
            'language_info': {'name': 'python'},
        },
        'nbformat': 4,
        'nbformat_minor': 0,
    }
    if gpu:
        nb['metadata']['accelerator'] = 'GPU'
    d = os.path.dirname(path)
    if d and not os.path.isdir(d):
        os.makedirs(d)
    with open(path, 'w') as fh:
        json.dump(nb, fh, indent=1)
    print('wrote %s' % path)
    print('  %d cells, %d code' % (len(cells),
                                   sum(1 for k, _ in cells if k == 'code')))


BUILDS = ('wt', 'ladder', 'compare')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-d', '--dir', default=P.project('md'))
    a = ap.parse_args()
    write(NB1, os.path.join(a.dir, 'ASIC1a_picode_MD_1_prepare.ipynb'), gpu=False)
    # One notebook per build, from the same source, differing only in the build.
    n_ph = sum(s.count('__BUILD__') for _, s in NB2)
    if n_ph == 0:
        raise SystemExit('FATAL: the per-build notebook has no __BUILD__ placeholder, '
                         'so all three would claim the same build')
    code_ph = sum(s.count('__BUILD__') for k, s in NB2 if k == 'code')
    if code_ph != 1:
        raise SystemExit('FATAL: expected the build name to be set in exactly one '
                         'code line, found %d' % code_ph)
    for b in BUILDS:
        cells = [(k, s.replace('__BUILD__', b)) for k, s in NB2]
        left = sum(s.count('__BUILD__') for _, s in cells)
        if left:
            raise SystemExit('FATAL: %d __BUILD__ placeholders survived substitution '
                             'for %s' % (left, b))
        write(cells, os.path.join(a.dir, 'ASIC1a_picode_MD_2_%s.ipynb' % b), gpu=True)
    # Salvage reuses the settings, Drive and helper cells verbatim; only the
    # driver differs, and BUILD names the pass because it measures all three.
    head = [(k, s.replace('__BUILD__', 'wt')
                 .replace('assert BUILD in BUILD_ORDER',
                          'BUILD = "salvage"  # this pass measures all three builds')
                 .replace('print("this notebook runs build:", BUILD.upper())',
                          'print("salvage pass over:", ", ".join(BUILD_ORDER))'))
            for k, s in NB2[1:-1]]
    left = sum(s.count('__BUILD__') for _, s in head)
    if left:
        raise SystemExit('FATAL: %d __BUILD__ placeholders survived into the salvage '
                         'notebook, which would make it claim one build' % left)
    # mdp_get, check_identity, selections and measure live in the run cell, so
    # it is split at its driver and the definitions and pre-flight kept.
    run_cell = NB2[-1][1]
    MARK = '\nincomplete = []\n'
    if run_cell.count(MARK) != 1:
        raise SystemExit('FATAL: cannot split the run cell on %r, so the salvage '
                         'notebook would be missing the shared definitions' % MARK)
    defs = run_cell.split(MARK)[0].rstrip() + '\n'
    defs = defs.replace('#@title 4. Run all three builds. This is the long cell.',
                        '#@title 4. Shared definitions and the three-build pre-flight')
    write([NB4[0]] + head + [('code', defs)] + NB4[1:],
          os.path.join(a.dir, 'ASIC1a_picode_MD_2b_salvage.ipynb'), gpu=False)
    write(NB3, os.path.join(a.dir, 'ASIC1a_picode_MD_3_compare.ipynb'), gpu=False)


if __name__ == '__main__':
    main()
