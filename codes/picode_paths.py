#!/usr/bin/env python3
"""
picode_paths.py - resolve the folders and constants this project reads from.

Environment variables, which win over codes/picode_local.json:

    PICODE_METHOD_DIR   the pI-code method folder (holds targets/ and tools/)
    PICODE_MEMCTRL_DIR  the folder holding data/humsavar.txt and the
                        AlphaMissense substitution table
    PICODE_LPO_DIR      the folder holding the burial script

The pI-code scoring scripts are not part of this project. Where they are
absent, load_ladder and load_lpo read the committed snapshot in data/, which
holds the ladder edges, the residue charges, the accessibility cutoffs and the
maximum-ASA table. The method folder takes precedence where it is present;
codes/export_ladder_definition.py --check compares the two.

Importing this module touches no files. Each resolver exits non-zero naming the
file it could not find.
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# Folder locations, kept out of the repository; see picode_local.json.example.
LOCAL_FILE = os.path.join(HERE, 'picode_local.json')


def _load_local():
    if not os.path.isfile(LOCAL_FILE):
        return {}
    try:
        with open(LOCAL_FILE) as fh:
            d = json.load(fh)
    except Exception as e:
        sys.exit('FATAL: cannot parse %s: %s' % (LOCAL_FILE, e))
    if not isinstance(d, dict):
        sys.exit('FATAL: %s must hold a JSON object' % LOCAL_FILE)
    return d


_LOCAL = _load_local()
_PARENTS = _LOCAL.get('search_parents') or ['']


def _search(key):
    """Candidate locations for a configured folder, nearest first.

    An absolute path in the configuration is taken as given. A relative one is
    looked for under each ancestor of the project root in turn, optionally
    beneath one of the search_parents, so the project can be moved without
    editing anything. marker is a file that must be present, so a folder of the
    right name but the wrong contents is not silently accepted.
    """
    spec = _LOCAL.get(key) or {}
    rel = spec.get('path')
    marker = spec.get('marker') or ''
    if not rel:
        return '', []
    if os.path.isabs(rel):
        return rel, [rel]
    out, d = [], ROOT
    for _ in range(4):
        d = os.path.dirname(d)
        if not d or d == os.sep:
            break
        for mid in _PARENTS:
            c = os.path.join(d, mid, rel) if mid else os.path.join(d, rel)
            c = os.path.normpath(c)
            if c not in out:
                out.append(c)
    if marker:
        hit = [c for c in out if os.path.isfile(os.path.join(c, marker))]
    else:
        hit = [c for c in out if os.path.isdir(c)]
    return (hit or out)[0], out


DEFAULT_METHOD, _CAND_METHOD = _search('method')
DEFAULT_MEMCTRL, _CAND_MEMCTRL = _search('memctrl')
DEFAULT_LPO, _CAND_LPO = _search('lpo')

_CANDIDATES = {'PICODE_METHOD_DIR': _CAND_METHOD,
               'PICODE_MEMCTRL_DIR': _CAND_MEMCTRL,
               'PICODE_LPO_DIR': _CAND_LPO}


def _dir(env, default, what):
    d = os.environ.get(env, default)
    if d and os.path.isdir(d):
        return d
    tried = _CANDIDATES.get(env, [])
    if not _LOCAL:
        sys.exit('FATAL: cannot find %s, and there is no local configuration.\n'
                 'This folder is not part of the repository. Only '
                 'extract_edge_costs.py, extract_interface_carboxylates.py and\n'
                 'asic_edge_vs_cost.py need it; their committed output is under '
                 'results/, which\nevery other script reads.\n'
                 'To point at a local copy, set %s, or copy\n'
                 '  %s.example\nto\n  %s\nand fill it in.'
                 % (what, env, LOCAL_FILE, LOCAL_FILE))
    sys.exit('FATAL: cannot find %s.\nlooked in:\n%s\nset %s to its location, '
             'or correct %s'
             % (what, '\n'.join('  ' + c for c in tried) or '  ' + str(d),
                env, LOCAL_FILE))


def _file(path, what, env):
    if not os.path.exists(path):
        sys.exit('FATAL: cannot find %s at %s\nset %s to the folder holding it'
                 % (what, path, env))
    return path


def method_dir():
    return _dir('PICODE_METHOD_DIR', DEFAULT_METHOD, 'the pi-code method folder')


def ladder_edge_cost():
    return _file(os.path.join(method_dir(), 'targets', 'ladder_edge_cost.py'),
                 'ladder_edge_cost.py', 'PICODE_METHOD_DIR')


def ladder_humsavar():
    return _file(os.path.join(method_dir(), 'targets', 'ladder_humsavar.py'),
                 'ladder_humsavar.py', 'PICODE_METHOD_DIR')


def interface_carboxylates():
    return _file(os.path.join(method_dir(), 'targets', 'interface_carboxylates.py'),
                 'interface_carboxylates.py', 'PICODE_METHOD_DIR')


def lpo_design():
    d = _dir('PICODE_LPO_DIR', DEFAULT_LPO, 'the folder holding the burial script')
    return _file(os.path.join(d, 'lpo_picode_design.py'),
                 'lpo_picode_design.py', 'PICODE_LPO_DIR')


LADDER_SNAPSHOT = os.path.join(ROOT, 'data', 'ladder_definition.json')
BURIAL_SNAPSHOT = os.path.join(ROOT, 'data', 'burial_thresholds.json')

_NOTED = set()


class _Snapshot(object):
    """Attribute holder presenting the same names as the method-folder module."""

    def __init__(self, source):
        self.SOURCE = source


def _note(what, path):
    """Say once per process that a snapshot stood in for the method folder."""
    if what in _NOTED:
        return
    _NOTED.add(what)
    sys.stderr.write('note: %s read from %s; the original scoring scripts are '
                     'not present\n' % (what, os.path.relpath(path, ROOT)))


def _snapshot_json(path, what, keys):
    import json
    if not os.path.isfile(path):
        sys.exit('FATAL: no method folder and no %s at %s\nthe repository should '
                 'ship this file; regenerate it with\n  python3 '
                 'codes/export_ladder_definition.py' % (what, path))
    try:
        with open(path) as fh:
            d = json.load(fh)
    except Exception as e:
        sys.exit('FATAL: cannot parse %s: %s' % (path, e))
    missing = [k for k in keys if k not in d]
    if missing:
        sys.exit('FATAL: %s is missing %s' % (path, ', '.join(missing)))
    _note(what, path)
    return d


def _exec_module(name, path):
    import importlib.util
    if not os.path.isfile(path):
        sys.exit('FATAL: no such file: %s' % path)
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        sys.exit('FATAL: cannot load %s' % path)
    m = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(m)
    except Exception as e:
        sys.exit('FATAL: cannot exec %s: %s' % (path, e))
    m.SOURCE = path
    return m


def _method_file(name):
    """Path to a method-folder script, or None where the folder is absent."""
    d = os.environ.get('PICODE_METHOD_DIR', DEFAULT_METHOD)
    f = os.path.join(d, 'targets', name)
    return f if os.path.isfile(f) else None


def _lpo_file():
    d = os.environ.get('PICODE_LPO_DIR', DEFAULT_LPO)
    f = os.path.join(d, 'lpo_picode_design.py')
    return f if os.path.isfile(f) else None


def load_lpo(allow_snapshot=True):
    """Burial thresholds and the maximum-ASA table, imported not restated.

    allow_snapshot=False is for the exporter, which must read the method folder
    or fail; everything else accepts the committed snapshot.
    """
    p = _lpo_file()
    if p is None and not allow_snapshot:
        p = lpo_design()          # exits naming the folder it looked in
    if p is None:
        d = _snapshot_json(BURIAL_SNAPSHOT, 'burial thresholds',
                           ('buried_rsa', 'partial_rsa', 'salt_bridge', 'maxasa'))
        m = _Snapshot(BURIAL_SNAPSHOT)
        m.BURIED_RSA = float(d['buried_rsa'])
        m.PARTIAL_RSA = float(d['partial_rsa'])
        m.SALT_BRIDGE = float(d['salt_bridge'])
        m.MAXASA = dict(d['maxasa'])
    else:
        m = _exec_module('lpo_picode_design', p)
    for n in ('BURIED_RSA', 'PARTIAL_RSA', 'SALT_BRIDGE', 'MAXASA'):
        if not hasattr(m, n):
            sys.exit('FATAL: %s does not define %s' % (m.SOURCE, n))
    if len(m.MAXASA) != 20:
        sys.exit('FATAL: expected 20 residues in the maximum-ASA table of %s, '
                 'found %d' % (m.SOURCE, len(m.MAXASA)))
    return m


def memctrl_dir():
    return _dir('PICODE_MEMCTRL_DIR', DEFAULT_MEMCTRL,
                'the folder holding humsavar.txt and the substitution table')


def nontm_substitutions():
    return _file(os.path.join(memctrl_dir(), 'data', 'nontm_substitutions.csv.gz'),
                 'nontm_substitutions.csv.gz', 'PICODE_MEMCTRL_DIR')


def humsavar():
    return _file(os.path.join(memctrl_dir(), 'data', 'humsavar.txt'),
                 'humsavar.txt', 'PICODE_MEMCTRL_DIR')


def project(*parts):
    return os.path.join(ROOT, *parts)


def load_ladder(allow_snapshot=True, path=None):
    """The ladder definitions, from the method folder or the committed snapshot.

    path names ladder_edge_cost.py explicitly and wins over both, which is what
    the --ladder option of the analysis scripts passes through.
    allow_snapshot=False is for the exporter, which must read the method folder
    or fail; everything else accepts the committed snapshot.
    """
    p = path or _method_file('ladder_edge_cost.py')
    if p is None and not allow_snapshot:
        p = ladder_edge_cost()    # exits naming the folder it looked in
    if p is None:
        d = _snapshot_json(LADDER_SNAPSHOT, 'ladder definition',
                           ('core', 'arg', 'ext', 'basicdir', 'charge'))
        m = _Snapshot(LADDER_SNAPSHOT)
        m.CORE = [tuple(e) for e in d['core']]
        m.ARG = [tuple(e) for e in d['arg']]
        m.EXT = [tuple(e) for e in d['ext']]
        m.BASICDIR = set(tuple(e) for e in d['basicdir'])
        m.CHARGE = dict(d['charge'])
    else:
        m = _exec_module('ladder_edge_cost', p)
    for n in ('CORE', 'ARG', 'EXT', 'BASICDIR', 'CHARGE'):
        if not hasattr(m, n):
            sys.exit('FATAL: %s does not define %s; the ladder definition has moved'
                     % (m.SOURCE, n))
    if len(m.CORE) != 12 or len(m.ARG) != 5:
        sys.exit('FATAL: expected 12 core and 5 arginine edges in %s, found %d and %d'
                 % (m.SOURCE, len(m.CORE), len(m.ARG)))
    return m


if __name__ == '__main__':
    print('project root   %s' % ROOT)
    if _method_file('ladder_edge_cost.py') is None:
        L, lpo = load_ladder(), load_lpo()
        print('method folder  not present; running from the committed snapshot')
        print('  ladder definition  %s' % os.path.relpath(L.SOURCE, ROOT))
        print('  burial thresholds  %s' % os.path.relpath(lpo.SOURCE, ROOT))
        print('ladder loaded: %d core, %d arginine, %d proposed extensions'
              % (len(L.CORE), len(L.ARG), len(L.EXT)))
        sys.exit(0)
    print('method folder  %s' % method_dir())
    print('  ladder_edge_cost.py       %s' % ladder_edge_cost())
    print('  ladder_humsavar.py        %s' % ladder_humsavar())
    print('  interface_carboxylates.py %s' % interface_carboxylates())
    print('membrane control folder %s' % memctrl_dir())
    print('  nontm_substitutions.csv.gz %s' % nontm_substitutions())
    print('  humsavar.txt               %s' % humsavar())
    L = load_ladder()
    print('ladder loaded: %d core, %d arginine, %d proposed extensions'
          % (len(L.CORE), len(L.ARG), len(L.EXT)))
