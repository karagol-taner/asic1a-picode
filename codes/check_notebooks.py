#!/usr/bin/env python3
"""
check_notebooks.py - does every name in these notebooks resolve on a Run All?

    python3 check_notebooks.py md/*.ipynb

Reports three cases, the second being legal Python and the other two not:

  - a name used at module level that no earlier cell has bound
  - a name called inside a function body and defined further down, reported
    only when it is defined nowhere in the notebook
  - a name called by a top-level statement but defined lower in the same cell

Reads the notebooks themselves. Colab line magics are neutralised for the parse
only. Exits non-zero on any unresolved name.
"""
import json, re, ast, builtins, sys, os

def demagic(s):
    s = re.sub(r"^(\s*)!.*$", r"\1pass", s, flags=re.M)
    return re.sub(r"^%.*$", "pass", s, flags=re.M)

def binds(tree):
    out = set()
    for n in ast.walk(tree):
        if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)): out.add(n.name)
        elif isinstance(n, ast.Name) and isinstance(n.ctx, (ast.Store, ast.Del)): out.add(n.id)
        elif isinstance(n, ast.alias): out.add((n.asname or n.name).split(".")[0])
        elif isinstance(n, ast.arg): out.add(n.arg)
        elif isinstance(n, ast.ExceptHandler) and n.name: out.add(n.name)
    return out

def loads(tree, inside_func):
    """Names loaded either at module level or inside a function body."""
    out = []
    def walk(node, depth):
        for ch in ast.iter_child_nodes(node):
            d = depth + 1 if isinstance(ch, (ast.FunctionDef, ast.AsyncFunctionDef,
                                             ast.ClassDef, ast.Lambda)) else depth
            if isinstance(ch, ast.Name) and isinstance(ch.ctx, ast.Load):
                if (depth > 0) == inside_func:
                    out.append((ch.id, ch.lineno))
            walk(ch, d)
    walk(tree, 0)
    return out

def check(path):
    nb = json.load(open(path))
    cells = [(i, demagic("".join(c["source"]))) for i, c in enumerate(nb["cells"])
             if c["cell_type"] == "code"]
    trees = {}
    for i, s in cells:
        try: trees[i] = ast.parse(s)
        except SyntaxError as e:
            print("  SYNTAX cell %d line %s: %s" % (i, e.lineno, e.msg)); return False
    everything = set(dir(builtins)) | {"get_ipython", "__name__", "In", "Out"}
    for i in trees: everything |= binds(trees[i])
    ok, seen = True, set(dir(builtins)) | {"get_ipython", "__name__", "In", "Out"}
    for i, s in cells:
        seen |= binds(trees[i])
        for nm, ln in loads(trees[i], inside_func=False):
            if nm not in seen:
                print("  FATAL  cell %d line %d: %s used at module level before it exists"
                      % (i, ln, nm)); ok = False
        for nm, ln in loads(trees[i], inside_func=True):
            if nm not in everything:
                print("  FATAL  cell %d line %d: %s called from a function but never "
                      "defined in this notebook" % (i, ln, nm)); ok = False
    # A def below an early form field is normal, so only a call that runs
    # too early counts.
    for i, s in cells:
        defs = {n.name: n.lineno for n in trees[i].body
                if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))}
        if not defs: continue
        for n in trees[i].body:
            if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef,
                              ast.Import, ast.ImportFrom)): continue
            for x in ast.walk(n):
                if isinstance(x, ast.Name) and isinstance(x.ctx, ast.Load):
                    if x.id in defs and defs[x.id] > x.lineno:
                        print("  FATAL  cell %d line %d: %s is called here but "
                              "defined at line %d of the same cell"
                              % (i, x.lineno, x.id, defs[x.id])); ok = False
    return ok

allok = True
for p in sys.argv[1:]:
    print("=" * 66); print(os.path.basename(p))
    r = check(p)
    allok &= r
    print("  " + ("every name resolves in Run All order" if r else "SEE ABOVE"))
print()
if allok:
    print("ALL CLEAN")
else:
    raise SystemExit("notebook name resolution failed, see above")
