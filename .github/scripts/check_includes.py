#!/usr/bin/env python3
"""Include-DAG ratchets for the G+Smo module graph.

1. Include ratchet. Scans #include directives in src/<module>/ (all of
   .h/.hpp/.cpp) and compares the resulting module->module edge set
   against the committed baseline cmake/dag_allowed.txt. The baseline
   may only shrink over time:

   - an edge in the scan but not in the baseline  -> ERROR (new coupling)
   - an edge in the baseline but not in the scan  -> NOTE  (run --update
     to shrink the baseline and lock in the improvement)

2. Declared-graph check. Every src/<module>/CMakeLists.txt declares its
   dependencies with gismo_add_module(<module> DEPENDS ...). The public
   headers (.h) of a module may only include modules it depends on,
   transitively. Interface edges outside the declared graph are the
   maintainers' work list and live in cmake/dag_undeclared.txt; that
   list may only shrink, new undeclared edges are an ERROR.

Usage:
    check_includes.py [--src SRC] [--baseline FILE] [--undeclared FILE]
                      [--update]
"""
import argparse
import os
import re
import sys

INC_RE = re.compile(r'^\s*#\s*include\s*[<"]([^">]+)[">]', re.M)
DECL_RE = re.compile(r'gismo_add_module\s*\(\s*(\w+)([^)]*)\)', re.S)
DECL_KEYWORDS = ('DEPENDS', 'SUBDIRS', 'NO_COMPONENT')


def scan(src, exts=('.h', '.hpp', '.cpp')):
    modules = sorted(d for d in os.listdir(src)
                     if os.path.isdir(os.path.join(src, d)))
    edges = set()
    origin = {}
    for mod in modules:
        for dirpath, _, files in os.walk(os.path.join(src, mod)):
            for fname in files:
                if not fname.endswith(exts):
                    continue
                fpath = os.path.join(dirpath, fname)
                with open(fpath, errors='ignore') as f:
                    text = f.read()
                # forwarding headers left behind by relocations only
                # redirect to the new location; they are not a dependency
                # of the old module on the new one
                if '@brief Forwarding header' in text:
                    continue
                for m in INC_RE.finditer(text):
                    parts = m.group(1).split('/')
                    if len(parts) >= 2 and parts[0] in modules \
                            and parts[0] != mod:
                        edge = (mod, parts[0])
                        edges.add(edge)
                        origin.setdefault(edge, os.path.relpath(fpath, src))
    return edges, origin


def declared(src):
    """module -> declared dependencies, from src/<module>/CMakeLists.txt"""
    graph = {}
    for mod in sorted(os.listdir(src)):
        path = os.path.join(src, mod, 'CMakeLists.txt')
        if not os.path.isfile(path):
            continue
        with open(path, errors='ignore') as f:
            text = '\n'.join(line.split('#')[0] for line in f)
        for m in DECL_RE.finditer(text):
            name, deps, key = m.group(1), [], None
            for tok in m.group(2).split():
                if tok in DECL_KEYWORDS:
                    key = tok
                elif key == 'DEPENDS':
                    deps.append(tok)
            graph[name] = deps
    return graph


def closure(graph):
    """module -> set of modules reachable through declared dependencies"""
    reach = {}

    def visit(m, seen):
        for d in graph.get(m, ()):
            if d not in seen:
                seen.add(d)
                visit(d, seen)
        return seen
    for m in graph:
        reach[m] = visit(m, set())
    return reach


def read_edges(path):
    edges = set()
    with open(path) as f:
        for line in f:
            line = line.split('#')[0].strip()
            if not line:
                continue
            a, _, b = line.partition('->')
            edges.add((a.strip(), b.strip()))
    return edges


def write_edges(path, edges, header):
    with open(path, 'w') as f:
        f.write(header)
        for a, b in sorted(edges):
            f.write("%s -> %s\n" % (a, b))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--src', default='src')
    ap.add_argument('--baseline', default='cmake/dag_allowed.txt')
    ap.add_argument('--undeclared', default='cmake/dag_undeclared.txt',
                    help='known interface edges outside the declared graph')
    ap.add_argument('--update', action='store_true',
                    help='rewrite both lists from the current scan')
    args = ap.parse_args()

    edges, origin = scan(args.src)

    # interface (.h) edges outside the declared module graph
    graph = declared(args.src)
    reach = closure(graph)
    iedges, iorigin = scan(args.src, exts=('.h',))
    undeclared = set()
    for a, b in iedges:
        if a in graph and b not in reach[a]:
            undeclared.add((a, b))
    # src/misc is build glue, not a module (see src/CMakeLists.txt)
    missing = sorted({a for a, _ in iedges
                      if a not in graph and a.startswith('gs')})
    if missing:
        print("WARNING: modules without gismo_add_module declaration: %s"
              % ', '.join(missing))

    if args.update:
        write_edges(args.baseline, edges,
                    "# Module include-dependency baseline (ratchet).\n"
                    "# Regenerate with: .github/scripts/check_includes.py"
                    " --update\n"
                    "# This list may only shrink; CI fails on edges not"
                    " listed here.\n")
        print("baseline written: %d edges" % len(edges))
        write_edges(args.undeclared, undeclared,
                    "# Interface (.h) include edges outside the declared"
                    " module graph\n"
                    "# (gismo_add_module DEPENDS in src/<module>/"
                    "CMakeLists.txt): the work list.\n"
                    "# Regenerate with: .github/scripts/check_includes.py"
                    " --update\n"
                    "# This list may only shrink; CI fails on edges not"
                    " listed here.\n")
        print("undeclared-edge list written: %d edges" % len(undeclared))
        return 0

    rc = 0
    allowed = read_edges(args.baseline)
    new = sorted(edges - allowed)
    gone = sorted(allowed - edges)
    for a, b in gone:
        print("NOTE: baseline edge no longer present: %s -> %s "
              "(shrink with --update)" % (a, b))
    if new:
        for a, b in new:
            print("ERROR: new module dependency %s -> %s (first seen in "
                  "src/%s)" % (a, b, origin[(a, b)]))
        print("\n%d new edge(s) violate the include ratchet. If a new "
              "dependency is intentional, discuss it in review and add it "
              "to %s explicitly." % (len(new), args.baseline))
        rc = 1
    else:
        print("include ratchet OK: %d edges within baseline (%d retired)"
              % (len(edges), len(gone)))

    known = read_edges(args.undeclared)
    new = sorted(undeclared - known)
    gone = sorted(known - undeclared)
    for a, b in gone:
        print("NOTE: undeclared interface edge resolved: %s -> %s "
              "(shrink with --update)" % (a, b))
    if new:
        for a, b in new:
            print("ERROR: public header of %s includes %s, which %s does "
                  "not depend on (first seen in src/%s); declare it with "
                  "DEPENDS in src/%s/CMakeLists.txt or do not include it"
                  % (a, b, a, iorigin[(a, b)], a))
        print("\n%d interface edge(s) outside the declared module graph."
              % len(new))
        rc = 1
    else:
        print("declared-graph check OK: %d interface edges, %d known "
              "undeclared (%d resolved)" % (len(iedges), len(undeclared),
                                            len(gone)))
    return rc


if __name__ == '__main__':
    sys.exit(main())
