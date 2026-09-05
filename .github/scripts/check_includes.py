#!/usr/bin/env python3
"""Include-DAG ratchet for the G+Smo module graph.

Scans #include directives in src/<module>/ and compares the resulting
module->module edge set against the committed baseline
(cmake/dag_allowed.txt). The baseline may only shrink over time:

- an edge in the scan but not in the baseline  -> ERROR (new coupling)
- an edge in the baseline but not in the scan  -> NOTE  (run --update to
  shrink the baseline and lock in the improvement)

Usage:
    check_includes.py [--src SRC] [--baseline FILE] [--update]
"""
import argparse
import os
import re
import sys

INC_RE = re.compile(r'^\s*#\s*include\s*[<"]([^">]+)[">]', re.M)


def scan(src):
    modules = sorted(d for d in os.listdir(src)
                     if os.path.isdir(os.path.join(src, d)))
    edges = set()
    origin = {}
    for mod in modules:
        for dirpath, _, files in os.walk(os.path.join(src, mod)):
            for fname in files:
                if not fname.endswith(('.h', '.hpp', '.cpp')):
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--src', default='src')
    ap.add_argument('--baseline', default='cmake/dag_allowed.txt')
    ap.add_argument('--update', action='store_true',
                    help='rewrite the baseline from the current scan')
    args = ap.parse_args()

    edges, origin = scan(args.src)

    if args.update:
        with open(args.baseline, 'w') as f:
            f.write("# Module include-dependency baseline (ratchet).\n"
                    "# Regenerate with: .github/scripts/check_includes.py"
                    " --update\n"
                    "# This list may only shrink; CI fails on edges not"
                    " listed here.\n")
            for a, b in sorted(edges):
                f.write("%s -> %s\n" % (a, b))
        print("baseline written: %d edges" % len(edges))
        return 0

    allowed = set()
    with open(args.baseline) as f:
        for line in f:
            line = line.split('#')[0].strip()
            if not line:
                continue
            a, _, b = line.partition('->')
            allowed.add((a.strip(), b.strip()))

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
        return 1
    print("include ratchet OK: %d edges within baseline (%d retired)"
          % (len(edges), len(gone)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
