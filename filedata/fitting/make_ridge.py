#!/usr/bin/env python3
"""Fabricate a RIDGE point cloud for validating anisotropic direction indicators.

A "ridge" is a tanh transition that varies in exactly ONE parametric direction
and is EXACTLY constant in the other. This makes the ground-truth refinement
direction known a priori and one-sided: unlike ``lshape_step`` (whose L-shaped
front mixes both directions ~50/50), a ridge gives an unambiguous pass/fail
target for a direction indicator under test.

    ridge_x: z = tanh(sigma*(u - 0.5))   -- varies in u ONLY, constant in v
    ridge_y: z = tanh(sigma*(v - 0.5))   -- the mirror, varies in v ONLY

``tanh`` is monotone, so the signed and unsigned (``|.|``-integrated) element
derivative integrals never disagree from cancellation -- that confound is
removed by construction.

Usage (all positional, with defaults matching the task's reference run)::

    make_ridge.py [axis] [N] [sigma] [OUT]

    axis   x | y            default x
    N      grid points PER DIRECTION (total N*N samples on [0,1]^2)   default 120
    sigma  tanh steepness; transition width ~ 4/sigma                 default 20
    OUT    output file, default "ridge_<axis>.xml", relative to the CWD

I/O format
----------
The output is the two-Matrix G+Smo XML consumed by the fitting examples, the
same layout as ``filedata/fitting/lshape_step.xml``::

    <Matrix id="0" rows="2" cols="N*N">   parameters (u; v), each row N*N long
    <Matrix id="1" rows="3" cols="N*N">   points     (x; y; z) = (u; v; z)

Row-major reading: row k of a Matrix is one coordinate of ALL samples.  The
samples run in the order  j (v, slow) outer, i (u, fast) inner, with
u = i/(N-1), v = j/(N-1).  Values are printed with 10 significant digits.
Stdlib only; no randomness -- the output is fully reproducible.
"""
from __future__ import annotations

import math
import sys
from typing import List, Sequence


def fmt(xs: Sequence[float]) -> str:
    """Format one Matrix row: 10 significant digits, single-space separated."""
    return ' '.join(f'{x:.10g}' for x in xs)


def usage() -> None:
    sys.exit(__doc__)


def main(argv: List[str]) -> None:
    axis: str = argv[1] if len(argv) > 1 else "x"
    N: int = int(argv[2]) if len(argv) > 2 else 120
    sigma: float = float(argv[3]) if len(argv) > 3 else 20.0
    out: str = argv[4] if len(argv) > 4 else f"ridge_{axis}.xml"

    if axis not in ("x", "y"):
        usage()

    us: List[float] = []
    vs: List[float] = []
    zs: List[float] = []
    for j in range(N):
        for i in range(N):
            u = i/(N-1); v = j/(N-1)
            t = u if axis == "x" else v
            z = math.tanh(sigma*(t-0.5))
            us.append(u); vs.append(v); zs.append(z)

    n = len(us)
    with open(out, 'w') as f:
        f.write('<?xml version="1.0" encoding="UTF-8"?>\n<xml>\n')
        f.write(f' <Matrix id="0" rows="2" cols="{n}">\n  {fmt(us)}\n  {fmt(vs)}\n </Matrix>\n')
        f.write(f' <Matrix id="1" rows="3" cols="{n}">\n  {fmt(us)}\n  {fmt(vs)}\n  {fmt(zs)}\n </Matrix>\n')
        f.write('</xml>\n')

    width = 4.0/sigma
    print(f"ridge_{axis} -> {out}")
    print(f"N={N}x{N} ({n} pts)  axis={axis}  sigma={sigma:g}  transition width~{width:.4f}"
          f"  grid spacing={1/(N-1):.5f}")
    print(f"z range: [{min(zs):.6g}, {max(zs):.6g}]")


if __name__ == "__main__":
    main(sys.argv)
