#!/usr/bin/env python3
"""Fabricate a BAND point cloud (circular or sinusoidal) for rh-adaptive fitting.

A "band" is a thin tube of half-width ``W`` around a decider curve.  The height
field is either a sharp step (z = 1 inside, 0 outside -- discontinuous) or a
smooth Gaussian ridge.  Both are deliberately NOT axis-alignable: the optimal
mesh is a thin *curved* strip of cells, which a tensor/THB mesh can only
approximate by a staircase.  This is the feature ``lshape_step`` lacks -- its
step follows the boundary of [0.5,1]^2 and is therefore perfectly axis-aligned.

This script is the generalisation of ``make_circle_band.py``; the
``circle step`` branch reproduces ``circle_band.xml`` byte-identically.

Usage (all positional)::

    make_band.py <shape> <profile> <N> <GEOM> <W> [OUT] [K]

    shape    circle | sine
    profile  step   | smooth
    N        grid points PER DIRECTION (total N*N samples on [0,1]^2)
    GEOM     circle: radius R about the hard-coded centre (0.5, 0.5)
             sine  : amplitude A of  f(x) = 0.5 + A*sin(2*pi*K*x)
    W        band HALF-width (also the Gaussian sigma of the smooth profile)
    OUT      output file, default "<shape>_band.xml", relative to the CWD
    K        sine periods over [0,1], default 1 (ignored for circle)

Distance to the decider curve
-----------------------------
circle:  d = | |(u,v) - c| - R |                              (exact)
sine  :  d = |v - f(u)| / sqrt(1 + f'(u)^2),  f'(x) = A*2*pi*K*cos(2*pi*K*x)
         This is the first-order (linearised) normal distance: it is exact to
         O(kappa*W) in relative terms.  For the datasets generated here
         kappa*W ~ 3.33*0.025 ~ 0.08, i.e. the band half-width is accurate to
         about 8 %, uniformly and smoothly along the curve -- irrelevant for a
         fitting benchmark, where only the geometry of the feature matters.

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

CX: float = 0.5  # decider-circle centre, hard-coded (as in make_circle_band.py)
CY: float = 0.5


def fmt(xs: Sequence[float]) -> str:
    """Format one Matrix row: 10 significant digits, single-space separated."""
    return ' '.join(f'{x:.10g}' for x in xs)


def usage() -> None:
    sys.exit(__doc__)


def main(argv: List[str]) -> None:
    if len(argv) < 6:
        usage()

    shape: str = argv[1]
    profile: str = argv[2]
    N: int = int(argv[3])
    GEOM: float = float(argv[4])
    W: float = float(argv[5])
    out: str = argv[6] if len(argv) > 6 else f"{shape}_band.xml"
    K: float = float(argv[7]) if len(argv) > 7 else 1.0

    if shape not in ("circle", "sine"):
        usage()
    if profile not in ("step", "smooth"):
        usage()

    us: List[float] = []
    vs: List[float] = []
    zs: List[float] = []
    for j in range(N):
        for i in range(N):
            u = i/(N-1); v = j/(N-1)
            if shape == "circle":
                r = ((u-CX)**2 + (v-CY)**2)**0.5
                d = abs(r-GEOM)
            else:
                fx = 0.5 + GEOM*math.sin(2*math.pi*K*u)
                fp = GEOM*2*math.pi*K*math.cos(2*math.pi*K*u)
                d = abs(v-fx)/(1.0+fp*fp)**0.5
            us.append(u); vs.append(v)
            if profile == "step":
                zs.append(1.0 if d <= W else 0.0)
            else:
                zs.append(math.exp(-d*d/(2.0*W*W)))

    n = len(us)
    with open(out, 'w') as f:
        f.write('<?xml version="1.0" encoding="UTF-8"?>\n<xml>\n')
        f.write(f' <Matrix id="0" rows="2" cols="{n}">\n  {fmt(us)}\n  {fmt(vs)}\n </Matrix>\n')
        f.write(f' <Matrix id="1" rows="3" cols="{n}">\n  {fmt(us)}\n  {fmt(vs)}\n  {fmt(zs)}\n </Matrix>\n')
        f.write('</xml>\n')

    geom_name = "R" if shape == "circle" else "A"
    per = "" if shape == "circle" else f"  K={K:g}"
    across = int(2*W*(N-1)) + 1
    print(f"{shape}/{profile} -> {out}")
    print(f"N={N}x{N} ({n} pts)  {geom_name}={GEOM}{per}  half-width={W}  grid spacing={1/(N-1):.5f}")
    if profile == "step":
        inside = int(sum(zs))
        print(f"samples across the band: ~{across}   points with z=1: {inside}"
              f" ({100.0*inside/n:.4f}%)")
    else:
        zmax = max(zs); zmin = min(zs)
        kmax = zs.index(zmax)
        exact = sum(1 for z in zs if z == 1.0 or z == 0.0)
        print(f"samples across one sigma: ~{across}   z in [{zmin:.6g}, {zmax:.10g}]"
              f"   argmax z at (u,v)=({us[kmax]:.6g}, {vs[kmax]:.6g})")
        print(f"samples with z exactly 1.0 or 0.0: {exact}")


if __name__ == "__main__":
    main(sys.argv)
