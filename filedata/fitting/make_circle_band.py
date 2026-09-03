#!/usr/bin/env python3
"""Fabricate a sharp CIRCULAR BAND point cloud for rh-adaptive fitting.

z = 1 inside the annulus | |x-c| - R | <= w,  0 outside  (discontinuous).
Deliberately NOT axis-alignable: the optimal mesh is a thin curved annulus of
cells, which a tensor/THB mesh can only approximate by a staircase. This is the
feature lshape_step lacks -- its step follows the boundary of [0.5,1]^2 and is
therefore perfectly axis-aligned.

Same XML layout as filedata/fitting/lshape_step.xml:
  id=0 : 2 x N parameters in [0,1]^2
  id=1 : 3 x N points (x, y, z)
"""
import sys

N   = int(sys.argv[1]) if len(sys.argv) > 1 else 160   # grid points per direction
R   = float(sys.argv[2]) if len(sys.argv) > 2 else 0.30   # band radius
W   = float(sys.argv[3]) if len(sys.argv) > 3 else 0.025  # band HALF-width
out = sys.argv[4] if len(sys.argv) > 4 else "circle_band.xml"
cx = cy = 0.5

us, vs, zs = [], [], []
for j in range(N):
    for i in range(N):
        u = i/(N-1); v = j/(N-1)
        r = ((u-cx)**2 + (v-cy)**2)**0.5
        us.append(u); vs.append(v)
        zs.append(1.0 if abs(r-R) <= W else 0.0)

n = len(us)
fmt = lambda xs: ' '.join(f'{x:.10g}' for x in xs)
with open(out,'w') as f:
    f.write('<?xml version="1.0" encoding="UTF-8"?>\n<xml>\n')
    f.write(f' <Matrix id="0" rows="2" cols="{n}">\n  {fmt(us)}\n  {fmt(vs)}\n </Matrix>\n')
    f.write(f' <Matrix id="1" rows="3" cols="{n}">\n  {fmt(us)}\n  {fmt(vs)}\n  {fmt(zs)}\n </Matrix>\n')
    f.write('</xml>\n')

across = int(2*W*(N-1)) + 1
print(f"N={N}x{N} ({n} pts)  R={R}  half-width={W}  grid spacing={1/(N-1):.5f}")
print(f"samples across the band: ~{across}   points with z=1: {int(sum(zs))} ({100*sum(zs)/n:.1f}%)")
