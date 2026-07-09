#!/usr/bin/env python3
"""
Geometry creation and evaluation using the nanobind pygismo bindings.

Run from the build directory:
    python3 geometry_nanobind_example.py

Or point to the build directory manually:
    PYTHONPATH=/tmp/gismo-nb-test python3 geometry_nanobind_example.py
"""

import sys, os
import numpy as np

build_dir = os.environ.get(
    "GISMO_BUILD_DIR",
    os.path.join(os.path.dirname(__file__), "..", "build"),
)
if os.path.isdir(os.path.join(build_dir, "pygismo")):
    sys.path.insert(0, build_dir)

import pygismo as gs

print("pygismo version:", gs.__version__)
print()

# ---------------------------------------------------------------------------
# 1. Knot vectors and B-spline bases
# ---------------------------------------------------------------------------
print("=== Knot vectors and B-spline bases ===")

kv = gs.gsKnotVector(0.0, 1.0, 3, 1, 1, 2)
print("Knot vector:", kv)
print(f"  size={kv.size()}, uSize={kv.uSize()}, first={kv.first()}, last={kv.last()}")

kv.uniformRefine(1)
print("After uniform refine:", kv)

basis = gs.gsBSplineBasis(0.0, 1.0, 3, 2, 1)
print("\nBasis:", basis)
print(
    f"  degree={basis.degree(0)}, ndof={basis.size()}, #elements={basis.numElements()}"
)

u = gs.gsMatrix(1, 1)
u.set(0, 0, 0.5)
vals = basis.eval(u)
print(f"  Basis functions at 0.5: {vals.numpy().flatten()}")
print()

# ---------------------------------------------------------------------------
# 2. Create geometries with gsNurbsCreator
# ---------------------------------------------------------------------------
print("=== Geometry creation ===")

square = gs.gsNurbsCreator.BSplineSquare(2.0, 0.0, 0.0)
print("Unit square [0,2]x[0,2]:", square)
print("  parDim:", square.parDim(), " geoDim:", square.geoDim())
print("  coefs:\n", square.coefs().numpy())

cube = gs.gsNurbsCreator.BSplineCube()
print("\nUnit cube:", cube)
print("  parDim:", cube.parDim(), " geoDim:", cube.geoDim())
print()

# ---------------------------------------------------------------------------
# 3. Evaluate geometry and derivatives
# ---------------------------------------------------------------------------
print("=== Evaluation ===")

pts = gs.gsMatrix(2, 3)
pts.set(0, 0, 0.0)
pts.set(1, 0, 0.0)
pts.set(0, 1, 0.5)
pts.set(1, 1, 0.5)
pts.set(0, 2, 1.0)
pts.set(1, 2, 1.0)

pts_np = pts.numpy()
vals = square.eval(pts)
jac = square.jacobian(pts)

print("Parameter points:\n", pts_np)
print("Physical points (eval):\n", vals.numpy())
print("Jacobian at last point:\n", jac.numpy())
print()

# ---------------------------------------------------------------------------
# 4. Refine the geometry
# ---------------------------------------------------------------------------
print("=== Refinement ===")

sq = gs.gsNurbsCreator.BSplineSquare()
print("Before refine:  ndof =", sq.basis().size(), " #el =", sq.basis().numElements())
sq.uniformRefine(2)
print("After 2x refine: ndof =", sq.basis().size(), " #el =", sq.basis().numElements())
sq.degreeElevate(1)
print("After p-elevate:  deg =", sq.basis().degree(0), " ndof =", sq.basis().size())
print()

# ---------------------------------------------------------------------------
# 5. Manipulate coefficients
# ---------------------------------------------------------------------------
print("=== Coefficient manipulation ===")

geo = gs.gsNurbsCreator.BSplineSquare()
print("Original coefs:\n", geo.coefs().numpy())

geo.coefs().set(2, 1, 0.8)
geo.coefs().set(3, 1, 1.2)
print("Modified coefs (via .set on reference):\n", geo.coefs().numpy())

pt = gs.gsMatrix(2, 1)
pt.set(0, 0, 0.5)
pt.set(1, 0, 1.0)
print(f"eval at (0.5, 1.0) = {geo.eval(pt).numpy().flatten()}")
print()

# ---------------------------------------------------------------------------
# 6. Assemble into a gsMultiPatch
# ---------------------------------------------------------------------------
print("=== MultiPatch ===")

mp = gs.gsMultiPatch()
mp.addPatch(gs.gsNurbsCreator.BSplineSquare())
mp.addPatch(gs.gsNurbsCreator.BSplineSquare(1.0, 1.0, 0.0))
print(f"MultiPatch: {mp.nPatches()} patches")
for i in range(mp.nPatches()):
    p = mp.patch(i)
    print(
        f"  Patch {i}: parDim={p.parDim()}, geoDim={p.geoDim()}, ndof={p.basis().size()}"
    )

mb = gs.gsMultiBasis(mp)
print("MultiBasis created from MultiPatch")
print()

# ---------------------------------------------------------------------------
# 7. Function expressions
# ---------------------------------------------------------------------------
print("=== Function expressions ===")

f = gs.gsFunctionExpr("x^2 + y^2", 2)
print("f = x^2 + y^2")

p1 = gs.gsMatrix(2, 1)
p1.set(0, 0, 1.0)
p1.set(1, 0, 1.0)
print(f"  f(1,1)   = {f.eval(p1).numpy().flatten()}")

p0 = gs.gsMatrix(2, 1)
p0.set(0, 0, 0.0)
p0.set(1, 0, 0.0)
print(f"  f(0,0)   = {f.eval(p0).numpy().flatten()}")
print(f"  grad(1,1)= {f.deriv(p1).numpy().flatten()}")
print()

# ---------------------------------------------------------------------------
# 8. Knot insertion on a B-spline curve
# ---------------------------------------------------------------------------
print("=== Knot insertion ===")

kv_curve = gs.gsKnotVector(0.0, 1.0, 0, 1, 1, 2)
print("Curve knot vector:", kv_curve)

b = gs.gsBSplineBasis(0.0, 1.0, 0, 2, 1)
print("Curve basis before insert:", b.knots())
b.uniformRefine(2)
print("Curve basis after 2x refine:", b.knots())
print(f"  degree={b.degree(0)}, ndof={b.size()}")
print()

print("Done!")
