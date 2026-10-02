/** @file immersed_gauss_green_quadrature_example.cpp

    @brief Stand-alone Gauss-Green cut-cell quadrature for a 2D domain
    described by closed loops of exact boundary curves (no level set), with
    method-independent oracle checks and an Algoim adaptive-rule baseline.

    Theory (Sommariva & Vianello 2009; Gunderman, Weiss & Evans 2021): let
    Omega be bounded by one outer loop and any number of non-overlapping,
    non-nested hole loops, oriented so that Omega always lies on the left of
    the direction of travel (outer loop CCW, hole loops CW). For a background
    cell C = [xL,xR] x [yB,yT] and F(x,y) = int_{xL}^{x} f(s,y) ds, the
    divergence theorem applied to the vector field (F,0) gives

      int_{C ^ Omega} f dA = int_{d(C ^ Omega)} F n_x ds
                           = sum_pieces int F(x(t),y(t)) y'(t) dt
                             + int_{right edge ^ Omega} F(xR,y) dy,

    because F == 0 identically on the left edge x = xL and n_x == 0 on the
    top/bottom edges, and the outward normal of a boundary piece traversed
    with Omega on the left is (y',-x')/|c'| so n_x ds = y' dt (Sommariva &
    Vianello's "dimension reduction": a 2D cut-cell integral becomes a sum of
    1D integrals over the trimming curve pieces plus one 1D integral over the
    cell's right edge). The piece-length (boundary) rule is the ordinary arc-
    length integral over the same pieces: ds = |c'(t)| dt, normal
    (y',-x')/|c'| (Omega on the left).

    Cell classification and cell assignment (no point-in-polygon; symbolic
    perturbation, no tolerance): geometry that aligns exactly with the grid
    -- a curve through a grid vertex, a break or polyline vertex on a knot
    line, a curve segment lying exactly along a knot line -- is resolved by
    treating every vertical line x = X_i as X_i + eps and every horizontal
    line y = Y_j as Y_j + eps (eps infinitesimal), so a point exactly ON a
    line belongs to the column/row to its LEFT/BELOW. The two exact
    predicates right(p,i) := p.x > X_i and above(p,j) := p.y > Y_j drive
    every topological decision:
    - a crossing of a curve with line X_i, between two consecutive samples
      ta < tb, exists iff right(c(ta),i) != right(c(tb),i); its parameter is
      found by bisection ON THE PREDICATE (not on x(t)-X_i with a
      tolerance), and its direction is the flip direction (left-to-right =
      +1, right-to-left = -1), not the sign of x'(t*) -- at a break or a
      polyline vertex the tangent is discontinuous and x'(t*) can be zero or
      carry the wrong one-sided sign;
    - collecting every vertical line's crossings (y = y(t*), direction),
      sorting by y, grouping exact ties and accumulating a winding number
      w(y) = sum of directions of crossings below y (Omega on the left of
      the oriented boundary => a +1 crossing enters, so w in {0,1}
      identically) yields that line's "inside" y-intervals, which double as
      the right-edge integration domain of the two cells adjacent to it;
    - a curve is split into pieces at every crossing parameter (vertical and
      horizontal) plus its own parametric breaks. Cell assignment is by a
      tracker, not by re-evaluating the curve at a piece's own midpoint (that
      second, independent evaluation is not exact -- x(t) == a can wiggle by
      an ULP away from a on a straight piece, let alone a rational one):
      each loop's boundary starts in the cell of its first curve's own start
      point (shared, by construction, with the previous curve's end point --
      every curve-to-curve junction is evaluated ONCE per loop and reused by
      both curves it borders); every detected crossing of an INTERIOR grid
      line then moves the running cell by its flip direction (never a box
      line x = x0 or x = x0+n*h, which only feeds the winding walk above),
      and split parameters that coincide exactly have their moves summed, so
      simultaneous vertical and horizontal crossings (a curve through a grid
      vertex) move the tracker diagonally in one step, not in an order-
      dependent pair of steps. A curve segment lying exactly along X_i is
      not guaranteed to be crossing-free: its evaluated x(t) can round an
      ULP to either side of X_i, and wherever that happens AT a sample the
      predicate flips and a crossing is recorded. Every recorded flip
      enters the tracker and the winding walk alike, and such flips come in
      pairs, so the two stay consistent. Between recorded flips the tracker
      keeps the column it carried onto the segment, which is the LEFT
      column where the samples give x == X_i, and the segment's F.y' term
      supplies (or cancels) that column's right-edge integral up to
      O(eps*|x|*h) rounding. The curve's arithmetic path can still round
      x(t) to a value an ULP off X_i BETWEEN samples, but that wiggle is
      invisible to both
      the crossing search and the winding walk, which only ever read the
      predicate AT samples and AT bisected crossing parameters -- so it can
      never make the tracker and the classification disagree with each
      other, even where it disagrees with the true continuous curve at the
      ULP level.
    - an uncut cell (no pieces) is Full iff the midpoint of its right edge
      lies in the inside set of that edge's line, else Empty -- unchanged by
      the above, since it never touches the predicate at a grid line itself.

    Exactness: on a straight piece the volume rule integrates x^a y^b exactly
    whenever a <= 2p+1 and a+b+1 <= 2*nq-1 (inner Gauss p+1 in x, outer Gauss
    nq in t), i.e. nq >= p+1 suffices for every moment with a+b <= 2p; on a
    curved (rational) piece the rule converges spectrally in nq once pieces
    are split at every parametric break of the curve.

    Traps documented, not handled:
    - exactly one outer loop plus any number of non-overlapping, non-nested
      hole loops (no support for touching or overlapping loops);
    - Omega must lie strictly inside the background box (no boundary piece
      may touch x0, x0+n*h, y0 or y0+n*h);
    - near-tangency of a loop's boundary with a grid line (detectNearTangencies:
      any x- or y-extremum of the loop, at a sample point, a break, or a
      junction between two consecutive curves, within 1e-8*h of a grid line)
      is only detected and reported (gsWarn "near-tangency"), never resolved
      -- the caller must change the background grid (--n0 or -r) so the
      extremum no longer sits that close to a line; a break or vertex that
      sits exactly ON a line with both neighbours on the same side of it is
      a legitimate touch, not a near-tangency, and is not reported;
    - a genuinely missed crossing pair inside one sample interval (the curve
      crosses the same grid line twice between two consecutive samples) aborts
      with GISMO_ENSURE instead -- the caller must raise --nsamp for that case;
    - quadrature weights can be negative (sign of y' on the curve pieces),
      which is correct, not a bug;
    - the Algoim baseline (V4) is skipped automatically (a gsWarn plus
      "Algoim baseline skipped" on stdout, exit code unaffected) whenever a
      square edge lies exactly on a background grid line (any theta with
      sin(theta)*cos(theta) ~ 0, and +-half within 1e-12 of some grid line):
      the level set then vanishes identically on a whole cell face, and
      gsAlgoimAdaptiveRule::mapTo enumerates every sub-box down to maxDepth
      without ever classifying one as uncut, exhausting memory -- a
      gsAlgoim adapter limitation, not a defect in this driver.

    Geometry sources (--source hand|xml|occ), all producing the same
    std::vector<Loop> the algorithm above consumes -- every topological
    decision (predicates, shared junctions, tracker, winding walk) stays
    inside gaussGreenQuadrature regardless of source:
    - hand (default): the DefaultGeometry square/disk built directly as
      G+Smo splines (handBuiltLoops), the reference geometry the V5 (occ)
      and V6 (xml) checks compare against.
    - xml (-f file.xml, required): a <PlanarDomain> read via gsFileData,
      every gsCurve of every gsCurveLoop wrapped as an exact spline
      evaluation (fromGeometry) -- no polyline conversion. --write-xml
      writes the hand-built DefaultGeometry as a <PlanarDomain>; the
      committed filedata/planar/square_minus_disk.xml was generated this
      way and is the file V6 reads back.
    - occ (-f file.step|.stp|.brep, optional): either the built-in
      square-minus-disk boolean (BRepAlgoAPI_Cut of a polygon face and a
      disk face, occDefaultShape) or a CAD file read by OCCT
      (STEPControl_Reader / BRepTools::Read). Either way the shape's single
      planar face's wires become loops of BoundaryCurves that evaluate the
      OCC edge geometry directly (BRepAdaptor_Curve::D1) -- again no
      conversion to a G+Smo spline. --write-cad writes the shape actually
      used (built-in or read) to a .brep/.step/.stp file.
    --assume-default asserts that an -f geometry IS DefaultGeometry's own
    numbers (the same --theta/--half the run was given): it enables V2 in
    spline mode, V5 (occ vs hand) / V6 (xml vs hand) -- the source's GG rule
    against the hand-built rule at nq = 12 -- and the Algoim baseline, all
    of which compare against DefaultGeometry's
    exact oracle or level set. Without it (the default for any -f run) V2
    and the Algoim baseline each print a "... skipped: geometry from -f is
    not known to be DefaultGeometry (pass --assume-default)" line instead
    of running, V5/V6 are not run at all (no line is printed), and V1/V3
    (which need no oracle beyond the loops' own moments) still run and must
    pass.

    Further traps, specific to xml/occ:
    - occ needs the shape to reduce to exactly one planar face in a plane
      parallel to XY (occLoops); a multi-face shape or a curved face aborts.
    - every loop, from any source, must be watertight: curve k's own fresh
      evaluation at t1 must match curve k+1's own start point to
      1e-10*(n*h) (the existing ENSURE in gaussGreenQuadrature). For xml
      this means the file's <CurveLoop> curves must already run head-to-
      tail; for occ this holds by construction once every TopAbs_REVERSED
      edge (BRepTools_WireExplorer order) is fed through reverseCurve.
      ensureInsideBox aborts ("outside the background box") before the GG
      algorithm runs at all if any sampled boundary point (the same nsamp
      samples the crossing search uses) lies outside the box -- the most
      common way to trip this is a foreign file: OCCT reads STEP files in
      millimetres by default, so a model authored in metres lands two to
      three orders of magnitude outside a box sized for --n0/-r.
    - gsPlanarDomain's constructor treats loop 0 as the outer loop and
      reorients any loop whose control-polygon shoelace area disagrees, via
      gsCurve::reverse() -- which for a NURBS curve reaches the
      unimplemented gsRationalBasis::reverse (GISMO_NO_IMPLEMENTATION) and
      aborts inside G+Smo's own reader, not this driver. writeDefaultXml
      sidesteps this by writing the hole loop already CW (explicit
      knot/weight/control-point reversal), so neither writing nor reading
      the committed file ever exercises it; a hand-authored file whose
      NURBS loop is misoriented the other way will still hit it.
    - an OCC gp_Circ edge carries no breaks (fromOccEdge only extracts knots
      from a GeomAbs_BSplineCurve), so its x/y extrema fall inside a sample
      interval rather than on a break, exercising the interior-extremum /
      missed-crossing-pair search in gaussGreenQuadrature. At the default
      --nsamp 64 this is not an issue for DefaultGeometry's disk; at a much
      finer grid the missed-pair ENSURE ("raise --nsamp") may legitimately
      ask for more samples.

    Mesh source (--source mesh): Omega is a planar triangle mesh (-f
    file.stl/.off/.obj, read via gsSurfMesh/gsReadSurfMesh) integrated
    EXACTLY, cell by cell, by clipping -- no Gauss-Green identity, and it
    does not go through gaussGreenQuadrature's own crossing search, tracker
    or winding walk (no in/out test). Every triangle is clipped against
    every background cell box it overlaps (Sutherland-Hodgman, the four
    half-planes x>=xL, x<=xR, y>=yB, y<=yT), fan-triangulated from the
    clipped polygon's own first vertex, and integrated with a collapsed
    (Duffy) Gauss rule (2p+1 points per collapsed direction, exact to
    degree 4p+1, so every checked moment a+b <= 2p is exact); every
    boundary edge is clipped the same way (Liang-Barsky) and integrated
    with a straight-segment Gauss rule (p+1 points). A Full cell KEEPS its
    clipped rule rather than being replaced by a tensor Gauss rule: the
    clipped rule is already exact to the stated degree, the Full label is
    purely informational, and replacing it would let a sliver below the
    classification tolerance silently change the integral -- the cost is a
    Full cell carrying (#fan triangles overlapping it)*(2p+1)^2 nodes
    instead of the (p+1)^2 a tensor rule would use. Complexity:
    O(#triangles x (n + overlapped cells per triangle)) for the volume,
    plus O(#boundary edges x (n + overlapped cells per edge)) for the
    boundary, n = cells per direction (the cell range is located by
    linear scans over the grid lines); each clip is O(1) (a triangle
    clipped by 4 half-planes has at most 7 vertices). Trap: the
    gsSurfMesh binary-STL reader is float/double-
    broken (it freads sizeof(Point) = 24 bytes per vertex where a binary
    STL stores 12 bytes of float32), so every shipped .stl here is ASCII
    (gmsh "Mesh.Binary = 0"), and a vertex line over the reader's 100-char
    fgets buffer is silently truncated. --polyline does not apply to
    --source mesh (rejected). The Algoim baseline is always skipped for
    --source mesh (gsWarn plus "Algoim baseline skipped" on stdout, exit
    code unaffected): the mesh is a faceted approximation of an arbitrary
    -f mesh, while the Algoim level set describes DefaultGeometry, so a
    per-cell diff would measure mesh discretisation error, not quadrature
    error.

    Example command lines:
      ./immersed_gauss_green_quadrature_example
      ./immersed_gauss_green_quadrature_example -r 3 --polyline 64 --no-algoim
      ./immersed_gauss_green_quadrature_example -r 1 --plot -o output_gg
      ./immersed_gauss_green_quadrature_example --theta 0 --half 0.5 -r 1 --no-algoim
        (square edges exactly on grid lines whenever --n0 is a multiple of 4;
        the Algoim baseline would be skipped automatically even without
        --no-algoim, see the trap above)
      ./immersed_gauss_green_quadrature_example --source xml -f filedata/planar/square_minus_disk.xml --assume-default --no-algoim
      ./immersed_gauss_green_quadrature_example --source occ --no-algoim
      ./immersed_gauss_green_quadrature_example --source mesh -r 2
      ./immersed_gauss_green_quadrature_example --source mesh -f planar/square_on_grid_minus_disk_mesh.stl --theta 0 --half 0.5 -r 1

    References:
    - A. Sommariva, M. Vianello, "Gauss-Green cubature and moment computation
      over arbitrary geometries", J. Comput. Appl. Math. 231 (2009) 886-896.
    - D. Gunderman, K. Weiss, J.A. Evans, "Spectral mesh-free quadrature for
      planar regions bounded by rational parametric curves", Comput.-Aided
      Des. 130 (2021) 102944.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <gismo.h>
#include <gsAlgoim/gsAlgoimAdaptiveRule.h>

#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepTools.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <GeomAbs_CurveType.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <Geom_BSplineCurve.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <STEPControl_Reader.hxx>
#include <STEPControl_Writer.hxx>
#include <TopAbs_Orientation.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Ax2.hxx>
#include <gp_Circ.hxx>
#include <gp_Dir.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace gismo;

namespace {

std::string fmtSci(real_t v, int prec = 6)
{
    std::ostringstream os;
    os << std::scientific << std::setprecision(prec) << v;
    return os.str();
}

/// One boundary curve c: [t0,t1] -> R^2. eval(t, x, dx): t is 1xN, x and dx
/// are 2xN (points and d/dt tangents). breaks: interior parameters where c
/// is only finitely smooth (spline knots, OCC continuity intervals); every
/// break is a mandatory split point for the crossing search below.
struct BoundaryCurve
{
    real_t t0, t1;
    std::vector<real_t> breaks;
    std::function<void(const gsMatrix<real_t>&, gsMatrix<real_t>&, gsMatrix<real_t>&)> eval;
};
typedef std::vector<BoundaryCurve> Loop;   // closed, consecutive curves, Omega on the left

/// Uniform background grid: n x n cells, cell (i,j) = [x0+i*h, x0+(i+1)*h] x
/// [y0+j*h, y0+(j+1)*h], linear id = i + n*j.
struct Grid { real_t x0, y0, h; index_t n; };

enum CellStatus { Empty = -1, Cut = 0, Full = 1 };

struct CellRule    { gsMatrix<real_t> nodes; gsVector<real_t> weights; };        // 2 x m, m
struct CellBdrRule { gsMatrix<real_t> nodes; gsVector<real_t> weights;
                     gsMatrix<real_t> normals; };                                // 2 x m, m, 2 x m (unit outward)

struct CutCellQuadrature
{
    Grid grid;
    std::vector<int>         status;   // CellStatus per cell id
    std::vector<CellRule>    vol;      // size n*n; empty rule for Empty cells
    std::vector<CellBdrRule> bdr;      // size n*n; empty rule for cells without boundary
};

/// Diagnostics filled by gaussGreenQuadrature.
struct GGStats
{
    index_t nVerticalCrossings   = 0;
    index_t nHorizontalCrossings = 0;
    index_t nPieces              = 0;
    index_t nCut                 = 0;
    index_t nFull                = 0;
    index_t nEmpty               = 0;
    index_t nNegativeWeights     = 0;
    index_t nNearTangency        = 0;
};

/// Default hand-built geometry: a square of half-side \a a rotated CCW by
/// \a thetaDeg about the origin, minus a disk of radius \a R centred at
/// (\a cx,\a cy). Shared by handBuiltLoops(), the exact moment oracle and
/// the Algoim level set (AlgoimLevelSet), so all three describe the same
/// numbers. \a a and \a thetaDeg are exposed on the command line (--half,
/// --theta) so an "edge-along-grid-line" configuration (--theta 0 --half
/// 0.5) can be requested without touching the disk; the default constructor
/// keeps the original defaults, 0.7 and 17 degrees.
struct DefaultGeometry
{
    real_t a        = 0.7;
    real_t thetaDeg = 17.0;
    real_t theta    = thetaDeg * (real_t)EIGEN_PI / (real_t)180.0;
    real_t ct       = math::cos(theta);
    real_t st       = math::sin(theta);
    real_t cx = 0.05, cy = -0.03, R = 0.35;
    gsMatrix<real_t> corners;   // 4 x 2, rows P0..P3, CCW

    DefaultGeometry() : corners(4,2) { buildCorners(); }

    DefaultGeometry(real_t aIn, real_t thetaDegIn)
    : a(aIn), thetaDeg(thetaDegIn),
      theta(thetaDeg * (real_t)EIGEN_PI / (real_t)180.0),
      ct(math::cos(theta)), st(math::sin(theta)), corners(4,2)
    { buildCorners(); }

    real_t exactArea()      const { return 4.0*a*a - (real_t)EIGEN_PI*R*R; }
    real_t exactPerimeter() const { return 8.0*a + 2.0*(real_t)EIGEN_PI*R; }

private:
    void buildCorners()
    {
        const real_t s[4][2] = { {-a,-a}, {a,-a}, {a,a}, {-a,a} };
        for (int k = 0; k != 4; ++k)
        {
            corners(k,0) = ct*s[k][0] - st*s[k][1];
            corners(k,1) = st*s[k][0] + ct*s[k][1];
        }
    }
};

/// Wraps a gsGeometry curve (kept alive via the shared_ptr capture) as a
/// BoundaryCurve. t0/t1 and the interior breaks come from \a kv
/// (gsKnotVector::breaks() returns the unique knots of the domain including
/// both ends; the interior entries are the mandatory split points).
BoundaryCurve fromGeometry(memory::shared_ptr<gsGeometry<real_t> > g,
                           const gsKnotVector<real_t> & kv)
{
    BoundaryCurve c;
    const std::vector<real_t> br = kv.breaks();
    c.t0 = br.front();
    c.t1 = br.back();
    if (br.size() > 2)
        c.breaks.assign(br.begin()+1, br.end()-1);
    c.eval = [g](const gsMatrix<real_t> & t, gsMatrix<real_t> & x, gsMatrix<real_t> & dx)
    {
        g->eval_into(t, x);
        g->deriv_into(t, dx);
    };
    return c;
}

/// Straight segment A -> B, t in [0,1], no breaks.
BoundaryCurve straightSegment(const gsVector<real_t,2> & A, const gsVector<real_t,2> & B)
{
    BoundaryCurve c;
    c.t0 = 0.0; c.t1 = 1.0;
    c.eval = [A,B](const gsMatrix<real_t> & t, gsMatrix<real_t> & x, gsMatrix<real_t> & dx)
    {
        const index_t n = t.cols();
        const gsVector<real_t,2> D = B - A;
        x.resize(2,n); dx.resize(2,n);
        for (index_t k = 0; k != n; ++k)
        {
            x(0,k) = A[0] + t(0,k)*D[0];
            x(1,k) = A[1] + t(0,k)*D[1];
            dx(0,k) = D[0]; dx(1,k) = D[1];
        }
    };
    return c;
}

/// The default hand-built geometry: a square (4 degree-1 B-spline segments,
/// CCW, outer loop) minus a disk (a single full-circle NURBS curve, CCW as
/// returned by gsNurbsCreator -- normalizeOrientation below must reverse it
/// to make it a CW hole; that reversal is the point of using it as-is here).
std::vector<Loop> handBuiltLoops(const DefaultGeometry & geo)
{
    std::vector<Loop> loops;

    Loop square;
    for (int k = 0; k != 4; ++k)
    {
        gsMatrix<real_t> C(2,2);
        C(0,0) = geo.corners(k,0);         C(0,1) = geo.corners(k,1);
        C(1,0) = geo.corners((k+1)%4,0);   C(1,1) = geo.corners((k+1)%4,1);
        gsBSpline<real_t> seg(0.0, 1.0, 0, 1, C);
        const gsKnotVector<real_t> kv = seg.knots();
        memory::shared_ptr<gsGeometry<real_t> > sp(new gsBSpline<real_t>(seg));
        square.push_back(fromGeometry(sp, kv));
    }
    loops.push_back(square);

    Loop circle;
    gsNurbsCreator<real_t>::NurbsPtr circ = gsNurbsCreator<real_t>::NurbsCircle(geo.R, geo.cx, geo.cy);
    const gsKnotVector<real_t> kv = circ->knots();
    memory::shared_ptr<gsGeometry<real_t> > sp(circ.release());
    circle.push_back(fromGeometry(sp, kv));
    loops.push_back(circle);

    return loops;
}

/// Replaces every curve of every loop by N straight segments sampled at
/// c(t0 + k(t1-t0)/N), k = 0..N-1. Each loop's vertex list (2 x M, M =
/// N*#curves) is returned via \a vertices (fed to the polygon oracle, never
/// to the GG rule); the new loops close exactly by construction, independent
/// of any floating-point mismatch between the original curve endpoints.
std::vector<Loop> toPolyline(const std::vector<Loop> & loops, index_t N,
                             std::vector<gsMatrix<real_t> > & vertices)
{
    std::vector<Loop> out;
    vertices.clear();
    for (const Loop & loop : loops)
    {
        const index_t M = N * (index_t)loop.size();
        gsMatrix<real_t> V(2, M);
        index_t col = 0;
        for (const BoundaryCurve & c : loop)
        {
            gsMatrix<real_t> t(1, N);
            for (index_t k = 0; k != N; ++k)
                t(0,k) = c.t0 + (real_t)k * (c.t1 - c.t0) / (real_t)N;
            gsMatrix<real_t> x, dx;
            c.eval(t, x, dx);
            V.block(0, col, 2, N) = x;
            col += N;
        }
        vertices.push_back(V);

        Loop newLoop;
        for (index_t m = 0; m != M; ++m)
        {
            const gsVector<real_t,2> A = V.col(m);
            const gsVector<real_t,2> B = V.col((m+1)%M);
            newLoop.push_back(straightSegment(A, B));
        }
        out.push_back(newLoop);
    }
    return out;
}

/// Sorted, deduplicated sample set for the crossing search AND for
/// globalMoment: uniform samples t0 + k(t1-t0)/nsamp, k=0..nsamp, union the
/// curve's own breaks. The same nsamp is used for every curve.
std::vector<real_t> sampleSet(const BoundaryCurve & c, index_t nsamp)
{
    std::vector<real_t> S;
    S.reserve(nsamp + 1 + c.breaks.size());
    for (index_t k = 0; k <= nsamp; ++k)
        S.push_back(c.t0 + (real_t)k * (c.t1 - c.t0) / (real_t)nsamp);
    S.insert(S.end(), c.breaks.begin(), c.breaks.end());
    std::sort(S.begin(), S.end());
    S.erase(std::unique(S.begin(), S.end()), S.end());
    return S;
}

/// Evaluates \a c at the single parameter \a t.
void evalCurve1(const BoundaryCurve & c, real_t t,
                real_t & x, real_t & y, real_t & dx, real_t & dy)
{
    gsMatrix<real_t> tm(1,1); tm(0,0) = t;
    gsMatrix<real_t> xm, dxm;
    c.eval(tm, xm, dxm);
    x = xm(0,0); y = xm(1,0); dx = dxm(0,0); dy = dxm(1,0);
}

/// Bisection ON THE PREDICATE \a pred (no tolerance, no value comparison),
/// given the CALLER'S OWN predicate values \a loVal, \a hiVal at the bracket
/// ends: \a pred is evaluated only at interior midpoints, never at \a lo or
/// \a hi again, since a caller whose bracket ends are shared junction points
/// (evaluated once and reused by two curves, see gaussGreenQuadrature) may
/// legitimately have a fresh re-evaluation there disagree at the ULP level
/// with the stored value the rest of the algorithm is using -- re-deriving
/// loVal/hiVal here would just reintroduce that inconsistency one call
/// later. The invariant loVal != hiVal is maintained every iteration by
/// moving whichever of lo/hi has the mid's predicate value to mid; the
/// returned parameter is always \a hi, so it consistently sits on the side
/// of the flip that the caller's \a hi endpoint started on. Terminates after
/// at most 64 iterations, or sooner once mid stops landing strictly inside
/// (lo,hi) -- the ULP limit of \a real_t.
real_t bisectPredicate(const std::function<bool(real_t)> & pred, real_t lo, real_t hi,
                       bool loVal, bool hiVal)
{
    GISMO_ENSURE(loVal != hiVal, "bisectPredicate: bracket has no predicate flip.");
    for (int it = 0; it != 64; ++it)
    {
        const real_t mid = 0.5*(lo+hi);
        if (!(mid > lo) || !(mid < hi)) break;
        if (pred(mid) == hiVal) hi = mid; else lo = mid;
    }
    return hi;
}

/// Column index of a point with x-coordinate \a x, given the grid's own
/// precomputed vertical line positions \a X (size n+1, X[i] = x0+i*h): the
/// number of interior lines X[1]..X[n-1] with x > X[i] -- the exact
/// predicate right(p,i), never a floor or a tolerance, so a point exactly ON
/// a line (x == X[i]) always belongs to the LEFT column (i-1). \a X must be
/// the SAME array the crossing search compares against (not a fresh
/// x0+i*h at each call site): grid.x0+i*h is a rounded sum, and two call
/// sites that each round it independently can disagree by an ULP at large i
/// or non-dyadic h, which would silently break the tracker's end-of-curve
/// and end-of-loop consistency ENSUREs below.
index_t colOf(real_t x, const std::vector<real_t> & X)
{
    index_t c = 0;
    for (size_t i = 1; i + 1 < X.size(); ++i)
        if (x > X[i]) ++c;
    return c;
}

/// Row index of a point with y-coordinate \a y, same convention as colOf.
index_t rowOf(real_t y, const std::vector<real_t> & Y)
{
    index_t c = 0;
    for (size_t j = 1; j + 1 < Y.size(); ++j)
        if (y > Y[j]) ++c;
    return c;
}

/// Plain bisection on f (no derivative available: f is itself a derivative
/// of the curve, and no second derivative is provided). Used only by the
/// interior-extremum / missed-crossing-pair search below, which is
/// unrelated to the predicate-based crossing search: it locates an x'(resp.
/// y') root, not a grid-line crossing, so it still needs a genuine value
/// (sign of the derivative), not a predicate.
real_t bisectRoot(const std::function<real_t(real_t)> & f, real_t ta, real_t tb)
{
    const real_t eps = std::numeric_limits<real_t>::epsilon();
    real_t va = f(ta), vb = f(tb);
    if (0.0 == va) return ta;
    if (0.0 == vb) return tb;
    GISMO_ENSURE((va > 0) != (vb > 0), "bisectRoot: bracket has no sign change.");

    real_t t = 0.5*(ta+tb);
    while (tb - ta > 4.0*eps*math::max((real_t)1.0, math::max(math::abs(ta), math::abs(tb))))
    {
        t = 0.5*(ta+tb);
        const real_t v = f(t);
        if (0.0 == v) return t;
        if ((v > 0) == (va > 0)) { ta = t; va = v; }
        else                     { tb = t; vb = v; }
    }
    return 0.5*(ta+tb);
}

bool strictlyBetween(real_t v, real_t p, real_t q)
{ return v > math::min(p,q) && v < math::max(p,q); }

/// Reversal of a curve: t -> -t (not t0+t1-t), tangent negated, breaks
/// negated and re-sorted ascending. Negation, unlike t0+t1-t, is exact in
/// floating point and injective on every representable double, so a
/// parameter found by bisection on the reversed curve maps back to exactly
/// the original parameter it was found at. t0+t1-t is not injective near
/// t = (t0+t1)/2 when doubles below the midpoint are twice as dense as
/// those above it (true whenever t0+t1 has this form, in particular for the
/// [0,1]-parametrised circle used by handBuiltLoops): two adjacent
/// candidates t and t+eps on one side of the midpoint can both round to the
/// same t0+t1-t on the other side, so a bisection converging to within 1
/// ULP of the true root, on the reversed curve, could evaluate the ORIGINAL
/// curve one or more ULP away from where the search actually terminated --
/// silently reintroducing exactly the precision loss the predicate-based
/// crossing search is designed to avoid. The lambda captures a copy of the
/// original eval function: capturing \a c by reference would dangle into
/// the loop that reverseLoop() overwrites with the reversed curves.
BoundaryCurve reverseCurve(const BoundaryCurve & c)
{
    BoundaryCurve r;
    r.t0 = -c.t1; r.t1 = -c.t0;
    r.breaks.resize(c.breaks.size());
    for (size_t i = 0; i != c.breaks.size(); ++i)
        r.breaks[c.breaks.size()-1-i] = -c.breaks[i];
    std::sort(r.breaks.begin(), r.breaks.end());

    const std::function<void(const gsMatrix<real_t>&, gsMatrix<real_t>&, gsMatrix<real_t>&)> oldEval = c.eval;
    r.eval = [oldEval](const gsMatrix<real_t> & t, gsMatrix<real_t> & x, gsMatrix<real_t> & dx)
    {
        gsMatrix<real_t> tt(t.rows(), t.cols());
        for (index_t k = 0; k != t.cols(); ++k) tt(0,k) = -t(0,k);
        oldEval(tt, x, dx);
        dx = -dx;
    };
    return r;
}

void reverseLoop(Loop & loop)
{
    Loop rev(loop.size());
    const size_t n = loop.size();
    for (size_t i = 0; i != n; ++i)
        rev[i] = reverseCurve(loop[n-1-i]);
    loop = rev;
}

/// Reads a <PlanarDomain> from \a path and turns every curve of every loop
/// into a BoundaryCurve via fromGeometry(), so V6 exercises the exact
/// spline evaluation path (no polyline conversion). Trap: the file's curves
/// must already be head-to-tail within each loop (the driver's own
/// watertightness ENSURE in gaussGreenQuadrature is what catches a violation).
std::vector<Loop> xmlLoops(const std::string & path)
{
    gsFileData<real_t> fd(path);
    memory::unique_ptr<gsPlanarDomain<real_t> > pd = fd.getFirst<gsPlanarDomain<real_t> >();
    GISMO_ENSURE(pd, "no <PlanarDomain> found in " << path);

    std::vector<Loop> loops;
    for (int i = 0; i != pd->numLoops(); ++i)
    {
        const gsCurveLoop<real_t> & cl = pd->loop(i);
        Loop loop;
        for (int j = 0; j != cl.numCurves(); ++j)
        {
            const gsCurve<real_t> & c = cl.curve(j);
            memory::shared_ptr<gsGeometry<real_t> > sp(c.clone().release());
            gsKnotVector<real_t> kv;
            if (const gsBSpline<real_t> * bs = dynamic_cast<const gsBSpline<real_t> *>(&c))
                kv = bs->knots();
            else if (const gsNurbs<real_t> * nb = dynamic_cast<const gsNurbs<real_t> *>(&c))
                kv = nb->knots();
            else
                GISMO_ERROR("xmlLoops: loop " << i << ", curve " << j << " in " << path
                            << " is neither a gsBSpline nor a gsNurbs.");
            loop.push_back(fromGeometry(sp, kv));
        }
        loops.push_back(loop);
    }
    return loops;
}

/// Writes the hand-built DefaultGeometry as a <PlanarDomain> to \a fn: the
/// outer loop as 4 degree-1 gsBSpline segments (identical construction to
/// handBuiltLoops), the hole loop as a single full-circle gsNurbs written
/// CW by explicit knot/weight/control-point reversal. Writing the hole CW
/// (rather than relying on gsPlanarDomain's own reorientation) matters
/// because gsPlanarDomain reorients loop i by calling gsCurve::reverse(),
/// which for a NURBS curve reaches the unimplemented gsRationalBasis::reverse
/// (GISMO_NO_IMPLEMENTATION) -- a G+Smo limitation this writer sidesteps so
/// that both writing and reading back never trigger it.
void writeDefaultXml(const DefaultGeometry & geo, const std::string & fn)
{
    std::vector<gsCurve<real_t> *> outerCurves;
    for (int k = 0; k != 4; ++k)
    {
        gsMatrix<real_t> C(2,2);
        C(0,0) = geo.corners(k,0);         C(0,1) = geo.corners(k,1);
        C(1,0) = geo.corners((k+1)%4,0);   C(1,1) = geo.corners((k+1)%4,1);
        outerCurves.push_back(new gsBSpline<real_t>(0.0, 1.0, 0, 1, C));
    }
    gsCurveLoop<real_t> * outer = new gsCurveLoop<real_t>(outerCurves);

    gsNurbsCreator<real_t>::NurbsPtr circ = gsNurbsCreator<real_t>::NurbsCircle(geo.R, geo.cx, geo.cy);
    gsKnotVector<real_t> kv = circ->knots();
    kv.reverse();
    gsMatrix<real_t> w = circ->weights().colwise().reverse();
    gsMatrix<real_t> C = circ->coefs().colwise().reverse();
    gsCurve<real_t> * holeCurve = new gsNurbs<real_t>(kv, w, C);
    gsCurveLoop<real_t> * hole = new gsCurveLoop<real_t>(holeCurve);

    gsPlanarDomain<real_t> pd(std::vector<gsCurveLoop<real_t>*>{outer, hole});

    gsFileData<real_t> fd;
    fd.setFloatPrecision(17);
    fd.add(pd);
    fd.save(fn);
    gsInfo << "wrote DefaultGeometry as <PlanarDomain> to " << fn << "\n";
}

/// Built-in OCC shape for the default geometry: a square face (planar wire
/// through the four DefaultGeometry corners) minus a disk face (circular
/// edge of radius geo.R centred at (geo.cx,geo.cy)), via BRepAlgoAPI_Cut.
TopoDS_Shape occDefaultShape(const DefaultGeometry & geo)
{
    BRepBuilderAPI_MakePolygon poly;
    for (int k = 0; k != 4; ++k)
        poly.Add(gp_Pnt(geo.corners(k,0), geo.corners(k,1), 0.0));
    poly.Close();
    gp_Pln plane(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(0.0, 0.0, 1.0));
    TopoDS_Face square = BRepBuilderAPI_MakeFace(plane, poly.Wire(), Standard_True).Face();

    gp_Circ circ(gp_Ax2(gp_Pnt(geo.cx, geo.cy, 0.0), gp_Dir(0.0, 0.0, 1.0)), geo.R);
    TopoDS_Edge ce = BRepBuilderAPI_MakeEdge(circ).Edge();
    TopoDS_Wire cw = BRepBuilderAPI_MakeWire(ce).Wire();
    TopoDS_Face disk = BRepBuilderAPI_MakeFace(plane, cw, Standard_True).Face();

    BRepAlgoAPI_Cut cut(square, disk);
    GISMO_ENSURE(cut.IsDone() && !cut.HasErrors(), "BRepAlgoAPI_Cut failed");
    return cut.Shape();
}

/// Reads a CAD shape from \a path, dispatching on the lower-cased extension
/// (brep, or step/stp via STEPControl_Reader). IFSelect_RetDone is the only
/// acceptable status -- IFSelect_ReturnStatus is RetVoid=0, RetDone=1,
/// RetError, RetFail, RetStop, so a truthiness test on the return value
/// silently accepts RetError/RetFail/RetStop as success.
TopoDS_Shape occReadShape(const std::string & path)
{
    std::string ext = gsFileManager::getExtension(path);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

    TopoDS_Shape s;
    if ("brep" == ext)
    {
        BRep_Builder b;
        GISMO_ENSURE(BRepTools::Read(s, path.c_str(), b), "occReadShape: failed to read " << path);
    }
    else if ("step" == ext || "stp" == ext)
    {
        STEPControl_Reader rd;
        GISMO_ENSURE(rd.ReadFile(path.c_str()) == IFSelect_RetDone,
                    "occReadShape: STEPControl_Reader::ReadFile failed on " << path);
        rd.TransferRoots();
        s = rd.OneShape();
    }
    else
        GISMO_ERROR("occReadShape: unsupported extension '" << ext << "' for " << path);

    GISMO_ENSURE(!s.IsNull(), "occReadShape: " << path << " produced a null shape");
    return s;
}

/// Wraps one OCC edge \a e (face-adjacent z0 already checked by the caller)
/// as a BoundaryCurve, evaluating the underlying curve directly via
/// BRepAdaptor_Curve::D1 -- no conversion to a G+Smo spline. The handle is
/// captured BY VALUE in the lambda (BRepAdaptor_Curve is a handle class in
/// OCCT 7.6, DEFINE_STANDARD_HANDLE), so the adaptor stays alive as long as
/// the BoundaryCurve does. BRepAdaptor_Curve applies the edge's placement
/// but ignores its TopAbs_Orientation; the caller reverses a REVERSED edge
/// with reverseCurve() instead.
BoundaryCurve fromOccEdge(const TopoDS_Edge & e, real_t z0)
{
    Handle(BRepAdaptor_Curve) ad = new BRepAdaptor_Curve(e);
    GISMO_ENSURE(ad->Is3DCurve(), "fromOccEdge: edge has no 3D curve");

    BoundaryCurve c;
    c.t0 = ad->FirstParameter();
    c.t1 = ad->LastParameter();

    if (GeomAbs_BSplineCurve == ad->GetType())
    {
        Handle(Geom_BSplineCurve) bs = ad->BSpline();
        std::vector<real_t> br;
        for (int i = 1; i <= bs->NbKnots(); ++i)
        {
            const real_t ki = bs->Knot(i);
            if (ki > c.t0 && ki < c.t1) br.push_back(ki);
        }
        std::sort(br.begin(), br.end());
        br.erase(std::unique(br.begin(), br.end()), br.end());
        c.breaks = br;
    }

    {
        gp_Pnt P0; gp_Vec V0; ad->D1(c.t0, P0, V0);
        gp_Pnt P1; gp_Vec V1; ad->D1(c.t1, P1, V1);
        GISMO_ENSURE(math::abs(P0.Z()-z0) <= 1e-9 && math::abs(P1.Z()-z0) <= 1e-9,
                    "fromOccEdge: edge is not at z = const (z0 = " << z0 << ")");
    }

    c.eval = [ad](const gsMatrix<real_t> & t, gsMatrix<real_t> & x, gsMatrix<real_t> & dx)
    {
        const index_t n = t.cols();
        x.resize(2,n); dx.resize(2,n);
        for (index_t k = 0; k != n; ++k)
        {
            gp_Pnt P; gp_Vec V;
            ad->D1(t(0,k), P, V);
            x(0,k) = P.X(); x(1,k) = P.Y();
            dx(0,k) = V.X(); dx(1,k) = V.Y();
        }
    };
    return c;
}

/// Turns the single planar face of \a shape into loops of BoundaryCurves
/// evaluating the OCC edge geometry directly. \a shape must contain exactly
/// one face (checked via TopExp::MapShapes' de-duplicating map), and that
/// face must lie in a plane parallel to the XY plane (checked via
/// BRepAdaptor_Surface::Plane()). Each wire's edges are visited in
/// BRepTools_WireExplorer order; a TopAbs_REVERSED edge is reversed with
/// reverseCurve() -- the driver's own orientation normalisation runs
/// afterwards regardless, so the face's own orientation is irrelevant here.
std::vector<Loop> occLoops(const TopoDS_Shape & shape)
{
    TopTools_IndexedMapOfShape faces;
    TopExp::MapShapes(shape, TopAbs_FACE, faces);
    GISMO_ENSURE(faces.Extent() == 1, "OCC source needs exactly one face, found " << faces.Extent());
    TopoDS_Face face = TopoDS::Face(faces(1));

    BRepAdaptor_Surface sf(face);
    GISMO_ENSURE(sf.GetType() == GeomAbs_Plane, "OCC source needs a planar face");
    gp_Pln pl = sf.Plane();
    GISMO_ENSURE(math::abs(pl.Axis().Direction().Z()) >= 1 - 1e-12,
                "OCC source: face is not parallel to the XY plane");
    const real_t z0 = pl.Location().Z();

    TopoDS_Wire outer = BRepTools::OuterWire(face);
    std::vector<TopoDS_Wire> wires;
    wires.push_back(outer);
    for (TopExp_Explorer ex(face, TopAbs_WIRE); ex.More(); ex.Next())
    {
        const TopoDS_Wire w = TopoDS::Wire(ex.Current());
        if (!w.IsSame(outer)) wires.push_back(w);
    }

    std::vector<Loop> loops;
    for (const TopoDS_Wire & wire : wires)
    {
        Loop loop;
        for (BRepTools_WireExplorer we(wire, face); we.More(); we.Next())
        {
            const TopoDS_Edge & e = we.Current();
            if (BRep_Tool::Degenerated(e)) continue;
            BoundaryCurve c = fromOccEdge(e, z0);
            if (TopAbs_REVERSED == we.Orientation()) c = reverseCurve(c);
            loop.push_back(c);
        }
        loops.push_back(loop);
    }
    return loops;
}

/// Writes \a s to \a fn, dispatching on the lower-cased extension: brep via
/// the 2-argument BRepTools::Write alias (full precision), step/stp via
/// STEPControl_Writer. Both write calls are checked against IFSelect_RetDone
/// explicitly -- see occReadShape for why a truthiness test is unsafe.
void writeCad(const TopoDS_Shape & s, const std::string & fn)
{
    std::string ext = gsFileManager::getExtension(fn);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

    if ("brep" == ext)
        GISMO_ENSURE(BRepTools::Write(s, fn.c_str()), "writeCad: failed to write " << fn);
    else if ("step" == ext || "stp" == ext)
    {
        STEPControl_Writer w;
        GISMO_ENSURE(w.Transfer(s, STEPControl_AsIs) == IFSelect_RetDone,
                    "writeCad: STEPControl_Writer::Transfer failed for " << fn);
        GISMO_ENSURE(w.Write(fn.c_str()) == IFSelect_RetDone,
                    "writeCad: STEPControl_Writer::Write failed for " << fn);
    }
    else
        GISMO_ERROR("writeCad: unsupported extension '" << ext << "' for " << fn);
    gsInfo << "wrote OCC shape to " << fn << "\n";
}

/// Single geometry switch point: every geometry source is dispatched here.
/// \a source has already been validated by main; GISMO_ERROR is a defensive
/// fallback, not a user-facing path. \a geo carries --theta/--half (and the
/// fixed disk parameters) into whichever source builds the loops. For
/// "occ", \a occShape receives the shape actually used (built-in or read
/// from \a file), so main can pass it on to --write-cad.
std::vector<Loop> makeLoops(const std::string & source, const DefaultGeometry & geo,
                            const std::string & file, TopoDS_Shape * occShape)
{
    if ("hand" == source) return handBuiltLoops(geo);
    if ("xml"  == source) return xmlLoops(file);
    if ("occ"  == source)
    {
        *occShape = file.empty() ? occDefaultShape(geo) : occReadShape(file);
        return occLoops(*occShape);
    }
    GISMO_ERROR("Unknown geometry source '" << source << "'.");
}

/// Neumaier's variant of Kahan compensated summation (the branch on
/// |m_sum| >= |x| chooses which of the two operands the rounding error is
/// extracted against, so the compensation term \a m_c is valid regardless of
/// which addend is larger). A plain sequential sum over the O(10^4-10^5)
/// quadrature nodes of a fine background mesh has rounding error that grows
/// with the term count (each addition rounds independently, and the error
/// bound scales with n*eps*sum(|x_i|)); Neumaier's compensation reduces that
/// to O(n*eps^2*sum(|x_i|)), which is what keeps the moment checks
/// (V1/V3a/V3d) inside their 1e-13 threshold at the finer refinement levels.
class KahanSum
{
public:
    void add(real_t x)
    {
        const real_t t = m_sum + x;
        m_c += (math::abs(m_sum) >= math::abs(x)) ? ((m_sum-t)+x) : ((x-t)+m_sum);
        m_sum = t;
    }
    real_t value() const { return m_sum + m_c; }
private:
    real_t m_sum = 0, m_c = 0;
};

/// Global Gauss-Green moment via the SAME divergence-theorem identity used
/// by the volume rule, applied to the whole loop (no knot-line splitting --
/// only \a c.breaks split a curve here), so V3a below is an independent
/// check of the crossing search and cell classification:
///   int_Omega x^a y^b dA = oint (x^(a+1)/(a+1)) y^b y'(t) dt.
real_t globalMoment(const std::vector<Loop> & loops, int a, int b, index_t nsamp)
{
    KahanSum total;
    gsGaussRule<real_t> gauss20(20);
    for (const Loop & loop : loops)
    for (const BoundaryCurve & c : loop)
    {
        const std::vector<real_t> S = sampleSet(c, nsamp);
        for (size_t i = 0; i+1 < S.size(); ++i)
        {
            gsMatrix<real_t> nodes; gsVector<real_t> w;
            gauss20.mapTo(S[i], S[i+1], nodes, w);
            gsMatrix<real_t> x, dx;
            c.eval(nodes, x, dx);
            for (index_t q = 0; q != nodes.cols(); ++q)
            {
                const real_t xa1 = std::pow(x(0,q), a+1) / (real_t)(a+1);
                const real_t yb  = std::pow(x(1,q), b);
                total.add(w[q] * xa1 * yb * dx(1,q));
            }
        }
    }
    return total.value();
}

/// Signed area of each loop via globalMoment(.,0,0,.) (A = oint x y' dt):
/// the loop with the largest |A| is the outer loop and must have A > 0; all
/// others are holes and must have A < 0. Returns the number of loops
/// reversed and prints it, unless \a verbose is false (used for the hidden
/// hand-built reference of the V5/V6 checks, so that a non-"hand" run still
/// prints exactly one "orientation:" line, its own).
index_t normalizeOrientation(std::vector<Loop> & loops, index_t nsamp, bool verbose = true)
{
    std::vector<real_t> A(loops.size());
    for (size_t i = 0; i != loops.size(); ++i)
    {
        std::vector<Loop> single(1, loops[i]);
        A[i] = globalMoment(single, 0, 0, nsamp);
    }
    size_t outer = 0;
    for (size_t i = 1; i != loops.size(); ++i)
        if (math::abs(A[i]) > math::abs(A[outer])) outer = i;

    index_t nReversed = 0;
    for (size_t i = 0; i != loops.size(); ++i)
    {
        const bool shouldBePositive = (i == outer);
        if ((shouldBePositive && A[i] < 0) || (!shouldBePositive && A[i] > 0))
        {
            reverseLoop(loops[i]);
            ++nReversed;
        }
    }
    if (verbose) gsInfo << "orientation: reversed " << nReversed << " loop(s)\n";
    return nReversed;
}

std::vector<std::pair<real_t,real_t> > clipIntervals(
    const std::vector<std::pair<real_t,real_t> > & ivals, real_t lo, real_t hi)
{
    std::vector<std::pair<real_t,real_t> > out;
    for (const auto & iv : ivals)
    {
        const real_t a = math::max(iv.first, lo);
        const real_t b = math::min(iv.second, hi);
        if (b > a) out.push_back(std::make_pair(a,b));
    }
    return out;
}

bool isInside(const std::vector<std::pair<real_t,real_t> > & ivals, real_t y)
{
    for (const auto & iv : ivals)
        if (y >= iv.first && y <= iv.second) return true;
    return false;
}

/// One split point of a curve (a break, a domain endpoint, or a grid-line
/// crossing), carrying the tracker move that happens there: dcol/drow are
/// nonzero only for a crossing of an INTERIOR grid line (a box line, x = x0
/// or x = x0+n*h, moves nothing -- see colOf/rowOf). Two records at exactly
/// the same t (e.g. a curve through a grid vertex, where a vertical and a
/// horizontal crossing coincide) are summed into one before the tracker
/// walks them, never merged by dropping either.
struct SplitRecord { real_t t; int dcol; int drow; };

/// A crossing of a horizontal background grid line Y_j, keyed by its own
/// split parameter t (see the closing-column comment at gaussGreenQuadrature).
struct HCrossing { real_t t, Yj, x, y; };

/// A closing-column correction attached to one endpoint of a piece: the
/// piece's curve point (x,y) at that endpoint overshoots the row line Yj it
/// was meant to reach by a few ULP of y, so the F.dy contribution of the
/// thin strip between y and Yj (width x-xL, the piece's own cell) must be
/// added at sb (sign +1) or removed at sa (sign -1) to keep the divergence
/// identity exact at the row boundary.
struct HCloseTerm { real_t Yj, x, y; int sign; };

/// One sub-interval [sa,sb] of a curve, assigned to background cell (ci,cj),
/// plus any closing-column corrections at its endpoints.
struct PieceRef { const BoundaryCurve * curve; real_t sa, sb; std::vector<HCloseTerm> close; };

/// A crossing of a vertical background grid line: y-coordinate and the flip
/// direction (+1 left-to-right, -1 right-to-left) used by the winding walk
/// below -- never the sign of x'(t*), which is discontinuous exactly at the
/// breaks and polyline vertices that most often coincide with a grid line.
struct VCrossing { real_t y; int dir; };

/// In/out walk on one vertical grid line's crossings (symbolic-perturbation
/// classification, see the file header): sorts by y, groups EXACT ties
/// (a curve through a grid vertex, or several curves meeting a line at the
/// same y) and applies each group's summed direction at once, so that
/// classification never depends on the order two simultaneous crossings
/// happen to be discovered in. The winding number w is asserted in {0,1}
/// after every group -- Omega is on the left of the oriented boundary, so a
/// properly closed loop can never wind twice over the same line -- and the
/// "inside" y-intervals (w == 1) double as the right-edge integration
/// domain of the two cells adjacent to this line.
std::vector<std::pair<real_t,real_t> > windingInsideIntervals(
    std::vector<VCrossing> & list, index_t lineIndex, real_t linePos)
{
    std::sort(list.begin(), list.end(),
             [](const VCrossing & a, const VCrossing & b){ return a.y < b.y; });
    std::vector<std::pair<real_t,real_t> > inside;
    int w = 0;
    bool haveStart = false;
    real_t startY = 0;
    size_t k = 0;
    const size_t m = list.size();
    while (k != m)
    {
        const real_t y = list[k].y;
        int delta = 0;
        size_t k2 = k;
        while (k2 != m && list[k2].y == y) { delta += list[k2].dir; ++k2; }
        const int wOld = w;
        w += delta;
        GISMO_ENSURE(0 == w || 1 == w,
                    "vertical line i="<<lineIndex<<" (x="<<linePos<<"): winding number "<<w
                    <<" out of {0,1} at y="<<y<<" -- Omega is not a simple, correctly oriented "
                      "region there.");
        if (0 == wOld && 1 == w) { startY = y; haveStart = true; }
        else if (1 == wOld && 0 == w)
        {
            GISMO_ENSURE(haveStart, "internal: winding walk closed an interval it never opened.");
            inside.push_back(std::make_pair(startY, y));
            haveStart = false;
        }
        k = k2;
    }
    GISMO_ENSURE(0 == w, "vertical line i="<<lineIndex<<" (x="<<linePos<<"): winding number "<<w
                <<" nonzero at the end of the line -- the boundary does not close on it.");
    return inside;
}

/// Flags every local x- or y-extremum of \a loop (including at a curve's own
/// break, at a junction between two consecutive curves -- a polyline vertex
/// or a square corner -- or exactly at a sample point) that lies within
/// 1e-8*h of a background grid line. Works directly on the loop's sampled
/// (x,y) values (never on the derivative), by concatenating every curve's
/// sampleSet() evaluation into one cyclic sequence (dropping the point each
/// curve shares with its neighbour, and the point the last curve shares with
/// the first, since a closed loop repeats it) and flagging position m
/// whenever the secant slopes (px[m]-px[m-1]) and (px[m+1]-px[m]) do not
/// share a strict sign -- which also catches a slope that is exactly zero at
/// m, the case a same-curve x'(t)==0 test at the extremum's own root cannot
/// see when the extremum sits exactly on a sample point or break (as every
/// extremum of the hand-built circle does).
void detectNearTangencies(const Loop & loop, const Grid & grid, index_t nsamp, GGStats & stats)
{
    std::vector<real_t> px, py;
    for (const BoundaryCurve & c : loop)
    {
        const std::vector<real_t> S = sampleSet(c, nsamp);
        gsMatrix<real_t> tS(1,(index_t)S.size());
        for (size_t k = 0; k != S.size(); ++k) tS(0,(index_t)k) = S[k];
        gsMatrix<real_t> Xs, DXs;
        c.eval(tS, Xs, DXs);
        const size_t start = px.empty() ? 0 : 1;
        for (size_t k = start; k != S.size(); ++k) { px.push_back(Xs(0,k)); py.push_back(Xs(1,k)); }
    }
    if (px.size() > 1) { px.pop_back(); py.pop_back(); }
    const size_t M = px.size();
    if (M < 3) return;

    const index_t n = grid.n;
    for (size_t m = 0; m != M; ++m)
    {
        const size_t mPrev = (m + M - 1) % M, mNext = (m + 1) % M;
        const real_t dx1 = px[m]-px[mPrev], dx2 = px[mNext]-px[m];
        // dx1 == dx2 == 0 identically is a sample strictly inside a run of
        // collinear points (a curve segment lying exactly along a vertical
        // grid line, or a polyline edge parallel to one): every interior
        // sample there has this property, none of them is a turning point,
        // so flagging them all would drown a legitimate edge-along-line
        // configuration in near-tangency noise without locating anything
        // new (the segment's own two endpoints, where the direction
        // actually changes, are still flagged normally).
        if (dx1*dx2 <= 0.0 && !(0.0 == dx1 && 0.0 == dx2))
            for (index_t i = 0; i <= n; ++i)
            {
                const real_t Xi = grid.x0 + i*grid.h;
                if (math::abs(px[m]-Xi) < 1e-8*grid.h)
                {
                    gsWarn << "near-tangency: loop extremum/vertex at ("<<px[m]<<", "<<py[m]
                           <<") with x = X_"<<i<<" ("<<Xi<<")\n";
                    ++stats.nNearTangency;
                }
            }
        const real_t dy1 = py[m]-py[mPrev], dy2 = py[mNext]-py[m];
        if (dy1*dy2 <= 0.0 && !(0.0 == dy1 && 0.0 == dy2))
            for (index_t j = 0; j <= n; ++j)
            {
                const real_t Yj = grid.y0 + j*grid.h;
                if (math::abs(py[m]-Yj) < 1e-8*grid.h)
                {
                    gsWarn << "near-tangency: loop extremum/vertex at ("<<px[m]<<", "<<py[m]
                           <<") with y = Y_"<<j<<" ("<<Yj<<")\n";
                    ++stats.nNearTangency;
                }
            }
    }
}

/// Builds the Gauss-Green volume + boundary rules for cell ^ Omega on every
/// cell of \a grid, from the exact boundary curves \a loops.
///
/// Volume rule, cell C = [xL,xR] x [yB,yT], F(x,y) = int_{xL}^x f(s,y) ds:
///   int_{C^Omega} f dA = sum_pieces int F(x(t),y(t)) y'(t) dt
///                       + int_{right edge ^ Omega} F(xR,y) dy.
/// Inner Gauss p+1 maps [-1,1] to [0,1] via xi_k=(1+z_k)/2, eta_k=w_k/2:
///   piece [ta,tb]: outer Gauss nq point (t_q,W_q); node
///     (xL + (x_q-xL)*xi_k, y_q), weight W_q * y'(t_q) * (x_q-xL) * eta_k;
///   right-edge interval [ya,yb]: outer Gauss nq point (y_q,W_q); node
///     (xL + (xR-xL)*xi_k, y_q), weight W_q * (xR-xL) * eta_k -- xR-xL, the
///     cell's own box width, not the nominal grid spacing, so this term is
///     self-consistent with the box V3b checks against.
/// Boundary rule, per piece: node c(t_q), weight W_q*|c'(t_q)|, normal
/// (y',-x')/|c'| at t_q (Omega on the left).
/// Every node lies in the cell's closed box; weights may be negative (sign
/// of y' on the piece contributions -- correct, not a bug).
/// Exactness: nq >= p+1 integrates every moment x^a y^b, a+b <= 2p, exactly
/// on straight pieces; on curved pieces the rule converges spectrally in nq
/// once pieces are split at every parametric break.
///
/// Algorithm (symbolic perturbation throughout -- see the file header for
/// the predicates right(p,i)/above(p,j), col(p)/row(p)):
///   1. near-tangency scan (detectNearTangencies), once per loop;
///   2. crossing search: for every curve and every background grid line,
///      a crossing is a flip of right(.,i) (resp. above(.,j)) between two
///      consecutive samples, refined by bisectPredicate (never by comparing
///      x(t)-X_i to a tolerance); its direction is the flip direction, not
///      the sign of x'(t*); vertical crossings (y, direction) are collected
///      per line for the winding walk, and every crossing parameter --
///      vertical, horizontal, or a curve's own break -- becomes a split
///      point of that curve;
///   3. cell assignment by a tracker, not by re-evaluating pm and asking
///      col(pm)/row(pm) (that second evaluation is not exact -- see the
///      file header): each loop starts at the cell of its first curve's
///      start point, shared by construction with the previous curve's end
///      point (every junction between two consecutive curves of a loop is
///      evaluated once and reused by both); every crossing of an interior
///      grid line then moves the running cell by its flip direction, split
///      parameters that coincide exactly have their moves summed, and each
///      piece between two consecutive split points is assigned the cell the
///      tracker holds right after the split point that opens it. A segment
///      lying exactly along X_i is not guaranteed to be crossing-free: its
///      evaluated x(t) can round an ULP to either side of X_i, and wherever
///      that happens AT a sample the predicate flips and a crossing is
///      recorded; every recorded flip enters the tracker and the winding
///      walk alike, and such flips come in pairs, so the two stay
///      consistent. Between recorded flips the tracker keeps the column it
///      carried onto the segment -- the LEFT column where the samples give
///      x == X_i -- and the segment's F.y' term supplies (or cancels) that
///      column's right-edge integral up to O(eps*|x|*h) rounding. A wiggle
///      of x(t) by an ULP BETWEEN samples is invisible to both the crossing
///      search and the winding walk, so it can never make the two disagree
///      with each other, even where it disagrees with the true curve;
///   4. cell classification from windingInsideIntervals, one call per
///      vertical line;
///   5. per cut cell: inner Gauss p+1 in x, outer Gauss nq on each piece
///      and on the right-edge-inside-Omega sub-intervals; full cells get a
///      plain tensor Gauss (p+1)^2 rule; empty cells get no rule.
CutCellQuadrature gaussGreenQuadrature(const std::vector<Loop> & loops, const Grid & grid,
                                       index_t p, index_t nq, index_t nsamp, GGStats & stats)
{
    stats = GGStats();
    const index_t n = grid.n;

    // Precomputed grid-line positions, read by EVERY call site that needs
    // X_i/Y_j (the crossing search, colOf/rowOf, windingInsideIntervals):
    // grid.x0+i*grid.h is a rounded sum, and two call sites that each round
    // it independently can disagree by an ULP at large i or non-dyadic h --
    // which would break the tracker's own end-of-curve/end-of-loop
    // consistency ENSUREs below, since those compare a tracker total built
    // from crossing tests against a fresh colOf/rowOf call.
    std::vector<real_t> X(n+1), Y(n+1);
    for (index_t i = 0; i <= n; ++i) X[i] = grid.x0 + i*grid.h;
    for (index_t j = 0; j <= n; ++j) Y[j] = grid.y0 + j*grid.h;

    CutCellQuadrature Q;
    Q.grid = grid;
    Q.status.assign((size_t)n*n, (int)Empty);
    Q.vol.assign((size_t)n*n, CellRule());
    Q.bdr.assign((size_t)n*n, CellBdrRule());

    for (const Loop & loop : loops)
        detectNearTangencies(loop, grid, nsamp, stats);

    // Inner rule: Gauss p+1 on [-1,1] mapped to [0,1] -- xi_k = (1+z_k)/2,
    // eta_k = w_k/2.
    gsGaussRule<real_t> innerRule(p+1);
    const gsMatrix<real_t> & zref = innerRule.referenceNodes();
    const gsVector<real_t> & wref = innerRule.referenceWeights();
    std::vector<real_t> xi(p+1), eta(p+1);
    for (index_t k = 0; k <= p; ++k) { xi[k] = 0.5*(1.0+zref(0,k)); eta[k] = 0.5*wref[k]; }

    gsGaussRule<real_t> outerRule(nq);
    gsGaussRule<real_t> fullRule(gsVector<index_t,2>::Constant(p+1));

    std::vector<std::vector<VCrossing> > vcross(n+1);
    std::vector<std::vector<PieceRef> > cellPieces((size_t)n*n);

    const real_t watertightScale = n * grid.h;

    for (size_t li = 0; li != loops.size(); ++li)
    {
        const Loop & loop = loops[li];
        const size_t K = loop.size();

        // Junction points: c_k(t0_k) for every curve of this loop, each
        // evaluated ONCE and shared by both curves it borders -- curve k's
        // own start, and curve k-1's own end. K == 1 (a single closed
        // curve, e.g. the hand-built circle) makes curve 0's end its own
        // start. Sharing these avoids two independent evaluations of what
        // should be the same point: a degree-1 gsBSpline edge x(t) == a is
        // not evaluated exactly at every t, so a fresh re-evaluation of a
        // junction can wiggle by an ULP and land on the wrong side of a
        // coincident grid line.
        std::vector<real_t> Jx(K), Jy(K);
        for (size_t k = 0; k != K; ++k)
        {
            real_t x,y,dx,dy;
            evalCurve1(loop[k], loop[k].t0, x, y, dx, dy);
            Jx[k] = x; Jy[k] = y;
        }

        // Tracker state, carried across every curve of this loop; must
        // return to this same value once the whole loop has been walked
        // (checked below).
        index_t col = colOf(Jx[0], X), row = rowOf(Jy[0], Y);

        for (size_t ci = 0; ci != K; ++ci)
        {
            const BoundaryCurve & curve = loop[ci];
            const size_t ciNext = (ci+1) % K;
            const std::vector<real_t> S = sampleSet(curve, nsamp);
            const index_t ns = (index_t)S.size();

            gsMatrix<real_t> tS(1,ns);
            for (index_t k = 0; k != ns; ++k) tS(0,k) = S[k];
            gsMatrix<real_t> Xs, DXs;
            curve.eval(tS, Xs, DXs);

            // Watertightness: the curve's own FRESH evaluation at t1 must
            // agree with the next curve's junction point to within a loose
            // absolute tolerance -- a real modelling bug (the loop does not
            // close), not a degeneracy, if this fires.
            const real_t dxw = Xs(0,ns-1) - Jx[ciNext], dyw = Xs(1,ns-1) - Jy[ciNext];
            GISMO_ENSURE(math::sqrt(dxw*dxw+dyw*dyw) <= 1e-10*watertightScale,
                        "loop "<<li<<", curve "<<ci<<": c(t1) = ("<<Xs(0,ns-1)<<", "<<Xs(1,ns-1)
                        <<") does not match the next curve's own start point ("<<Jx[ciNext]<<", "
                        <<Jy[ciNext]<<") -- the loop is not watertight.");

            // Every curve's own two endpoint SAMPLES are overwritten by the
            // loop's shared junction points (positions only; derivatives
            // stay the curve's own, since the boundary rule and the
            // interior-extremum search are curve-local): this is what makes
            // the crossing test at ci's own t0/t1 agree, bit-for-bit, with
            // the crossing test the neighbouring curve performs at its
            // matching endpoint.
            Xs(0,0) = Jx[ci]; Xs(1,0) = Jy[ci];
            Xs(0,ns-1) = Jx[ciNext]; Xs(1,ns-1) = Jy[ciNext];

            // Split records = this curve's own breaks and two domain
            // endpoints (delta (0,0) always) plus every crossing parameter
            // found below, each carrying the tracker move that happens
            // there.
            std::vector<SplitRecord> records;
            records.reserve(curve.breaks.size() + 2);
            records.push_back(SplitRecord{curve.t0, 0, 0});
            records.push_back(SplitRecord{curve.t1, 0, 0});
            for (real_t b : curve.breaks) records.push_back(SplitRecord{b, 0, 0});

            // Every horizontal crossing found below, keyed by its own split
            // parameter: bisectPredicate returns t* = hi, one specific side
            // of the predicate flip, so c(t*) generally overshoots the true
            // row boundary y = Y_j by a few ULP of y rather than landing on
            // it exactly. A piece ending or starting at such a t* is closed
            // against the row line it was meant to reach -- see the
            // closing-column comment below.
            std::vector<HCrossing> hcross;

            for (index_t k = 0; k+1 < ns; ++k)
            {
                const real_t ta = S[k], tb = S[k+1];
                const real_t xa = Xs(0,k), xb = Xs(0,k+1);
                const real_t ya = Xs(1,k), yb = Xs(1,k+1);
                const real_t dxa = DXs(0,k), dxb = DXs(0,k+1);
                const real_t dya = DXs(1,k), dyb = DXs(1,k+1);

                for (index_t i = 0; i <= n; ++i)
                {
                    const real_t Xi = X[i];
                    const bool ra = xa > Xi, rb = xb > Xi;
                    if (ra != rb)
                    {
                        const real_t tstar = bisectPredicate(
                            [&curve,Xi](real_t t)
                            { real_t x,y,dxv,dyv; evalCurve1(curve,t,x,y,dxv,dyv); return x > Xi; },
                            ta, tb, ra, rb);
                        real_t x,y,dxv,dyv; evalCurve1(curve, tstar, x, y, dxv, dyv);
                        // Direction = the flip direction (left-to-right =
                        // +1, right-to-left = -1), not the sign of x'(t*):
                        // the two agree away from a break or polyline
                        // vertex, and at one the flip direction is the only
                        // notion of "entering"/"leaving" that is well
                        // defined.
                        const int dir = (!ra && rb) ? +1 : -1;
                        vcross[i].push_back(VCrossing{y, dir});
                        ++stats.nVerticalCrossings;
                        // A box line (i == 0 or i == n) moves the tracker
                        // nowhere: colOf/rowOf only count INTERIOR lines,
                        // and Omega must lie strictly inside the box (file
                        // header trap) -- but the crossing still becomes a
                        // split point, and still feeds the winding walk
                        // above, in case that trap is ever violated.
                        const int dcol = (i >= 1 && i <= n-1) ? dir : 0;
                        records.push_back(SplitRecord{tstar, dcol, 0});
                        if (0.0 == dxv)
                        {
                            gsWarn << "near-tangency: curve (loop "<<li<<", curve "<<ci
                                   <<") with x = X_"<<i<<" ("<<Xi<<") near ("<<x<<", "<<y<<")\n";
                            ++stats.nNearTangency;
                        }
                    }
                }

                for (index_t j = 0; j <= n; ++j)
                {
                    const real_t Yj = Y[j];
                    const bool ra = ya > Yj, rb = yb > Yj;
                    if (ra != rb)
                    {
                        const real_t tstar = bisectPredicate(
                            [&curve,Yj](real_t t)
                            { real_t x,y,dxv,dyv; evalCurve1(curve,t,x,y,dxv,dyv); return y > Yj; },
                            ta, tb, ra, rb);
                        real_t x,y,dxv,dyv; evalCurve1(curve, tstar, x, y, dxv, dyv);
                        ++stats.nHorizontalCrossings;
                        const int dir = (!ra && rb) ? +1 : -1;
                        const int drow = (j >= 1 && j <= n-1) ? dir : 0;
                        records.push_back(SplitRecord{tstar, 0, drow});
                        hcross.push_back(HCrossing{tstar, Yj, x, y});
                        if (0.0 == dyv)
                        {
                            gsWarn << "near-tangency: curve (loop "<<li<<", curve "<<ci
                                   <<") with y = Y_"<<j<<" ("<<Yj<<") near ("<<x<<", "<<y<<")\n";
                            ++stats.nNearTangency;
                        }
                    }
                }

                if (dxa*dxb < 0.0)
                {
                    const real_t text = bisectRoot(
                        [&curve](real_t t)
                        { real_t x,y,dxv,dyv; evalCurve1(curve,t,x,y,dxv,dyv); return dxv; },
                        ta, tb);
                    real_t xext,yext,dxvE,dyvE; evalCurve1(curve, text, xext, yext, dxvE, dyvE);
                    for (index_t i = 0; i <= n; ++i)
                    {
                        const real_t Xi = X[i];
                        if (strictlyBetween(Xi,xext,xa) && strictlyBetween(Xi,xext,xb))
                            GISMO_ENSURE(false, "missed crossing pair of curve (loop "<<li<<", curve "<<ci
                                        <<") with x = X_"<<i<<" ("<<Xi<<") near ("<<xext<<", "<<yext
                                        <<"): raise --nsamp");
                        else if (math::abs(xext-Xi) < 1e-8*grid.h)
                        {
                            gsWarn << "near-tangency: curve (loop "<<li<<", curve "<<ci
                                   <<") with x = X_"<<i<<" near ("<<xext<<", "<<yext<<")\n";
                            ++stats.nNearTangency;
                        }
                    }
                }
                if (dya*dyb < 0.0)
                {
                    const real_t text = bisectRoot(
                        [&curve](real_t t)
                        { real_t x,y,dxv,dyv; evalCurve1(curve,t,x,y,dxv,dyv); return dyv; },
                        ta, tb);
                    real_t xext,yext,dxvE,dyvE; evalCurve1(curve, text, xext, yext, dxvE, dyvE);
                    for (index_t j = 0; j <= n; ++j)
                    {
                        const real_t Yj = Y[j];
                        if (strictlyBetween(Yj,yext,ya) && strictlyBetween(Yj,yext,yb))
                            GISMO_ENSURE(false, "missed crossing pair of curve (loop "<<li<<", curve "<<ci
                                        <<") with y = Y_"<<j<<" ("<<Yj<<") near ("<<xext<<", "<<yext
                                        <<"): raise --nsamp");
                        else if (math::abs(yext-Yj) < 1e-8*grid.h)
                        {
                            gsWarn << "near-tangency: curve (loop "<<li<<", curve "<<ci
                                   <<") with y = Y_"<<j<<" near ("<<xext<<", "<<yext<<")\n";
                            ++stats.nNearTangency;
                        }
                    }
                }
            }

            // Sort by t, then merge EXACTLY equal t by SUMMING dcol/drow:
            // a curve through a grid vertex crosses a vertical and a
            // horizontal line at (numerically) the same t, and std::unique
            // would silently drop one of the two simultaneous increments,
            // leaving the tracker's cell one column or row short of where
            // the crossing search itself says the piece boundary lies.
            std::sort(records.begin(), records.end(),
                     [](const SplitRecord & a, const SplitRecord & b){ return a.t < b.t; });
            std::vector<SplitRecord> merged;
            merged.reserve(records.size());
            for (const SplitRecord & rec : records)
            {
                if (!merged.empty() && merged.back().t == rec.t)
                { merged.back().dcol += rec.dcol; merged.back().drow += rec.drow; }
                else merged.push_back(rec);
            }

            // Tracker walk: apply record k's own move, THEN assign the
            // piece it opens (up to record k+1) the resulting cell -- so a
            // piece never straddles the split point that bounds it, and
            // the state right after this curve's own last record is, by
            // construction, the cell of the NEXT curve's start (checked
            // below). A segment lying exactly along X_i is not guaranteed
            // to be crossing-free: its evaluated x(t) can round an ULP to
            // either side of X_i, and wherever that happens AT a sample the
            // predicate flips and a crossing is recorded. Every recorded
            // flip enters this walk and the winding walk above alike, and
            // such flips come in pairs, so the two stay consistent. Between
            // recorded flips the tracker keeps the column it carried onto
            // the segment, which is the LEFT column where the samples give
            // x == X_i, and the segment's F.y' term supplies (or cancels)
            // that column's right-edge integral up to O(eps*|x|*h)
            // rounding. An ULP wiggle of x(t) BETWEEN samples is invisible
            // to both this walk and the winding walk above, so it can never
            // make the two disagree with each other.
            for (size_t k = 0; k != merged.size(); ++k)
            {
                col += merged[k].dcol; row += merged[k].drow;
                GISMO_ENSURE(0 <= col && col < n && 0 <= row && row < n,
                            "loop "<<li<<", curve "<<ci<<": tracker left the grid at t="
                            <<merged[k].t<<" (col="<<col<<", row="<<row<<").");
                if (k+1 != merged.size() && merged[k+1].t > merged[k].t)
                {
                    const real_t sa = merged[k].t, sb = merged[k+1].t;
                    PieceRef pc{&curve, sa, sb};
                    for (const HCrossing & hc : hcross)
                    {
                        if (hc.t == sb) pc.close.push_back(HCloseTerm{ hc.Yj, hc.x, hc.y, +1 });
                        if (hc.t == sa) pc.close.push_back(HCloseTerm{ hc.Yj, hc.x, hc.y, -1 });
                    }
                    cellPieces[col + n*row].push_back(pc);
                    ++stats.nPieces;
                }
            }
            GISMO_ENSURE(col == colOf(Jx[ciNext], X) && row == rowOf(Jy[ciNext], Y),
                        "loop "<<li<<", curve "<<ci<<": tracker ends at (col="<<col<<", row="<<row
                        <<") but the next curve's own start point is in cell (col="
                        <<colOf(Jx[ciNext],X)<<", row="<<rowOf(Jy[ciNext],Y)<<").");
        }
        GISMO_ENSURE(col == colOf(Jx[0], X) && row == rowOf(Jy[0], Y),
                    "loop "<<li<<": tracker did not return to its own start cell after "
                      "walking the whole loop (col="<<col<<", row="<<row<<", start cell col="
                    <<colOf(Jx[0],X)<<", row="<<rowOf(Jy[0],Y)<<").");
    }

    std::vector<std::vector<std::pair<real_t,real_t> > > insideIntervals(n+1);
    for (index_t i = 0; i <= n; ++i)
        insideIntervals[i] = windingInsideIntervals(vcross[i], i, X[i]);

    for (index_t j = 0; j != n; ++j)
    for (index_t i = 0; i != n; ++i)
    {
        const index_t id = i + n*j;
        const real_t xL = X[i], xR = X[i+1];
        const real_t yB = Y[j], yT = Y[j+1];

        if (!cellPieces[id].empty())
        {
            Q.status[id] = Cut; ++stats.nCut;
            std::vector<real_t> vx,vy,vw, bx,by,bw,bnx,bny;

            for (const PieceRef & pc : cellPieces[id])
            {
                gsMatrix<real_t> touter; gsVector<real_t> Wouter;
                outerRule.mapTo(pc.sa, pc.sb, touter, Wouter);
                gsMatrix<real_t> Xc, DXc;
                pc.curve->eval(touter, Xc, DXc);
                for (index_t q = 0; q != touter.cols(); ++q)
                {
                    const real_t xq = Xc(0,q), yq = Xc(1,q), dxq = DXc(0,q), dyq = DXc(1,q);
                    const real_t speed = math::sqrt(dxq*dxq + dyq*dyq);
                    bx.push_back(xq); by.push_back(yq);
                    bw.push_back(Wouter[q]*speed);
                    bnx.push_back(speed > 0 ? dyq/speed : 0.0);
                    bny.push_back(speed > 0 ? -dxq/speed : 0.0);

                    for (index_t k = 0; k <= p; ++k)
                    {
                        const real_t xNode = xL + (xq-xL)*xi[k];
                        const real_t wgt = Wouter[q]*dyq*(xq-xL)*eta[k];
                        vx.push_back(xNode); vy.push_back(yq); vw.push_back(wgt);
                        if (wgt < 0.0) ++stats.nNegativeWeights;
                    }
                }

                // Closing-column correction (volume rule only, no boundary
                // node -- the thin strip it represents is not part of the
                // actual curve): a rectangle of width (x-xL) and signed
                // height (Yj-y) at the piece endpoint that overshot its row
                // line, restoring exactness of the divergence identity there.
                for (const HCloseTerm & hc : pc.close)
                    for (index_t k = 0; k <= p; ++k)
                    {
                        const real_t xNode = xL + (hc.x-xL)*xi[k];
                        const real_t wgt = (real_t)hc.sign*(hc.Yj-hc.y)*(hc.x-xL)*eta[k];
                        vx.push_back(xNode); vy.push_back(0.5*(hc.y+hc.Yj)); vw.push_back(wgt);
                        if (wgt < 0.0) ++stats.nNegativeWeights;
                    }
            }

            const std::vector<std::pair<real_t,real_t> > edgeIv = clipIntervals(insideIntervals[i+1], yB, yT);
            for (const auto & seg : edgeIv)
            {
                gsMatrix<real_t> touter; gsVector<real_t> Wouter;
                outerRule.mapTo(seg.first, seg.second, touter, Wouter);
                for (index_t q = 0; q != touter.cols(); ++q)
                {
                    const real_t yq = touter(0,q);
                    // Width of the right-edge contribution = xR-xL, the
                    // cell's OWN box width (matching the header's C =
                    // [xL,xR] x [yB,yT]), not the nominal grid spacing h:
                    // xL, xR are independently rounded sums grid.x0+i*grid.h,
                    // and for large i or non-dyadic h their difference can be
                    // a few ULP away from h itself, which is exactly the
                    // discrepancy V3b bounds against via this same xR-xL.
                    for (index_t k = 0; k <= p; ++k)
                    {
                        const real_t xNode = xL + (xR-xL)*xi[k];
                        const real_t wgt = Wouter[q]*(xR-xL)*eta[k];
                        vx.push_back(xNode); vy.push_back(yq); vw.push_back(wgt);
                    }
                }
            }

            CellRule vr;
            vr.nodes.resize(2, (index_t)vx.size()); vr.weights.resize((index_t)vx.size());
            for (size_t k = 0; k != vx.size(); ++k)
            { vr.nodes(0,(index_t)k) = vx[k]; vr.nodes(1,(index_t)k) = vy[k]; vr.weights[(index_t)k] = vw[k]; }
            Q.vol[id] = vr;

            CellBdrRule br;
            br.nodes.resize(2, (index_t)bx.size()); br.weights.resize((index_t)bx.size());
            br.normals.resize(2, (index_t)bx.size());
            for (size_t k = 0; k != bx.size(); ++k)
            {
                br.nodes(0,(index_t)k) = bx[k]; br.nodes(1,(index_t)k) = by[k];
                br.weights[(index_t)k] = bw[k];
                br.normals(0,(index_t)k) = bnx[k]; br.normals(1,(index_t)k) = bny[k];
            }
            Q.bdr[id] = br;
        }
        else
        {
            const real_t ymid = 0.5*(yB+yT);
            if (isInside(insideIntervals[i+1], ymid))
            {
                Q.status[id] = Full; ++stats.nFull;
                gsVector<real_t> lower(2), upper(2);
                lower << xL, yB; upper << xR, yT;
                gsMatrix<real_t> nodes; gsVector<real_t> weights;
                fullRule.mapTo(lower, upper, nodes, weights);
                Q.vol[id].nodes = nodes; Q.vol[id].weights = weights;
            }
            else { Q.status[id] = Empty; ++stats.nEmpty; }
        }
    }

    return Q;
}

std::vector<std::pair<int,int> > momentSet(index_t p)
{
    std::vector<std::pair<int,int> > s;
    for (int a = 0; a <= 2*p; ++a)
        for (int b = 0; b <= 2*p-a; ++b)
            s.push_back(std::make_pair(a,b));
    return s;
}

/// Area of one cell's volume rule, Kahan-summed for consistency with the
/// global moments below (ggVolumeMoment etc.): this cell-local sum is over
/// only O(10-100) terms, so a plain Eigen::sum() is normally accurate enough
/// on its own -- the same compensation is used here regardless, so that no
/// V3b comparison can ever depend on which of the two summation orders was
/// picked.
real_t cellArea(const CellRule & r)
{
    KahanSum s;
    for (index_t k = 0; k != r.weights.size(); ++k) s.add(r.weights[k]);
    return s.value();
}

real_t ggVolumeMoment(const CutCellQuadrature & Q, int a, int b)
{
    KahanSum total;
    for (const CellRule & r : Q.vol)
        for (index_t k = 0; k != r.weights.size(); ++k)
            total.add(r.weights[k] * std::pow(r.nodes(0,k),a) * std::pow(r.nodes(1,k),b));
    return total.value();
}

real_t ggBoundaryMoment(const CutCellQuadrature & Q, int a, int b)
{
    KahanSum total;
    for (const CellBdrRule & r : Q.bdr)
        for (index_t k = 0; k != r.weights.size(); ++k)
            total.add(r.weights[k] * std::pow(r.nodes(0,k),a) * std::pow(r.nodes(1,k),b));
    return total.value();
}

/// e = |I_h - I| / max(|I|, refScale); refScale is the oracle area
/// (volume checks) or perimeter (boundary checks) -- plain relative error is
/// undefined for moments that vanish exactly (e.g. int x*y over the rotated
/// square, whose second-moment tensor is isotropic).
real_t scaledErr(real_t computed, real_t exact, real_t refScale)
{ return math::abs(computed - exact) / math::max(math::abs(exact), refScale); }

real_t shoelaceArea(const gsMatrix<real_t> & V)
{
    real_t A = 0; const index_t M = V.cols();
    for (index_t m = 0; m != M; ++m)
    {
        const index_t m1 = (m+1)%M;
        A += V(0,m)*V(1,m1) - V(0,m1)*V(1,m);
    }
    return 0.5*A;
}

void reverseVertexLoop(gsMatrix<real_t> & V)
{ V = V.rowwise().reverse().eval(); }

/// The polygon oracle's OWN orientation normalisation (shoelace area,
/// largest |A| positive, others negative), independent of
/// normalizeOrientation() above: this is a deliberately separate code path,
/// so V1 exercises the GG rule's own orientation handling against a
/// reference that could not share a bug with it.
std::vector<gsMatrix<real_t> > normalizeVertexOrientation(std::vector<gsMatrix<real_t> > V)
{
    std::vector<real_t> A(V.size());
    for (size_t i = 0; i != V.size(); ++i) A[i] = shoelaceArea(V[i]);
    size_t outer = 0;
    for (size_t i = 1; i != V.size(); ++i)
        if (math::abs(A[i]) > math::abs(A[outer])) outer = i;
    for (size_t i = 0; i != V.size(); ++i)
    {
        const bool shouldBePositive = (i == outer);
        if ((shouldBePositive && A[i] < 0) || (!shouldBePositive && A[i] > 0))
            reverseVertexLoop(V[i]);
    }
    return V;
}

/// Polygon moment via a signed fan of Duffy-collapsed triangles from the
/// origin: for edge (A,B), int_T g = int_0^1 int_0^1 g(u*A + u*v*(B-A)) * u *
/// det(A,B-A) du dv, Gauss p+2 in u and v (exact to degree 2p+3 >= a+b+1).
real_t polygonMoment(const std::vector<gsMatrix<real_t> > & loops, int a, int b, index_t p)
{
    gsGaussRule<real_t> gauss(p+2);
    const gsMatrix<real_t> & z = gauss.referenceNodes();
    const gsVector<real_t> & w = gauss.referenceWeights();
    const index_t m = p+2;
    std::vector<real_t> u(m), wu(m);
    for (index_t k = 0; k != m; ++k) { u[k] = 0.5*(1.0+z(0,k)); wu[k] = 0.5*w[k]; }

    real_t total = 0;
    for (const gsMatrix<real_t> & V : loops)
    {
        const index_t M = V.cols();
        for (index_t e = 0; e != M; ++e)
        {
            const gsVector<real_t,2> A = V.col(e);
            const gsVector<real_t,2> B = V.col((e+1)%M);
            const real_t det = A[0]*(B[1]-A[1]) - A[1]*(B[0]-A[0]);
            for (index_t iu = 0; iu != m; ++iu)
            for (index_t iv = 0; iv != m; ++iv)
            {
                const real_t uu = u[iu], vv = u[iv];
                const real_t x = uu*A[0] + uu*vv*(B[0]-A[0]);
                const real_t y = uu*A[1] + uu*vv*(B[1]-A[1]);
                const real_t gval = std::pow(x,a)*std::pow(y,b);
                total += wu[iu]*wu[iv] * gval * uu * det;
            }
        }
    }
    return total;
}

/// Boundary oracle: per edge, Gauss p+2 in the edge parameter, weight x |B-A|.
real_t polygonBoundaryMoment(const std::vector<gsMatrix<real_t> > & loops, int a, int b, index_t p)
{
    gsGaussRule<real_t> gauss(p+2);
    const gsMatrix<real_t> & z = gauss.referenceNodes();
    const gsVector<real_t> & w = gauss.referenceWeights();
    const index_t m = p+2;

    real_t total = 0;
    for (const gsMatrix<real_t> & V : loops)
    {
        const index_t M = V.cols();
        for (index_t e = 0; e != M; ++e)
        {
            const gsVector<real_t,2> A = V.col(e);
            const gsVector<real_t,2> B = V.col((e+1)%M);
            const real_t len = (B-A).norm();
            for (index_t k = 0; k != m; ++k)
            {
                const real_t s = 0.5*(1.0+z(0,k));
                const real_t x = A[0] + s*(B[0]-A[0]), y = A[1] + s*(B[1]-A[1]);
                total += 0.5*w[k] * std::pow(x,a)*std::pow(y,b) * len;
            }
        }
    }
    return total;
}

/// Disk moment by polar coordinates x = cx+rho*cos(phi), y = cy+rho*sin(phi):
/// Gauss p+2 in rho in [0,R] (integrand degree a+b+1 in rho, including the
/// rho Jacobian) times a trapezoid rule with M = 4p+4 equispaced phi (exact
/// for trigonometric degree < M).
real_t diskMoment(real_t cx, real_t cy, real_t R, int a, int b, index_t p)
{
    gsGaussRule<real_t> gauss(p+2);
    const gsMatrix<real_t> & z = gauss.referenceNodes();
    const gsVector<real_t> & w = gauss.referenceWeights();
    const index_t mr = p+2;
    const index_t M = 4*p+4;
    const real_t wphi = 2.0*(real_t)EIGEN_PI/M;

    real_t total = 0;
    for (index_t k = 0; k != mr; ++k)
    {
        const real_t rho = R*0.5*(1.0+z(0,k));
        const real_t wr  = 0.5*w[k]*R;
        real_t sumPhi = 0;
        for (index_t mIdx = 0; mIdx != M; ++mIdx)
        {
            const real_t phi = 2.0*(real_t)EIGEN_PI*mIdx/M;
            const real_t x = cx + rho*math::cos(phi), y = cy + rho*math::sin(phi);
            sumPhi += std::pow(x,a)*std::pow(y,b);
        }
        total += wr * rho * sumPhi * wphi;
    }
    return total;
}

/// Circle boundary moment, same trapezoid rule, ds = R dphi.
real_t circleBoundaryMoment(real_t cx, real_t cy, real_t R, int a, int b, index_t p)
{
    const index_t M = 4*p+4;
    real_t total = 0;
    for (index_t m = 0; m != M; ++m)
    {
        const real_t phi = 2.0*(real_t)EIGEN_PI*m/M;
        const real_t x = cx + R*math::cos(phi), y = cy + R*math::sin(phi);
        total += std::pow(x,a)*std::pow(y,b);
    }
    return total * R * (2.0*(real_t)EIGEN_PI/M);
}

real_t exactMoment(const DefaultGeometry & geo, int a, int b, index_t p)
{
    std::vector<gsMatrix<real_t> > sq(1, geo.corners.transpose());
    return polygonMoment(sq, a, b, p) - diskMoment(geo.cx, geo.cy, geo.R, a, b, p);
}

real_t exactBoundaryMoment(const DefaultGeometry & geo, int a, int b, index_t p)
{
    std::vector<gsMatrix<real_t> > sq(1, geo.corners.transpose());
    return polygonBoundaryMoment(sq, a, b, p) + circleBoundaryMoment(geo.cx, geo.cy, geo.R, a, b, p);
}

/// V1 (polyline mode only): GG rule at --nq vs the polygon oracle.
bool checkV1(const CutCellQuadrature & Q, const std::vector<gsMatrix<real_t> > & rawVertices, index_t p)
{
    const std::vector<gsMatrix<real_t> > V = normalizeVertexOrientation(rawVertices);
    const real_t polyArea  = polygonMoment(V, 0, 0, p);
    const real_t polyPerim = polygonBoundaryMoment(V, 0, 0, p);
    const std::vector<std::pair<int,int> > moments = momentSet(p);

    real_t volErr = 0;
    for (const auto & ab : moments)
        volErr = math::max(volErr, scaledErr(ggVolumeMoment(Q,ab.first,ab.second),
                                             polygonMoment(V,ab.first,ab.second,p), polyArea));

    const real_t perimErr = scaledErr(ggBoundaryMoment(Q,0,0), polyPerim, polyPerim);

    real_t bdrErr = 0;
    for (const auto & ab : moments)
        bdrErr = math::max(bdrErr, scaledErr(ggBoundaryMoment(Q,ab.first,ab.second),
                                             polygonBoundaryMoment(V,ab.first,ab.second,p), polyPerim));

    const bool pass = (volErr <= 1e-13 && perimErr <= 1e-13 && bdrErr <= 1e-13);
    gsInfo << "V1 polyline: vol " << fmtSci(volErr) << " perim " << fmtSci(perimErr)
           << " bdr " << fmtSci(bdrErr) << "  " << (pass ? "PASS" : "FAIL") << "\n";
    return pass;
}

/// V2 (spline mode only): nq ladder 2..12 vs the exact oracle, using \a
/// producer to rebuild the rule at each nq (never touches the loops
/// directly: any rule producer with this signature can be substituted here
/// unchanged).
bool checkV2(const std::function<CutCellQuadrature(index_t)> & producer,
            const DefaultGeometry & geo, index_t p)
{
    const real_t exactArea = geo.exactArea(), exactPerim = geo.exactPerimeter();
    const std::vector<std::pair<int,int> > moments = momentSet(p);

    gsInfo << std::right << std::setw(4) << "nq" << std::setw(18) << "vol scaled err" << std::setw(18) << "bdr scaled err\n";
    real_t lastVol = 0, lastBdr = 0;
    for (index_t nqTest = 2; nqTest <= 12; ++nqTest)
    {
        const CutCellQuadrature Qt = producer(nqTest);
        real_t volErr = 0, bdrErr = 0;
        for (const auto & ab : moments)
        {
            volErr = math::max(volErr, scaledErr(ggVolumeMoment(Qt,ab.first,ab.second),
                                                 exactMoment(geo,ab.first,ab.second,p), exactArea));
            bdrErr = math::max(bdrErr, scaledErr(ggBoundaryMoment(Qt,ab.first,ab.second),
                                                 exactBoundaryMoment(geo,ab.first,ab.second,p), exactPerim));
        }
        gsInfo << std::setw(4) << nqTest << std::setw(18) << fmtSci(volErr) << std::setw(18) << fmtSci(bdrErr) << "\n";
        if (12 == nqTest) { lastVol = volErr; lastBdr = bdrErr; }
    }
    const bool pass = (lastVol <= 1e-13 && lastBdr <= 1e-13);
    gsInfo << "V2 spline: nq=12 vol " << fmtSci(lastVol) << " bdr " << fmtSci(lastBdr)
           << "  " << (pass ? "PASS" : "FAIL") << "\n";
    return pass;
}

/// V3a-d: self-consistency checks of the GG rule alone (no oracle other
/// than the identities themselves). Every argument here is either \a Q or a
/// precomputed reference -- never the loops -- so any rule producer's
/// output can be checked by the same function unchanged. \a areaScale is
/// the V3a scale: geo.exactArea() for the default geometry, or the global
/// GG area |refMoments[0]| (momentSet's first entry is (0,0)) when the
/// geometry is not known to be DefaultGeometry.
bool checkV3(const CutCellQuadrature & Q, const std::vector<real_t> & refMoments,
            const std::vector<std::pair<int,int> > & moments, real_t areaScale)
{
    real_t v3a = 0;
    for (size_t k = 0; k != moments.size(); ++k)
        v3a = math::max(v3a, scaledErr(ggVolumeMoment(Q,moments[k].first,moments[k].second),
                                       refMoments[k], areaScale));
    const bool passA = v3a <= 1e-13;
    gsInfo << "V3a: " << fmtSci(v3a) << "  " << (passA ? "PASS" : "FAIL") << "\n";

    const real_t hh = Q.grid.h*Q.grid.h;
    const index_t n = Q.grid.n;
    // A boundary piece lying exactly on, or within an ULP of, a knot line
    // carries an evaluation rounding of x of order ulp(L), L = the largest
    // background-box coordinate magnitude; the volume rule multiplies that
    // rounding by the cell height through the (x-xL) factor of the inner
    // Gauss node/weight (see gaussGreenQuadrature), so the irreducible
    // per-cell area error at such a cell is O(eps*L*h), not O(eps*h^2). A
    // purely h^2-relative bound (the box*(1+-1e-14) form below, box ~ h^2) therefore
    // rejects a correct rule once h drops below about L*eps/1e-14 ~ 0.02 for
    // L ~ 1. The symmetric absolute slack tolAbs adds an O(eps*L*h) margin
    // (measured maximum 1.004*eps*L*h on an exactly grid-aligned edge; the
    // factor 16 below is a 16x margin over that measurement) to both sides
    // of the bound, so a cell whose area is wrong by more than about
    // 1e-12*h^2 at h = 0.005 -- an order of magnitude above the measured
    // rounding -- still fails.
    const real_t L = math::max(math::max(math::abs(Q.grid.x0), math::abs(Q.grid.x0 + n*Q.grid.h)),
                               math::max(math::abs(Q.grid.y0), math::abs(Q.grid.y0 + n*Q.grid.h)));
    const real_t tolAbs = 16*std::numeric_limits<real_t>::epsilon()*L*Q.grid.h;
    real_t cutMin = std::numeric_limits<real_t>::infinity();
    real_t cutMax = -std::numeric_limits<real_t>::infinity();
    bool boundsOk = true;
    for (size_t id = 0; id != Q.vol.size(); ++id)
    {
        const index_t i = (index_t)id % n, j = (index_t)id / n;
        const real_t xL = Q.grid.x0 + i*Q.grid.h, xR = Q.grid.x0 + (i+1)*Q.grid.h;
        const real_t yB = Q.grid.y0 + j*Q.grid.h, yT = Q.grid.y0 + (j+1)*Q.grid.h;
        // The true bound of a cell's area is its own box (xR-xL)*(yT-yB), not
        // the nominal h^2: for non-dyadic h these differ by a few ULP
        // (grid.x0 + i*grid.h is rounded independently per i), and the
        // full-cell tensor Gauss rule integrates the box exactly, not h^2.
        const real_t box = (xR-xL)*(yT-yB);
        const real_t area = cellArea(Q.vol[id]);
        if (area < -(1e-14*box + tolAbs) || area > box*(1.0+1e-14) + tolAbs) boundsOk = false;
        if (Cut == Q.status[id])
        {
            cutMin = math::min(cutMin, area/hh);
            cutMax = math::max(cutMax, area/hh);
        }
    }
    const bool passB = boundsOk;
    gsInfo << "V3b: cut-cell min " << fmtSci(cutMin) << " max " << fmtSci(cutMax) << "  " << (passB ? "PASS" : "FAIL") << "\n";

    KahanSum sx, sy;
    for (const CellBdrRule & r : Q.bdr)
        for (index_t k = 0; k != r.weights.size(); ++k)
        { sx.add(r.weights[k]*r.normals(0,k)); sy.add(r.weights[k]*r.normals(1,k)); }
    const real_t v3c = math::max(math::abs(sx.value()), math::abs(sy.value()));
    const bool passC = v3c <= 1e-13;
    gsInfo << "V3c: " << fmtSci(v3c) << "  " << (passC ? "PASS" : "FAIL") << "\n";

    KahanSum xdotn;
    for (const CellBdrRule & r : Q.bdr)
        for (index_t k = 0; k != r.weights.size(); ++k)
            xdotn.add(r.weights[k]*(r.nodes(0,k)*r.normals(0,k) + r.nodes(1,k)*r.normals(1,k)));
    const real_t areaGG = ggVolumeMoment(Q,0,0);
    const real_t v3d = math::abs(xdotn.value() - 2.0*areaGG) / math::max((real_t)1.0, areaGG);
    const bool passD = v3d <= 1e-13;
    gsInfo << "V3d: " << fmtSci(v3d) << "  " << (passD ? "PASS" : "FAIL") << "\n";

    return passA && passB && passC && passD;
}

void printHeader(index_t p, index_t r, const Grid & grid, index_t nq, index_t nsamp, index_t polyN,
                 const std::string & source, const CutCellQuadrature & Q, const GGStats & stats, real_t ggTime)
{
    index_t nVolNodes = 0, nBdrNodes = 0;
    for (const CellRule & rl : Q.vol) nVolNodes += rl.weights.size();
    for (const CellBdrRule & rl : Q.bdr) nBdrNodes += rl.weights.size();

    gsInfo << "p=" << p << " r=" << r << " n=" << grid.n << " h=" << fmtSci(grid.h)
           << " nq=" << nq << " nsamp=" << nsamp << " polyline=" << polyN << " source=" << source << "\n";
    gsInfo << "cells: cut=" << stats.nCut << " full=" << stats.nFull << " empty=" << stats.nEmpty
           << "  crossings: vert=" << stats.nVerticalCrossings << " horiz=" << stats.nHorizontalCrossings
           << "  pieces=" << stats.nPieces << "\n";
    gsInfo << "GG volume nodes=" << nVolNodes << " (negative weights=" << stats.nNegativeWeights << ")"
           << "  GG boundary nodes=" << nBdrNodes << "  near-tangency warnings=" << stats.nNearTangency << "\n";
    gsInfo << "GG build time = " << fmtSci(ggTime,2) << " s\n";
}

/// R-conjunction a^b = a + b + sqrt(a^2+b^2) and its chain rule through two
/// coordinates (da1 = da/dx1 etc.). At s == 0 exactly both weights a/s, b/s
/// are set to 0: phi is not differentiable there (the square corners, where
/// both conjoined factors vanish simultaneously), and any bounded value is a
/// legitimate weak sub-gradient.
void conjunction(real_t a, real_t b, real_t da1, real_t da2, real_t db1, real_t db2,
                 real_t & val, real_t & dval1, real_t & dval2)
{
    const real_t s  = math::sqrt(a*a + b*b);
    const real_t wa = (0.0 == s) ? 0.0 : a/s;
    const real_t wb = (0.0 == s) ? 0.0 : b/s;
    val   = a + b + s;
    dval1 = (1.0+wa)*da1 + (1.0+wb)*db1;
    dval2 = (1.0+wa)*da2 + (1.0+wb)*db2;
}

/// R-function conjunction level set of the same square-minus-disk geometry
/// as DefaultGeometry, negative inside: xi = R(-theta) x (centre 0),
/// f = (xi1^2-a^2)/a^2, g = (xi2^2-a^2)/a^2 (negative inside the square),
/// h = (R^2-|x-c|^2)/R^2 (positive inside the disk), phi = (f^g)^h. Since
/// sign(a^b) = sign(max(a,b)), the zero set of phi is exactly the boundary
/// of the square-minus-disk domain -- unlike max(phi_square,-phi_disk),
/// whose gradient also jumps on the interior medial curve
/// {phi_square = -phi_disk}, which runs through cut cells near the disk.
/// deriv_into is load-bearing: the Algoim wrapper's interval (Taylor) bounds
/// and surface normal both come from it.
class AlgoimLevelSet : public gsFunction<real_t>
{
public:
    GISMO_CLONE_FUNCTION(AlgoimLevelSet)

    AlgoimLevelSet(const DefaultGeometry & geo, const gsMatrix<real_t> & box)
    : m_geo(geo), m_box(box) { }

    short_t domainDim() const override { return 2; }
    short_t targetDim() const override { return 1; }
    gsMatrix<real_t> support() const override { return m_box; }

    void eval_into(const gsMatrix<real_t> & u, gsMatrix<real_t> & result) const override
    {
        result.resize(1, u.cols());
        for (index_t k = 0; k != u.cols(); ++k)
            evalPhi(u(0,k), u(1,k), result(0,k));
    }

    /// Row 0 = d(phi)/dx, row 1 = d(phi)/dy.
    void deriv_into(const gsMatrix<real_t> & u, gsMatrix<real_t> & result) const override
    {
        result.resize(2, u.cols());
        real_t phi;
        for (index_t k = 0; k != u.cols(); ++k)
            evalPhiGrad(u(0,k), u(1,k), phi, result(0,k), result(1,k));
    }

private:
    void toXi(real_t x, real_t y, real_t & xi1, real_t & xi2) const
    { xi1 = m_geo.ct*x + m_geo.st*y; xi2 = -m_geo.st*x + m_geo.ct*y; }

    void rotateToPhysical(real_t g1, real_t g2, real_t & gx, real_t & gy) const
    { gx = m_geo.ct*g1 - m_geo.st*g2; gy = m_geo.st*g1 + m_geo.ct*g2; }

    void evalPhi(real_t x, real_t y, real_t & phi) const
    {
        real_t xi1, xi2; toXi(x, y, xi1, xi2);
        const real_t aa = m_geo.a*m_geo.a;
        const real_t f = (xi1*xi1 - aa)/aa;
        const real_t g = (xi2*xi2 - aa)/aa;
        const real_t phi1 = f + g + math::sqrt(f*f+g*g);
        const real_t dx = x-m_geo.cx, dy = y-m_geo.cy;
        const real_t RR = m_geo.R*m_geo.R;
        const real_t h = (RR - (dx*dx+dy*dy))/RR;
        phi = phi1 + h + math::sqrt(phi1*phi1+h*h);
    }

    void evalPhiGrad(real_t x, real_t y, real_t & phi, real_t & gx, real_t & gy) const
    {
        real_t xi1, xi2; toXi(x, y, xi1, xi2);
        const real_t aa = m_geo.a*m_geo.a;
        const real_t f   = (xi1*xi1 - aa)/aa;
        const real_t df1 = 2.0*xi1/aa;
        const real_t g   = (xi2*xi2 - aa)/aa;
        const real_t dg2 = 2.0*xi2/aa;

        real_t phi1, dphi1_1, dphi1_2;
        conjunction(f, g, df1, 0.0, 0.0, dg2, phi1, dphi1_1, dphi1_2);

        real_t g1x, g1y;
        rotateToPhysical(dphi1_1, dphi1_2, g1x, g1y);

        const real_t dx = x-m_geo.cx, dy = y-m_geo.cy;
        const real_t RR = m_geo.R*m_geo.R;
        const real_t h   = (RR - (dx*dx+dy*dy))/RR;
        const real_t dh1 = -2.0*dx/RR;
        const real_t dh2 = -2.0*dy/RR;

        conjunction(phi1, h, g1x, g1y, dh1, dh2, phi, gx, gy);
    }

    DefaultGeometry  m_geo;
    gsMatrix<real_t> m_box;
};

/// V4 baseline (INFO, no threshold): direct construction of
/// gsAlgoimAdaptiveRule (rule 15), mirroring the gsQuadrature::getPtr recipe
/// (src/gsAssembler/gsQuadrature.h, makeAlgoimAdaptivePtr) line by line, but
/// built here directly because the per-cell comparison needs the rule on
/// our own cell indexing and an independent classification. Returns the
/// concatenated Algoim nodes/weights over all cells (for --plot).
void runAlgoimBaseline(const DefaultGeometry & geo, const Grid & grid, index_t p, index_t maxDepth,
                       real_t indicatorTol, real_t lipschitz, const CutCellQuadrature & Q,
                       bool polylineMode, gsMatrix<real_t> & algoimNodes, gsVector<real_t> & algoimWeights)
{
    gsMatrix<real_t> box(2,2);
    box << grid.x0, grid.x0 + grid.n*grid.h, grid.y0, grid.y0 + grid.n*grid.h;
    AlgoimLevelSet phi(geo, box);   // named local: gsAlgoimAdaptiveRule stores phi non-owned

    gsVector<index_t> nn(2); nn.setConstant(5);
    gsLobattoRule<real_t> QR(nn);
    real_t maxGrad = 0;
    for (index_t j = 0; j != grid.n; ++j)
    for (index_t i = 0; i != grid.n; ++i)
    {
        gsVector<real_t> lower(2), upper(2);
        lower << grid.x0+i*grid.h,     grid.y0+j*grid.h;
        upper << grid.x0+(i+1)*grid.h, grid.y0+(j+1)*grid.h;
        gsMatrix<real_t> pts, grad; gsVector<real_t> wts;
        QR.mapTo(lower, upper, pts, wts);
        phi.deriv_into(pts, grad);
        for (index_t c = 0; c != grad.cols(); ++c)
            maxGrad = math::max(maxGrad, math::sqrt(grad(0,c)*grad(0,c) + grad(1,c)*grad(1,c)));
    }
    GISMO_ENSURE(lipschitz >= 1.5*maxGrad,
                "--lipschitz "<<lipschitz<<" is below 1.5*maxGrad = "<<1.5*maxGrad
                <<": raise --lipschitz, or gsAlgoimAdaptiveRule::classify may misclassify a "
                  "genuinely cut sub-box as uncut and silently drop area.");
    gsInfo << "Algoim: maxGrad = " << fmtSci(maxGrad) << " (Lipschitz = " << lipschitz << ")\n";

    gsOptionList o = gsAlgoimAdaptiveRule<real_t>::defaultOptions();
    o.setInt   ("dim", -1);
    o.setReal  ("quA", 1.0);
    o.setInt   ("quB", 1);
    o.setInt   ("maxDepth", maxDepth);
    o.setString("indicator", "integralChange");
    o.setReal  ("indicatorTol", indicatorTol);
    o.setReal  ("LipschitzConstant", lipschitz);
    gsAlgoimAdaptiveRule<real_t> rule(phi, (short_t)p, o);

    gsStopwatch clk;
    std::vector<real_t> areaAlgoim((size_t)grid.n*grid.n, 0.0);
    std::vector<gsMatrix<real_t> > cellNodes((size_t)grid.n*grid.n);
    std::vector<gsVector<real_t> > cellWeights((size_t)grid.n*grid.n);
    index_t totalAlgoimNodes = 0;
    for (index_t j = 0; j != grid.n; ++j)
    for (index_t i = 0; i != grid.n; ++i)
    {
        const index_t id = i + grid.n*j;
        gsVector<real_t> lower(2), upper(2);
        lower << grid.x0+i*grid.h,     grid.y0+j*grid.h;
        upper << grid.x0+(i+1)*grid.h, grid.y0+(j+1)*grid.h;
        gsMatrix<real_t> nodes; gsVector<real_t> weights;
        rule.mapTo(lower, upper, nodes, weights);
        areaAlgoim[id] = weights.size() ? weights.sum() : 0.0;
        cellNodes[id] = nodes; cellWeights[id] = weights;
        totalAlgoimNodes += nodes.cols();
    }
    const real_t algoimTime = clk.stop();

    real_t areaGGtotal = 0, areaAlgoimTotal = 0;
    for (index_t id = 0; id != grid.n*grid.n; ++id)
    {
        areaGGtotal     += cellArea(Q.vol[id]);
        areaAlgoimTotal += areaAlgoim[id];
    }
    const real_t exactArea = geo.exactArea();
    const real_t ggAreaErr     = scaledErr(areaGGtotal, exactArea, exactArea);
    const real_t algoimAreaErr = scaledErr(areaAlgoimTotal, exactArea, exactArea);

    // Max scaled err over M_vol = {x^a y^b : a+b <= 2p} for both methods
    // against the exact oracle (the GG rule's own totals above only cover
    // a = b = 0). The Algoim moment is summed with the same compensated
    // accumulator as ggVolumeMoment, over the concatenated per-cell nodes.
    const std::vector<std::pair<int,int> > mvol = momentSet(p);
    real_t ggVolErr = 0, algoimVolErr = 0;
    for (const auto & ab : mvol)
    {
        const real_t exact = exactMoment(geo, ab.first, ab.second, p);
        ggVolErr = math::max(ggVolErr, scaledErr(ggVolumeMoment(Q,ab.first,ab.second), exact, exactArea));

        KahanSum algMoment;
        for (size_t id = 0; id != cellNodes.size(); ++id)
        {
            const gsMatrix<real_t> & nd = cellNodes[id];
            const gsVector<real_t> & wt = cellWeights[id];
            for (index_t k = 0; k != wt.size(); ++k)
                algMoment.add(wt[k] * std::pow(nd(0,k),ab.first) * std::pow(nd(1,k),ab.second));
        }
        algoimVolErr = math::max(algoimVolErr, scaledErr(algMoment.value(), exact, exactArea));
    }

    if (polylineMode)
        gsInfo << "V4 note: --polyline is active, so the GG rule describes the polygon while "
                  "Algoim always describes the exact geometry; the per-cell diff below mixes "
                  "geometry and quadrature error.\n";

    gsInfo << "V4 totals: exact area = " << fmtSci(exactArea)
           << "  GG area = " << fmtSci(areaGGtotal) << " (scaled err " << fmtSci(ggAreaErr)
           << ", max M_vol scaled err " << fmtSci(ggVolErr) << ")"
           << "  Algoim area = " << fmtSci(areaAlgoimTotal) << " (scaled err " << fmtSci(algoimAreaErr)
           << ", max M_vol scaled err " << fmtSci(algoimVolErr) << ")  INFO\n";

    struct Worst { index_t i,j; real_t diff; };
    std::vector<Worst> worst;
    real_t maxDiff = 0;
    const real_t hh = grid.h*grid.h;
    for (index_t j = 0; j != grid.n; ++j)
    for (index_t i = 0; i != grid.n; ++i)
    {
        const index_t id = i + grid.n*j;
        const real_t aGG = cellArea(Q.vol[id]);
        const real_t diff = math::abs(aGG - areaAlgoim[id]);
        maxDiff = math::max(maxDiff, diff/hh);
        worst.push_back(Worst{i,j,diff});
    }
    std::sort(worst.begin(), worst.end(), [](const Worst & lhs, const Worst & rhs){ return lhs.diff > rhs.diff; });

    gsInfo << "V4 max |area_GG - area_Algoim| / h^2 = " << fmtSci(maxDiff) << "  INFO\n";
    gsInfo << "V4 worst cells:\n";
    for (index_t k = 0; k != 5 && k != (index_t)worst.size(); ++k)
    {
        const index_t i = worst[k].i, j = worst[k].j, id = i + grid.n*j;
        const real_t cx = grid.x0+(i+0.5)*grid.h, cy = grid.y0+(j+0.5)*grid.h;
        const real_t aGG = cellArea(Q.vol[id]);
        const real_t xL = grid.x0+i*grid.h, xR = grid.x0+(i+1)*grid.h;
        const real_t yB = grid.y0+j*grid.h, yT = grid.y0+(j+1)*grid.h;
        bool corner = false;
        for (index_t c = 0; c != 4; ++c)
        {
            const real_t px = geo.corners(c,0), py = geo.corners(c,1);
            if (px >= xL && px <= xR && py >= yB && py <= yT) { corner = true; break; }
        }
        gsInfo << "  (" << i << "," << j << ") centre=(" << cx << "," << cy << ")"
               << " GG/h^2=" << fmtSci(aGG/hh) << " Algoim/h^2=" << fmtSci(areaAlgoim[id]/hh)
               << " diff/h^2=" << fmtSci(worst[k].diff/hh) << " corner=" << (corner ? "yes" : "no") << "\n";
    }

    gsInfo << "V4 Algoim: nodes=" << totalAlgoimNodes
           << " nFallbackLeaves=" << rule.stats().nFallbackLeaves
           << " nSubBoxes=" << rule.stats().nSubBoxes
           << " build time = " << fmtSci(algoimTime,2) << " s  INFO\n";

    index_t total = 0;
    for (const auto & w : cellWeights) total += w.size();
    algoimNodes.resize(2, total);
    algoimWeights.resize(total);
    index_t col = 0;
    for (size_t id = 0; id != cellNodes.size(); ++id)
    {
        const index_t m = cellWeights[id].size();
        if (m > 0)
        {
            algoimNodes.block(0,col,2,m) = cellNodes[id];
            algoimWeights.segment(col,m) = cellWeights[id];
            col += m;
        }
    }
}

/// False iff a square edge lies exactly on a background grid line: theta = 0
/// or 180 degrees makes the square axis-aligned (sin(theta)*cos(theta) ~ 0,
/// tested to 1e-14 since the product itself, not theta, is what determines
/// whether the edges are axis-aligned), and \a half is then compared against
/// every grid line to 1e-12 (the square can be axis-aligned without any
/// coincident line). In that configuration the level set vanishes
/// identically on a whole cell face, and gsAlgoimAdaptiveRule::mapTo
/// enumerates every sub-box down to maxDepth without ever classifying one as
/// uncut, exhausting memory (a gsAlgoim adapter limitation -- see the file
/// header trap). \a geo.a is the SAME half-side used to check +-half, since
/// half == geo.a here.
bool algoimBaselineUsable(const DefaultGeometry & geo, const Grid & grid)
{
    if (math::abs(geo.st*geo.ct) >= 1e-14) return true;
    for (index_t i = 0; i <= grid.n; ++i)
    {
        const real_t Xi = grid.x0 + i*grid.h;
        if (math::abs(geo.a-Xi) <= 1e-12 || math::abs(-geo.a-Xi) <= 1e-12) return false;
    }
    for (index_t j = 0; j <= grid.n; ++j)
    {
        const real_t Yj = grid.y0 + j*grid.h;
        if (math::abs(geo.a-Yj) <= 1e-12 || math::abs(-geo.a-Yj) <= 1e-12) return false;
    }
    return true;
}

/// Writes the GG volume/boundary nodes and (unless \a algoimNodes is null)
/// the Algoim nodes as ParaView point sets, weight carried as the point
/// value V.
void writePlot(const CutCellQuadrature & Q, const gsMatrix<real_t> * algoimNodes,
              const gsVector<real_t> * algoimWeights, const std::string & out)
{
    index_t totalVol = 0, totalBdr = 0;
    for (const CellRule & r : Q.vol) totalVol += r.weights.size();
    for (const CellBdrRule & r : Q.bdr) totalBdr += r.weights.size();

    gsMatrix<real_t> Xv(1,totalVol), Yv(1,totalVol), Zv(1,totalVol), Vv(1,totalVol);
    index_t col = 0;
    for (const CellRule & r : Q.vol)
        for (index_t k = 0; k != r.weights.size(); ++k)
        { Xv(0,col)=r.nodes(0,k); Yv(0,col)=r.nodes(1,k); Zv(0,col)=0.0; Vv(0,col)=r.weights[k]; ++col; }
    gsWriteParaviewPoints(Xv, Yv, Zv, Vv, out + "/gg_volume_nodes", 12);

    gsMatrix<real_t> Xb(1,totalBdr), Yb(1,totalBdr), Zb(1,totalBdr), Vb(1,totalBdr);
    col = 0;
    for (const CellBdrRule & r : Q.bdr)
        for (index_t k = 0; k != r.weights.size(); ++k)
        { Xb(0,col)=r.nodes(0,k); Yb(0,col)=r.nodes(1,k); Zb(0,col)=0.0; Vb(0,col)=r.weights[k]; ++col; }
    gsWriteParaviewPoints(Xb, Yb, Zb, Vb, out + "/gg_boundary_nodes", 12);

    if (algoimNodes)
    {
        const index_t m = algoimWeights->size();
        gsMatrix<real_t> Xa(1,m), Ya(1,m), Za(1,m), Va(1,m);
        for (index_t k = 0; k != m; ++k)
        { Xa(0,k)=(*algoimNodes)(0,k); Ya(0,k)=(*algoimNodes)(1,k); Za(0,k)=0.0; Va(0,k)=(*algoimWeights)[k]; }
        gsWriteParaviewPoints(Xa, Ya, Za, Va, out + "/algoim_nodes", 12);
    }

    gsInfo << "ParaView output written to " << out << "\n";
}

/// V5/V6: compares a non-"hand" source's GG rule (\a Qsrc, built at the
/// fixed nq = 12 so the check never false-alarms at a low --nq, where two
/// legitimately different parametrisations of the same circle can still
/// disagree) against the hand-built rule at the same nq (\a Qhand) -- rules
/// only, never loops, so this is the same producer-only contract as V2.
/// \a label is "V5 occ vs hand" or "V6 xml vs hand"; \a tol is 1e-12 resp.
/// 1e-14.
bool checkVsHand(const CutCellQuadrature & Qsrc, const CutCellQuadrature & Qhand,
                 const DefaultGeometry & geo, index_t p, const std::string & label, real_t tol)
{
    const std::vector<std::pair<int,int> > moments = momentSet(p);
    real_t volErr = 0, bdrErr = 0;
    for (const auto & ab : moments)
    {
        volErr = math::max(volErr, scaledErr(ggVolumeMoment(Qsrc,ab.first,ab.second),
                                             ggVolumeMoment(Qhand,ab.first,ab.second), geo.exactArea()));
        bdrErr = math::max(bdrErr, scaledErr(ggBoundaryMoment(Qsrc,ab.first,ab.second),
                                             ggBoundaryMoment(Qhand,ab.first,ab.second), geo.exactPerimeter()));
    }
    const bool pass = (volErr <= tol && bdrErr <= tol);
    gsInfo << label << " (nq=12): vol " << fmtSci(volErr) << " bdr " << fmtSci(bdrErr)
           << "  " << (pass ? "PASS" : "FAIL") << "\n";
    return pass;
}

/// Aborts (GISMO_ENSURE) unless every sampled point of every curve of every
/// loop lies STRICTLY inside the background box: the file header's "Omega
/// strictly inside the box" trap, turned into an early, specific abort
/// instead of a confusing failure deep inside the tracker. The main
/// offender is a foreign CAD file: OCCT reads STEP in millimetres by
/// default, so a model authored in metres lands two-three orders of
/// magnitude outside the box built from --n0/-r.
void ensureInsideBox(const std::vector<Loop> & loops, const Grid & grid, index_t nsamp)
{
    const real_t xlo = grid.x0, xhi = grid.x0 + grid.n*grid.h;
    const real_t ylo = grid.y0, yhi = grid.y0 + grid.n*grid.h;
    for (const Loop & loop : loops)
        for (const BoundaryCurve & c : loop)
        {
            const std::vector<real_t> S = sampleSet(c, nsamp);
            gsMatrix<real_t> tS(1,(index_t)S.size());
            for (size_t k = 0; k != S.size(); ++k) tS(0,(index_t)k) = S[k];
            gsMatrix<real_t> Xs, DXs;
            c.eval(tS, Xs, DXs);
            for (index_t k = 0; k != Xs.cols(); ++k)
            {
                const real_t x = Xs(0,k), y = Xs(1,k);
                GISMO_ENSURE(x > xlo && x < xhi && y > ylo && y < yhi,
                            "point (" << x << ", " << y << ") lies outside the background box "
                              "[" << xlo << ", " << xhi << "] x [" << ylo << ", " << yhi << "]");
            }
        }
}

/// A planar triangle mesh, flattened out of gsSurfMesh into plain arrays so
/// that meshClipQuadrature() and the oracles below never touch gsSurfMesh
/// itself. \a nV/\a nE/\a nF are the reader's own counts, used only by V7's
/// Euler-characteristic line; a face silently dropped by
/// gsSurfMesh::add_face() (a complex edge/vertex) shows up there and
/// nowhere else, since the reader's own return value stays true.
struct PlanarMesh
{
    std::vector<gsMatrix<real_t> > tri;       // 2 x 3 per triangle
    std::vector<gsMatrix<real_t> > bdrEdge;   // 2 x 3 per boundary edge: A (from), B (to), C (third vertex of the adjacent face)
    index_t nV = 0, nE = 0, nF = 0;
    index_t nDegenerate = 0;
    real_t xmin = 0, xmax = 0, ymin = 0, ymax = 0;
};

/// Reads \a file via gsReadSurfMesh (STL/OFF/OBJ, resolved through the
/// search paths) and flattens it to a PlanarMesh. Every vertex must be
/// planar to 1e-12 times the larger xy-bbox extent; z is then dropped.
/// Degenerate triangles (zero 2D edge-vector determinant) are skipped and
/// counted, never pushed into \a tri. Boundary edges are found by walking
/// every boundary halfedge once (a boundary halfedge has no face; its
/// OPPOSITE halfedge does, and that face is the triangle the edge borders),
/// which is why the boundary rule needs no orientation convention of its
/// own -- the outward normal is derived per edge from that triangle's third
/// vertex C in meshClipQuadrature().
PlanarMesh readPlanarMesh(const std::string & file)
{
    gsSurfMesh<real_t> mesh;
    GISMO_ENSURE(gsReadSurfMesh(file, mesh), "readPlanarMesh: cannot read mesh '" << file << "'");

    PlanarMesh M;
    M.nV = (index_t)mesh.n_vertices();
    M.nE = (index_t)mesh.n_edges();
    M.nF = (index_t)mesh.n_faces();

    real_t xmin =  std::numeric_limits<real_t>::infinity(), xmax = -xmin;
    real_t ymin =  std::numeric_limits<real_t>::infinity(), ymax = -ymin;
    for (gsSurfMesh<real_t>::Vertex v : mesh.vertices())
    {
        const gsSurfMesh<real_t>::Point & pt = mesh.position(v);
        xmin = math::min(xmin, pt[0]); xmax = math::max(xmax, pt[0]);
        ymin = math::min(ymin, pt[1]); ymax = math::max(ymax, pt[1]);
    }
    M.xmin = xmin; M.xmax = xmax; M.ymin = ymin; M.ymax = ymax;

    const real_t scale = math::max(xmax-xmin, ymax-ymin);
    real_t z0 = 0; bool haveZ0 = false;
    for (gsSurfMesh<real_t>::Vertex v : mesh.vertices())
    {
        const gsSurfMesh<real_t>::Point & pt = mesh.position(v);
        if (!haveZ0) { z0 = pt[2]; haveZ0 = true; }
        GISMO_ENSURE(math::abs(pt[2]-z0) <= 1e-12*scale,
                    "readPlanarMesh: '" << file << "' is not planar (z=" << pt[2]
                      << ", z0=" << z0 << ")");
    }

    for (gsSurfMesh<real_t>::Face f : mesh.faces())
    {
        gsMatrix<real_t> T(2,3);
        index_t col = 0;
        for (gsSurfMesh<real_t>::Vertex v : mesh.vertices(f))
        {
            GISMO_ENSURE(col < 3, "readPlanarMesh: '" << file << "' has a non-triangular face");
            const gsSurfMesh<real_t>::Point & pt = mesh.position(v);
            T(0,col) = pt[0]; T(1,col) = pt[1]; ++col;
        }
        const real_t det = (T(0,1)-T(0,0))*(T(1,2)-T(1,0)) - (T(1,1)-T(1,0))*(T(0,2)-T(0,0));
        if (0.0 == det) { ++M.nDegenerate; continue; }
        M.tri.push_back(T);
    }
    if (M.nDegenerate)
        gsInfo << "readPlanarMesh: skipped " << M.nDegenerate << " degenerate triangle(s)\n";

    for (gsSurfMesh<real_t>::Halfedge h : mesh.halfedges())
    {
        if (!mesh.is_boundary(h)) continue;
        const gsSurfMesh<real_t>::Vertex vA = mesh.from_vertex(h);
        const gsSurfMesh<real_t>::Vertex vB = mesh.to_vertex(h);
        const gsSurfMesh<real_t>::Face f = mesh.face(mesh.opposite_halfedge(h));
        gsSurfMesh<real_t>::Vertex vC;
        bool haveC = false;
        for (gsSurfMesh<real_t>::Vertex v : mesh.vertices(f))
            if (v != vA && v != vB) { vC = v; haveC = true; break; }
        GISMO_ENSURE(haveC, "readPlanarMesh: boundary edge's adjacent face is not a triangle");

        const gsSurfMesh<real_t>::Point & A = mesh.position(vA);
        const gsSurfMesh<real_t>::Point & B = mesh.position(vB);
        const gsSurfMesh<real_t>::Point & C = mesh.position(vC);
        gsMatrix<real_t> E(2,3);
        E(0,0)=A[0]; E(1,0)=A[1];
        E(0,1)=B[0]; E(1,1)=B[1];
        E(0,2)=C[0]; E(1,2)=C[1];
        M.bdrEdge.push_back(E);
    }

    return M;
}

/// A polygon as a plain list of (x,y) pairs.
typedef std::vector<std::pair<real_t,real_t> > Poly2;

/// x-clip intersection, exact at x = \a xc, interpolated from the endpoint
/// with the smaller x, so the result does not depend on the edge's
/// direction. Two cells sharing the line x = xc get bitwise-identical points
/// only when both clip the same segment; when an earlier clip has already
/// shortened the edge in one of them, the two points agree to rounding only
/// (an O(eps) gap/overlap, absorbed by the V7a partition tolerance).
std::pair<real_t,real_t> intersectX(const std::pair<real_t,real_t> & S, const std::pair<real_t,real_t> & E, real_t xc)
{
    if (S.first <= E.first)
        return std::make_pair(xc, S.second + (xc-S.first)/(E.first-S.first)*(E.second-S.second));
    return std::make_pair(xc, E.second + (xc-E.first)/(S.first-E.first)*(S.second-E.second));
}

/// y-clip intersection, same convention as intersectX with the roles of x
/// and y exchanged.
std::pair<real_t,real_t> intersectY(const std::pair<real_t,real_t> & S, const std::pair<real_t,real_t> & E, real_t yc)
{
    if (S.second <= E.second)
        return std::make_pair(S.first + (yc-S.second)/(E.second-S.second)*(E.first-S.first), yc);
    return std::make_pair(E.first + (yc-E.second)/(S.second-E.second)*(S.first-E.first), yc);
}

/// One Sutherland-Hodgman half-plane clip: keeps every vertex satisfying
/// \a inside, plus the exact intersection point (\a isect) at every edge
/// where \a inside flips.
Poly2 clipHalfPlane(const Poly2 & in, const std::function<bool(real_t,real_t)> & inside,
                    const std::function<std::pair<real_t,real_t>(const std::pair<real_t,real_t> &,
                                                                  const std::pair<real_t,real_t> &)> & isect)
{
    const size_t m = in.size();
    Poly2 out;
    if (0 == m) return out;
    out.reserve(m+1);
    for (size_t k = 0; k != m; ++k)
    {
        const std::pair<real_t,real_t> & S = in[k];
        const std::pair<real_t,real_t> & E = in[(k+1)%m];
        const bool sIn = inside(S.first, S.second), eIn = inside(E.first, E.second);
        if (sIn) out.push_back(S);
        if (sIn != eIn) out.push_back(isect(S,E));
    }
    return out;
}

/// Clips triangle (x0,y0)-(x1,y1)-(x2,y2) against the CLOSED box
/// [xL,xR] x [yB,yT] by Sutherland-Hodgman, in the fixed order x>=xL,
/// x<=xR, y>=yB, y<=yT (every inside test non-strict). A triangle clipped
/// by 4 half-planes has at most 7 vertices; the result may still carry
/// duplicate or collinear vertices, left for the caller's fan triangulation
/// to skip (zero-determinant fan triangles).
Poly2 clipTriangleToBox(real_t x0, real_t y0, real_t x1, real_t y1, real_t x2, real_t y2,
                        real_t xL, real_t xR, real_t yB, real_t yT)
{
    Poly2 poly;
    poly.push_back(std::make_pair(x0,y0));
    poly.push_back(std::make_pair(x1,y1));
    poly.push_back(std::make_pair(x2,y2));

    poly = clipHalfPlane(poly,
        [xL](real_t x, real_t) { return x >= xL; },
        [xL](const std::pair<real_t,real_t> & S, const std::pair<real_t,real_t> & E) { return intersectX(S,E,xL); });
    poly = clipHalfPlane(poly,
        [xR](real_t x, real_t) { return x <= xR; },
        [xR](const std::pair<real_t,real_t> & S, const std::pair<real_t,real_t> & E) { return intersectX(S,E,xR); });
    poly = clipHalfPlane(poly,
        [yB](real_t, real_t y) { return y >= yB; },
        [yB](const std::pair<real_t,real_t> & S, const std::pair<real_t,real_t> & E) { return intersectY(S,E,yB); });
    poly = clipHalfPlane(poly,
        [yT](real_t, real_t y) { return y <= yT; },
        [yT](const std::pair<real_t,real_t> & S, const std::pair<real_t,real_t> & E) { return intersectY(S,E,yT); });
    return poly;
}

/// Diagnostics filled by meshClipQuadrature(), and the per-triangle
/// partition data V7a needs: \a triClippedArea is the Kahan-summed total
/// clipped fan area of triangle \a t over EVERY cell it overlaps,
/// \a triArea its own unclipped area, \a triDiam its longest edge.
struct MeshClipStats
{
    index_t nTriangles = 0, nDegenerate = 0, nBoundaryEdges = 0;
    index_t nClipPolys = 0, nFanTris = 0, nBdrPieces = 0;
    index_t nCut = 0, nFull = 0, nEmpty = 0;
    std::vector<real_t> triClippedArea, triArea, triDiam;
};

/// Exact cell-by-cell integration of Omega = the planar triangle mesh \a M
/// by clipping: no Gauss-Green identity, and it does not go through
/// gaussGreenQuadrature's own crossing search, tracker or winding walk (no
/// in/out test).
///
/// Complexity: O(#triangles x (n + overlapped cells per triangle)) for the
/// volume loop below, plus O(#boundary edges x (n + overlapped cells per
/// edge)) for the boundary loop, n = cells per direction: locating the cell
/// range costs O(n) per triangle/edge (colOf/rowOf are linear scans over the
/// grid lines), and each clip is O(1) (clipTriangleToBox: at most 7 vertices).
///
/// Volume, per triangle and per overlapped cell: Sutherland-Hodgman clip,
/// fan-triangulate the (convex) result from its own first vertex, and
/// integrate each fan triangle (Q0,Q1,Q2) with the Duffy collapse
/// x = Q0 + u(Q1-Q0) + u*v(Q2-Q1), weight = wu*wv*u*s*d (s = the triangle's
/// own orientation sign, d = det(Q1-Q0,Q2-Q0)), \a gv = 2p+1 points in each
/// of u,v. A polynomial of total degree d in x,y becomes degree <= d in v
/// and <= d+1 in u (the extra +1 from the Jacobian u), so Gauss with 2p+1
/// points (exact to degree 4p+1) integrates every x^a y^b with a+b <= 4p
/// exactly -- V7 only ever checks a+b <= 2p, as elsewhere in this file.
///
/// A Full cell (clipped area equal to its box within the V3b tolerance)
/// KEEPS its clipped rule rather than being replaced by a tensor Gauss rule
/// of (p+1)^2 points: the clipped rule is already exact to the degree
/// above, the Full label is purely
/// informational, and swapping in a tensor rule would let a sliver below
/// the classification tolerance silently change the integral. The cost is
/// that a Full cell carries (#fan triangles overlapping it)*(2p+1)^2 nodes
/// instead of (p+1)^2.
///
/// Boundary, per edge and per overlapped cell: Liang-Barsky clip of the
/// parametrised segment A+t*e, t in [0,1], against the box; the kept
/// sub-segment is integrated with a straight-segment Gauss rule, \a gb =
/// p+1 points (a product of two degree-p basis functions restricted to a
/// straight segment has degree 2p <= 2p+1). A boundary piece may land in an
/// Empty cell (an edge lying exactly on a grid line belongs to the LEFT/
/// BELOW column/row by the colOf/rowOf convention even when Omega lies on
/// the other side); this is still correct for every global check, since
/// none of them assume a boundary piece implies a Cut or Full cell.
CutCellQuadrature meshClipQuadrature(const PlanarMesh & M, const Grid & grid, index_t p,
                                     MeshClipStats & stats)
{
    stats = MeshClipStats();
    const index_t n = grid.n;
    std::vector<real_t> X(n+1), Y(n+1);
    for (index_t i = 0; i <= n; ++i) X[i] = grid.x0 + i*grid.h;
    for (index_t j = 0; j <= n; ++j) Y[j] = grid.y0 + j*grid.h;

    CutCellQuadrature Q;
    Q.grid = grid;
    Q.status.assign((size_t)n*n, (int)Empty);
    Q.vol.assign((size_t)n*n, CellRule());
    Q.bdr.assign((size_t)n*n, CellBdrRule());

    std::vector<std::vector<real_t> > vx((size_t)n*n), vy((size_t)n*n), vw((size_t)n*n);
    std::vector<std::vector<real_t> > bx((size_t)n*n), by((size_t)n*n), bw((size_t)n*n),
                                       bnx((size_t)n*n), bny((size_t)n*n);
    std::vector<char> touched((size_t)n*n, 0);

    gsGaussRule<real_t> gv(2*p+1);
    const gsMatrix<real_t> & zv = gv.referenceNodes();
    const gsVector<real_t> & wv = gv.referenceWeights();
    const index_t mv = 2*p+1;
    std::vector<real_t> uv(mv), wuv(mv);
    for (index_t k = 0; k != mv; ++k) { uv[k] = 0.5*(1.0+zv(0,k)); wuv[k] = 0.5*wv[k]; }

    gsGaussRule<real_t> gb(p+1);
    const gsMatrix<real_t> & zb = gb.referenceNodes();
    const gsVector<real_t> & wb = gb.referenceWeights();
    const index_t mb = p+1;
    std::vector<real_t> ub(mb), wub(mb);
    for (index_t k = 0; k != mb; ++k) { ub[k] = 0.5*(1.0+zb(0,k)); wub[k] = 0.5*wb[k]; }

    stats.nTriangles    = (index_t)M.tri.size();
    stats.nDegenerate   = M.nDegenerate;
    stats.nBoundaryEdges = (index_t)M.bdrEdge.size();
    stats.triClippedArea.assign(M.tri.size(), 0.0);
    stats.triArea.assign(M.tri.size(), 0.0);
    stats.triDiam.assign(M.tri.size(), 0.0);

    for (size_t t = 0; t != M.tri.size(); ++t)
    {
        const gsMatrix<real_t> & T = M.tri[t];
        const real_t x0 = T(0,0), y0 = T(1,0);
        const real_t x1 = T(0,1), y1 = T(1,1);
        const real_t x2 = T(0,2), y2 = T(1,2);
        const real_t det0 = (x1-x0)*(y2-y0) - (y1-y0)*(x2-x0);
        const real_t s = (det0 > 0.0) ? 1.0 : -1.0;
        stats.triArea[t] = 0.5*math::abs(det0);
        const real_t e01 = math::sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0));
        const real_t e12 = math::sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
        const real_t e20 = math::sqrt((x0-x2)*(x0-x2)+(y0-y2)*(y0-y2));
        stats.triDiam[t] = math::max(e01, math::max(e12,e20));

        const real_t xmin = math::min(x0, math::min(x1,x2)), xmax = math::max(x0, math::max(x1,x2));
        const real_t ymin = math::min(y0, math::min(y1,y2)), ymax = math::max(y0, math::max(y1,y2));
        const index_t i0 = colOf(xmin, X), i1 = colOf(xmax, X);
        const index_t j0 = rowOf(ymin, Y), j1 = rowOf(ymax, Y);

        KahanSum triClipped;
        for (index_t j = j0; j <= j1; ++j)
        for (index_t i = i0; i <= i1; ++i)
        {
            const Poly2 poly = clipTriangleToBox(x0,y0,x1,y1,x2,y2, X[i],X[i+1], Y[j],Y[j+1]);
            if (poly.size() < 3) continue;
            ++stats.nClipPolys;
            const index_t id = i + n*j;
            const real_t V0x = poly[0].first, V0y = poly[0].second;

            for (size_t k = 1; k+1 < poly.size(); ++k)
            {
                const real_t Vkx = poly[k].first,   Vky = poly[k].second;
                const real_t Vk1x = poly[k+1].first, Vk1y = poly[k+1].second;
                const real_t dk = (Vkx-V0x)*(Vk1y-V0y) - (Vky-V0y)*(Vk1x-V0x);
                if (0.0 == dk) continue;
                ++stats.nFanTris;
                touched[id] = 1;
                triClipped.add(s*dk*0.5);

                const real_t dvx = Vk1x-Vkx, dvy = Vk1y-Vky;   // Q2-Q1 direction of the u*v term
                for (index_t iu = 0; iu != mv; ++iu)
                {
                    const real_t u = uv[iu];
                    const real_t Qux = V0x + u*(Vkx-V0x), Quy = V0y + u*(Vky-V0y);
                    for (index_t iv = 0; iv != mv; ++iv)
                    {
                        const real_t v = uv[iv];
                        const real_t nx = Qux + u*v*dvx, ny = Quy + u*v*dvy;
                        const real_t wgt = wuv[iu]*wuv[iv]*u*s*dk;
                        vx[id].push_back(nx); vy[id].push_back(ny); vw[id].push_back(wgt);
                    }
                }
            }
        }
        stats.triClippedArea[t] = triClipped.value();
    }

    for (size_t e = 0; e != M.bdrEdge.size(); ++e)
    {
        const gsMatrix<real_t> & E = M.bdrEdge[e];
        const real_t Ax = E(0,0), Ay = E(1,0);
        const real_t Bx = E(0,1), By = E(1,1);
        const real_t Cx = E(0,2), Cy = E(1,2);
        const real_t ex = Bx-Ax, ey = By-Ay;
        const real_t len = math::sqrt(ex*ex+ey*ey);
        real_t nrmx = ey/len, nrmy = -ex/len;
        // Omega lies on the side of C: flip so the normal points away from it.
        if (nrmx*(Cx-Ax) + nrmy*(Cy-Ay) > 0.0) { nrmx = -nrmx; nrmy = -nrmy; }

        const real_t xmin = math::min(Ax,Bx), xmax = math::max(Ax,Bx);
        const real_t ymin = math::min(Ay,By), ymax = math::max(Ay,By);
        const index_t i0 = colOf(xmin, X), i1 = colOf(xmax, X);
        const index_t j0 = rowOf(ymin, Y), j1 = rowOf(ymax, Y);

        for (index_t j = j0; j <= j1; ++j)
        for (index_t i = i0; i <= i1; ++i)
        {
            real_t t0 = 0.0, t1 = 1.0;
            bool reject = false;
            const real_t pq[4][2] = {
                { -ex, Ax-X[i]   },
                {  ex, X[i+1]-Ax },
                { -ey, Ay-Y[j]   },
                {  ey, Y[j+1]-Ay }
            };
            for (int k = 0; k != 4 && !reject; ++k)
            {
                const real_t pk = pq[k][0], qk = pq[k][1];
                if (0.0 == pk) { if (qk < 0.0) reject = true; }
                else if (pk < 0.0) t0 = math::max(t0, qk/pk);
                else                t1 = math::min(t1, qk/pk);
            }
            if (reject || !(t1 > t0)) continue;

            const index_t id = i + n*j;
            ++stats.nBdrPieces;
            for (index_t k = 0; k != mb; ++k)
            {
                const real_t tq = t0 + ub[k]*(t1-t0);
                bx[id].push_back(Ax + tq*ex); by[id].push_back(Ay + tq*ey);
                bw[id].push_back(wub[k]*(t1-t0)*len);
                bnx[id].push_back(nrmx); bny[id].push_back(nrmy);
            }
        }
    }

    const real_t L = math::max(math::max(math::abs(grid.x0), math::abs(grid.x0+n*grid.h)),
                               math::max(math::abs(grid.y0), math::abs(grid.y0+n*grid.h)));
    const real_t tolAbs = 16*std::numeric_limits<real_t>::epsilon()*L*grid.h;

    for (index_t j = 0; j != n; ++j)
    for (index_t i = 0; i != n; ++i)
    {
        const index_t id = i + n*j;
        CellRule vr;
        vr.nodes.resize(2, (index_t)vx[id].size());
        vr.weights.resize((index_t)vx[id].size());
        for (size_t k = 0; k != vx[id].size(); ++k)
        { vr.nodes(0,(index_t)k) = vx[id][k]; vr.nodes(1,(index_t)k) = vy[id][k]; vr.weights[(index_t)k] = vw[id][k]; }
        Q.vol[id] = vr;

        CellBdrRule br;
        br.nodes.resize(2, (index_t)bx[id].size());
        br.weights.resize((index_t)bx[id].size());
        br.normals.resize(2, (index_t)bx[id].size());
        for (size_t k = 0; k != bx[id].size(); ++k)
        {
            br.nodes(0,(index_t)k) = bx[id][k]; br.nodes(1,(index_t)k) = by[id][k];
            br.weights[(index_t)k] = bw[id][k];
            br.normals(0,(index_t)k) = bnx[id][k]; br.normals(1,(index_t)k) = bny[id][k];
        }
        Q.bdr[id] = br;

        if (!touched[id]) { Q.status[id] = Empty; ++stats.nEmpty; continue; }

        const real_t box = (X[i+1]-X[i])*(Y[j+1]-Y[j]);
        const real_t area = cellArea(Q.vol[id]);
        if (area >= box-(1e-14*box+tolAbs) && area <= box*(1.0+1e-14)+tolAbs)
        { Q.status[id] = Full; ++stats.nFull; }
        else
        { Q.status[id] = Cut; ++stats.nCut; }
    }

    return Q;
}

/// Volume oracle: per triangle, polygonMoment() (its own independent code
/// path -- a signed fan from the ORIGIN, Gauss p+2, unrelated to the
/// clipped-and-collapsed rule above) over the triangle re-oriented CCW,
/// Kahan-summed over triangles. A_mesh = meshVolumeMoment(M,0,0,p).
real_t meshVolumeMoment(const PlanarMesh & M, int a, int b, index_t p)
{
    KahanSum total;
    for (const gsMatrix<real_t> & T : M.tri)
    {
        const real_t det = (T(0,1)-T(0,0))*(T(1,2)-T(1,0)) - (T(1,1)-T(1,0))*(T(0,2)-T(0,0));
        gsMatrix<real_t> Tccw(2,3);
        if (det >= 0.0) Tccw = T;
        else { Tccw.col(0) = T.col(0); Tccw.col(1) = T.col(2); Tccw.col(2) = T.col(1); }
        std::vector<gsMatrix<real_t> > loop(1, Tccw);
        total.add(polygonMoment(loop, a, b, p));
    }
    return total.value();
}

/// Boundary oracle: per boundary edge, straight-segment Gauss p+2 over the
/// WHOLE edge (no clipping), Kahan-summed -- the same per-edge formula as
/// polygonBoundaryMoment(), applied edge by edge since M.bdrEdge is not a
/// closed vertex loop (its third column is the adjacent triangle's apex C,
/// not a boundary vertex). P_mesh = meshBoundaryMoment(M,0,0,p).
real_t meshBoundaryMoment(const PlanarMesh & M, int a, int b, index_t p)
{
    gsGaussRule<real_t> gauss(p+2);
    const gsMatrix<real_t> & z = gauss.referenceNodes();
    const gsVector<real_t> & w = gauss.referenceWeights();
    const index_t m = p+2;

    KahanSum total;
    for (const gsMatrix<real_t> & E : M.bdrEdge)
    {
        const real_t Ax = E(0,0), Ay = E(1,0), Bx = E(0,1), By = E(1,1);
        const real_t len = math::sqrt((Bx-Ax)*(Bx-Ax)+(By-Ay)*(By-Ay));
        for (index_t k = 0; k != m; ++k)
        {
            const real_t s = 0.5*(1.0+z(0,k));
            const real_t x = Ax + s*(Bx-Ax), y = Ay + s*(By-Ay);
            total.add(0.5*w[k] * std::pow(x,a)*std::pow(y,b) * len);
        }
    }
    return total.value();
}

/// V7: checks the mesh-clipping rule \a Q against the mesh \a M's own
/// oracles (never against loops -- the mesh source builds none) and its own
/// per-triangle partition. Prints one line per item, PASS/FAIL/INFO as V3
/// does, and returns the AND of every PASS item.
bool checkV7(const CutCellQuadrature & Q, const PlanarMesh & M, const MeshClipStats & st,
            const DefaultGeometry & geo, index_t p)
{
    const index_t chi = M.nV - M.nE + M.nF;
    gsInfo << "V7 mesh: V=" << M.nV << " E=" << M.nE << " F=" << M.nF
           << " boundary edges=" << st.nBoundaryEdges << " chi=V-E+F=" << chi << "  INFO\n";

    const real_t L = math::max(math::max(math::abs(Q.grid.x0), math::abs(Q.grid.x0+Q.grid.n*Q.grid.h)),
                               math::max(math::abs(Q.grid.y0), math::abs(Q.grid.y0+Q.grid.n*Q.grid.h)));
    const real_t eps = std::numeric_limits<real_t>::epsilon();

    real_t v7a = 0; bool passA = true;
    for (size_t t = 0; t != st.triArea.size(); ++t)
    {
        const real_t diff  = math::abs(st.triClippedArea[t] - st.triArea[t]);
        const real_t bound = 1e-13*st.triArea[t] + 16*eps*L*st.triDiam[t];
        if (diff > bound) passA = false;
        v7a = math::max(v7a, diff / st.triArea[t]);
    }
    gsInfo << "V7a: " << fmtSci(v7a) << "  " << (passA ? "PASS" : "FAIL") << "\n";

    const std::vector<std::pair<int,int> > moments = momentSet(p);
    const real_t A_mesh = meshVolumeMoment(M, 0, 0, p);
    const real_t P_mesh = meshBoundaryMoment(M, 0, 0, p);

    real_t v7b = 0;
    for (const auto & ab : moments)
        v7b = math::max(v7b, scaledErr(ggVolumeMoment(Q,ab.first,ab.second),
                                       meshVolumeMoment(M,ab.first,ab.second,p), A_mesh));
    const bool passB = v7b <= 1e-13;
    gsInfo << "V7b: " << fmtSci(v7b) << "  " << (passB ? "PASS" : "FAIL") << "\n";

    real_t v7c = scaledErr(ggBoundaryMoment(Q,0,0), P_mesh, P_mesh);
    for (const auto & ab : moments)
        v7c = math::max(v7c, scaledErr(ggBoundaryMoment(Q,ab.first,ab.second),
                                       meshBoundaryMoment(M,ab.first,ab.second,p), P_mesh));
    const bool passC = v7c <= 1e-13;
    gsInfo << "V7c: " << fmtSci(v7c) << "  " << (passC ? "PASS" : "FAIL") << "\n";

    const index_t n = Q.grid.n;
    std::vector<real_t> X(n+1), Y(n+1);
    for (index_t i = 0; i <= n; ++i) X[i] = Q.grid.x0 + i*Q.grid.h;
    for (index_t j = 0; j <= n; ++j) Y[j] = Q.grid.y0 + j*Q.grid.h;
    const real_t tolAbs = 16*eps*L*Q.grid.h;
    const real_t hh = Q.grid.h*Q.grid.h;
    real_t cutMin = std::numeric_limits<real_t>::infinity();
    real_t cutMax = -std::numeric_limits<real_t>::infinity();
    bool passD = true;
    for (size_t id = 0; id != Q.vol.size(); ++id)
    {
        const index_t i = (index_t)id % n, j = (index_t)id / n;
        const real_t box = (X[i+1]-X[i])*(Y[j+1]-Y[j]);
        const real_t area = cellArea(Q.vol[id]);
        if (area < -(1e-14*box+tolAbs) || area > box*(1.0+1e-14)+tolAbs) passD = false;
        if (Cut == Q.status[id]) { cutMin = math::min(cutMin, area/hh); cutMax = math::max(cutMax, area/hh); }
    }
    gsInfo << "V7d: cut-cell min " << fmtSci(cutMin) << " max " << fmtSci(cutMax)
           << "  " << (passD ? "PASS" : "FAIL") << "\n";

    KahanSum sx, sy;
    for (const CellBdrRule & r : Q.bdr)
        for (index_t k = 0; k != r.weights.size(); ++k)
        { sx.add(r.weights[k]*r.normals(0,k)); sy.add(r.weights[k]*r.normals(1,k)); }
    const real_t v7e = math::max(math::abs(sx.value()), math::abs(sy.value()));
    const bool passE = v7e <= 1e-13;
    gsInfo << "V7e: " << fmtSci(v7e) << "  " << (passE ? "PASS" : "FAIL") << "\n";

    KahanSum xdotn;
    for (const CellBdrRule & r : Q.bdr)
        for (index_t k = 0; k != r.weights.size(); ++k)
            xdotn.add(r.weights[k]*(r.nodes(0,k)*r.normals(0,k) + r.nodes(1,k)*r.normals(1,k)));
    const real_t areaGG = ggVolumeMoment(Q,0,0);
    const real_t v7f = math::abs(xdotn.value() - 2.0*areaGG) / math::max((real_t)1.0, areaGG);
    const bool passF = v7f <= 1e-13;
    gsInfo << "V7f: " << fmtSci(v7f) << "  " << (passF ? "PASS" : "FAIL") << "\n";

    const real_t exactArea = geo.exactArea();
    const real_t infoDiff  = scaledErr(A_mesh, exactArea, exactArea);
    gsInfo << "V7 info: mesh area=" << fmtSci(A_mesh) << " exact geometry area=" << fmtSci(exactArea)
           << " scaled diff=" << fmtSci(infoDiff)
           << "  INFO (meaningful only when the mesh approximates DefaultGeometry)\n";

    return passA && passB && passC && passD && passE && passF;
}

/// --source mesh entry point. Builds the grid the same way main() does,
/// reads and clip-integrates the mesh, prints the header/stats/checks, and
/// always skips the Algoim baseline -- see the file header for why a
/// per-cell diff against DefaultGeometry's level set would not measure
/// quadrature error here.
int runMeshSource(const std::string & file, index_t p, index_t r, index_t n0,
                  const DefaultGeometry & geo, bool noAlgoim, bool plot, const std::string & out)
{
    const index_t n = n0 * (index_t(1) << r);
    Grid grid{ -1.0, -1.0, 2.0/(real_t)n, n };

    const PlanarMesh M = readPlanarMesh(file);
    GISMO_ENSURE(M.xmin > grid.x0 && M.xmax < grid.x0 + n*grid.h &&
                M.ymin > grid.y0 && M.ymax < grid.y0 + n*grid.h,
                "runMeshSource: mesh bbox [" << M.xmin << "," << M.xmax << "] x ["
                  << M.ymin << "," << M.ymax << "] does not lie strictly inside the "
                  "background box [" << grid.x0 << "," << grid.x0+n*grid.h << "] x ["
                  << grid.y0 << "," << grid.y0+n*grid.h << "]");

    gsStopwatch clk;
    MeshClipStats stats;
    CutCellQuadrature Q = meshClipQuadrature(M, grid, p, stats);
    const real_t ggTime = clk.stop();

    GGStats ggStats;
    ggStats.nCut = stats.nCut; ggStats.nFull = stats.nFull; ggStats.nEmpty = stats.nEmpty;
    ggStats.nPieces = stats.nClipPolys;
    printHeader(p, r, grid, 0, 0, 0, "mesh", Q, ggStats, ggTime);
    gsInfo << "mesh: file=" << file << " triangles=" << stats.nTriangles
           << " degenerate=" << stats.nDegenerate << " boundary edges=" << stats.nBoundaryEdges
           << " clipped polys=" << stats.nClipPolys << " fan triangles=" << stats.nFanTris
           << " boundary pieces=" << stats.nBdrPieces << "\n";

    const bool pass = checkV7(Q, M, stats, geo, p);

    if (!noAlgoim)
    {
        gsWarn << "Algoim baseline skipped: --source mesh is a faceted approximation of an "
                  "arbitrary -f mesh; phi describes DefaultGeometry, so a per-cell diff would "
                  "measure mesh discretisation error, not quadrature error.\n";
        gsInfo << "Algoim baseline skipped\n";
    }

    if (plot) writePlot(Q, nullptr, nullptr, out);

    return pass ? EXIT_SUCCESS : EXIT_FAILURE;
}

} // anonymous namespace

int main(int argc, char *argv[])
{
    index_t p             = 2;
    index_t r             = 2;
    index_t n0            = 8;
    index_t nq            = 12;
    index_t nsamp         = 64;
    index_t polyN         = 0;
    std::string source    = "hand";
    real_t  thetaDeg      = 17.0;
    real_t  halfA         = 0.7;
    index_t maxDepth      = 6;
    real_t  indicatorTol  = 1e-9;
    real_t  lipschitz     = 200;
    bool    noAlgoim      = false;
    bool    plot          = false;
    std::string outFolder = "output_immersed_gauss_green";
    std::string fileName, writeXml, writeCadFile;
    bool    assumeDefault = false;

    gsCmdLine cmd("Gauss-Green cut-cell quadrature from boundary curves (2D), verified against "
                 "method-independent oracles and an Algoim adaptive-rule baseline.");
    cmd.addInt   ("k", "degree",  "Background spline degree p", p);
    cmd.addInt   ("r", "refine",  "Refinement level: n0*2^r cells per direction", r);
    cmd.addInt   ("",  "n0",      "Cells per direction at r = 0", n0);
    cmd.addInt   ("",  "nq",      "Gauss points per curve piece / right-edge interval", nq);
    cmd.addInt   ("",  "nsamp",   "Uniform samples per curve for the crossing search", nsamp);
    cmd.addInt   ("",  "polyline","Replace every curve by N straight segments (0 = exact curves)", polyN);
    cmd.addString("",  "source",  "Geometry source: hand | xml | occ | mesh", source);
    cmd.addString("f", "file",    "Geometry file: .xml (--source xml), .step/.stp/.brep (--source occ), "
                                 "or .stl/.off/.obj (--source mesh)", fileName);
    cmd.addString("",  "write-xml", "Write the hand-built DefaultGeometry as a <PlanarDomain> to this file", writeXml);
    cmd.addString("",  "write-cad", "Write the OCC shape (--source occ) to this .brep/.step/.stp file", writeCadFile);
    cmd.addSwitch("assume-default", "The -f geometry IS DefaultGeometry (--theta/--half): enables V2, V5/V6 and the Algoim baseline", assumeDefault);
    cmd.addReal  ("",  "theta",   "Square rotation in degrees (DefaultGeometry)", thetaDeg);
    cmd.addReal  ("",  "half",    "Square half-side length a (DefaultGeometry)", halfA);
    cmd.addInt   ("",  "maxDepth",     "Algoim adaptive maximum depth", maxDepth);
    cmd.addReal  ("",  "indicatorTol", "Algoim integralChange tolerance", indicatorTol);
    cmd.addReal  ("",  "lipschitz",    "Lipschitz constant for Algoim box classification", lipschitz);
    cmd.addSwitch("no-algoim", "Skip the Algoim baseline", noAlgoim);
    cmd.addSwitch("plot", "Write ParaView point sets of the quadrature nodes", plot);
    cmd.addString("o", "output", "Output folder for --plot", outFolder);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    if (p < 1)      { gsWarn << "-k/--degree must be >= 1\n"; return EXIT_FAILURE; }
    if (r < 0)      { gsWarn << "-r/--refine must be >= 0\n"; return EXIT_FAILURE; }
    if (n0 < 1)     { gsWarn << "--n0 must be >= 1\n"; return EXIT_FAILURE; }
    if (nq < 1)     { gsWarn << "--nq must be >= 1\n"; return EXIT_FAILURE; }
    if (nsamp < 2)  { gsWarn << "--nsamp must be >= 2\n"; return EXIT_FAILURE; }
    if (polyN < 0)  { gsWarn << "--polyline must be >= 0\n"; return EXIT_FAILURE; }
    if (maxDepth < 0) { gsWarn << "--maxDepth must be >= 0\n"; return EXIT_FAILURE; }
    if (halfA <= 0.0) { gsWarn << "--half must be > 0\n"; return EXIT_FAILURE; }
    if ("hand" != source && "xml" != source && "occ" != source && "mesh" != source)
    { gsWarn << "--source must be one of hand|xml|occ|mesh\n"; return EXIT_FAILURE; }
    if ("mesh" == source && polyN > 0)
    { gsWarn << "--polyline does not apply to --source mesh\n"; return EXIT_FAILURE; }
    if ("xml" == source && fileName.empty())
    { gsWarn << "--source xml requires -f <file.xml>\n"; return EXIT_FAILURE; }
    if ("hand" == source && !fileName.empty())
    { gsWarn << "--source hand: -f is not used\n"; return EXIT_FAILURE; }
    if (!writeCadFile.empty() && "occ" != source)
    { gsWarn << "--write-cad requires --source occ\n"; return EXIT_FAILURE; }
    if (!writeCadFile.empty())
    {
        std::string wext = gsFileManager::getExtension(writeCadFile);
        std::transform(wext.begin(), wext.end(), wext.begin(), ::tolower);
        if ("brep" != wext && "step" != wext && "stp" != wext)
        { gsWarn << "--write-cad extension must be brep, step or stp\n"; return EXIT_FAILURE; }
    }

    std::string resolvedFile;
    if (!fileName.empty())
    {
        resolvedFile = gsFileManager::find(fileName);
        if (resolvedFile.empty())
        { gsWarn << "-f " << fileName << " not found\n"; return EXIT_FAILURE; }
        if ("occ" == source)
        {
            std::string fext = gsFileManager::getExtension(resolvedFile);
            std::transform(fext.begin(), fext.end(), fext.begin(), ::tolower);
            if ("brep" != fext && "step" != fext && "stp" != fext)
            { gsWarn << "-f extension for --source occ must be brep, step or stp\n"; return EXIT_FAILURE; }
        }
    }

    std::string outPath = outFolder;
    if (outPath.empty()) outPath = "output_immersed_gauss_green";
    const bool isAbsolutePath = (!outPath.empty() && outPath[0] == '/');
    if (!isAbsolutePath) outPath = gsFileManager::getCurrentPath() + "/" + outPath;
    const std::string out = gsFileManager::getCanonicRepresentation(outPath);
    if (plot) gsFileManager::mkdir(out);

    DefaultGeometry geo(halfA, thetaDeg);

    if ("mesh" == source)
        return runMeshSource(fileName.empty() ? "planar/square_minus_disk_mesh.stl" : fileName,
                             p, r, n0, geo, noAlgoim, plot, out);

    if (!writeXml.empty()) writeDefaultXml(geo, writeXml);

    TopoDS_Shape occShape;
    std::vector<Loop> loops = makeLoops(source, geo, resolvedFile, &occShape);

    if (!writeCadFile.empty()) writeCad(occShape, writeCadFile);

    std::vector<gsMatrix<real_t> > polyVertices;
    const bool polylineMode = (polyN > 0);
    if (polylineMode)
        loops = toPolyline(loops, polyN, polyVertices);

    normalizeOrientation(loops, nsamp);

    const index_t n = n0 * (index_t(1) << r);
    Grid grid{ -1.0, -1.0, 2.0/(real_t)n, n };

    ensureInsideBox(loops, grid, nsamp);

    gsStopwatch clk;
    GGStats stats;
    CutCellQuadrature Q = gaussGreenQuadrature(loops, grid, p, nq, nsamp, stats);
    const real_t ggTime = clk.stop();

    printHeader(p, r, grid, nq, nsamp, polyN, source, Q, stats, ggTime);
    if (!fileName.empty()) gsInfo << "file=" << resolvedFile << "\n";

    // True for "hand" (no -f) and for the built-in "occ" default (no -f);
    // false for an -f file unless the caller vouches for it with
    // --assume-default. Gates every check whose oracle is DefaultGeometry's
    // own numbers (V2, V5/V6, the Algoim baseline) -- a foreign file is not
    // known to be that geometry.
    const bool geometryIsDefault = fileName.empty() || assumeDefault;

    bool allPass = true;

    if (polylineMode)
        allPass = checkV1(Q, polyVertices, p) && allPass;
    else
    {
        const std::vector<Loop> & loopsRef = loops;
        auto producer = [&loopsRef, &grid, p, nsamp](index_t nqTest)
        {
            GGStats dummy;
            return gaussGreenQuadrature(loopsRef, grid, p, nqTest, nsamp, dummy);
        };
        if (geometryIsDefault)
            allPass = checkV2(producer, geo, p) && allPass;
        else
            gsInfo << "V2 skipped: geometry from -f is not known to be DefaultGeometry "
                      "(pass --assume-default)\n";

        if ("hand" != source && geometryIsDefault)
        {
            std::vector<Loop> handLoops = handBuiltLoops(geo);
            normalizeOrientation(handLoops, nsamp, false);
            GGStats dummyHand;
            const CutCellQuadrature Qhand = gaussGreenQuadrature(handLoops, grid, p, 12, nsamp, dummyHand);
            const CutCellQuadrature Qsrc  = producer(12);
            const std::string label = ("occ" == source) ? "V5 occ vs hand" : "V6 xml vs hand";
            const real_t tol        = ("occ" == source) ? 1e-12 : 1e-14;
            allPass = checkVsHand(Qsrc, Qhand, geo, p, label, tol) && allPass;
        }
    }

    const std::vector<std::pair<int,int> > mvol = momentSet(p);
    std::vector<real_t> refMoments(mvol.size());
    for (size_t k = 0; k != mvol.size(); ++k)
        refMoments[k] = globalMoment(loops, mvol[k].first, mvol[k].second, nsamp);
    const real_t areaScale = geometryIsDefault ? geo.exactArea() : math::abs(refMoments[0]);
    allPass = checkV3(Q, refMoments, mvol, areaScale) && allPass;

    gsMatrix<real_t> algoimNodes;
    gsVector<real_t> algoimWeights;
    bool algoimSkipped = false;
    if (!noAlgoim)
    {
        if (!geometryIsDefault)
        {
            gsInfo << "Algoim baseline skipped: geometry from -f is not known to be DefaultGeometry "
                      "(pass --assume-default)\n";
            algoimSkipped = true;
        }
        else if (algoimBaselineUsable(geo, grid))
            runAlgoimBaseline(geo, grid, p, maxDepth, indicatorTol, lipschitz, Q, polylineMode,
                              algoimNodes, algoimWeights);
        else
        {
            gsWarn << "Algoim baseline skipped: a square edge lies exactly on a background "
                      "grid line; gsAlgoimAdaptiveRule::mapTo exhausts memory when phi "
                      "vanishes on a whole cell face.\n";
            gsInfo << "Algoim baseline skipped\n";
            algoimSkipped = true;
        }
    }

    const bool haveAlgoim = !noAlgoim && !algoimSkipped;
    if (plot)
        writePlot(Q, haveAlgoim ? &algoimNodes : nullptr, haveAlgoim ? &algoimWeights : nullptr, out);

    return allPass ? EXIT_SUCCESS : EXIT_FAILURE;
}
