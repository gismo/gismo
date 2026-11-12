#!/usr/bin/env python3
"""
gs_domain_generator.py

Parametric Planar Boundary Domain Generator for G+Smo and IGA Multi-Patch Splitting.

Generates 2D planar boundary curves and domains based on parameters:
- num_corners: Number of C^0 continuity reduction points (0, 1, 2, 3, 4, 5, 6, ...)
- convex: Whether the domain is convex (True / False)
- star_shaped: Whether the domain is star-shaped with respect to a kernel (True / False)
- num_holes: Number of interior holes (0, 1, 2, ...)
- degree: B-spline degree p (default 3)
- curved: Whether boundary edges between corners are curved or straight chords
- randomness: Stochastic perturbation to corner positions and curvatures
- seed: Random seed for reproducible benchmarks

Export formats:
1. G+Smo native <PlanarDomain> XML (compatible with gsReadFile, gsPlanarDomain, gsWriteParaview)
2. ParaView VTK XML PolyData (.vtp) and Collection (.pvd)
3. Mathematica script (.m) for interactive visualization

Author: G+Smo Multi-Patch & Domain Splitting Toolset
Compatible with Python 3.6+ without external dependencies.
"""

import math
import random
import os
import xml.etree.ElementTree as ET
from typing import List, Tuple, Dict, Any, Optional


# =============================================================================
# Mathematical Helpers
# =============================================================================

def comb(n: int, k: int) -> int:
    """Binomial coefficient C(n, k), Python 3.6 compatible."""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))


def bernstein(p: int, j: int, t: float) -> float:
    """Evaluate j-th Bernstein polynomial of degree p at t in [0, 1]."""
    if t <= 0.0:
        return 1.0 if j == 0 else 0.0
    if t >= 1.0:
        return 1.0 if j == p else 0.0
    return comb(p, j) * (t ** j) * ((1.0 - t) ** (p - j))


def bernstein_deriv(p: int, j: int, t: float) -> float:
    """Derivative of j-th Bernstein polynomial of degree p at t in [0, 1]."""
    term1 = bernstein(p - 1, j - 1, t) if j > 0 else 0.0
    term2 = bernstein(p - 1, j, t) if j < p else 0.0
    return p * (term1 - term2)


# =============================================================================
# Geometric Primitives: Point2D
# =============================================================================

class Point2D:
    """2D Euclidean point / vector."""

    def __init__(self, x: float = 0.0, y: float = 0.0):
        self.x = float(x)
        self.y = float(y)

    def __repr__(self) -> str:
        return "({:.6f}, {:.6f})".format(self.x, self.y)

    def __add__(self, other: 'Point2D') -> 'Point2D':
        return Point2D(self.x + other.x, self.y + other.y)

    def __sub__(self, other: 'Point2D') -> 'Point2D':
        return Point2D(self.x - other.x, self.y - other.y)

    def __mul__(self, scalar: float) -> 'Point2D':
        return Point2D(self.x * scalar, self.y * scalar)

    def __rmul__(self, scalar: float) -> 'Point2D':
        return Point2D(self.x * scalar, self.y * scalar)

    def __truediv__(self, scalar: float) -> 'Point2D':
        return Point2D(self.x / scalar, self.y / scalar)

    def __neg__(self) -> 'Point2D':
        return Point2D(-self.x, -self.y)

    def dot(self, other: 'Point2D') -> float:
        return self.x * other.x + self.y * other.y

    def cross(self, other: 'Point2D') -> float:
        """2D cross product: x1*y2 - y1*x2."""
        return self.x * other.y - self.y * other.x

    def norm(self) -> float:
        return math.hypot(self.x, self.y)

    def norm_sq(self) -> float:
        return self.x * self.x + self.y * self.y

    def dist(self, other: 'Point2D') -> float:
        return (self - other).norm()

    def normalized(self) -> 'Point2D':
        n = self.norm()
        if n < 1e-14:
            return Point2D(0.0, 0.0)
        return self / n

    def normal(self) -> 'Point2D':
        """Outward normal (rotated 90 deg clockwise: (y, -x) for CCW tangent)."""
        return Point2D(self.y, -self.x).normalized()

    def normal_ccw(self) -> 'Point2D':
        """Rotated 90 deg counter-clockwise: (-y, x)."""
        return Point2D(-self.y, self.x).normalized()

    def rotated(self, angle_rad: float) -> 'Point2D':
        c = math.cos(angle_rad)
        s = math.sin(angle_rad)
        return Point2D(self.x * c - self.y * s, self.x * s + self.y * c)

    def to_list(self) -> List[float]:
        return [self.x, self.y]

    def to_tuple(self) -> Tuple[float, float]:
        return (self.x, self.y)

    def __getitem__(self, idx: int) -> float:
        if idx == 0:
            return self.x
        elif idx == 1:
            return self.y
        raise IndexError("Point2D index out of range (must be 0 or 1)")

    def __len__(self) -> int:
        return 2


# =============================================================================
# B-Spline Curve Segment
# =============================================================================

class BSplineCurve:
    """
    1D to 2D B-spline curve segment.
    Uses clamped open knot vector [0]*(p+1) + [1]*(p+1) with p+1 control points.
    """

    def __init__(self, control_points: List[Point2D], degree: int = 3, knot_vector: Optional[List[float]] = None):
        self.control_points = [Point2D(p.x, p.y) for p in control_points]
        self.degree = int(degree)
        if knot_vector is None:
            self.knot_vector = [0.0] * (self.degree + 1) + [1.0] * (self.degree + 1)
        else:
            self.knot_vector = list(knot_vector)

    def start_point(self) -> Point2D:
        return self.control_points[0]

    def end_point(self) -> Point2D:
        return self.control_points[-1]

    def eval(self, t: float) -> Point2D:
        """
        Evaluate B-spline curve at parameter t.
        Strictly respects the valid knot span domain [u_front, u_back].
        """
        u_min = self.knot_vector[self.degree]
        u_max = self.knot_vector[-self.degree - 1]
        t_clamped = max(u_min, min(u_max, float(t)))
        t_norm = (t_clamped - u_min) / (u_max - u_min) if (u_max - u_min) > 1e-14 else 0.0
        p = self.degree
        res = Point2D(0.0, 0.0)
        for j, cp in enumerate(self.control_points):
            b = bernstein(p, j, t_norm)
            res = res + cp * b
        return res

    def deriv(self, t: float) -> Point2D:
        """Evaluate first derivative C'(t) strictly within valid knot span."""
        u_min = self.knot_vector[self.degree]
        u_max = self.knot_vector[-self.degree - 1]
        t_clamped = max(u_min, min(u_max, float(t)))
        t_norm = (t_clamped - u_min) / (u_max - u_min) if (u_max - u_min) > 1e-14 else 0.0
        p = self.degree
        res = Point2D(0.0, 0.0)
        for j, cp in enumerate(self.control_points):
            db = bernstein_deriv(p, j, t_norm)
            res = res + cp * db
        scaling = 1.0 / (u_max - u_min) if (u_max - u_min) > 1e-14 else 1.0
        return res * scaling

    def sample(self, num_points: int = 40) -> List[Point2D]:
        """Sample points along the curve strictly within [u_min, u_max]."""
        u_min = self.knot_vector[self.degree]
        u_max = self.knot_vector[-self.degree - 1]
        pts = []
        for i in range(num_points):
            t = u_min + (u_max - u_min) * (float(i) / float(num_points - 1)) if num_points > 1 else u_min
            pts.append(self.eval(t))
        return pts

    def reverse(self) -> 'BSplineCurve':
        """
        Return a reversed B-spline curve.
        Reverses the control points and transforms knot vector: u'_i = u_min + u_max - u_{m-i}.
        """
        rev_cps = [Point2D(p.x, p.y) for p in reversed(self.control_points)]
        u_min = self.knot_vector[0]
        u_max = self.knot_vector[-1]
        rev_knots = [u_min + u_max - k for k in reversed(self.knot_vector)]
        return BSplineCurve(rev_cps, self.degree, rev_knots)

    def signed_area_control_polygon(self) -> float:
        """Signed area contribution of control polygon."""
        acc = 0.0
        n = len(self.control_points)
        for i in range(n - 1):
            acc += self.control_points[i].x * self.control_points[i + 1].y - \
                   self.control_points[i].y * self.control_points[i + 1].x
        return 0.5 * acc


# =============================================================================
# CurveLoop (Loop of B-Spline Curves)
# =============================================================================

class CurveLoop:
    """
    Closed loop of B-spline curves.
    In G+Smo gsPlanarDomain:
    - Outer loop is oriented CCW (counter-clockwise, signed area > 0).
    - Inner loops (holes) are oriented CW (clockwise, signed area < 0).
    """

    def __init__(self, curves: Optional[List[BSplineCurve]] = None, is_outer: bool = True):
        self.curves = curves if curves is not None else []
        self.is_outer = is_outer
        self._sample_cache = {}

    def num_curves(self) -> int:
        return len(self.curves)

    def add_curve(self, curve: BSplineCurve):
        self.curves.append(curve)
        if hasattr(self, '_sample_cache'):
            self._sample_cache.clear()

    def sample_loop(self, num_points_per_curve: int = 40) -> List[Point2D]:
        """
        Sample closed polyline of points around loop.
        Sequential, closed, strictly avoids duplicate vertices at junctions and closures.
        Guarantees zero degenerate zero-length segments.
        Cached for high-performance point-in-polygon queries.
        """
        if not hasattr(self, '_sample_cache'):
            self._sample_cache = {}
        if num_points_per_curve in self._sample_cache:
            return self._sample_cache[num_points_per_curve]

        poly = []
        for i, curve in enumerate(self.curves):
            pts = curve.sample(num_points_per_curve)
            if len(pts) > 1:
                poly.extend(pts[:-1])
            elif len(pts) == 1:
                poly.extend(pts)

        # Deduplicate consecutive close points
        cleaned = []
        for pt in poly:
            if not cleaned or (pt - cleaned[-1]).norm() > 1e-8:
                cleaned.append(pt)
        if len(cleaned) > 1 and (cleaned[0] - cleaned[-1]).norm() < 1e-8:
            cleaned.pop()

        self._sample_cache[num_points_per_curve] = cleaned
        return cleaned

    def sample(self, num_points_per_curve: int = 40) -> List[Point2D]:
        """Sample closed polyline around loop."""
        return self.sample_loop(num_points_per_curve)

    def signed_area(self, num_points_per_curve: int = 50) -> float:
        """
        Compute the signed polygon area for discretized loop:
        A = 0.5 * sum(x_i * y_{i+1} - x_{i+1} * y_i).
        Returns > 0 for Counter-Clockwise (CCW), < 0 for Clockwise (CW).
        """
        poly = self.sample_loop(num_points_per_curve)
        n = len(poly)
        if n < 3:
            return 0.0
        acc = 0.0
        for i in range(n):
            p1 = poly[i]
            p2 = poly[(i + 1) % n]
            acc += p1.x * p2.y - p2.x * p1.y
        return 0.5 * acc

    def is_ccw(self, num_points_per_curve: int = 50) -> bool:
        """
        Strict B-Rep check: Signed Area > 0.
        """
        return self.signed_area(num_points_per_curve) > 0.0

    def reverse(self):
        """
        Reverse orientation of the loop:
        Reverse order of curves and reverse parameter direction of each B-spline curve.
        """
        rev_curves = [c.reverse() for c in reversed(self.curves)]
        self.curves = rev_curves
        if hasattr(self, '_sample_cache'):
            self._sample_cache.clear()

    def ensure_orientation(self, should_be_ccw: bool):
        """
        Ensure loop orientation matches required B-Rep standard:
        - Outer loop: Counter-Clockwise (CCW, signed area > 0)
        - Hole loop: Clockwise (CW, signed area < 0)
        """
        if self.is_ccw() != should_be_ccw:
            self.reverse()

    def corners(self, tol: float = 1e-4) -> List[Point2D]:
        """
        Identify corner points (junction points between curves where continuity is reduced to C^0).
        """
        corners_list = []
        n = len(self.curves)
        if n == 0:
            return corners_list

        for i in range(n):
            c_prev = self.curves[(i - 1) % n]
            c_curr = self.curves[i]
            # Junction point
            pt = c_curr.start_point()
            # Check incoming and outgoing tangent directions
            v_in = (c_prev.control_points[-1] - c_prev.control_points[-2]).normalized()
            v_out = (c_curr.control_points[1] - c_curr.control_points[0]).normalized()

            # Cross product or dot product to determine if there is a kink / direction change
            cross_val = abs(v_in.cross(v_out))
            dot_val = v_in.dot(v_out)

            # If tangents are not parallel and in the same direction, it's a C^0 corner
            if cross_val > tol or dot_val < 1.0 - tol:
                corners_list.append(pt)

        return corners_list

    def bounding_box(self) -> Tuple[float, float, float, float]:
        """Returns (xmin, ymin, xmax, ymax)."""
        xmin = 1e30
        ymin = 1e30
        xmax = -1e30
        ymax = -1e30
        for curve in self.curves:
            for cp in curve.control_points:
                xmin = min(xmin, cp.x)
                ymin = min(ymin, cp.y)
                xmax = max(xmax, cp.x)
                ymax = max(ymax, cp.y)
        return (xmin, ymin, xmax, ymax)

    def contains_point(self, pt: Point2D) -> bool:
        """Winding number test for point inside closed loop."""
        poly = self.sample_loop(40)
        n = len(poly)
        if n < 3:
            return False
        wn = 0
        for i in range(n):
            p1 = poly[i]
            p2 = poly[(i + 1) % n]
            if p1.y <= pt.y:
                if p2.y > pt.y:
                    if (p2.x - p1.x) * (pt.y - p1.y) - (pt.x - p1.x) * (p2.y - p1.y) > 1e-12:
                        wn += 1
            else:
                if p2.y <= pt.y:
                    if (p2.x - p1.x) * (pt.y - p1.y) - (pt.x - p1.x) * (p2.y - p1.y) < -1e-12:
                        wn -= 1
        return wn != 0


# =============================================================================
# PlanarDomain Class
# =============================================================================

class PlanarDomain:
    """
    Representation of a Planar Domain with boundary curves and holes.
    Fully compatible with G+Smo gsPlanarDomain.
    """

    def __init__(self,
                 outer_loop: CurveLoop,
                 hole_loops: Optional[List[CurveLoop]] = None,
                 metadata: Optional[Dict[str, Any]] = None):
        self.outer_loop = outer_loop
        self.outer_loop.is_outer = True
        self.outer_loop.ensure_orientation(should_be_ccw=True)

        self.hole_loops = hole_loops if hole_loops is not None else []
        for hl in self.hole_loops:
            hl.is_outer = False
            hl.ensure_orientation(should_be_ccw=False)

        self.metadata = metadata if metadata is not None else {}
        self._cache_corners = self.outer_loop.corners()

    @property
    def num_corners(self) -> int:
        if "num_corners" in self.metadata:
            return int(self.metadata["num_corners"])
        return len(self._cache_corners)

    @property
    def num_holes(self) -> int:
        return len(self.hole_loops)

    @property
    def is_convex(self) -> bool:
        return bool(self.metadata.get("convex", False))

    @property
    def is_star_shaped(self) -> bool:
        return bool(self.metadata.get("star_shaped", True))

    def get_corners(self) -> List[Point2D]:
        """Returns corner points where continuity is reduced to C^0."""
        return self._cache_corners

    def bounding_box(self) -> Tuple[float, float, float, float]:
        return self.outer_loop.bounding_box()

    # -------------------------------------------------------------------------
    # G+Smo XML Export
    # -------------------------------------------------------------------------
    def to_gismo_xml(self, filepath: str):
        """
        Export domain to G+Smo XML format (<PlanarDomain>).
        Can be read directly in C++ via:
            gsPlanarDomain<real_t>::uPtr pd;
            gsReadFile<gsPlanarDomain<real_t>>(filename, pd);
        """
        root = ET.Element("xml")
        pd_node = ET.SubElement(root, "PlanarDomain", {"id": "0"})

        # Loop 0: Outer boundary (CCW)
        self._write_loop_xml(pd_node, self.outer_loop, index=0)

        # Loops 1..H: Holes (CW)
        for i, hole in enumerate(self.hole_loops):
            self._write_loop_xml(pd_node, hole, index=i + 1)

        tree = ET.ElementTree(root)
        self._indent_xml(root)

        # Ensure directory exists
        dirname = os.path.dirname(filepath)
        if dirname and not os.path.exists(dirname):
            os.makedirs(dirname, exist_ok=True)

        with open(filepath, "w", encoding="utf-8") as f:
            f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
            tree.write(f, encoding="unicode")

    def _write_loop_xml(self, parent_node: ET.Element, loop: CurveLoop, index: int):
        loop_node = ET.SubElement(parent_node, "CurveLoop", {"index": str(index)})
        for c_idx, curve in enumerate(loop.curves):
            geom_node = ET.SubElement(loop_node, "Geometry", {"type": "BSpline", "id": str(c_idx)})
            basis_node = ET.SubElement(geom_node, "Basis", {"type": "BSplineBasis"})

            kv_str = " ".join("{:.15g}".format(k) for k in curve.knot_vector) + " "
            kv_node = ET.SubElement(basis_node, "KnotVector", {"degree": str(curve.degree)})
            kv_node.text = kv_str

            coefs_lines = []
            for cp in curve.control_points:
                coefs_lines.append("{:.15g} {:.15g}".format(cp.x, cp.y))
            coefs_str = "\n" + "\n".join(coefs_lines) + "\n"

            coefs_node = ET.SubElement(geom_node, "coefs", {"geoDim": "2"})
            coefs_node.text = coefs_str

    def _indent_xml(self, elem: ET.Element, level: int = 0):
        indent_str = " " * (level * 2)
        if len(elem):
            if not elem.text or not elem.text.strip():
                elem.text = "\n" + indent_str + "  "
            if not elem.tail or not elem.tail.strip():
                elem.tail = "\n" + indent_str
            for sub_elem in elem:
                self._indent_xml(sub_elem, level + 1)
            if not elem.tail or not elem.tail.strip():
                elem.tail = "\n" + indent_str
        else:
            if level and (not elem.tail or not elem.tail.strip()):
                elem.tail = "\n" + indent_str

    # -------------------------------------------------------------------------
    # ParaView VTK PolyData (.vtp), UnstructuredGrid (.vtu) and Collection (.pvd)
    # -------------------------------------------------------------------------
    def to_paraview(self, filepath_base: str, num_samples_per_curve: int = 50):
        """
        Export domain to ParaView formats:
        - <base>_surface.vtp: 2D planar triangulated mesh (VTK_TRIANGLE cells, holes cleanly excluded)
        - <base>_curves.vtp: Sampled boundary and hole curves (sequentially ordered, no duplicate closures)
        - <base>_corners.vtp: Corner vertices (points of reduced C^0 continuity)
        - <base>_cnet.vtp: B-spline control polygons
        - <base>.vtu: Native VTK UnstructuredGrid with explicit VTK_TRIANGLE (type 5) cells
        - <base>.pvd: Master collection file for 1-click ParaView inspection
        """
        dirname = os.path.dirname(filepath_base)
        if dirname and not os.path.exists(dirname):
            os.makedirs(dirname, exist_ok=True)

        basename = os.path.basename(filepath_base)

        surface_file = filepath_base + "_surface.vtp"
        curves_file = filepath_base + "_curves.vtp"
        corners_file = filepath_base + "_corners.vtp"
        cnet_file = filepath_base + "_cnet.vtp"
        vtu_file = filepath_base + ".vtu"
        pvd_file = filepath_base + ".pvd"

        # 1. Compute robust triangulation once
        mesh = triangulate_polygon_with_holes(self.outer_loop, self.hole_loops, num_samples_per_curve=25)

        # 2. Export Triangulated Planar Surface (.vtp)
        self._export_vtp_surface(surface_file, num_samples_per_curve=25, precomputed_mesh=mesh)

        # 3. Export Boundary & Hole Curves (.vtp)
        self._export_vtp_curves(curves_file, num_samples_per_curve)

        # 4. Export Corners (.vtp)
        self._export_vtp_corners(corners_file)

        # 5. Export Control Net (.vtp)
        self._export_vtp_control_net(cnet_file)

        # 6. Export Native UnstructuredGrid (.vtu)
        self.to_vtu(vtu_file, num_samples_per_curve=25, precomputed_mesh=mesh)

        # 6. Export Master PVD
        with open(pvd_file, "w", encoding="utf-8") as f:
            f.write('<?xml version="1.0"?>\n')
            f.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n')
            f.write('  <Collection>\n')
            f.write(f'    <DataSet part="0" file="{basename}_surface.vtp" name="PlanarSurface"/>\n')
            f.write(f'    <DataSet part="1" file="{basename}_curves.vtp" name="BoundaryCurves"/>\n')
            f.write(f'    <DataSet part="2" file="{basename}_corners.vtp" name="Corners"/>\n')
            f.write(f'    <DataSet part="3" file="{basename}_cnet.vtp" name="ControlPolygons"/>\n')
            f.write('  </Collection>\n')
            f.write('</VTKFile>\n')

    def _export_vtp_curves(self, filename: str, num_samples: int):
        """
        Export closed boundary and hole polylines to VTK PolyData.
        Guarantees:
        - Only sampled boundary points along curves are exported (NO control points).
        - Loop vertices are ordered sequentially without repeating closure vertex twice.
        - Zero zero-length or degenerate edges.
        """
        points = []
        lines_conn = []
        lines_offsets = []
        loop_ids = []
        point_ids = []

        curr_pt_idx = 0
        all_loops = [self.outer_loop] + self.hole_loops

        for l_idx, loop in enumerate(all_loops):
            loop_pts = loop.sample_loop(num_samples)
            L = len(loop_pts)
            if L < 2:
                continue

            start_idx = curr_pt_idx
            for p_idx, pt in enumerate(loop_pts):
                points.append((pt.x, pt.y, 0.0))
                loop_ids.append(l_idx)
                point_ids.append(p_idx)
                curr_pt_idx += 1

            # Connect points sequentially along the closed loop:
            # (0, 1), (1, 2), ..., (L-2, L-1), (L-1, 0)
            for k in range(L):
                lines_conn.extend([start_idx + k, start_idx + ((k + 1) % L)])
                lines_offsets.append(len(lines_conn))

        num_points = len(points)
        num_lines = len(lines_offsets)

        pts_str = " ".join(f"{p[0]:.6g} {p[1]:.6g} {p[2]:.6g}" for p in points)
        conn_str = " ".join(str(c) for c in lines_conn)
        offs_str = " ".join(str(o) for o in lines_offsets)
        loop_str = " ".join(str(lid) for lid in loop_ids)
        pid_str = " ".join(str(pid) for pid in point_ids)

        vtp_content = f'''<?xml version="1.0"?>
<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">
  <PolyData>
    <Piece NumberOfPoints="{num_points}" NumberOfVerts="0" NumberOfLines="{num_lines}" NumberOfStrips="0" NumberOfPolys="0">
      <Points>
        <DataArray type="Float32" NumberOfComponents="3" format="ascii">
          {pts_str}
        </DataArray>
      </Points>
      <Lines>
        <DataArray type="Int32" Name="connectivity" format="ascii">
          {conn_str}
        </DataArray>
        <DataArray type="Int32" Name="offsets" format="ascii">
          {offs_str}
        </DataArray>
      </Lines>
      <PointData Scalars="LoopID">
        <DataArray type="Int32" Name="LoopID" format="ascii">
          {loop_str}
        </DataArray>
        <DataArray type="Int32" Name="PointID" format="ascii">
          {pid_str}
        </DataArray>
      </PointData>
    </Piece>
  </PolyData>
</VTKFile>'''
        with open(filename, "w", encoding="utf-8") as f:
            f.write(vtp_content)

    def _export_vtp_corners(self, filename: str):
        corners = self.get_corners()
        num_corners = len(corners)

        if num_corners == 0:
            pts_str = ""
            conn_str = ""
            offs_str = ""
            id_str = ""
        else:
            pts_str = " ".join(f"{c.x:.6g} {c.y:.6g} 0.0" for c in corners)
            conn_str = " ".join(str(i) for i in range(num_corners))
            offs_str = " ".join(str(i + 1) for i in range(num_corners))
            id_str = " ".join(str(i) for i in range(num_corners))

        vtp_content = f'''<?xml version="1.0"?>
<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">
  <PolyData>
    <Piece NumberOfPoints="{num_corners}" NumberOfVerts="{num_corners}" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="0">
      <Points>
        <DataArray type="Float32" NumberOfComponents="3" format="ascii">
          {pts_str}
        </DataArray>
      </Points>
      <Verts>
        <DataArray type="Int32" Name="connectivity" format="ascii">
          {conn_str}
        </DataArray>
        <DataArray type="Int32" Name="offsets" format="ascii">
          {offs_str}
        </DataArray>
      </Verts>
      <PointData Scalars="CornerID">
        <DataArray type="Int32" Name="CornerID" format="ascii">
          {id_str}
        </DataArray>
      </PointData>
    </Piece>
  </PolyData>
</VTKFile>'''
        with open(filename, "w", encoding="utf-8") as f:
            f.write(vtp_content)

    def _export_vtp_control_net(self, filename: str):
        points = []
        lines_conn = []
        lines_offsets = []

        curr_pt_idx = 0
        all_loops = [self.outer_loop] + self.hole_loops

        for loop in all_loops:
            for curve in loop.curves:
                cps = curve.control_points
                start_idx = curr_pt_idx
                for cp in cps:
                    points.append((cp.x, cp.y, 0.0))
                    curr_pt_idx += 1
                for k in range(len(cps) - 1):
                    lines_conn.extend([start_idx + k, start_idx + k + 1])
                    lines_offsets.append(len(lines_conn))

        num_points = len(points)
        num_lines = len(lines_offsets)

        pts_str = " ".join(f"{p[0]:.6g} {p[1]:.6g} {p[2]:.6g}" for p in points)
        conn_str = " ".join(str(c) for c in lines_conn)
        offs_str = " ".join(str(o) for o in lines_offsets)

        vtp_content = f'''<?xml version="1.0"?>
<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">
  <PolyData>
    <Piece NumberOfPoints="{num_points}" NumberOfVerts="0" NumberOfLines="{num_lines}" NumberOfStrips="0" NumberOfPolys="0">
      <Points>
        <DataArray type="Float32" NumberOfComponents="3" format="ascii">
          {pts_str}
        </DataArray>
      </Points>
      <Lines>
        <DataArray type="Int32" Name="connectivity" format="ascii">
          {conn_str}
        </DataArray>
        <DataArray type="Int32" Name="offsets" format="ascii">
          {offs_str}
        </DataArray>
      </Lines>
    </Piece>
  </PolyData>
</VTKFile>'''
        with open(filename, "w", encoding="utf-8") as f:
            f.write(vtp_content)

    def _export_vtp_surface(self, filename: str, num_samples_per_curve: int = 25, precomputed_mesh=None):
        """
        Triangulate the planar domain with holes into explicit VTK_TRIANGLE cells.
        Outputs 2D polygon triangles in VTK PolyData with strict boundary and hole conformance.
        """
        if precomputed_mesh is not None:
            points, triangles = precomputed_mesh
        else:
            points, triangles = triangulate_polygon_with_holes(self.outer_loop, self.hole_loops, num_samples_per_curve)

        num_points = len(points)
        num_triangles = len(triangles)

        pts_str = " ".join(f"{p.x:.6g} {p.y:.6g} 0.0" for p in points)
        conn_str = " ".join(f"{t[0]} {t[1]} {t[2]}" for t in triangles)
        offs_str = " ".join(str((i + 1) * 3) for i in range(num_triangles))

        vtp_content = f'''<?xml version="1.0"?>
<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">
  <PolyData>
    <Piece NumberOfPoints="{num_points}" NumberOfVerts="0" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="{num_triangles}">
      <Points>
        <DataArray type="Float32" NumberOfComponents="3" format="ascii">
          {pts_str}
        </DataArray>
      </Points>
      <Polys>
        <DataArray type="Int32" Name="connectivity" format="ascii">
          {conn_str}
        </DataArray>
        <DataArray type="Int32" Name="offsets" format="ascii">
          {offs_str}
        </DataArray>
      </Polys>
    </Piece>
  </PolyData>
</VTKFile>'''
        with open(filename, "w", encoding="utf-8") as f:
            f.write(vtp_content)

    def to_vtu(self, filename: str, num_samples_per_curve: int = 25, precomputed_mesh=None):
        """
        Export domain as native VTK UnstructuredGrid (.vtu) with explicit VTK_TRIANGLE (cell type 5) cells.
        """
        if precomputed_mesh is not None:
            points, triangles = precomputed_mesh
        else:
            points, triangles = triangulate_polygon_with_holes(self.outer_loop, self.hole_loops, num_samples_per_curve)

        num_points = len(points)
        num_triangles = len(triangles)

        pts_str = " ".join(f"{p.x:.6g} {p.y:.6g} 0.0" for p in points)
        conn_str = " ".join(f"{t[0]} {t[1]} {t[2]}" for t in triangles)
        offs_str = " ".join(str((i + 1) * 3) for i in range(num_triangles))
        types_str = " ".join("5" for _ in range(num_triangles))

        dirname = os.path.dirname(filename)
        if dirname and not os.path.exists(dirname):
            os.makedirs(dirname, exist_ok=True)

        vtu_content = f'''<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
  <UnstructuredGrid>
    <Piece NumberOfPoints="{num_points}" NumberOfCells="{num_triangles}">
      <Points>
        <DataArray type="Float32" NumberOfComponents="3" format="ascii">
          {pts_str}
        </DataArray>
      </Points>
      <Cells>
        <DataArray type="Int32" Name="connectivity" format="ascii">
          {conn_str}
        </DataArray>
        <DataArray type="Int32" Name="offsets" format="ascii">
          {offs_str}
        </DataArray>
        <DataArray type="UInt8" Name="types" format="ascii">
          {types_str}
        </DataArray>
      </Cells>
    </Piece>
  </UnstructuredGrid>
</VTKFile>'''
        with open(filename, "w", encoding="utf-8") as f:
            f.write(vtu_content)

    # -------------------------------------------------------------------------
    # Mathematica (.m) Export
    # -------------------------------------------------------------------------
    def to_mathematica(self, filepath: str):
        """
        Export domain to a Mathematica script (.m) for interactive visualization.
        """
        dirname = os.path.dirname(filepath)
        if dirname and not os.path.exists(dirname):
            os.makedirs(dirname, exist_ok=True)

        lines = [
            "(* ========================================================================= *)",
            "(* G+Smo Planar Domain Definition & Mathematica Visualization Script        *)",
            f"(* Domain: corners={self.num_corners}, convex={self.is_convex}, star={self.is_star_shaped}, holes={self.num_holes} *)",
            "(* ========================================================================= *)",
            "",
            "ClearAll[gismoDomain, plotGismoDomain];",
            "",
            "(* 1. Domain Specification *)",
            "gismoDomain = <|",
            f'  "NumCorners" -> {self.num_corners},',
            f'  "IsConvex" -> {str(self.is_convex).lower()},',
            f'  "IsStarShaped" -> {str(self.is_star_shaped).lower()},',
            f'  "NumHoles" -> {self.num_holes},',
        ]

        # Corners
        corners = self.get_corners()
        c_strs = [f"{{{c.x:.8g}, {c.y:.8g}}}" for c in corners]
        lines.append(f'  "Corners" -> {{{", ".join(c_strs)}}},')

        # Outer loop curves
        lines.append('  "OuterLoop" -> {')
        for i, curve in enumerate(self.outer_loop.curves):
            cps_str = ", ".join(f"{{{cp.x:.8g}, {cp.y:.8g}}}" for cp in curve.control_points)
            kv_str = ", ".join(f"{k:.8g}" for k in curve.knot_vector)
            comma = "," if i < len(self.outer_loop.curves) - 1 else ""
            lines.append(f'    <|"Degree" -> {curve.degree}, "ControlPoints" -> {{{cps_str}}}, "Knots" -> {{{kv_str}}}|>{comma}')
        lines.append('  },')

        # Hole loops
        lines.append('  "Holes" -> {')
        for h_idx, hole in enumerate(self.hole_loops):
            lines.append('    {')
            for i, curve in enumerate(hole.curves):
                cps_str = ", ".join(f"{{{cp.x:.8g}, {cp.y:.8g}}}" for cp in curve.control_points)
                kv_str = ", ".join(f"{k:.8g}" for k in curve.knot_vector)
                comma = "," if i < len(hole.curves) - 1 else ""
                lines.append(f'      <|"Degree" -> {curve.degree}, "ControlPoints" -> {{{cps_str}}}, "Knots" -> {{{kv_str}}}|>{comma}')
            h_comma = "," if h_idx < len(self.hole_loops) - 1 else ""
            lines.append(f'    }}{h_comma}')
        lines.append('  }')
        lines.append('|>;')
        lines.append('')

        # Mathematica Visualization Function
        lines.extend([
            "(* 2. Visualization Function *)",
            "bernsteinPoly[p_, j_, t_] := Binomial[p, j] * t^j * (1 - t)^(p - j);",
            "evalBSpline[curve_, t_] := Module[{p = curve[\"Degree\"], cps = curve[\"ControlPoints\"]},",
            "  Sum[cps[[j + 1]] * bernsteinPoly[p, j, t], {j, 0, p}]",
            "];",
            "",
            "plotGismoDomain[domain_] := Module[{",
            "  outerCurves, holeCurves, cornerPts, cnetLines, curvePlots, plotRange, bbox",
            "},",
            "  outerCurves = domain[\"OuterLoop\"];",
            "  holeCurves = Flatten[domain[\"Holes\"]];",
            "  cornerPts = domain[\"Corners\"];",
            "",
            "  (* Control polygons *)",
            "  cnetLines = Join[",
            "    Line /@ (outerCurves[[All, \"ControlPoints\"]]),",
            "    Line /@ (holeCurves[[All, \"ControlPoints\"]])",
            "  ];",
            "",
            "  (* Parametric boundary curves *)",
            "  curvePlots = Join[",
            "    Table[ParametricPlot[evalBSpline[c, t], {t, 0, 1}, PlotStyle -> Directive[Thick, Blue]], {c, outerCurves}],",
            "    Table[ParametricPlot[evalBSpline[c, t], {t, 0, 1}, PlotStyle -> Directive[Thick, Red]], {c, holeCurves}]",
            "  ];",
            "",
            "  Show[",
            "    Graphics[{",
            "      {Directive[LightGray, Opacity[0.3]], EdgeForm[None]},",
            "      {Directive[Dashed, Gray, Opacity[0.6]], cnetLines},",
            "      {Directive[PointSize[0.015], Darker[Green]], Point /@ (outerCurves[[All, \"ControlPoints\"]][[All, 1]])},",
            "      {Directive[PointSize[0.025], Red], Point[cornerPts]},",
            "      MapIndexed[Text[Style[Row[{\"C\", #2[[1]] - 1}], 12, Bold, Darker[Red]], #1 + {0.05, 0.05}] &, cornerPts]",
            "    }],",
            "    curvePlots,",
            "    Axes -> True, GridLines -> Automatic, AspectRatio -> Automatic,",
            "    PlotLabel -> Style[StringTemplate[\"Planar Domain: `` Corners, Convex: ``, Holes: ``\"][",
            "      domain[\"NumCorners\"], domain[\"IsConvex\"], domain[\"NumHoles\"]], 14, Bold]",
            "  ]",
            "];",
            "",
            "(* 3. Render Plot *)",
            "domainVisualization = plotGismoDomain[gismoDomain];",
            "Print[\"Domain loaded successfully. Displaying plot:\"];",
            "domainVisualization"
        ])

        with open(filepath, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")

    # -------------------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------------------
    def summary(self) -> str:
        bbox = self.bounding_box()
        corners = self.get_corners()
        corner_strs = [f"C{i}:{c}" for i, c in enumerate(corners)]
        res = [
            "=" * 60,
            "G+SMO PLANAR DOMAIN SUMMARY",
            "=" * 60,
            f"Corners (C^0 reduction points): {self.num_corners}",
            f"Convex:                        {self.is_convex}",
            f"Star-shaped:                   {self.is_star_shaped}",
            f"Holes:                         {self.num_holes}",
            f"Outer loop curves:             {self.outer_loop.num_curves()}",
            f"Outer loop orientation:        {'CCW (valid)' if self.outer_loop.is_ccw() else 'CW (INVALID)'}",
            f"Bounding Box:                  [{bbox[0]:.4f}, {bbox[1]:.4f}] x [{bbox[2]:.4f}, {bbox[3]:.4f}]",
            f"Corners:                       {', '.join(corner_strs) if corner_strs else 'None (smooth closed boundary)'}",
            "=" * 60
        ]
        return "\n".join(res)


# =============================================================================
# Planar Domain Triangulation (Ear-Clipping with Holes)
# =============================================================================

# =============================================================================
# Robust Constrained Planar Domain Triangulation (Earcut with Bridge-Cuts)
# =============================================================================

class _EarcutNode:
    __slots__ = ('i', 'x', 'y', 'prev', 'next', 'z', 'prevZ', 'nextZ', 'steiner')
    def __init__(self, i: int, x: float, y: float):
        self.i = i
        self.x = float(x)
        self.y = float(y)
        self.prev = None
        self.next = None
        self.z = None
        self.prevZ = None
        self.nextZ = None
        self.steiner = False


def _earcut_signed_area(coords: List[Tuple[float, float]]) -> float:
    acc = 0.0
    n = len(coords)
    for i in range(n):
        p1 = coords[i]
        p2 = coords[(i + 1) % n]
        acc += p1[0] * p2[1] - p2[0] * p1[1]
    return 0.5 * acc


def _earcut_linked_list(coords: List[Tuple[float, float]], point_offset: int) -> Optional[_EarcutNode]:
    last = None
    n = len(coords)
    if n < 3:
        return None
    for i in range(n):
        node = _EarcutNode(point_offset + i, coords[i][0], coords[i][1])
        if last is None:
            node.prev = node
            node.next = node
        else:
            node.next = last.next
            node.prev = last
            last.next.prev = node
            last.next = node
        last = node
    return _earcut_filter_points(last)


def _earcut_filter_points(start: Optional[_EarcutNode], end: Optional[_EarcutNode] = None) -> Optional[_EarcutNode]:
    if not start:
        return start
    if end is None:
        end = start

    p = start
    again = True
    while again or p != end:
        again = False
        if not p.steiner and (_earcut_equals(p, p.next) or _earcut_area(p.prev, p, p.next) == 0):
            _earcut_remove_node(p)
            p = end = p.prev
            if p == p.next:
                break
            again = True
        else:
            p = p.next
    return end


def _earcut_remove_node(p: _EarcutNode):
    p.next.prev = p.prev
    p.prev.next = p.next
    if p.prevZ:
        p.prevZ.nextZ = p.nextZ
    if p.nextZ:
        p.nextZ.prevZ = p.prevZ


def _earcut_equals(p1: _EarcutNode, p2: _EarcutNode) -> bool:
    return abs(p1.x - p2.x) < 1e-11 and abs(p1.y - p2.y) < 1e-11


def _earcut_area(p: _EarcutNode, q: _EarcutNode, r: _EarcutNode) -> float:
    return (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y)


def _earcut_split_polygon(a: _EarcutNode, b: _EarcutNode) -> _EarcutNode:
    a2 = _EarcutNode(a.i, a.x, a.y)
    b2 = _EarcutNode(b.i, b.x, b.y)
    an = a.next
    bp = b.prev

    a.next = b
    b.prev = a

    a2.next = an
    an.prev = a2

    b2.next = a2
    a2.prev = b2

    bp.next = b2
    b2.prev = bp

    return b2


def _earcut_point_in_triangle(ax: float, ay: float, bx: float, by: float, cx: float, cy: float, px: float, py: float) -> bool:
    return (cx - px) * (ay - py) - (ax - px) * (cy - py) >= 0 and \
           (ax - px) * (by - py) - (bx - px) * (ay - py) >= 0 and \
           (bx - px) * (cy - py) - (cx - px) * (by - py) >= 0


def _earcut_is_ear(ear: _EarcutNode) -> bool:
    a = ear.prev
    b = ear
    c = ear.next

    if _earcut_area(a, b, c) >= 0:
        return False

    p = ear.next.next
    while p != ear.prev:
        if _earcut_point_in_triangle(a.x, a.y, b.x, b.y, c.x, c.y, p.x, p.y) and _earcut_area(p.prev, p, p.next) >= 0:
            return False
        p = p.next
    return True


def _earcut_locally_inside(a: _EarcutNode, b: _EarcutNode) -> bool:
    return _earcut_area(a.prev, a, a.next) < 0 if _earcut_area(a, b, a.next) >= 0 and _earcut_area(a, a.prev, b) >= 0 else \
           _earcut_area(a, b, a.prev) < 0 or _earcut_area(a, a.next, b) < 0


def _earcut_intersects(p1: _EarcutNode, q1: _EarcutNode, p2: _EarcutNode, q2: _EarcutNode) -> bool:
    def sgn(val):
        return 1 if val > 1e-11 else (-1 if val < -1e-11 else 0)
    def on_seg(p, q, r):
        return q.x <= max(p.x, r.x) and q.x >= min(p.x, r.x) and q.y <= max(p.y, r.y) and q.y >= min(p.y, r.y)

    o1 = sgn(_earcut_area(p1, q1, p2))
    o2 = sgn(_earcut_area(p1, q1, q2))
    o3 = sgn(_earcut_area(p2, q2, p1))
    o4 = sgn(_earcut_area(p2, q2, q1))
    if o1 != o2 and o3 != o4:
        return True
    if o1 == 0 and on_seg(p1, p2, q1): return True
    if o2 == 0 and on_seg(p1, q2, q1): return True
    if o3 == 0 and on_seg(p2, p1, q2): return True
    if o4 == 0 and on_seg(p2, q1, q2): return True
    return False


def _earcut_cure_local_intersections(start: _EarcutNode, triangles: List[Tuple[int, int, int]]) -> _EarcutNode:
    p = start
    while True:
        a = p.prev
        b = p.next.next
        if not _earcut_equals(a, b) and _earcut_intersects(a, p, p.next, b) and _earcut_locally_inside(a, b) and _earcut_locally_inside(b, a):
            triangles.append((a.i, p.i, b.i))
            _earcut_remove_node(p)
            _earcut_remove_node(p.next)
            p = start = b
        p = p.next
        if p == start:
            break
    return _earcut_filter_points(p)


def _earcut_linked(ear: Optional[_EarcutNode], triangles: List[Tuple[int, int, int]], pass_num: int = 0):
    if not ear:
        return

    stop = ear
    curr = ear
    while curr.prev != curr.next:
        prev_n = curr.prev
        next_n = curr.next

        if _earcut_is_ear(curr):
            triangles.append((prev_n.i, curr.i, next_n.i))
            _earcut_remove_node(curr)
            curr = next_n.next
            stop = next_n.next
            continue

        curr = next_n
        if curr == stop:
            if pass_num == 0:
                _earcut_linked(_earcut_filter_points(curr), triangles, 1)
            elif pass_num == 1:
                curr = _earcut_cure_local_intersections(_earcut_filter_points(curr), triangles)
                _earcut_linked(curr, triangles, 2)
            break


def _segments_intersect_coords(p1: Tuple[float, float], p2: Tuple[float, float],
                              p3: Tuple[float, float], p4: Tuple[float, float]) -> bool:
    def ccw(a, b, c):
        return (c[1] - a[1]) * (b[0] - a[0]) > (b[1] - a[1]) * (c[0] - a[0])
    if (abs(p1[0] - p3[0]) < 1e-7 and abs(p1[1] - p3[1]) < 1e-7) or \
       (abs(p1[0] - p4[0]) < 1e-7 and abs(p1[1] - p4[1]) < 1e-7) or \
       (abs(p2[0] - p3[0]) < 1e-7 and abs(p2[1] - p3[1]) < 1e-7) or \
       (abs(p2[0] - p4[0]) < 1e-7 and abs(p2[1] - p4[1]) < 1e-7):
        return False
    return (ccw(p1, p3, p4) != ccw(p2, p3, p4)) and (ccw(p1, p2, p3) != ccw(p1, p2, p4))


def triangulate_polygon_with_holes(outer_input, holes_input, num_samples_per_curve: int = 25) -> Tuple[List[Point2D], List[Tuple[int, int, int]]]:
    """
    Robust constrained triangulation of a 2D planar domain with holes.
    Guarantees:
    - Accurate shortest-distance, mutually visible bridge cuts without edge intersections
    - Linked-list ear clipping with local intersection resolution and degenerate handling
    - Point-In-Domain filtering: centroid must be strictly INSIDE outer boundary AND OUTSIDE all interior hole loops
    - 0 triangles inside holes, 0 triangles outside domain boundary, and 100% domain area coverage.
    """
    # 1. Sample outer loop
    if hasattr(outer_input, 'sample_loop'):
        outer_poly_pts = outer_input.sample_loop(num_samples_per_curve)
        outer_obj = outer_input
    elif hasattr(outer_input, 'sample'):
        outer_poly_pts = outer_input.sample(num_samples_per_curve)
        outer_obj = outer_input
    else:
        outer_poly_pts = [Point2D(p[0], p[1]) for p in outer_input]
        outer_obj = None

    out_coords = [(p.x, p.y) for p in outer_poly_pts]
    if _earcut_signed_area(out_coords) < 0:
        out_coords.reverse()
        outer_poly_pts.reverse()

    # 2. Sample hole loops
    hole_coords_list = []
    hole_objs = []
    if holes_input:
        for h in holes_input:
            if hasattr(h, 'sample_loop'):
                h_poly_pts = h.sample_loop(num_samples_per_curve)
                hole_objs.append(h)
            elif hasattr(h, 'sample'):
                h_poly_pts = h.sample(num_samples_per_curve)
                hole_objs.append(h)
            else:
                h_poly_pts = [Point2D(p[0], p[1]) for p in h]
                hole_objs.append(None)
            h_coords = [(p.x, p.y) for p in h_poly_pts]
            if _earcut_signed_area(h_coords) > 0:
                h_coords.reverse()
                h_poly_pts.reverse()
            hole_coords_list.append(h_coords)

    # All initial boundary edges for visibility/intersection testing
    edges = []
    n_out = len(out_coords)
    for i in range(n_out):
        edges.append((out_coords[i], out_coords[(i + 1) % n_out]))
    for h_coords in hole_coords_list:
        nh = len(h_coords)
        for i in range(nh):
            edges.append((h_coords[i], h_coords[(i + 1) % nh]))

    all_pts = []
    all_pts.extend(out_coords)
    outer_node = _earcut_linked_list(out_coords, 0)

    hole_nodes = []
    for h_coords in hole_coords_list:
        offset = len(all_pts)
        all_pts.extend(h_coords)
        h_node = _earcut_linked_list(h_coords, offset)
        if h_node:
            hole_nodes.append(h_node)

    # 3. Connect each hole via shortest mutually visible bridge
    for h_node in hole_nodes:
        # Fast candidate sorting: collect nearest outer vertices for each hole vertex
        candidates = []
        curr_h = h_node
        while True:
            curr_p = outer_node
            d_list = []
            while True:
                d_sq = (curr_h.x - curr_p.x) ** 2 + (curr_h.y - curr_p.y) ** 2
                d_list.append((d_sq, curr_h, curr_p))
                curr_p = curr_p.next
                if curr_p == outer_node:
                    break
            d_list.sort(key=lambda item: item[0])
            candidates.extend(d_list[:3])
            curr_h = curr_h.next
            if curr_h == h_node:
                break

        candidates.sort(key=lambda item: item[0])

        best_h = None
        best_p = None
        for d_sq, cand_h, cand_p in candidates:
            p1 = (cand_h.x, cand_h.y)
            p2 = (cand_p.x, cand_p.y)
            hit = False
            for e in edges:
                if _segments_intersect_coords(p1, p2, e[0], e[1]):
                    hit = True
                    break
            if not hit:
                mid = Point2D((p1[0] + p2[0]) * 0.5, (p1[1] + p2[1]) * 0.5)
                in_domain = True
                if outer_obj and hasattr(outer_obj, 'contains_point'):
                    if not outer_obj.contains_point(mid):
                        in_domain = False
                if in_domain and hole_objs:
                    for hl in hole_objs:
                        if hl and hasattr(hl, 'contains_point') and hl.contains_point(mid):
                            in_domain = False
                            break
                if in_domain:
                    best_h = cand_h
                    best_p = cand_p
                    break

        if best_p is None:
            # Fallback to exhaustive search if top candidates were somehow blocked
            best_dist = float('inf')
            curr_h = h_node
            while True:
                curr_p = outer_node
                while True:
                    d_sq = (curr_h.x - curr_p.x) ** 2 + (curr_h.y - curr_p.y) ** 2
                    if d_sq < best_dist:
                        p1 = (curr_h.x, curr_h.y)
                        p2 = (curr_p.x, curr_p.y)
                        hit = False
                        for e in edges:
                            if _segments_intersect_coords(p1, p2, e[0], e[1]):
                                hit = True
                                break
                        if not hit:
                            mid = Point2D((p1[0] + p2[0]) * 0.5, (p1[1] + p2[1]) * 0.5)
                            in_domain = True
                            if outer_obj and hasattr(outer_obj, 'contains_point'):
                                if not outer_obj.contains_point(mid):
                                    in_domain = False
                            if in_domain and hole_objs:
                                for hl in hole_objs:
                                    if hl and hasattr(hl, 'contains_point') and hl.contains_point(mid):
                                        in_domain = False
                                        break
                            if in_domain:
                                best_dist = d_sq
                                best_h = curr_h
                                best_p = curr_p
                    curr_p = curr_p.next
                    if curr_p == outer_node:
                        break
                curr_h = curr_h.next
                if curr_h == h_node:
                    break

        if best_p is not None and best_h is not None:
            outer_node = _earcut_filter_points(_earcut_split_polygon(best_p, best_h))
            edges.append(((best_p.x, best_p.y), (best_h.x, best_h.y)))

    # 4. Triangulate merged boundary
    triangles = []
    _earcut_linked(outer_node, triangles)

    # 5. Point-In-Domain filtering
    valid_tris = []
    for t in triangles:
        p0, p1, p2 = all_pts[t[0]], all_pts[t[1]], all_pts[t[2]]
        c = Point2D((p0[0] + p1[0] + p2[0]) / 3.0, (p0[1] + p1[1] + p2[1]) / 3.0)

        # Check centroid inside outer boundary
        if outer_obj and hasattr(outer_obj, 'contains_point'):
            if not outer_obj.contains_point(c):
                continue

        # Check centroid outside all hole boundaries
        if hole_objs:
            in_hole = False
            for hl in hole_objs:
                if hl and hasattr(hl, 'contains_point') and hl.contains_point(c):
                    in_hole = True
                    break
            if in_hole:
                continue

        # Discard degenerate zero-area triangles
        tri_a = 0.5 * abs((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p1[1] - p0[1]) * (p2[0] - p0[0]))
        if tri_a < 1e-12:
            continue

        valid_tris.append(t)

    res_pts = [Point2D(p[0], p[1]) for p in all_pts]
    return res_pts, valid_tris


# =============================================================================
# Planar Domain Generator
# =============================================================================

class PlanarDomainGenerator:
    """
    Parametric Boundary & Planar Domain Generator.
    Produces B-spline boundary loops based on:
    - num_corners (0, 1, 2, 3, 4, 5, 6, ...)
    - convex (bool)
    - star_shaped (bool)
    - num_holes (int)
    - degree (int)
    - curved (bool/float)
    - randomness / seed
    """

    @classmethod
    def generate(cls,
                 num_corners: int = 4,
                 convex: bool = True,
                 star_shaped: bool = True,
                 num_holes: int = 0,
                 degree: int = 3,
                 curved: bool = True,
                 curvature_amount: float = 0.2,
                 scale: float = 1.0,
                 center: Tuple[float, float] = (0.0, 0.0),
                 randomness: float = 0.0,
                 seed: Optional[int] = None,
                 hole_shape: str = "circle",
                 hole_radius_ratio: float = 0.22) -> PlanarDomain:
        """
        High-level parametric generator.
        """
        if seed is not None:
            random.seed(seed)

        center_pt = Point2D(center[0], center[1])

        # Validate logic: convex domain is necessarily star-shaped
        if convex and not star_shaped:
            star_shaped = True

        # Generate outer loop based on num_corners
        if num_corners == 0:
            outer_loop = cls._generate_0_corner(convex, star_shaped, degree, curved, curvature_amount, scale, randomness)
        elif num_corners == 1:
            outer_loop = cls._generate_1_corner(convex, star_shaped, degree, curved, curvature_amount, scale, randomness)
        elif num_corners == 2:
            outer_loop = cls._generate_2_corner(convex, star_shaped, degree, curved, curvature_amount, scale, randomness)
        else:
            outer_loop = cls._generate_n_corner(num_corners, convex, star_shaped, degree, curved, curvature_amount, scale, randomness)

        # Translate to desired center
        if center_pt.norm() > 1e-9:
            for c in outer_loop.curves:
                for i in range(len(c.control_points)):
                    c.control_points[i] = c.control_points[i] + center_pt

        outer_loop.ensure_orientation(should_be_ccw=True)

        # Generate holes if requested
        hole_loops = []
        if num_holes > 0:
            hole_loops = cls._generate_holes(outer_loop, num_holes, hole_shape, hole_radius_ratio, degree, scale, center_pt, num_corners=num_corners, star_shaped=star_shaped, convex=convex)

        metadata = {
            "num_corners": num_corners,
            "convex": convex,
            "star_shaped": star_shaped,
            "num_holes": num_holes,
            "degree": degree,
            "curved": curved,
            "curvature_amount": curvature_amount,
            "scale": scale,
            "randomness": randomness,
            "seed": seed
        }

        return PlanarDomain(outer_loop=outer_loop, hole_loops=hole_loops, metadata=metadata)

    # -------------------------------------------------------------------------
    # 0 Corners: Smooth Boundary (C^1 continuous everywhere)
    # -------------------------------------------------------------------------
    @classmethod
    def _generate_0_corner(cls, convex: bool, star_shaped: bool, degree: int,
                           curved: bool, curvature: float, scale: float,
                           randomness: float) -> CurveLoop:
        """
        Smooth closed boundary with 0 corners (C^1 everywhere).
        Divided into 4 smooth cubic B-spline quadrants meeting with matched tangents.
        """
        if convex:
            # Ellipse / Circle
            a = scale * (1.0 + random.uniform(-0.2, 0.2) * randomness)
            b = scale * 0.8 * (1.0 + random.uniform(-0.2, 0.2) * randomness)
            # Cubic Bezier quarter-arc circle constant
            k = 4.0 / 3.0 * (math.sqrt(2.0) - 1.0)
            quads = [
                # Q0: (a, 0) -> (0, b)
                [Point2D(a, 0.0), Point2D(a, b * k), Point2D(a * k, b), Point2D(0.0, b)],
                # Q1: (0, b) -> (-a, 0)
                [Point2D(0.0, b), Point2D(-a * k, b), Point2D(-a, b * k), Point2D(-a, 0.0)],
                # Q2: (-a, 0) -> (0, -b)
                [Point2D(-a, 0.0), Point2D(-a, -b * k), Point2D(-a * k, -b), Point2D(0.0, -b)],
                # Q3: (0, -b) -> (a, 0)
                [Point2D(0.0, -b), Point2D(a * k, -b), Point2D(a, -b * k), Point2D(a, 0.0)]
            ]
        elif star_shaped:
            # Non-convex, star-shaped (Smooth bean / kidney shape)
            # Radial function r(theta) = R0 * (1 - indent * cos(theta))
            r0 = scale
            delta = 0.35 + 0.15 * curvature
            # Sample 4 quadrants from smooth analytic kidney curve and fit C^1 cubic arcs
            quads = cls._fit_radial_closed_curve(
                lambda th: r0 * (1.0 - delta * math.sin(th) + 0.1 * math.cos(2.0 * th)),
                num_segments=4
            )
        else:
            # Non-convex, NOT star-shaped (Smooth horseshoe / C-shape)
            quads = cls._generate_smooth_horseshoe(scale)

        curves = [BSplineCurve(cps, degree=degree) for cps in quads]
        return CurveLoop(curves, is_outer=True)

    # -------------------------------------------------------------------------
    # 1 Corner: Teardrop / Cardioid (C^0 only at apex)
    # -------------------------------------------------------------------------
    @classmethod
    def _generate_1_corner(cls, convex: bool, star_shaped: bool, degree: int,
                           curved: bool, curvature: float, scale: float,
                           randomness: float) -> CurveLoop:
        """
        Boundary with exactly 1 corner (point where continuity is reduced to C^0).
        Sharp apex at (+L, 0), smooth rounded back at (-R, 0).
        Consists of 2 curves meeting with C^1 continuity at the back and C^0 at the apex.
        """
        L = scale * 1.5
        R = scale * 0.9
        H = scale * 0.85

        if convex:
            # Top half from apex (L, 0) to back (-R, 0)
            # Control points engineered so tangent at (-R, 0) is strictly vertical (0, -1)
            cps_top = [
                Point2D(L, 0.0),
                Point2D(L * 0.5, H * 0.8),
                Point2D(-R, H * 0.5),
                Point2D(-R, 0.0)
            ]
            # Bottom half from back (-R, 0) to apex (L, 0)
            # Tangent at (-R, 0) is strictly vertical (0, -1) -> exact C^1 continuity at back!
            cps_bot = [
                Point2D(-R, 0.0),
                Point2D(-R, -H * 0.5),
                Point2D(L * 0.5, -H * 0.8),
                Point2D(L, 0.0)
            ]
        else:
            # Non-convex 1-corner domain: waist indentation (negative curvature)
            # Tangent at (-R, 0) remains strictly vertical (0, -1) -> C^1 at back, C^0 only at apex!
            w = 0.35 * (1.0 + curvature)
            cps_top = [
                Point2D(L, 0.0),
                Point2D(L * 0.5, H * (0.8 - w)),
                Point2D(-R, H * 0.5),
                Point2D(-R, 0.0)
            ]
            cps_bot = [
                Point2D(-R, 0.0),
                Point2D(-R, -H * 0.5),
                Point2D(L * 0.5, -H * (0.8 - w)),
                Point2D(L, 0.0)
            ]

        curves = [BSplineCurve(cps_top, degree=degree), BSplineCurve(cps_bot, degree=degree)]
        return CurveLoop(curves, is_outer=True)

    # -------------------------------------------------------------------------
    # 2 Corners: Lens / Crescent (C^0 at both tips)
    # -------------------------------------------------------------------------
    @classmethod
    def _generate_2_corner(cls, convex: bool, star_shaped: bool, degree: int,
                           curved: bool, curvature: float, scale: float,
                           randomness: float) -> CurveLoop:
        """
        Boundary with exactly 2 corners (e.g. lens, crescent).
        Corners at (-L, 0) and (+L, 0).
        """
        L = scale * 1.2
        H_top = scale * 0.8

        if convex:
            # Symmetric lens / vesica piscis
            H_bot = -H_top
            cps_top = [
                Point2D(-L, 0.0),
                Point2D(-L * 0.4, H_top),
                Point2D(L * 0.4, H_top),
                Point2D(L, 0.0)
            ]
            cps_bot = [
                Point2D(L, 0.0),
                Point2D(L * 0.4, H_bot),
                Point2D(-L * 0.4, H_bot),
                Point2D(-L, 0.0)
            ]
        else:
            # Crescent moon (inner curve arches in the same direction)
            if star_shaped:
                H_inner = H_top * 0.45  # Stays above origin
            else:
                H_inner = H_top * 0.85  # Deep crescent wrapping around

            cps_top = [
                Point2D(-L, 0.0),
                Point2D(-L * 0.5, H_top * 1.3),
                Point2D(L * 0.5, H_top * 1.3),
                Point2D(L, 0.0)
            ]
            cps_bot = [
                Point2D(L, 0.0),
                Point2D(L * 0.4, H_inner),
                Point2D(-L * 0.4, H_inner),
                Point2D(-L, 0.0)
            ]

        curves = [BSplineCurve(cps_top, degree=degree), BSplineCurve(cps_bot, degree=degree)]
        return CurveLoop(curves, is_outer=True)

    # -------------------------------------------------------------------------
    # N Corners (K >= 3): General Polygon / Curved Boundary
    # -------------------------------------------------------------------------
    @classmethod
    def _generate_n_corner(cls, num_corners: int, convex: bool, star_shaped: bool,
                           degree: int, curved: bool, curvature: float,
                           scale: float, randomness: float) -> CurveLoop:
        """
        General domain with K corners (K >= 3).
        """
        K = num_corners

        if not star_shaped:
            # Non-star-shaped geometry (L-shape, U-shape, etc.)
            return cls._generate_non_star_domain(K, degree, curved, curvature, scale, randomness)

        # Star-shaped domain with K corners
        # Compute angles theta_0, ..., theta_{K-1}
        angles = []
        d_theta = 2.0 * math.pi / float(K)
        for i in range(K):
            base_angle = i * d_theta
            if randomness > 1e-4:
                base_angle += random.uniform(-0.35, 0.35) * d_theta * randomness
            angles.append(base_angle)
        angles.sort()

        # Compute radii for each corner
        corners = []
        if convex:
            # Strictly convex corner distribution
            for th in angles:
                r = scale * (1.0 + random.uniform(-0.15, 0.15) * randomness)
                corners.append(Point2D(r * math.cos(th), r * math.sin(th)))
            # Ensure convexity of corner polygon
            corners = cls._ensure_convex_polygon(corners)
        else:
            # Non-convex star-shaped: alternating radii (star polygon)
            for i, th in enumerate(angles):
                if i % 2 == 1:
                    r = scale * (0.45 - 0.2 * curvature)  # Valley / reflex corner
                else:
                    r = scale * (1.1 + 0.2 * curvature)   # Outer tip
                if randomness > 1e-4:
                    r *= (1.0 + random.uniform(-0.15, 0.15) * randomness)
                corners.append(Point2D(r * math.cos(th), r * math.sin(th)))

        # Build B-spline curves connecting consecutive corners
        curves = []
        n_c = len(corners)
        for i in range(n_c):
            p_start = corners[i]
            p_end = corners[(i + 1) % n_c]
            edge_vec = p_end - p_start
            edge_len = edge_vec.norm()

            if not curved:
                # Straight chord (polygonal edge)
                cps = [
                    p_start,
                    p_start + edge_vec * (1.0 / 3.0),
                    p_start + edge_vec * (2.0 / 3.0),
                    p_end
                ]
            else:
                # Curved edge with controlled bulge
                # Normal vector pointing outward (CW 90 deg from edge vector in CCW polygon)
                norm_vec = edge_vec.normal()

                # Curvature bulge height
                if convex:
                    h = edge_len * curvature * 0.25  # Positive outward bulge preserves convexity
                else:
                    # Alternating outward / inward or slight outward
                    h = edge_len * curvature * (0.2 if i % 2 == 0 else -0.15)

                if randomness > 1e-4:
                    h *= (1.0 + random.uniform(-0.2, 0.2) * randomness)

                cps = [
                    p_start,
                    p_start + edge_vec * (1.0 / 3.0) + norm_vec * h,
                    p_start + edge_vec * (2.0 / 3.0) + norm_vec * h,
                    p_end
                ]

            curves.append(BSplineCurve(cps, degree=degree))

        return CurveLoop(curves, is_outer=True)

    # -------------------------------------------------------------------------
    # Non-Star-Shaped Domains (L-Shape, Horseshoe, Serpentine)
    # -------------------------------------------------------------------------
    @classmethod
    def _generate_non_star_domain(cls, num_corners: int, degree: int,
                                  curved: bool, curvature: float, scale: float,
                                  randomness: float) -> CurveLoop:
        """
        Generate non-star-shaped domains (L-shape, U-shape, etc.).
        """
        if num_corners == 6:
            # Canonical L-shaped domain with 6 corners:
            # (0,0) -> (2,0) -> (2,1) -> (1,1) -> (1,2) -> (0,2)
            pts = [
                Point2D(-scale, -scale),
                Point2D(scale, -scale),
                Point2D(scale, 0.0),
                Point2D(0.0, 0.0),          # Reflex corner (270 deg)
                Point2D(0.0, scale),
                Point2D(-scale, scale)
            ]
        elif num_corners == 8:
            # U-shaped / Horseshoe channel with 8 corners:
            w = scale * 0.4
            pts = [
                Point2D(-scale, -scale),
                Point2D(scale, -scale),
                Point2D(scale, scale),
                Point2D(scale - w, scale),
                Point2D(scale - w, -scale + w),   # Inner corner 1
                Point2D(-scale + w, -scale + w),  # Inner corner 2
                Point2D(-scale + w, scale),
                Point2D(-scale, scale)
            ]
        elif num_corners == 12:
            # Canonical 12-corner Swiss Cross (4 arms, 8 reflex 270 deg corners)
            w = scale * 0.38
            L = scale
            pts = [
                Point2D(-w, -L), Point2D(w, -L), Point2D(w, -w),
                Point2D(L, -w), Point2D(L, w), Point2D(w, w),
                Point2D(w, L), Point2D(-w, L), Point2D(-w, w),
                Point2D(-L, w), Point2D(-L, -w), Point2D(-w, -w)
            ]
        elif num_corners == 16:
            # 16-corner 4-finger comb domain (8 convex tips + 8 reflex corners)
            y_bot = -scale
            y_slot = -scale * 0.4
            y_top = scale
            pts = [
                Point2D(-scale, y_bot),
                Point2D(scale, y_bot),
                Point2D(scale, y_top),
                Point2D(scale * 0.6, y_top),
                Point2D(scale * 0.6, y_slot),      # Reflex 1
                Point2D(scale * 0.45, y_slot),     # Reflex 2
                Point2D(scale * 0.45, y_top),
                Point2D(scale * 0.05, y_top),
                Point2D(scale * 0.05, y_slot),     # Reflex 3
                Point2D(-scale * 0.05, y_slot),    # Reflex 4
                Point2D(-scale * 0.05, y_top),
                Point2D(-scale * 0.45, y_top),
                Point2D(-scale * 0.45, y_slot),    # Reflex 5
                Point2D(-scale * 0.6, y_slot),     # Reflex 6
                Point2D(-scale * 0.6, y_top),
                Point2D(-scale, y_top)
            ]
        else:
            # Generalized stepped / serpentine non-star polygon with K corners
            pts = cls._generate_serpentine_corners(num_corners, scale)

        curves = []
        n_pts = len(pts)
        for i in range(n_pts):
            p_start = pts[i]
            p_end = pts[(i + 1) % n_pts]
            edge_vec = p_end - p_start
            edge_len = edge_vec.norm()

            if not curved:
                cps = [
                    p_start,
                    p_start + edge_vec * (1.0 / 3.0),
                    p_start + edge_vec * (2.0 / 3.0),
                    p_end
                ]
            else:
                norm_vec = edge_vec.normal()
                h = edge_len * curvature * 0.15
                cps = [
                    p_start,
                    p_start + edge_vec * (1.0 / 3.0) + norm_vec * h,
                    p_start + edge_vec * (2.0 / 3.0) + norm_vec * h,
                    p_end
                ]
            curves.append(BSplineCurve(cps, degree=degree))

        return CurveLoop(curves, is_outer=True)

    @classmethod
    def _generate_serpentine_corners(cls, K: int, scale: float) -> List[Point2D]:
        """Generate non-star-shaped corner polygon for arbitrary K >= 4."""
        # Create an elbow with K-4 intermediate steps along the concave re-entrant cut
        pts = [
            Point2D(-scale, -scale),
            Point2D(scale, -scale),
            Point2D(scale, scale * 0.2)
        ]
        num_inner = K - 4
        for j in range(num_inner):
            frac = float(j + 1) / float(num_inner + 1)
            # Step in and out to produce re-entrant pockets
            if j % 2 == 0:
                pts.append(Point2D(scale * 0.1, scale * (0.2 + frac * 0.6)))
            else:
                pts.append(Point2D(scale * 0.6, scale * (0.2 + frac * 0.6)))
        pts.append(Point2D(-scale * 0.2, scale))
        pts.append(Point2D(-scale, scale))
        return pts[:K]

    # -------------------------------------------------------------------------
    # Hole Generation (CW oriented loops)
    # -------------------------------------------------------------------------
    @classmethod
    def _generate_holes(cls, outer_loop: CurveLoop, num_holes: int,
                        shape: str, radius_ratio: float, degree: int,
                        scale: float, domain_center: Point2D,
                        num_corners: int = 4, star_shaped: bool = True,
                        convex: bool = True) -> List[CurveLoop]:
        """
        Generate H non-overlapping holes positioned strictly inside the domain interior.
        Oriented CLOCKWISE (CW) as required by G+Smo gsPlanarDomain.
        """
        holes = []
        r_base = scale * radius_ratio

        shapes_pool = ["cross", "star", "triangle", "diamond", "hexagon", "square"]
        def get_hole_shape(idx: int) -> str:
            if isinstance(shape, (list, tuple)):
                return shape[idx % len(shape)]
            if str(shape).lower() == "mixed":
                return shapes_pool[idx % len(shapes_pool)]
            return shape

        raw_specs = []

        # Case 1: U-channel domain (8 corners, non-star)
        if num_corners == 8 and not star_shaped:
            r_hole = min(r_base, scale * 0.08 if num_holes >= 4 else scale * 0.11)
            if num_holes == 1:
                centers = [domain_center + Point2D(0.0, -scale * 0.8)]
            elif num_holes == 2:
                centers = [
                    domain_center + Point2D(-scale * 0.8, scale * 0.2),
                    domain_center + Point2D(scale * 0.8, scale * 0.2)
                ]
            elif num_holes == 3:
                centers = [
                    domain_center + Point2D(-scale * 0.8, scale * 0.3),
                    domain_center + Point2D(scale * 0.8, scale * 0.3),
                    domain_center + Point2D(0.0, -scale * 0.8)
                ]
            elif num_holes == 4:
                centers = [
                    domain_center + Point2D(-scale * 0.8, scale * 0.45),
                    domain_center + Point2D(-scale * 0.8, -scale * 0.05),
                    domain_center + Point2D(scale * 0.8, scale * 0.45),
                    domain_center + Point2D(scale * 0.8, -scale * 0.05)
                ]
            else:
                centers = [
                    domain_center + Point2D(-scale * 0.8, scale * 0.45),
                    domain_center + Point2D(-scale * 0.8, -scale * 0.05),
                    domain_center + Point2D(scale * 0.8, scale * 0.45),
                    domain_center + Point2D(scale * 0.8, -scale * 0.05),
                    domain_center + Point2D(0.0, -scale * 0.8)
                ]
            for idx, c_h in enumerate(centers[:num_holes]):
                raw_specs.append((c_h, r_hole, get_hole_shape(idx)))

        # Case 2: L-shaped domain (6 corners, non-star)
        elif num_corners == 6 and not star_shaped:
            r_hole = min(r_base, scale * 0.12 if num_holes >= 3 else scale * 0.16)
            if num_holes == 1:
                centers = [domain_center + Point2D(scale * 0.5, -scale * 0.5)]
            elif num_holes == 2:
                centers = [
                    domain_center + Point2D(scale * 0.5, -scale * 0.5),
                    domain_center + Point2D(-scale * 0.5, scale * 0.5)
                ]
            elif num_holes == 3:
                centers = [
                    domain_center + Point2D(scale * 0.55, -scale * 0.5),
                    domain_center + Point2D(-scale * 0.5, scale * 0.55),
                    domain_center + Point2D(-scale * 0.5, -scale * 0.5)
                ]
            else:
                centers = [
                    domain_center + Point2D(scale * 0.65, -scale * 0.5),
                    domain_center + Point2D(scale * 0.20, -scale * 0.5),
                    domain_center + Point2D(-scale * 0.5, scale * 0.65),
                    domain_center + Point2D(-scale * 0.5, -scale * 0.5)
                ]
            for idx, c_h in enumerate(centers[:num_holes]):
                raw_specs.append((c_h, r_hole, get_hole_shape(idx)))

        # Case 3: Swiss Cross domain (12 corners, non-star)
        elif num_corners == 12 and not star_shaped:
            r_hole = min(r_base, scale * 0.14)
            centers = [
                domain_center + Point2D(0.0, 0.0),
                domain_center + Point2D(0.0, scale * 0.65),
                domain_center + Point2D(0.0, -scale * 0.65),
                domain_center + Point2D(scale * 0.65, 0.0),
                domain_center + Point2D(-scale * 0.65, 0.0)
            ]
            for idx, c_h in enumerate(centers[:num_holes]):
                raw_specs.append((c_h, r_hole, get_hole_shape(idx)))

        # Case 4: 4-finger Comb domain (16 corners, non-star): placed in solid rectangular backbone
        elif num_corners == 16 and not star_shaped:
            r_hole = min(r_base, scale * 0.09)
            centers = [
                domain_center + Point2D(-scale * 0.75, -scale * 0.42),
                domain_center + Point2D(-scale * 0.25, -scale * 0.42),
                domain_center + Point2D(scale * 0.25, -scale * 0.42),
                domain_center + Point2D(scale * 0.75, -scale * 0.42),
                domain_center + Point2D(0.0, -scale * 0.42)
            ]
            for idx, c_h in enumerate(centers[:num_holes]):
                raw_specs.append((c_h, r_hole, get_hole_shape(idx)))

        # Case 4b: Crescent domain (2 corners, non-convex): placed along arched crescent centerline
        elif num_corners == 2 and not convex:
            r_hole = min(r_base, scale * 0.10)
            if num_holes == 1:
                centers = [domain_center + Point2D(0.0, scale * 0.52)]
            elif num_holes == 2:
                centers = [
                    domain_center + Point2D(-scale * 0.35, scale * 0.48),
                    domain_center + Point2D(scale * 0.35, scale * 0.48)
                ]
            elif num_holes == 3:
                centers = [
                    domain_center + Point2D(-scale * 0.45, scale * 0.46),
                    domain_center + Point2D(0.0, scale * 0.52),
                    domain_center + Point2D(scale * 0.45, scale * 0.46)
                ]
            else:
                centers = []
                for idx in range(num_holes):
                    frac = (idx + 0.5) / float(num_holes)
                    x = (-0.6 + frac * 1.2) * scale
                    y = (0.52 - 0.18 * (x / scale) ** 2) * scale
                    centers.append(domain_center + Point2D(x, y))
            for idx, c_h in enumerate(centers[:num_holes]):
                raw_specs.append((c_h, r_hole, get_hole_shape(idx)))

        # Case 5: Non-convex star domains (starburst, star5, etc. with reflex valleys)
        elif not convex and star_shaped and num_corners >= 4:
            r_hole = min(r_base, scale * 0.075)
            ring_r = scale * 0.18
            for idx in range(num_holes):
                th = idx * 2.0 * math.pi / float(num_holes) + math.pi / float(num_holes)
                off = Point2D(ring_r * math.cos(th), ring_r * math.sin(th))
                c_h = domain_center + off
                raw_specs.append((c_h, r_hole, get_hole_shape(idx)))

        # Case 6: General Star / Convex / Polygon / Bean domains
        else:
            outer_sample_pts = outer_loop.sample_loop(60)
            r_bnd_min = min((p - domain_center).norm() for p in outer_sample_pts)

            if num_holes == 1:
                center_h = domain_center
                if not outer_loop.contains_point(center_h):
                    center_h = domain_center + Point2D(scale * 0.3, 0.0)
                r_hole = min(r_base, r_bnd_min * 0.35)
                raw_specs.append((center_h, r_hole, get_hole_shape(0)))

            elif num_holes == 2:
                d_off = r_bnd_min * 0.45
                r_hole = min(r_base, r_bnd_min * 0.22)
                offsets = [Point2D(-d_off, 0.0), Point2D(d_off, 0.0)]
                for idx, off in enumerate(offsets):
                    c_h = domain_center + off
                    if not outer_loop.contains_point(c_h):
                        c_h = domain_center + off * 0.6
                    raw_specs.append((c_h, r_hole, get_hole_shape(idx)))

            elif num_holes == 3:
                ring_r = r_bnd_min * 0.42
                r_hole = min(r_base, r_bnd_min * 0.18)
                for idx in range(3):
                    th = idx * 2.0 * math.pi / 3.0 + math.pi / 6.0
                    off = Point2D(ring_r * math.cos(th), ring_r * math.sin(th))
                    c_h = domain_center + off
                    raw_specs.append((c_h, r_hole, get_hole_shape(idx)))

            elif num_holes == 4:
                d_off = r_bnd_min * 0.38
                r_hole = min(r_base, r_bnd_min * 0.16)
                offsets = [
                    Point2D(-d_off, -d_off),
                    Point2D(d_off, -d_off),
                    Point2D(d_off, d_off),
                    Point2D(-d_off, d_off)
                ]
                for idx, off in enumerate(offsets):
                    c_h = domain_center + off
                    if not outer_loop.contains_point(c_h):
                        c_h = domain_center + off * 0.6
                    raw_specs.append((c_h, r_hole, get_hole_shape(idx)))

            elif num_holes == 5:
                d_off = r_bnd_min * 0.38
                r_hole = min(r_base, r_bnd_min * 0.14)
                offsets = [
                    Point2D(0.0, 0.0),
                    Point2D(-d_off, -d_off),
                    Point2D(d_off, -d_off),
                    Point2D(d_off, d_off),
                    Point2D(-d_off, d_off)
                ]
                for idx, off in enumerate(offsets):
                    c_h = domain_center + off
                    raw_specs.append((c_h, r_hole, get_hole_shape(idx)))

            elif num_holes == 8:
                d_card = r_bnd_min * 0.55
                d_diag = r_bnd_min * 0.38
                r_hole = min(r_base, r_bnd_min * 0.11)
                offsets = [
                    Point2D(0.0, d_card),
                    Point2D(d_card, 0.0),
                    Point2D(0.0, -d_card),
                    Point2D(-d_card, 0.0),
                    Point2D(-d_diag, -d_diag),
                    Point2D(d_diag, -d_diag),
                    Point2D(d_diag, d_diag),
                    Point2D(-d_diag, d_diag)
                ]
                for idx, off in enumerate(offsets):
                    c_h = domain_center + off
                    raw_specs.append((c_h, r_hole, get_hole_shape(idx)))

            else:
                ring_r = r_bnd_min * 0.45
                r_hole = min(r_base, (r_bnd_min * 0.45) / math.sqrt(num_holes))
                for idx in range(num_holes):
                    th = idx * 2.0 * math.pi / float(num_holes)
                    off = Point2D(ring_r * math.cos(th), ring_r * math.sin(th))
                    c_h = domain_center + off
                    raw_specs.append((c_h, r_hole, get_hole_shape(idx)))

        # Active safety verification loop for every hole:
        # Guarantees that every hole is strictly inside the domain with >= 4% scale boundary clearance
        outer_sample = outer_loop.sample_loop(60)
        verified_holes = []
        for idx, (curr_c, curr_r, s_h) in enumerate(raw_specs):
            best_hole = None
            for attempt in range(15):
                test_h = cls._make_single_hole(curr_c, curr_r, s_h, degree)
                h_pts = test_h.sample_loop(20)
                all_inside = all(outer_loop.contains_point(p) for p in h_pts)
                min_bnd_dist = min(min((p - op).norm() for op in outer_sample) for p in h_pts) if all_inside else -1.0
                if all_inside and min_bnd_dist >= scale * 0.04:
                    best_hole = test_h
                    break
                else:
                    curr_c = domain_center + (curr_c - domain_center) * 0.88
                    curr_r *= 0.88
            if best_hole is None:
                best_hole = cls._make_single_hole(curr_c, curr_r, s_h, degree)
            verified_holes.append(best_hole)

        return verified_holes

    @classmethod
    def _make_single_hole(cls, center: Point2D, radius: float, shape: str, degree: int) -> CurveLoop:
        """
        Create a single hole loop oriented CLOCKWISE (CW).
        Supported shapes:
        - circle: smooth 4-quadrant cubic B-spline circle
        - square: 4-corner square hole
        - triangle / tri: 3-corner equilateral triangle hole
        - diamond / rhombus: 4-corner diamond hole
        - cross / plus / swiss_cross: 12-corner Swiss cross hole (8 reflex corners)
        - hexagon / hex: 6-corner regular hexagon hole
        - star: 10-corner 5-pointed star hole (5 reflex corners)
        - l_hole: 6-corner L-shaped hole (1 reflex corner)
        """
        s_lower = str(shape).lower().replace("-", "_").replace(" ", "_")

        if s_lower == "square":
            # 4-corner square hole (CW order: Top-Left -> Top-Right -> Bot-Right -> Bot-Left)
            corners = [
                Point2D(center.x - radius, center.y - radius),
                Point2D(center.x - radius, center.y + radius),
                Point2D(center.x + radius, center.y + radius),
                Point2D(center.x + radius, center.y - radius)
            ]
            curves = []
            for i in range(4):
                p_start = corners[i]
                p_end = corners[(i + 1) % 4]
                e = p_end - p_start
                cps = [p_start + e * (j / float(degree)) for j in range(degree + 1)]
                curves.append(BSplineCurve(cps, degree=degree))
            h_loop = CurveLoop(curves, is_outer=False)
            h_loop.ensure_orientation(should_be_ccw=False)
            return h_loop

        elif s_lower in ("triangle", "tri"):
            # Equilateral triangle in CW order (top -> bottom-right -> bottom-left)
            pts = [
                center + Point2D(0.0, radius),
                center + Point2D(radius * math.cos(-math.pi / 6.0), radius * math.sin(-math.pi / 6.0)),
                center + Point2D(radius * math.cos(-5.0 * math.pi / 6.0), radius * math.sin(-5.0 * math.pi / 6.0))
            ]
            curves = []
            for i in range(3):
                p_start = pts[i]
                p_end = pts[(i + 1) % 3]
                e = p_end - p_start
                cps = [p_start + e * (j / float(degree)) for j in range(degree + 1)]
                curves.append(BSplineCurve(cps, degree=degree))
            h_loop = CurveLoop(curves, is_outer=False)
            h_loop.ensure_orientation(should_be_ccw=False)
            return h_loop

        elif s_lower in ("cross", "plus", "swiss_cross"):
            # 12-corner Swiss cross hole (4 outer tips, 8 reflex re-entrant corners at 270 deg)
            w = radius * 0.35
            L = radius
            raw_pts = [
                Point2D(-w, L), Point2D(w, L), Point2D(w, w),
                Point2D(L, w), Point2D(L, -w), Point2D(w, -w),
                Point2D(w, -L), Point2D(-w, -L), Point2D(-w, -w),
                Point2D(-L, -w), Point2D(-L, w), Point2D(-w, w)
            ]
            curves = []
            for i in range(12):
                p_start = center + raw_pts[i]
                p_end = center + raw_pts[(i + 1) % 12]
                e = p_end - p_start
                cps = [p_start + e * (j / float(degree)) for j in range(degree + 1)]
                curves.append(BSplineCurve(cps, degree=degree))
            h_loop = CurveLoop(curves, is_outer=False)
            h_loop.ensure_orientation(should_be_ccw=False)
            return h_loop

        elif s_lower in ("diamond", "rhombus"):
            # 4-corner diamond hole (CW order: top -> right -> bottom -> left)
            pts = [
                center + Point2D(0.0, radius),
                center + Point2D(radius * 0.75, 0.0),
                center + Point2D(0.0, -radius),
                center + Point2D(-radius * 0.75, 0.0)
            ]
            curves = []
            for i in range(4):
                p_start = pts[i]
                p_end = pts[(i + 1) % 4]
                e = p_end - p_start
                cps = [p_start + e * (j / float(degree)) for j in range(degree + 1)]
                curves.append(BSplineCurve(cps, degree=degree))
            h_loop = CurveLoop(curves, is_outer=False)
            h_loop.ensure_orientation(should_be_ccw=False)
            return h_loop

        elif s_lower in ("hexagon", "hex"):
            # 6-corner regular hexagon in CW order
            pts = [
                center + Point2D(radius * math.cos(math.pi / 2.0 - i * math.pi / 3.0),
                                 radius * math.sin(math.pi / 2.0 - i * math.pi / 3.0))
                for i in range(6)
            ]
            curves = []
            for i in range(6):
                p_start = pts[i]
                p_end = pts[(i + 1) % 6]
                e = p_end - p_start
                cps = [p_start + e * (j / float(degree)) for j in range(degree + 1)]
                curves.append(BSplineCurve(cps, degree=degree))
            h_loop = CurveLoop(curves, is_outer=False)
            h_loop.ensure_orientation(should_be_ccw=False)
            return h_loop

        elif s_lower in ("l_hole", "l_shape"):
            # 6-corner L-shaped hole with reflex corner in CW order
            s = radius * 0.8
            raw_pts = [
                Point2D(-s, s), Point2D(s, s), Point2D(s, 0.0),
                Point2D(0.0, 0.0), Point2D(0.0, -s), Point2D(-s, -s)
            ]
            curves = []
            for i in range(6):
                p_start = center + raw_pts[i]
                p_end = center + raw_pts[(i + 1) % 6]
                e = p_end - p_start
                cps = [p_start + e * (j / float(degree)) for j in range(degree + 1)]
                curves.append(BSplineCurve(cps, degree=degree))
            h_loop = CurveLoop(curves, is_outer=False)
            h_loop.ensure_orientation(should_be_ccw=False)
            return h_loop

        elif s_lower == "star":
            # 5-pointed star hole: 10 vertices in CW order
            K = 5
            pts = []
            for i in range(2 * K):
                th = -i * math.pi / float(K)
                r = radius if i % 2 == 0 else radius * 0.5
                pts.append(center + Point2D(r * math.cos(th), r * math.sin(th)))
            curves = []
            for i in range(len(pts)):
                p_start = pts[i]
                p_end = pts[(i + 1) % len(pts)]
                e = p_end - p_start
                cps = [p_start + e * (j / float(degree)) for j in range(degree + 1)]
                curves.append(BSplineCurve(cps, degree=degree))
            h_loop = CurveLoop(curves, is_outer=False)
            h_loop.ensure_orientation(should_be_ccw=False)
            return h_loop

        else:
            # Circular hole: 4 smooth cubic B-spline quadrants in CW order
            k = radius * 4.0 / 3.0 * (math.sqrt(2.0) - 1.0)
            r = radius
            quads = [
                [center + Point2D(r, 0.0), center + Point2D(r, -k), center + Point2D(k, -r), center + Point2D(0.0, -r)],
                [center + Point2D(0.0, -r), center + Point2D(-k, -r), center + Point2D(-r, -k), center + Point2D(-r, 0.0)],
                [center + Point2D(-r, 0.0), center + Point2D(-r, k), center + Point2D(-k, r), center + Point2D(0.0, r)],
                [center + Point2D(0.0, r), center + Point2D(k, r), center + Point2D(r, k), center + Point2D(r, 0.0)]
            ]
            curves = [BSplineCurve(cps, degree=degree) for cps in quads]
            h_loop = CurveLoop(curves, is_outer=False)
            h_loop.ensure_orientation(should_be_ccw=False)
            return h_loop

    # -------------------------------------------------------------------------
    # Geometric Utilities
    # -------------------------------------------------------------------------
    @classmethod
    def _ensure_convex_polygon(cls, pts: List[Point2D]) -> List[Point2D]:
        """Ensure all interior cross products are strictly positive."""
        n = len(pts)
        for _ in range(5):
            for i in range(n):
                p0 = pts[(i - 1) % n]
                p1 = pts[i]
                p2 = pts[(i + 1) % n]
                cross = (p1 - p0).cross(p2 - p1)
                if cross <= 1e-4:
                    # Move p1 outwards along angle bisector
                    bisector = ((p1 - p0).normalized() - (p2 - p1).normalized()).normal().normalized()
                    pts[i] = pts[i] + bisector * 0.15
        return pts

    @classmethod
    def _fit_radial_closed_curve(cls, r_func, num_segments: int = 4) -> List[List[Point2D]]:
        """Fit C^1 cubic Bezier segments to radial curve r(theta)."""
        quads = []
        d_th = 2.0 * math.pi / float(num_segments)
        for s in range(num_segments):
            th0 = s * d_th
            th1 = (s + 1) * d_th
            # Sample start and end
            r0 = r_func(th0)
            r1 = r_func(th1)
            p0 = Point2D(r0 * math.cos(th0), r0 * math.sin(th0))
            p3 = Point2D(r1 * math.cos(th1), r1 * math.sin(th1))

            # Numerical tangent at endpoints
            eps = 1e-4
            dr0 = (r_func(th0 + eps) - r_func(th0 - eps)) / (2.0 * eps)
            dr1 = (r_func(th1 + eps) - r_func(th1 - eps)) / (2.0 * eps)

            t0 = Point2D(dr0 * math.cos(th0) - r0 * math.sin(th0), dr0 * math.sin(th0) + r0 * math.cos(th0)).normalized()
            t1 = Point2D(dr1 * math.cos(th1) - r1 * math.sin(th1), dr1 * math.sin(th1) + r1 * math.cos(th1)).normalized()

            chord = (p3 - p0).norm()
            p1 = p0 + t0 * (chord / 3.0)
            p2 = p3 - t1 * (chord / 3.0)
            quads.append([p0, p1, p2, p3])
        return quads

    @classmethod
    def _generate_smooth_horseshoe(cls, scale: float) -> List[List[Point2D]]:
        """Smooth non-convex, non-star-shaped horseshoe with 4 C^1 matching arcs."""
        # A curved tubular arc with rounded caps
        r_out = scale * 1.2
        r_in = scale * 0.5
        r_cap = (r_out - r_in) * 0.5
        y_mid = (r_out + r_in) * 0.5

        # Outer arc (half circle): (0, r_out) -> (-r_out, 0) -> (0, -r_out)
        k_out = r_out * 4.0 / 3.0 * (math.sqrt(2.0) - 1.0)
        c_out1 = [Point2D(0.0, r_out), Point2D(-k_out, r_out), Point2D(-r_out, k_out), Point2D(-r_out, 0.0)]
        c_out2 = [Point2D(-r_out, 0.0), Point2D(-r_out, -k_out), Point2D(-k_out, -r_out), Point2D(0.0, -r_out)]

        # Bottom rounded cap
        c_cap_bot = [Point2D(0.0, -r_out), Point2D(r_cap * 0.8, -r_out), Point2D(r_cap * 0.8, -r_in), Point2D(0.0, -r_in)]

        # Inner arc: (0, -r_in) -> (-r_in, 0) -> (0, r_in)
        k_in = r_in * 4.0 / 3.0 * (math.sqrt(2.0) - 1.0)
        c_in2 = [Point2D(0.0, -r_in), Point2D(-k_in, -r_in), Point2D(-r_in, -k_in), Point2D(-r_in, 0.0)]
        c_in1 = [Point2D(-r_in, 0.0), Point2D(-r_in, k_in), Point2D(-k_in, r_in), Point2D(0.0, r_in)]

        # Top rounded cap
        c_cap_top = [Point2D(0.0, r_in), Point2D(r_cap * 0.8, r_in), Point2D(r_cap * 0.8, r_out), Point2D(0.0, r_out)]

        return [c_out1, c_out2, c_cap_bot, c_in2, c_in1, c_cap_top]

    # -------------------------------------------------------------------------
    # Presets Library
    # -------------------------------------------------------------------------
    @classmethod
    def preset(cls, name: str, scale: float = 1.0, degree: int = 3, **kwargs) -> PlanarDomain:
        """
        Create domain from a curated library of standard benchmark presets.
        Presets:
        - circle, ellipse, bean, horseshoe_smooth (0 corners)
        - teardrop, heart (1 corner)
        - lens, crescent (2 corners)
        - triangle, curved_triangle (3 corners)
        - quad, curved_quad, indented_quad, rectangle (4 corners)
        - pentagon (5 corners)
        - star5 (10 corners)
        - lshape (6 corners, non-star)
        - u_channel (8 corners, non-star)
        - plate_with_hole, plate_with_2holes, lshape_with_hole
        """
        n = name.lower().replace("-", "_").replace(" ", "_")

        if n in ("circle", "disk"):
            return cls.generate(num_corners=0, convex=True, star_shaped=True, degree=degree, scale=scale, **kwargs)
        elif n == "ellipse":
            return cls.generate(num_corners=0, convex=True, star_shaped=True, degree=degree, scale=scale, **kwargs)
        elif n in ("bean", "kidney"):
            return cls.generate(num_corners=0, convex=False, star_shaped=True, degree=degree, scale=scale, **kwargs)
        elif n in ("horseshoe_smooth", "c_shape_smooth"):
            return cls.generate(num_corners=0, convex=False, star_shaped=False, degree=degree, scale=scale, **kwargs)

        elif n == "teardrop":
            return cls.generate(num_corners=1, convex=True, star_shaped=True, degree=degree, scale=scale, **kwargs)
        elif n == "heart":
            return cls.generate(num_corners=1, convex=False, star_shaped=True, degree=degree, scale=scale, **kwargs)

        elif n in ("lens", "football"):
            return cls.generate(num_corners=2, convex=True, star_shaped=True, degree=degree, scale=scale, **kwargs)
        elif n == "crescent":
            return cls.generate(num_corners=2, convex=False, star_shaped=True, degree=degree, scale=scale, **kwargs)

        elif n in ("triangle", "tri"):
            return cls.generate(num_corners=3, convex=True, star_shaped=True, curved=kwargs.pop("curved", False), degree=degree, scale=scale, **kwargs)
        elif n == "curved_triangle":
            kwargs.pop("curved", None)
            return cls.generate(num_corners=3, convex=True, star_shaped=True, curved=True, degree=degree, scale=scale, **kwargs)

        elif n in ("quad", "rectangle", "square"):
            return cls.generate(num_corners=4, convex=True, star_shaped=True, curved=kwargs.pop("curved", False), degree=degree, scale=scale, **kwargs)
        elif n == "curved_quad":
            kwargs.pop("curved", None)
            return cls.generate(num_corners=4, convex=True, star_shaped=True, curved=True, degree=degree, scale=scale, **kwargs)
        elif n == "indented_quad":
            kwargs.pop("curved", None)
            return cls.generate(num_corners=4, convex=False, star_shaped=True, curved=True, degree=degree, scale=scale, **kwargs)

        elif n == "pentagon":
            return cls.generate(num_corners=5, convex=True, star_shaped=True, curved=kwargs.pop("curved", False), degree=degree, scale=scale, **kwargs)
        elif n in ("star5", "star"):
            return cls.generate(num_corners=10, convex=False, star_shaped=True, curved=kwargs.pop("curved", False), degree=degree, scale=scale, **kwargs)

        elif n in ("lshape", "l_shape"):
            return cls.generate(num_corners=6, convex=False, star_shaped=False, curved=kwargs.pop("curved", False), degree=degree, scale=scale, **kwargs)
        elif n in ("curved_lshape", "curved_l_shape"):
            kwargs.pop("curved", None)
            return cls.generate(num_corners=6, convex=False, star_shaped=False, curved=True, degree=degree, scale=scale, **kwargs)
        elif n in ("u_channel", "horseshoe"):
            return cls.generate(num_corners=8, convex=False, star_shaped=False, curved=kwargs.pop("curved", False), degree=degree, scale=scale, **kwargs)

        elif n in ("annulus", "donut"):
            return cls.generate(num_corners=0, convex=True, star_shaped=True, num_holes=1, degree=degree, scale=scale, **kwargs)
        elif n == "plate_with_hole":
            return cls.generate(num_corners=4, convex=True, star_shaped=True, num_holes=1, curved=False, degree=degree, scale=scale, **kwargs)
        elif n == "plate_with_2holes":
            return cls.generate(num_corners=4, convex=True, star_shaped=True, num_holes=2, curved=kwargs.pop("curved", False), degree=degree, scale=scale, **kwargs)
        elif n in ("plate_5holes", "plate_with_5holes"):
            return cls.generate(num_corners=4, convex=True, star_shaped=True, num_holes=kwargs.pop("num_holes", 5), hole_shape=kwargs.pop("hole_shape", "mixed"), curved=kwargs.pop("curved", False), degree=degree, scale=scale, **kwargs)
        elif n in ("plate_8holes", "plate_with_8holes", "plate_8holes_weird"):
            return cls.generate(num_corners=4, convex=True, star_shaped=True, num_holes=kwargs.pop("num_holes", 8), hole_shape=kwargs.pop("hole_shape", "mixed"), curved=kwargs.pop("curved", False), degree=degree, scale=scale, **kwargs)
        elif n == "lshape_with_hole":
            return cls.generate(num_corners=6, convex=False, star_shaped=False, num_holes=1, curved=kwargs.pop("curved", False), degree=degree, scale=scale, **kwargs)
        elif n in ("lshape_3holes", "curved_lshape_3holes"):
            c_val = kwargs.pop("curved", ("curved" in n))
            return cls.generate(num_corners=6, convex=False, star_shaped=False, num_holes=kwargs.pop("num_holes", 3), hole_shape=kwargs.pop("hole_shape", "mixed"), curved=c_val, degree=degree, scale=scale, **kwargs)
        elif n in ("uchannel_4holes", "u_channel_4holes"):
            return cls.generate(num_corners=8, convex=False, star_shaped=False, num_holes=kwargs.pop("num_holes", 4), hole_shape=kwargs.pop("hole_shape", "mixed"), curved=kwargs.pop("curved", False), degree=degree, scale=scale, **kwargs)
        elif n in ("cross", "swiss_cross"):
            return cls.generate(num_corners=12, convex=False, star_shaped=False, curved=kwargs.pop("curved", False), degree=degree, scale=scale, **kwargs)
        elif n in ("cross_5holes", "swiss_cross_5holes"):
            return cls.generate(num_corners=12, convex=False, star_shaped=False, num_holes=kwargs.pop("num_holes", 5), hole_shape=kwargs.pop("hole_shape", "mixed"), curved=kwargs.pop("curved", False), degree=degree, scale=scale, **kwargs)
        elif n in ("comb", "comb_16corner"):
            return cls.generate(num_corners=16, convex=False, star_shaped=False, curved=kwargs.pop("curved", False), degree=degree, scale=scale, **kwargs)
        elif n in ("comb_4holes", "comb_16corner_4holes", "comb_16corner_4holes_curved", "curved_comb_4holes"):
            c_val = kwargs.pop("curved", ("curved" in n))
            return cls.generate(num_corners=16, convex=False, star_shaped=False, num_holes=kwargs.pop("num_holes", 4), hole_shape=kwargs.pop("hole_shape", "mixed"), curved=c_val, degree=degree, scale=scale, **kwargs)
        elif n in ("starburst", "star12"):
            return cls.generate(num_corners=12, convex=False, star_shaped=True, curved=kwargs.pop("curved", False), degree=degree, scale=scale, **kwargs)
        elif n in ("starburst_3holes", "star12_3holes", "starburst_12corner_3holes"):
            return cls.generate(num_corners=12, convex=False, star_shaped=True, num_holes=kwargs.pop("num_holes", 3), hole_shape=kwargs.pop("hole_shape", "mixed"), curved=kwargs.pop("curved", False), degree=degree, scale=scale, **kwargs)
        elif n in ("crescent_3holes", "crescent_with_3holes", "crescent_2corner_3holes"):
            return cls.generate(num_corners=2, convex=False, star_shaped=True, num_holes=kwargs.pop("num_holes", 3), hole_shape=kwargs.pop("hole_shape", "mixed"), degree=degree, scale=scale, **kwargs)

        else:
            raise ValueError(f"Unknown preset name '{name}'. Available presets: circle, ellipse, bean, teardrop, lens, crescent, triangle, curved_triangle, quad, curved_quad, indented_quad, pentagon, star5, starburst, lshape, u_channel, swiss_cross, comb, plate_with_hole, plate_with_2holes, plate_8holes, lshape_3holes, uchannel_4holes, cross_5holes, comb_4holes, starburst_3holes.")
