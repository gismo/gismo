// Axis-aligned square of half-side a = 0.5 (no rotation), minus the disk centred at
// (0.05,-0.03), R = 0.35, same disk as square_minus_disk_mesh.geo. The disk's x-range
// [-0.3, 0.4] and y-range [-0.38, 0.32] lie inside the square. With the default --n0 8,
// h = 0.25*2^-r, so x, y = +-0.5 are knot lines for every r.
// Generated with: gmsh square_on_grid_minus_disk_mesh.geo -2 -format stl -o square_on_grid_minus_disk_mesh.stl
Mesh.Binary = 0;          // ASCII STL: the gsSurfMesh binary-STL reader is float/double-broken
a = 0.5; th = 0*Pi/180; cx = 0.05; cy = -0.03; R = 0.35; lc = 0.08;
Point(1) = {Cos(th)*(-a) - Sin(th)*(-a), Sin(th)*(-a) + Cos(th)*(-a), 0, lc};
Point(2) = {Cos(th)*( a) - Sin(th)*(-a), Sin(th)*( a) + Cos(th)*(-a), 0, lc};
Point(3) = {Cos(th)*( a) - Sin(th)*( a), Sin(th)*( a) + Cos(th)*( a), 0, lc};
Point(4) = {Cos(th)*(-a) - Sin(th)*( a), Sin(th)*(-a) + Cos(th)*( a), 0, lc};
Point(5) = {cx, cy, 0, lc};  Point(6) = {cx+R, cy, 0, lc};  Point(7) = {cx, cy+R, 0, lc};
Point(8) = {cx-R, cy, 0, lc};  Point(9) = {cx, cy-R, 0, lc};
Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
Circle(5) = {6,5,7}; Circle(6) = {7,5,8}; Circle(7) = {8,5,9}; Circle(8) = {9,5,6};
Curve Loop(1) = {1,2,3,4}; Curve Loop(2) = {5,6,7,8};
Plane Surface(1) = {1,2};
