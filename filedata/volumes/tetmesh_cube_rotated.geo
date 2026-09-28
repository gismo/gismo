// Cube of side s = 0.8, centred at the origin, rotated about z by 17 deg then
// about x by 29 deg (both about the origin), then translated by
// c = (0.03,-0.02,0.01). Test geometry for immersed_tetmesh_clip_quadrature_example
// (--case rotcube).
// Generated with: gmsh tetmesh_cube_rotated.geo -3 -format msh41 -o tetmesh_cube_rotated.msh
SetFactory("OpenCASCADE");
Mesh.Binary = 0;
Mesh.MshFileVersion = 4.1;
Mesh.MeshSizeMax = 0.1;

s = 0.8;
Box(1) = {-s/2, -s/2, -s/2, s, s, s};
Rotate {{0,0,1}, {0,0,0}, 17*Pi/180} { Volume{1}; }
Rotate {{1,0,0}, {0,0,0}, 29*Pi/180} { Volume{1}; }
Translate {0.03,-0.02,0.01} { Volume{1}; }
