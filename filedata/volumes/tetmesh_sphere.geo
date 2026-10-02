// Sphere centred at (0.03,-0.02,0.01), radius 0.55. Test geometry for
// immersed_tetmesh_clip_quadrature_example (--case sphere): a curved
// boundary with no exact analytic moments, only the INFO volume check
// against 4/3 pi r^3.
// Generated with: gmsh tetmesh_sphere.geo -3 -format msh41 -o tetmesh_sphere.msh
SetFactory("OpenCASCADE");
Mesh.Binary = 0;
Mesh.MshFileVersion = 4.1;
Mesh.MeshSizeMax = 0.08;

Sphere(1) = {0.03,-0.02,0.01, 0.55};
