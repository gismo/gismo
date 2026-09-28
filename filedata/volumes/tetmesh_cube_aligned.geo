// Cube [-0.5,0.5]^3, axis-aligned. Test geometry for
// immersed_tetmesh_clip_quadrature_example (--case cube): with the default
// background box [-1,1]^3 and n = 4*2^r, the faces +-0.5 fall exactly on the
// grid lines X[n/4] and X[3n/4] (h dyadic), so this mesh exercises the
// face-on-knot-plane degeneracy the driver's T9 check looks for.
// Generated with: gmsh tetmesh_cube_aligned.geo -3 -format msh41 -o tetmesh_cube_aligned.msh
SetFactory("OpenCASCADE");
Mesh.Binary = 0;
Mesh.MshFileVersion = 4.1;
Mesh.MeshSizeMax = 0.125;

Box(1) = {-0.5, -0.5, -0.5, 1, 1, 1};
