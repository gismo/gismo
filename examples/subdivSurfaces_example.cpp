/** @file subdivSurfaces_example.cpp

    @brief Tests different subdivision schemes

    Author(s): A. Mantzaflaris, M.Marsala, D.Tolis
*/

#include <gismo.h>

using namespace gismo;

int main(int argc, char** argv)
{
    std::string fn("off/cube.off");
    bool plot = false;
    bool dm = false;
    bool cc = false;
    bool ds = false;
    bool loop = false;
    index_t r(1);
    index_t dsopt(0);
    index_t loopopt(1);

    gsCmdLine cmd("Hi, give me a mesh");
    cmd.addPlainString("filename", "File containing mesh", fn);
    cmd.addSwitch("plot", "Plot the results", plot);
    cmd.addSwitch("dual", "Create the dual mesh (graph)", dm);
    cmd.addSwitch("cc", "Catmull-Clark Subdivision Scheme", cc);
    cmd.addSwitch("ds", "Doo-Sabin Subdivision Scheme", ds);
    cmd.addSwitch("loop", "Loop Subdivision Scheme", loop);
    cmd.addInt("d", "ds_opt", "Option for mask in Doo-Sabin subdivision scheme", dsopt);
    cmd.addInt("l", "loop_opt", "Option for mask in Loop subdivision scheme", loopopt);
    cmd.addInt("r", "ref", "Number of refinement steps", r);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

       
    if (!(ds || cc || loop))
        cc =true;

    // Read the mesh
    gsSurfMesh mesh;
    gsReadFile<>(fn,mesh);

    gsInfo << "Input: " << mesh.n_vertices() << " vertices, "
        << mesh.n_edges() << " edges, " << mesh.n_faces() << " faces. \n";

    gsSubdivScheme smesh(mesh);
    smesh.options().setInt("ds_opt", dsopt);
    smesh.options().setInt("loop_opt", loopopt);

    if (cc) // Catmull-Clark
    {
        smesh.options().setInt("scheme", 0);
        gsInfo << "Catmull-Clark subdivision "<<r <<" times.\n";
    }
    if (ds) // Doo-Sabin
    {
        smesh.options().setInt("scheme", 1);
        gsInfo << "Doo-Sabin subdivision "<<r <<" times.\n";
    }
    if (loop) // Loop
    {
        smesh.options().setInt("scheme", 2);
        gsInfo << "Loop subdivision "<<r <<" times.\n";
    }

    smesh.subdivide(r);
   
    mesh.write("mesh_in.off");
    if (dm) // Dual mesh
        mesh.dual_mesh(1);

    if (plot)
    {
        gsWriteParaview(mesh,"mesh_in", { });
        gsFileManager::open("mesh_in.vtk");
    }

    return EXIT_SUCCESS;
}
