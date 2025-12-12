/** @file 

    @brief 

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

    // Create mesh
    

    gsSurfMesh mesh;
    gsReadFile<>(fn,mesh);

    gsInfo << "Input: " << mesh.n_vertices() << " vertices, "
        << mesh.n_edges() << " edges, " << mesh.n_faces() << " faces. \n";

    gsSubdivScheme smesh(mesh);
    smesh.options().setInt("ds_opt", dsopt);
    smesh.options().setInt("loop_opt", loopopt);

    if (ds==true) { // Doo-Sabin
        for (index_t i = 0; i < r; ++i)
        {
            smesh.ds_subdivide();
        }
    }
    else if (cc==true) {  // Catmull-Clark
        for (index_t i = 0; i < r; ++i)
        {
            smesh.cc_subdivide();
        }
    }
    else if (loop==true) { // Loop
        for (index_t i = 0; i < r; ++i)
        {
            smesh.loop_subdivide();
        }
    }
    

    mesh.write("mesh_in.off");
    if (dm) {
        mesh.dual_mesh(1);

    }
    if (plot)
    {
        gsWriteParaview(mesh,"mesh_in", { });
        gsFileManager::open("mesh_in.vtk");
    }

    


    return EXIT_SUCCESS;
}
