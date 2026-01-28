/** @file subdivisionSurfaces_example.cpp

    @brief Tests different subdivision schemes

    Author(s): A. Mantzaflaris, M.Marsala, D.Tolis, L. Mussmaecher
*/

#include "gsMesh2/gsCatmullClark.h"
#include "gsMesh2/gsDooSabin.h"
#include "gsMesh2/gsSubdivisionScheme.h"
#include <gismo.h>

using namespace gismo;

int main(int argc, char** argv)
{
    // Set-Up options for the command-line options.
    std::string fn("off/cube.off");
    bool plot = false;
    bool dm = false;
    std::string scheme_name = "Catmull-Clark";
    index_t r(1);
    index_t dsopt(0);
    index_t loopopt(1);

    // Read command-line options
    gsCmdLine cmd("Hi, give me a mesh");
    cmd.addPlainString("filename", "File containing mesh", fn);
    cmd.addSwitch("plot", "Plot the results", plot);
    cmd.addSwitch("dual", "Create the dual mesh (graph)", dm);
    cmd.addString("s", "scheme", "Choice of Subdivision Scheme", scheme_name);
    cmd.addInt("d", "ds.boundaryMask", "Option for mask in Doo-Sabin subdivision scheme", dsopt);
    cmd.addInt("l", "loop.maskType", "Option for mask in Loop subdivision scheme", loopopt);
    cmd.addInt("r", "ref", "Number of refinement steps", r);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
      

    // Read the mesh
    gsSurfMesh mesh;
    auto _readFile = gsReadFile<>(fn,mesh);

    gsInfo << "Input: " << mesh.n_vertices() << " vertices, "
        << mesh.n_edges() << " edges, " << mesh.n_faces() << " faces. \n";

    // Currently, Catmull-Clark is always chosen.
    // TODO: Fit other schemes into this class hierarchy and implement the example properly.
    gsSubdivisionScheme* scheme;

    if (scheme_name == "Catmull-Clark"){
            scheme = new gsCatmullClark();
    } else if (scheme_name == "Doo-Sabin"){
            scheme = new gsDooSabin();
    } else {
            scheme = new gsCatmullClark();
    }

    scheme->subdivide(mesh, r);  
   
    mesh.write("mesh_out.off");
    if (dm) // Dual mesh
        mesh.dual_mesh();

    if (plot)
    {
        gsWriteParaview(mesh,"mesh_out", { });
        gsFileManager::open("mesh_out.vtk");
    }

    delete scheme;

    return EXIT_SUCCESS;
}
