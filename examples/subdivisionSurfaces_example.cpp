/** @file subdivisionSurfaces_example.cpp

    @brief Tests different subdivision schemes

    Option flag gives different result depending on the scheme:

    Doo-Sabin options:
        0 - Interpolatory boundary using Chaikin scheme.
        1 - Trimmed boundary using Doo-Sabin scheme.
    
    Loop options:
       0 - Simplified Loop's scheme. (cf. book Warren, Weimer 2002) 
       1 - Original Loop's scheme.  (cf. book Loop 1987) 

    Author(s): A. Mantzaflaris, M.Marsala, D.Tolis, L. Mussmaecher
*/

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
    index_t option(0);
    bool normalize(true);

    // Read command-line options
    gsCmdLine cmd("Hi, give me a mesh");
    cmd.addPlainString("filename", "File containing mesh", fn);
    cmd.addSwitch("plot", "Plot the results", plot);
    cmd.addSwitch("dual", "Create the dual mesh (graph)", dm);
    cmd.addString("s", "scheme", "Choice of Subdivision Scheme", scheme_name);
    cmd.addInt("r", "ref", "Number of refinement steps", r);
    cmd.addInt("o","option","Option on subdivision scheme", option);
    cmd.addSwitch("normalize", "Normalize limit normals and tangents", normalize);
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
            gsInfo << "Catmull-Clark subdivision "<< r <<" times.\n";
    } else if (scheme_name == "Doo-Sabin"){
            scheme = new gsDooSabin();
            scheme->options().setInt("ds.boundaryMask", option); // option on doo-sabin boundary treatment
            gsInfo << "Doo-Sabin subdivision "<< r <<" times.\n";
        }
    else if (scheme_name == "Loop") {
            scheme = new gsLoop();
            scheme->options().setInt("loop.maskType", option); // option on loop mask type
            gsInfo << "Loop subdivision "<< r <<" times.\n";
    }
    else {
            scheme = new gsCatmullClark();
            gsInfo << "Catmull-Clark subdivision "<< r <<" times.\n";
    }

    scheme->assign(&mesh);
    scheme->subdivide(r);  
   
    mesh.write("mesh_out.off");
    if (dm) // Dual mesh
        mesh.dual_mesh();

    gsInfo << "Output: " << mesh.n_vertices() << " vertices, "
        << mesh.n_edges() << " edges, " << mesh.n_faces() << " faces. \n";

    if (plot)
    {
        gsWriteParaview(mesh,"mesh_out", { });
        gsFileManager::open("mesh_out.vtk");
    }

    delete scheme;

    return EXIT_SUCCESS;
}
