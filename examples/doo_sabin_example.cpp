/** @file 

    @brief 

    Author(s): A. Mantzaflaris
*/

#include <gismo.h>
//#include <gsG1ACC/gsSurfScheme.h>
#ifdef gsParasolid_ENABLED
#include <gsParasolid/gsReadParasolid.h>
#include <gsParasolid/gsWriteParasolid.h>
#endif

using namespace gismo;

// Marks an edge
void markEdge(gsSurfMesh & mesh, gsSurfMesh::Halfedge_property<bool> & sharp,
              unsigned i, unsigned j)
{
    typedef gsSurfMesh::Vertex Vertex;
    typedef gsSurfMesh::Halfedge Halfedge;
    Halfedge he = mesh.find_halfedge(Vertex(i), Vertex(j));
    sharp[he] = true;
    he = mesh.opposite_halfedge(he);
    sharp[he] = true;
}
//class readOffClass : public gsFileData<> {
//
//    public:
//        void readOFF(std::string filename)
//        {
//
//            readOffFile(filename);
//        }
//};
int main(int argc, char** argv)
{
    std::string fn("off/cube_with_boundary.off");
    bool plot = false;
    bool dm = false;
    index_t r(1);

    gsCmdLine cmd("Hi, give me a mesh");
    cmd.addPlainString("filename", "File containing mesh", fn);
    cmd.addSwitch("plot", "Plot the results", plot);
    cmd.addSwitch("dual", "Create the dual mesh (graph)", dm);
    cmd.addInt   ("r", "ref", "Number of refinement steps", r);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Create mesh
    
    gsSurfMesh mesh;
    gsReadFile<>(fn,mesh);

    typedef gsEigen::Vector<real_t,3> Point;
    typedef gsSurfMesh::Vertex Vertex;

    gsMatrix<real_t> M;
    M = mesh.get_image_vertex_coeffs(4);

    gsDebug << M << std::endl;

    gsInfo << "Input: " << mesh.n_vertices() << " vertices, "
        << mesh.n_edges() << " edges, " << mesh.n_faces() << " faces. \n";
    // r steps of Doo-Sabin subdivisions

    // Check if we have boundaries in our mesh
    //mesh.boundary_reconstruction();

    for( index_t i = 0; i<r; ++i)
    {
        mesh.ds_subdivide_robust();
      
        /*mesh = new_mesh;*/ // TODO: implement inPlace
        //mesh.ds_subdivide();
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

    //mesh.ds_subdivide();
    //mesh.write("out_ds_1.off");

    //gsSurfScheme sf(mesh);
    //const gsMultiPatch<real_t> & mp = sf.acc5();
    //mp.closeGaps(1e-5);
    //gsInfo << "Output: "<< mp <<"\n";

    // Export as an XML file, that can be used with gsView to visualize on ParaView
    /*
    * 
    *mesh.ds_subdivide();
    mesh.write("out_ds_1.off");
    mesh.ds_subdivide();
    mesh.write("out_ds_2.off");
    mesh.ds_subdivide();
    mesh.write("out_ds_3.off");
    mesh.ds_subdivide();
    mesh.write("out_ds_4.off");
    mesh.ds_subdivide();
    mesh.write("out_ds_5.off");
    mesh.ds_subdivide();
    mesh.write("out_ds_6.off");
    
    
    
    
    
    
    */


    //gsWrite(mesh, "out_acc5");

    // Export as ParaSolid file
    //extensions::gsWriteParasolid(mp, "out_acc5");

    return EXIT_SUCCESS;
}
