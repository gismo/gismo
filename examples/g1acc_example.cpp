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
    std::string fn("off/octtorus.off");
    bool plot = false;
    index_t r(1);

    gsCmdLine cmd("Hi, give me a mesh");
    cmd.addPlainString("filename", "File containing mesh", fn);
    cmd.addSwitch("plot", "Plot the results", plot);
    cmd.addInt   ("r", "ref", "Number of refinement steps", r);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Create mesh
    
    gsSurfMesh mesh;
    gsReadFile<>(fn,mesh);

    typedef gsEigen::Vector<real_t,3> Point;
    typedef gsSurfMesh::Vertex Vertex;

/*
    mesh.add_vertex(Point(-1,-1,-1)); // Vertices
    mesh.add_vertex(Point( 1,-1,-1));
    mesh.add_vertex(Point(-1, 1,-1));
    mesh.add_vertex(Point( 1, 1,-1));
    mesh.add_vertex(Point(-1,-1, 1));
    mesh.add_vertex(Point( 1,-1, 1));
    mesh.add_vertex(Point(-1, 1, 1));
    mesh.add_vertex(Point( 1, 1, 1));
    mesh.add_quad( Vertex(0), Vertex(4), Vertex(6), Vertex(2) ); // Faces
    mesh.add_quad( Vertex(1), Vertex(3), Vertex(7), Vertex(5) );
    mesh.add_quad( Vertex(0), Vertex(1), Vertex(5), Vertex(4) );
    mesh.add_quad( Vertex(2), Vertex(6), Vertex(7), Vertex(3) );
    mesh.add_quad( Vertex(0), Vertex(2), Vertex(3), Vertex(1) );
    mesh.add_quad( Vertex(4), Vertex(5), Vertex(7), Vertex(6) );
    //auto sharp = mesh.add_halfedge_property<bool>("h:sharp"); // Sharp edges
    //markEdge(mesh, sharp, 0, 1);
    //markEdge(mesh, sharp, 2, 3);
    //markEdge(mesh, sharp, 0, 2);
    //markEdge(mesh, sharp, 1, 3);
    mesh.write("out_ds_0.off");
    */
    gsInfo << "Input: " << mesh.n_vertices() << " vertices, "
        << mesh.n_edges() << " edges, " << mesh.n_faces() << " faces. \n";
    // One step of Catmull-Clark subdivisions
    for( index_t i = 0; i<r; ++i)
    {
        gsSurfMesh new_mesh = mesh.ds_subdivide_robust();
        mesh = new_mesh; // to do: implement inPlace
        //mesh.ds_subdivide();
    }

    mesh.write("out_ds_1.off");

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
