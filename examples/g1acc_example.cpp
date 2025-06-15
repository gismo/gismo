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

int main(int argc, char** argv)
{
    // Create mesh
    gsSurfMesh mesh;
    typedef gsEigen::Vector<real_t,3> Point;
    typedef gsSurfMesh::Vertex Vertex;

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
    mesh.write("initial_cube.off");
    // One step of Catmull-Clark subdivisions
    mesh.ds_subdivide();

    gsInfo << "Input: "<<mesh.n_vertices()<< " vertices, "
           << mesh.n_edges() << " edges, " << mesh.n_faces() << " faces. \n";

    //gsSurfScheme sf(mesh);
    //const gsMultiPatch<real_t> & mp = sf.acc5();
    //mp.closeGaps(1e-5);
    //gsInfo << "Output: "<< mp <<"\n";

    // Export as an XML file, that can be used with gsView to visualize on ParaView

    mesh.write("out_ds.off");
    //gsWrite(mesh, "out_acc5");

    // Export as ParaSolid file
    //extensions::gsWriteParasolid(mp, "out_acc5");

    return EXIT_SUCCESS;
}
