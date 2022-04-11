
#include <gismo.h>

#include <gsMesh2/gsSurfMesh.cpp>

using namespace gismo;



void cc_subdivide(gsSurfMesh & mesh)
{
    gsSurfMesh::Vertex v;
    gsSurfMesh::Halfedge he;

    auto points = mesh.get_vertex_property<Point>("v:point");

    index_t env = mesh.n_vertices(); // edge vertices start here
        
    // loop over all edges, add edge points
    for (auto eit : mesh.edges())
    {
        he = mesh.halfedge(eit,0);
        Point pt( (points[mesh.from_vertex(he)]+points[mesh.to_vertex(he)]) / 2 );
        v = mesh.add_vertex(pt);
        mesh.insert_vertex(he,v);
    }

    index_t fnv = mesh.n_vertices(); // face vertices start here

    //mesh.write("out0.off");

    // loop over all faces, add face points
    for (auto fit : mesh.faces())
    {
        Point pt(0,0,0);
        auto fv = mesh.vertices(fit);
        for (auto vc = fv.begin(); vc!=fv.end(); ++vc, ++vc)
            pt += points[*vc];
        pt /= 4;
        mesh.add_vertex(pt);  // vertex gets shifted face id
    }

    //mesh.write("out1.off");

    int i = 0;
    for (auto fit : mesh.faces())
    {
        
        //gsInfo << "face "<< fit <<"\n"; for (auto fv : mesh.vertices(fit)) gsInfo <<" "<< fv.idx(); gsInfo <<"\n";

        v = gsSurfMesh::Vertex(fnv+(i++));//face vertex id ?
        //Start from an original vertex
        auto fv = mesh.vertices(fit).begin();
        //gsInfo <<"before "<< (*fv).idx() <<"\n";
        // [! points to to_vertex().idx()]
        if ( (*fv).idx() >= env) ++fv; //todo: add -> operator
        //++fv INF LOOP
        //assert ( (*fv).idx() < nv )
        //gsInfo <<"after "<< (*fv).idx() <<"\n";
        
        mesh.quad_split(fit,v,fv.he());
        //mesh.split(fit,v);
    }

    //mesh.write("out2.off");

    for (i = env; i!=fnv;++i)
    {
        v = gsSurfMesh::Vertex(i); //edge points
        auto vit = mesh.vertices(v);
        auto vcp = vit;
        Point pt(0,0,0);
        //std::cout << "e-Vertex "<< v <<"\n";
        if (vit) do
                 {
                     //std::cout << "   v: "<< *vit <<"\n";
                     pt += points[*vit];
                 } while (++vit != vcp);
        pt /= 4 ; // ==mesh.valence(v);
        points[v] = pt;
    }
    mesh.write("out3.off");

    // avg(F) + avg( P)
    for (i = 0; i!=env;++i)
    {
        v = gsSurfMesh::Vertex(i); // original vertices
        auto vit = mesh.halfedges(v);
        auto vcp = vit;
        Point E(0,0,0);// 2*E = sum(avg(f) + avg(e) = F + R [=> R=2*E-F ]
        Point F(0,0,0);// F
        //Point R(0,0,0);// R
        //std::cout << "o-Vertex "<< v <<"\n";
        if (vit) do
                 {
                     //std::cout << "   Edge pt "<< mesh.to_vertex(*vit) <<"\n";
                     E += points[ mesh.to_vertex(*vit) ];
                     F += points[ mesh.to_vertex(mesh.next_halfedge(*vit)) ];
                     //R += points[ mesh.to_vertex( mesh.next_halfedge(mesh.opposite_halfedge(mesh.next_halfedge(*vit))) ) ];
                     //std::cout << "   Face pt "<< mesh.to_vertex(mesh.next_halfedge(*vit)) <<"\n";
                 } while (++vit != vcp);
        auto n = mesh.valence(v);
        points[v] = ( (F/n) + (2.0*(2*E-F)/n) + (n-3)* points[v] ) / n ;
    }

    //mesh.write("out4.off");

}

gsMatrix<> mesh_property_to_matrix(const gsSurfMesh & mesh, const std::string & prop)
{
    auto & pv = mesh.get_vertex_property<Point>(prop).vector();
    gsMatrix<> res(3,pv.size());
    index_t i = 0;
    for( auto & pt : pv)
    {
        res(0,i)   = pt[0];
        res(1,i)   = pt[1];
        res(2,i++) = pt[2];
    }
    return res;
}


void cc_limit_points(gsSurfMesh & mesh)
{
    gsSurfMesh::Vertex v;
    gsSurfMesh::Halfedge he;

    auto points = mesh.get_vertex_property<Point>("v:point");
    auto limits = mesh.add_vertex_property<Point>("v:limit");
    Point pt;
    for (auto vit : mesh.vertices())
    {
        const int n = mesh.valence(vit);
        pt = n*n*points[vit];
        for ( auto he : mesh.halfedges(vit) )
        {
            pt += 4 * points[ mesh.to_vertex(he) ] + 
                points[ mesh.to_vertex(mesh.next_halfedge(he)) ];
        }
        limits[vit] = pt /( n*(n+5) );
    }
}

void cc_limit_unormal(gsSurfMesh & mesh)
{
    gsSurfMesh::Vertex v;
    gsSurfMesh::Halfedge he;

    auto points = mesh.get_vertex_property<Point>("v:point");
    auto limits = mesh.add_vertex_property<Point>("v:unormal");
    Point pt;
    for (auto vit : mesh.vertices())
    {
        const int n = mesh.valence(vit);
        pt = n*n*points[vit];
        for ( auto he : mesh.halfedges(vit) )
        {// TODO
            pt += 4 * points[ mesh.to_vertex(he) ] + 
                points[ mesh.to_vertex(mesh.next_halfedge(he)) ];
        }
        limits[vit] = pt /( n*(n+5) );
    }
}


int main(int argc, char** argv)
{
    index_t r(0);
    std::string fn("");

    //! [Parse Command line]
    gsCmdLine cmd("Hi, give me a mesh");
    cmd.addInt   ("r", "ref", "Number of refinement steps", r);
    cmd.addPlainString("filename", "File containing data to draw (.xml or third-party)", fn);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsSurfMesh mesh;
    mesh.read(fn);
    std::cout << "Input:" << std::endl;
    std::cout << "  vertices: " << mesh.n_vertices() << std::endl;
    std::cout << "  edges:    " << mesh.n_edges()    << std::endl;
    std::cout << "  faces:    " << mesh.n_faces()    << std::endl;

    for( index_t i = 0; i<r; ++i)
        cc_subdivide(mesh);

    std::cout << "Output:" << std::endl;
    std::cout << "  vertices: " << mesh.n_vertices() << std::endl;
    std::cout << "  edges:    " << mesh.n_edges()    << std::endl;
    std::cout << "  faces:    " << mesh.n_faces()    << std::endl;

    cc_limit_points(mesh);

    gsMatrix<> lp = mesh_property_to_matrix(mesh,"v:limit");
    gsWrite(lp,"out_values");
    
    mesh.write("out_mesh.off");

    return 0;
}


//=============================================================================
