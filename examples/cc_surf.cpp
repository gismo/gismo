
#include <gismo.h>

#include <gsMesh2/gsSurfMesh.cpp>

using namespace gismo;



void cc_subdivide(gsSurfMesh & mesh)
{
    gsSurfMesh::Vertex v;
    gsSurfMesh::Halfedge he;

    // reserve vertices, edges, faces
    mesh.reserve( mesh.n_vertices()+mesh.n_edges()+mesh.n_faces(),
                  2*mesh.n_edges(), 4*mesh.n_faces() );
                  
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
        v = gsSurfMesh::Vertex(fnv+(i++));//face vertex id ?
        //Start from an original vertex
        auto fv = mesh.vertices(fit).begin();
        if ( (*fv).idx() >= env) ++fv; //todo: add -> operator
        //assert ( (*fv).idx() < nv )
        mesh.quad_split(fit,v,fv.he());
    }

    //mesh.write("out2.off");

    //parallel
#   pragma omp parallel for private(v)
    for (i = env; i!=fnv;++i)
    {
        v = gsSurfMesh::Vertex(i); //edge points
        auto vit = mesh.vertices(v);
        auto vcp = vit;
        Point pt(0,0,0);
        if (vit) do
                 {
                     pt += points[*vit];
                 } while (++vit != vcp);
        pt /= 4 ; // =mesh.valence(v);
        points[v] = pt;
    }

    //mesh.write("out3.off");

#   pragma omp parallel for private(v)
    for (i = 0; i!=env;++i)
    {
        v = gsSurfMesh::Vertex(i); // original vertices
        auto vit = mesh.halfedges(v);
        auto vcp = vit;
        Point E(0,0,0);// 2*E = avg(f) + avg(e) = F + R [=> R=2*E-F ]
        Point F(0,0,0);// F
        if (vit)
            do
            {
                E += points[ mesh.to_vertex(*vit) ];
                F += points[ mesh.to_vertex(mesh.next_halfedge(*vit)) ];
            } while (++vit != vcp);
        auto n = mesh.valence(v);
        //points[v] = ( (F/n) + (2*(2*E-F)/n) + (n-3) * points[v] ) / n ;
        points[v] = ( (n*(n-3))*points[v] + 4*E - F ) / (n*n);
    }

    //mesh.write("out4.off");
}

gsMatrix<> mesh_property_to_matrix(const gsSurfMesh & mesh, const std::string & prop)
{
    auto & pv = mesh.get_vertex_property<Point>(prop).vector();
    gsMatrix<> res(3,pv.size());
    index_t i = 0;
    for( auto & pt : pv)
        res.col(i++) = pt;
    return res;
}


void cc_limit_points(gsSurfMesh & mesh)
{
    auto points = mesh.get_vertex_property<Point>("v:point");
    auto limits = mesh.add_vertex_property<Point>("v:limit");
    Point pt;
    real_t n;
#   pragma omp parallel for private(pt,n)
    for (auto vit = mesh.vertices_begin(); vit!= mesh.vertices_end(); ++vit)
    {
        n = mesh.valence(*vit);
        pt = n*n*points[*vit];
        for ( auto he : mesh.halfedges(*vit) )
        {
            pt += 4 * points[ mesh.to_vertex(he) ] + 
                points[ mesh.to_vertex(mesh.next_halfedge(he)) ];
        }
        limits[*vit] = pt / (n*(n+5));
    }
}

void cc_limit_unormals(gsSurfMesh & mesh)
{
    gsSurfMesh::Vertex v;
    gsSurfMesh::Halfedge he, h2;

    auto points = mesh.get_vertex_property<Point>("v:point");
    auto limits = mesh.add_vertex_property<Point>("v:unormal");
    Point t1, t2;
    real_t c1, c2, cc1, cc2;
    index_t i;
#   pragma omp parallel for private(t1,t2,c1,c2,cc1,cc2,i)
    for (auto vit = mesh.vertices_begin(); vit!= mesh.vertices_end(); ++vit)
    {
        const real_t n = mesh.valence(*vit);
        const real_t cospin = math::cos(EIGEN_PI/n);
        cc2 = 1 / ( n * math::sqrt(4+cospin*cospin) );
        cc1 = 1/n + cospin*cc2;
        t1.setZero();
        t2.setZero();
        i = 0;
        for ( auto he : mesh.halfedges(*vit) )
        {
            h2 = mesh.ccw_rotated_halfedge(he);
            c1  = math::cos( 2*i   *EIGEN_PI/n)*cc1;
            c2  = math::cos((2*i+1)*EIGEN_PI/n)*cc2;
            t1 += c1 * points[ mesh.to_vertex(he ) ]
                + c2 * points[ mesh.to_vertex(mesh.next_halfedge(he)) ];
            t2 += c1 * points[ mesh.to_vertex(h2 ) ]
                + c2 * points[ mesh.to_vertex(mesh.next_halfedge(h2)) ];
            ++i;
        }
        limits[*vit] = t1.cross(t2).normalized();
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
    gsInfo << "Input:" << "\n";
    gsInfo << "  vertices: " << mesh.n_vertices() << "\n";
    gsInfo << "  edges:    " << mesh.n_edges()    << "\n";
    gsInfo << "  faces:    " << mesh.n_faces()    << "\n";

    for( index_t i = 0; i<r; ++i)
        cc_subdivide(mesh);

    gsInfo << "Output:" << "\n";
    gsInfo << "  vertices: " << mesh.n_vertices() << "\n";
    gsInfo << "  edges:    " << mesh.n_edges()    << "\n";
    gsInfo << "  faces:    " << mesh.n_faces()    << "\n";

    gsInfo << "Getting limit points..\n";
    cc_limit_points  (mesh);
    gsInfo << "Getting limit normals..\n";
    cc_limit_unormals(mesh);

    gsInfo << "Writing to files..\n";

    mesh.write("out_mesh.off");

    gsMatrix<> lp = mesh_property_to_matrix(mesh,"v:limit");
    gsWrite(lp,"out_values");
    gsMatrix<> ln = mesh_property_to_matrix(mesh,"v:unormal");
    gsWrite(ln,"out_normals");


    return 0;
}


//=============================================================================
