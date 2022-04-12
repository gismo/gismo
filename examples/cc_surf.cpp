
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

inline index_t idx(index_t i, index_t j, index_t s, index_t sz)
{
    switch (s)
    {
    case 0:
        return sz * j + i;
    case 1:
        return sz * i + (sz-1-j);
    case 2:
        return sz * (sz-1-j) + (sz-1-i);
    case 3:
        return sz * (sz-1-i) + j;
    default:
        GISMO_ERROR("idx error");
    }
}

gsMultiPatch<real_t> cc_acc3(gsSurfMesh & mesh)
{
    auto points = mesh.get_vertex_property<Point>("v:point");
    gsMultiPatch<real_t> mp;
    gsKnotVector<> kv(0,1,0,4);
    gsTensorBSplineBasis<2> bb(kv,kv);
    gsMatrix<> coefs(16,3); coefs.setZero();
    gsSurfMesh::Halfedge h2;
    gsSurfMesh::Vertex v;
    real_t n;
//#   pragma omp parallel for private(pt,n,he)
    for (auto fit = mesh.faces_begin(); fit!= mesh.faces_end(); ++fit)
    {
        index_t s = 0;
        for ( auto he : mesh.halfedges(*fit) )
        {
            auto c00 = coefs.row(idx(0,0,s,4)).transpose();
            auto c10 = coefs.row(idx(1,0,s,4)).transpose();
            auto c01 = coefs.row(idx(0,1,s,4)).transpose();
            auto c11 = coefs.row(idx(1,1,s,4)).transpose();
            v = mesh.from_vertex(he);
            n = mesh.valence(v);
            c00 = n*n*points[v];
            c10 = c01 = c11 = n * points[v];
            for (auto h : mesh.halfedges(v))
                c00 += 4 * points[ mesh.to_vertex(h) ] +
                    points[ mesh.to_vertex(mesh.next_halfedge(h)) ] ;

            c11 += 2 * points[mesh.to_vertex(he)] +
                points[mesh.to_vertex(mesh.next_halfedge(he))] +
                2 * points[mesh.to_vertex(mesh.ccw_rotated_halfedge(he))];

            c10 += 2 * points[mesh.to_vertex(he)] +
                0.5 * points[mesh.to_vertex(mesh.next_halfedge(he))] +
                points[mesh.to_vertex(mesh.ccw_rotated_halfedge(he))]+
                points[mesh.to_vertex(mesh.cw_rotated_halfedge(he)) ]+
                0.5 * points[mesh.to_vertex(mesh.next_halfedge(mesh.cw_rotated_halfedge(he))) ];

            h2 = mesh.ccw_rotated_halfedge(he);
            c01 += 2 * points[mesh.to_vertex(h2)] +
                0.5 * points[mesh.to_vertex(mesh.next_halfedge(h2))] +
                points[mesh.to_vertex(mesh.ccw_rotated_halfedge(h2))]+
                points[mesh.to_vertex(mesh.cw_rotated_halfedge(h2)) ]+
                0.5 * points[mesh.to_vertex(mesh.next_halfedge(mesh.cw_rotated_halfedge(h2))) ];

            c00 /= n*(n+5);
            c10 /= n+5;
            c01 /= n+5;
            c11 /= n+5;
            ++s;
        }
        mp.addPatch( bb.makeGeometry(coefs) );
    }
    mp.computeTopology();
    return mp;
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

void error_mesh_multipatch(const gsSurfMesh & mesh,
                           const gsMultiPatch<> & mp)
{
    index_t np = mesh.n_vertices();
    auto limits = mesh.get_vertex_property<Point>("v:limit");
    auto unv    = mesh.get_vertex_property<Point>("v:unormal");

    gsMapData<> gdata;
    gdata.addFlags( NEED_VALUE|NEED_NORMAL );
    std::pair<index_t,gsVector<> > cp;
    gsMatrix<> pt, err(1,np), nerr(1,np);

    index_t i = 0;
    for (auto v : mesh.vertices())
    {
        gsInfo << "\r"<< i << std::flush;
        pt = limits[v];
        cp = mp.closestPointTo(pt, 1e-15 );
        gdata.points = cp.second;
        gdata.patchId = cp.first;
        mp.patch(cp.first).computeMap(gdata);
        err.at(i) = (pt-gdata.eval(0)).norm();
        nerr.at(i) = (unv[v]-gdata.normal(0).normalized()).norm();
        ++i;
    }
    gsInfo <<"\nGeometry error   : "<< err.norm() <<"\n";
    gsInfo <<  "Unit normal error: "<< nerr.norm() <<"\n";
    //gsWriteParaview( *msh, "ccmesh", err);
}

int main(int argc, char** argv)
{
    index_t r(0), s(0);
    std::string fn("");

    //! [Parse Command line]
    gsCmdLine cmd("Hi, give me a mesh");
    cmd.addInt   ("r", "ref", "Number of refinement steps", r);
    cmd.addInt   ("s", "sam", "Number of sampling steps",   s);
    cmd.addPlainString("filename", "File containing data to draw (.xml or third-party)", fn);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsSurfMesh mesh;
    mesh.read(fn);
    gsInfo << "Input: "<<mesh.n_vertices()<< " vertices, "
           << mesh.n_edges() << " edges, " << mesh.n_faces() << " faces. \n";

    for( index_t i = 0; i<r; ++i)
        cc_subdivide(mesh);

    gsInfo << "Extracting "<<mesh.n_faces()<<" Bezier patches..\n";
    gsMultiPatch<> mp = cc_acc3(mesh);
    gsWrite(mp,"out_acc3");

    for( index_t i = 0; i<s; ++i)
        cc_subdivide(mesh);

    gsInfo << "Sampled at "<<mesh.n_vertices()<< " points. \n";

    gsInfo << "Getting limit points..\n";
    cc_limit_points  (mesh);
    gsInfo << "Getting limit normals..\n";
    cc_limit_unormals(mesh);

    /*
    gsInfo << "Writing to files..\n";

    mesh.write("out_mesh.off");

    gsMatrix<> lp = mesh_property_to_matrix(mesh,"v:limit");
    gsWrite(lp,"out_values");
    gsMatrix<> ln = mesh_property_to_matrix(mesh,"v:unormal");
    gsWrite(ln,"out_normals");
    */

    gsInfo << "Computing errors..\n";
    error_mesh_multipatch(mesh,mp);

    return 0;
}


//=============================================================================
