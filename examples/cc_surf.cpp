/** @file cc_surf

    @brief CC surface

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gismo.h>

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
    Point tmp;
    for (auto eit : mesh.edges())
    {
        he  = mesh.halfedge(eit,0);
        tmp = (points[mesh.from_vertex(he)]+points[mesh.to_vertex(he)]) / 2;
        v   = mesh.add_vertex(tmp);
        mesh.insert_vertex(he,v);
    }

    index_t fnv = mesh.n_vertices(); // face vertices start here

    // loop over all faces, add face points
    for (auto fit : mesh.faces())
    {
        auto fv = mesh.vertices(fit);
        tmp.setZero();
        for (auto vc = fv.begin(); vc!=fv.end(); ++vc, ++vc)
            tmp += points[*vc];
        tmp /= 4;
        mesh.add_vertex(tmp);  // vertex gets shifted face id
    }

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

#   pragma omp parallel for private(v,tmp)
    for (i = env; i!=fnv;++i)
    {
        v = gsSurfMesh::Vertex(i); //edge points
        auto vit = mesh.vertices(v);
        auto vcp = vit;
        tmp.setZero();
        if (vit) do
                 {
                     tmp += points[*vit];
                 } while (++vit != vcp);
        tmp /= 4 ; // =mesh.valence(v);
        points[v] = tmp;
    }

#   pragma omp parallel for private(v)
    for (i = 0; i!=env;++i)
    {
        v = gsSurfMesh::Vertex(i); // original vertices
        auto vit = mesh.halfedges(v);
        auto vcp = vit;
        auto & pt = points[v];//original vertex positions are computed using new edge/face points only
        auto n = mesh.valence(v);
        //formula: pt = ( (n*(n-3))*points[v] + 4*E - F ) / (n*n);
        pt *= n*(n-3);
        if (vit)
            do
            { //pt += 4*E-F
                pt += 4*points[ mesh.to_vertex(*vit) ]
                    - points[ mesh.to_vertex(mesh.next_halfedge(*vit)) ];
            } while (++vit != vcp);
        pt /= n*n;
    }
}

void cc_limit_points(gsSurfMesh & mesh, bool swap_pts = false)
{
    auto points = mesh.get_vertex_property<Point>("v:point");
    auto limits = mesh.add_vertex_property<Point>("v:limit");
    real_t n;
#   pragma omp parallel for private(n)
    for (auto vit = mesh.vertices_begin(); vit!= mesh.vertices_end(); ++vit)
    {
        n = mesh.valence(*vit);
        auto & pt = limits[*vit];
        pt = n*n*points[*vit];
        for ( auto he : mesh.halfedges(*vit) )
        {
            pt += 4 * points[ mesh.to_vertex(he) ] +
                points[ mesh.to_vertex(mesh.next_halfedge(he)) ];
        }
        pt /= (n*(n+5));
    }

    if (swap_pts)//vertices are moved to their limit positions
    {
        mesh.swap_vertex_property("v:point","v:limit");
        mesh.rename_vertex_property(points,"v:original_point");
    }
}

void cc_limit_normals(gsSurfMesh & mesh)
{
    gsSurfMesh::Vertex v;
    gsSurfMesh::Halfedge he, h2;

    auto points = mesh.get_vertex_property<Point>("v:point");
    auto limits = mesh.add_vertex_property<Point>("v:normal");
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

inline index_t face_pt_idx(index_t i, index_t j, index_t s, index_t sz)
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
    gsKnotVector<> kv(0,1,0,4);//cubic degree
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
            auto c00 = coefs.row(face_pt_idx(0,0,s,4)).transpose();
            auto c10 = coefs.row(face_pt_idx(1,0,s,4)).transpose();
            auto c01 = coefs.row(face_pt_idx(0,1,s,4)).transpose();
            auto c11 = coefs.row(face_pt_idx(1,1,s,4)).transpose();
            v = mesh.from_vertex(he);
            n = mesh.valence(v);
            c00 = n*n*points[v];
            c10 = c01 = c11 = n * points[v];
            for (auto h : mesh.halfedges(v))
                c00 += 4 * points[ mesh.to_vertex(h) ] +
                    points[ mesh.to_vertex(mesh.next_halfedge(h)) ] ;
            c00 /= n*(n+5);

            c11 += 2 * points[mesh.to_vertex(he)] +
                points[mesh.to_vertex(mesh.next_halfedge(he))] +
                2 * points[mesh.to_vertex(mesh.ccw_rotated_halfedge(he))];
            c11 /= n+5;

            c10 += 2 * points[mesh.to_vertex(he)] +
                0.5 * points[mesh.to_vertex(mesh.next_halfedge(he))] +
                points[mesh.to_vertex(mesh.ccw_rotated_halfedge(he))]+
                points[mesh.to_vertex(mesh.cw_rotated_halfedge(he)) ]+
                0.5 * points[mesh.to_vertex(mesh.next_halfedge(mesh.cw_rotated_halfedge(he))) ];
            c10 /= n+5;

            h2 = mesh.ccw_rotated_halfedge(he);
            c01 += 2 * points[mesh.to_vertex(h2)] +
                0.5 * points[mesh.to_vertex(mesh.next_halfedge(h2))] +
                points[mesh.to_vertex(mesh.ccw_rotated_halfedge(h2))]+
                points[mesh.to_vertex(mesh.cw_rotated_halfedge(h2)) ]+
                0.5 * points[mesh.to_vertex(mesh.next_halfedge(mesh.cw_rotated_halfedge(h2))) ];
            c01 /= n+5;

            ++s;
        }
        mp.addPatch( bb.makeGeometry(coefs) );
    }
    mp.computeTopology();
    return mp;
}

gsMatrix<> mesh_property_to_matrix(const gsSurfMesh & mesh,
                                   const std::string & prop)
{
    auto & pv = mesh.get_vertex_property<Point>(prop).vector();
    gsMatrix<> res(3,pv.size());
    index_t i = 0;
    for( auto & pt : pv)
        res.col(i++) = pt;
    return res;
}

void error_mesh_multipatch(gsSurfMesh & mesh,
                           const gsMultiPatch<> & mp)
{
    auto limit = mesh.get_vertex_property<Point>("v:limit");
    auto unv    = mesh.get_vertex_property<Point>("v:normal");
    auto gerr    = mesh.add_vertex_property<real_t>("v:geometry_error");
    auto nerr    = mesh.add_vertex_property<real_t>("v:normal_error");

    gsMapData<> gdata;
    gdata.addFlags( NEED_VALUE|NEED_NORMAL );
    std::pair<index_t,gsVector<> > cp;
    gsMatrix<> pt;
    real_t tge = 0, tne = 0;
    index_t i = 0;
    for (auto v : mesh.vertices())
    {
        gsInfo << "\r"<< i << std::flush;
        pt = limit[v];
        cp = mp.closestPointTo(pt, 1e-15 );
        gdata.points = cp.second;
        gdata.patchId = cp.first;
        mp.patch(cp.first).computeMap(gdata);
        tge = math::max(tge, gerr[v] = (pt-gdata.eval(0)).norm() );
        tne = math::max(tne, nerr[v] = (unv[v]-gdata.normal(0).normalized()).norm() );
        ++i;
    }
    gsInfo <<"\nMax geometry error   : "<< tge <<"\n";
    gsInfo <<  "Max unit-normal error: "<< tne <<"\n";
}

int main(int argc, char** argv)
{
    index_t r(0), s(0);
    std::string fn("off/octtorus.off");
    std::string fn_patch("");
    bool plot = false;
    
    //! [Parse Command line]
    gsCmdLine cmd("Hi, give me a CC mesh");
    cmd.addInt   ("c", "ref", "Number of refinement steps", r);
    cmd.addInt   ("s", "sam", "Number of sampling refinement steps",   s);
    cmd.addPlainString("filename", "File containing mesh", fn);
    cmd.addString("g","geometry", "File containing multipatch geometry to compare)", fn_patch);
    cmd.addSwitch("plot", "Plot the results", plot);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsSurfMesh mesh;
    //mesh.read(fn);
    gsReadFile<>(fn,mesh);

    gsInfo << "Input: "<<mesh.n_vertices()<< " vertices, "
           << mesh.n_edges() << " edges, " << mesh.n_faces() << " faces. \n";

    //subdivision before creating ACC3 patches
    for( index_t i = 0; i<r; ++i)
        cc_subdivide(mesh);

    gsMultiPatch<> mp;
    if (!fn_patch.empty())
    {   //read patches from file
        gsReadFile<>(fn_patch,mp);
        gsInfo << "Input: "<<mp.nPatches()<<" spline patches.\n";
    }
    else //if geometry is not given then use ACC3
    {
        gsInfo << "Extracting "<<mesh.n_faces()<<" Bezier patches (ACC3)..\n";
        mp = cc_acc3(mesh);
    }

    for( index_t i = 0; i<s; ++i)
        cc_subdivide(mesh);
    gsInfo << "Sampled at "<<mesh.n_vertices()<< " points. \n";

    gsInfo << "Getting limit points..\n";
    cc_limit_points(mesh, false); //(!)overwriting points spoils error comp.

    gsInfo << "Getting limit normals..\n";
    cc_limit_normals(mesh);

    gsInfo << "Computing errors..\n";
    error_mesh_multipatch(mesh,mp);
    //mesh.property_stats();

    if (plot)
    {
        //replace vertices by their limit positions
        mesh.swap_vertex_property("v:point","v:limit");
        gsWriteParaview(mesh,"mesh", {"v:geometry_error","v:normal_error","v:normal"});
        gsFileManager::open("mesh.vtk");
    }

    return EXIT_SUCCESS;
}


//=============================================================================