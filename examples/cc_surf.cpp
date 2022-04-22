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
        mesh.cc_subdivide();

    gsMultiPatch<> mp;
    if (!fn_patch.empty())
    {   //read patches from file
        gsReadFile<>(fn_patch,mp);
        gsInfo << "Input: "<<mp.nPatches()<<" spline patches.\n";
    }
    else //if geometry is not given then use ACC3
    {
        gsInfo << "Extracting "<<mesh.n_faces()<<" Bezier patches (ACC3)..\n";
        mp = mesh.cc_acc3();
    }

    for( index_t i = 0; i<s; ++i)
        mesh.cc_subdivide();
    gsInfo << "Sampled at "<<mesh.n_vertices()<< " points. \n";

    gsInfo << "Getting limit points..\n";
    mesh.cc_limit_points("v:limit", false); //(!)overwriting points spoils error comp.

    gsInfo << "Getting limit normals..\n";
    mesh.cc_limit_normals("v:normal");

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
