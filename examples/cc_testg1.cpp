/** @file 

    @brief 

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <iostream>

#include <gismo.h>

#include <cfenv> //for float exceptions

using namespace gismo;

int main(int argc, char *argv[])
{
    // Enable floating point exceptions
    //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

    gsDebugVar( omp_get_max_threads() );
    
    std::string fn("cube");
    index_t numSamples(1000);

    //! [Parse Command line]
    gsCmdLine cmd("Hi, give me a file (eg: .xml) and I will try to draw it!");
    cmd.addPlainString("filename", "File containing data to draw (.xml or third-party)", fn);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse Command line]

    gsFileData<>  fd;
    gsMultiPatch<> mp;
    gsMatrix<> values, normals, vproj, vnormals;
    gsMesh<> ccmesh;
    fd.read(fn+".xml");
    //fd.read(fn+"_acc3.xml");
    fd.getFirst(mp);
    gsInfo<< "Got "<< mp <<"\n";
    fd.read(fn+"_values.xml");
    fd.getFirst(values);
    fd.read(fn+"_normals.xml");
    //fd.read(fn+"_t2.xml");
    fd.getFirst(normals);
    index_t np = values.cols(); //bitorus: 48

    values.conservativeResize(3,np);
    gsInfo<<std::setprecision(15)<<"limit var : "<< values.col(0).transpose()<<"\n";
    gsInfo<<std::setprecision(15)<<"bezier pt : "<< mp.patch(0).coef(30)<<"\n" ;
    normals.conservativeResize(3,np);
    vproj.conservativeResize(3,np);
    vnormals.conservativeResize(3,np);
    gsMapData<> gdata;
    gdata.addFlags( NEED_VALUE|NEED_NORMAL );
    std::pair<index_t,gsVector<> > cp;

//    /*
    for (index_t i = 0; i!= np; ++i)
    {
        gsInfo << "\r"<< i << std::flush;
        cp = mp.closestPointTo(values.col(i), 1e-15 );
        gdata.points = cp.second;
        gdata.patchId = cp.first;
        mp.patch(cp.first).computeMap(gdata);
        vproj.col(i) = gdata.eval(0);
        vnormals.col(i) = gdata.normal(0).normalized();

        //real_t dist = (values.col(i)-mp.patch(cp.first).eval(cp.second)).norm();
        //gsInfo <<"**Pid="<<cp.first<<", Dist("<<i<<"): "<< dist  <<"\n";
    }
    gsInfo << "\n";

    gsInfo<<std::setprecision(15)<<"project pt: "<< vproj.col(0).transpose()<<"\n" ;

    //*/
    //fd.read(fn+"_refined_mesh_cc.axl");
    fd.read(fn+"_refined_mesh_cc.off");
    //fd.getAnyFirst(ccmesh); //error
    gsMesh<>::uPtr msh = fd.getAnyFirst< gsMesh<> >();
    gsInfo<<std::setprecision(15)<<"mesh vert : "<< msh->vertex(0).transpose()<<"\n" ;

    gsWriteParaview( *msh, "ccmesh", normals);
    /*
    gsInfo << "NV: "<< msh->numVertices() <<"\n";
    gsInfo << "NN: "<< normals.cols() <<"\n";
    gsInfo << "m: "<< msh->vertex(0).transpose() <<"\n";
    gsInfo << "v: "<< values.col(0).transpose() <<"\n";
    gsInfo << "n: "<< normals.col(0).transpose() <<"\n";
    return 0;
    //*/

    gsMatrix<> err(1,np), nerr(1,np);
    for (index_t i = 0; i!=np; ++i)
    {
        //err(i) = (msh->vertex(i)-values.col(i)).norm(); //diff between vertex and limit pt
        err(i) = (values.col(i)-vproj.col(i)).norm(); //diff CC and G1-ACC5
        //err(i) = (normals.col(i)-vnormals.col(i)).norm(); //diff unv(CC) and unv(G1-ACC5)
        //GISMO_ENSURE( !math::isinf(err(i)), "Inf: ("<< values.col(i).transpose()<<") - ("<< vproj.col(i).transpose()<<")\n");

        // msh has the refined CC mesh
//        msh->vertex(i).setCoords( values.col(i) ); //use CC values as the plotted geometry
        //msh->vertex(i).setCoords( vproj.col(i) ); //use projection on G1 as the plotted geometry
    }

    gsWriteParaview( *msh, "ccmesh", err);

    //ccmesh = give(*msh);
    gsWriteParaviewPoints(values, "values");
    //gsWriteParaviewPoints(vproj, "vproj");
    //gsWriteParaview(mp, "patches", numSamples);
    //gsWriteParaview( ccmesh, "ccmesh");

    gsFileManager::open("ccmesh.vtk");
    return EXIT_SUCCESS;
}
