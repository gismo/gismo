/** @file gsPatchFromBoundary.cpp

    @brief Constructs a patch (surface, volume,..) from a set of
    boundaries.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gismo.h>
#include <iostream>
#include <string>
#include <unistd.h>

using std::cout; using std::cin;
using std::endl; using std::string;


using namespace gismo;



template <typename T = real_t>
std::vector<gsMultiPatch<>> loadSplineFromObj(std::string fn)
{
    std::vector<gsMultiPatch<>> result(15);
    std::vector<gsVector3d<T>> verts;
    std::vector<int> coefNumbers;
    std::vector<T> knotsU;
    short_t degree;
    short_t patchID;

    std::ifstream in("../filedata/" + fn, std::ios::in);
    if(!in)
    {
        gsWarn << "Could not open file.\n";
        exit(1);
    }
    std::string line;
    while(std::getline(in, line))
    {
        gsInfo << line << "\n";
        //Get Object
        if(line.substr(0,2) == "o ")
        {
            verts.clear();
            coefNumbers.clear();
            knotsU.clear();

            std::istringstream v(line.substr(7));
            v >> patchID;
        }

        //Get vertices
        if(line.substr(0,2) == "v ")
        {
            T x,y,z;
            std::istringstream v(line.substr(2));
            v >> x;
            v >> y;
            v >> z;
            verts.emplace_back(x,y,z);
        }
        if(line.substr(0,4) == "deg ")
        {
            std::istringstream v(line.substr(4));
            v >> degree;
        }
        if(line.substr(0,7) == "parm u ")
        {
            std::istringstream v(line.substr(7));
            T a;
            while(v >> a)
            {
                knotsU.push_back(a);
            }
        }

        if(line.substr(0,3) == "end")
        {
            // Reparameterize to [0,1]
            const double lengthU = knotsU.back() - knotsU.front();
            //const double lengthV = knotsV.back() - knotsV.front();
            for(auto & knot : knotsU)
            {
                knot = knot / lengthU;
            }

            //for (index_t i = 0; i < knotsU.size(); i++)
            //    gsInfo << knotsU[i] << " ";

            // Make Basis
            gsKnotVector<T> kv0(knotsU, degree);

            gsMatrix<T> coefs(verts.size(), 3);
            for(typename std::vector<gsVector3d<T>>::size_type i = 0; i < verts.size(); i++)
            {
                coefs.row(i) = verts[i].transpose();
            }

            gsMatrix<T> coefs2d = 0.02*coefs.block(0,0,verts.size(),2);

            gsBSplineBasis<T> basis(kv0);
            gsBSpline<T> curve(basis, coefs2d);

            result[patchID-1].addPatch(curve);
        }
    }
    return result;
}





int main(int argc, char* argv[])
{
    bool save      = false;
    bool plot      = false;
    index_t method = 0;
    real_t tol = 1e-4;
    //std::string fn = "rhino/test2.obj";
    std::string fn = "rhino/turtle_15patches_g4.obj";


    // Read input from command line arguments
    gsCmdLine cmd("Constructs a patch given a domain boundary.");
    cmd.addPlainString("filename", "File containing boundary data", fn);
    cmd.addInt("m","method" ,"Method: 0 Coons' patch (default), 1 Spring patch, 2: Cross-Ap. patch", method);
    cmd.addReal  ("t","tolerance","Tolerance for identifing patch interfaces", tol);
    cmd.addSwitch("save", "Save result in XML format", save);
    cmd.addSwitch("plot", "Paraview", plot);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
/*
    // Load XML file
    gsMultiPatch<> boundary;
    gsReadFile<>(fn, boundary);
    GISMO_ENSURE(!boundary.empty(), "The gsMultiPatch is empty - maybe file is missing or corrupt.");
    gsInfo<<"Got "<< boundary <<"\n";
    boundary.computeTopology(tol);
    GISMO_ENSURE( boundary.isClosed(), "The boundary is not closed, adjust tolerance.");
    boundary.closeGaps(tol);
*/

    std::vector<gsMultiPatch<>> boundary = loadSplineFromObj(fn);
    for (size_t np = 0; np < boundary.size(); np++)
    {
        boundary[np].computeTopology(tol);
        GISMO_ENSURE( boundary[np].isClosed(), "The boundary is not closed, adjust tolerance.");
        boundary[np].closeGaps(tol);

    }

    switch (method)
    {
        case 1:
        {
            /*
            gsInfo<<"Using spring patch construction.\n";
            gsSpringPatch<real_t> spring(boundary);
            gsInfo<<"Created a " << spring.compute() <<"\n";
            if (save) gsWrite(spring.result(), "result_patch");
             */
            break;
        }
        case 2:
        {
            /*
            gsInfo<<"Using cross approximation construction.\n";
            gsCrossApPatch<real_t> cross(boundary);
            gsInfo<<"Created a " << cross.compute() <<"\n";
            if (save) gsWrite(cross.result(), "result_patch");
             */
            break;
        }
        case 0:
        default:
            gsMultiPatch<> multiPatch;
            gsInfo<<"Using Coons' patch construction.\n";
            for (size_t np = 0; np < boundary.size(); np++)
            {
                gsCoonsPatch<real_t> coons(boundary[np]);
                gsInfo<<"Created a " << coons.compute() <<"\n";
                multiPatch.addPatch(coons.result());

            }
            multiPatch.computeTopology();
            if (save) gsWrite(multiPatch, "result_multiPatch");
            if (plot) gsWriteParaview(multiPatch, "result_multiPatch");
            break;
    }

    if (save)
        gsInfo << "Result saved to result_multiPatch.xml\n";
    else
        gsInfo << "Done. No output created, re-run with --save to get xml "
                  "file containing the data.\n";
    return 0;
}
