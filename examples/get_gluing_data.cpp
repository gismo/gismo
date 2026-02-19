/** @file get_gluing_data.cpp

    @brief Tutorial on gsBasis class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Hasanova, S. Takacs
*/

#include <iostream>
#include <set>
#include <map>
#include <algorithm>
#include <numeric>
#include <gismo.h>

using namespace gismo;

template<class T>
std::pair< gsMatrix<T>, gsMatrix<T> >
evalGeoAtSide(const gsGeometry<T>& geo, patchSide ps, index_t samplingPoints)
{
    gsMatrix<T> support = geo.support();
    //gsInfo << "support:\n" << support << "\n";
    
    const index_t normDir = ps.direction();
    const index_t ambDir = 1-normDir;
    const index_t par = ps.parameter();
    
    //gsInfo << "normDir:" << normDir << "\n";
    //gsInfo << "par:" << par << "\n";
    
    gsMatrix<T> points(2,samplingPoints);
    for (index_t i=0; i<samplingPoints; ++i)
    {
        points(ambDir, i) = support( ambDir, 0 ) * (1. - i/(samplingPoints-1.)) 
                          + support( ambDir, 1 ) * i / (samplingPoints-1.);
        points(normDir,i) = support( normDir, par );
    }
    //gsInfo << "points:\n" << points << "\n";
    
    std::pair< gsMatrix<T>, gsMatrix<T> > result;
    
    geo.eval_into(points, result.first);
    geo.deriv_into(points, result.second);
    
    return result;
}




int main(int argc, char* argv[])
{
    using T = double;

    std::string geometry("domain2d/simple_two_patch_geometry.xml");
   
    gsCmdLine cmd("Example for get gluing data.");
    cmd.addString("f", "file", "G+Smo input multi patch file.", geometry);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }


    
    // Read geometry
    gsMultiPatch<>::uPtr mpPtr = gsReadFile<>(geometry);
    if (!mpPtr)
    {
        gsInfo << "No geometry found in file " << geometry << ".\n";
        return -1;
    }
    gsMultiPatch<>& mp = *mpPtr;
    
    mp.computeTopology(); // if the xml file does not have a topolgoy, we need to compute it first
    
    if (mp.nPatches() != 2)
    {
        gsInfo << "The file contains " << mp.nPatches() << " patches. This is not ok!\n";
        return -1;
    }
    
    if (mp.nInterfaces() != 1)
    {
        gsInfo << "The mp contains " << mp.nInterfaces() << " interfaces.  This is not ok!\n";
        return -1;
    }

    const index_t samplingPoints = 5;
    
    for (auto it = mp.iBegin(); it<mp.iEnd(); ++it)
    {
        const boundaryInterface& interf = *it;
        gsInfo << "The one and only interfaces is:\n" << interf << "\n";
        
        patchSide firstPS = interf.first();
        patchSide secondPS = interf.second();
        
        
        // A. First patch
        const gsGeometry<>& firstGeo = mp.patch(firstPS.patch);
        gsInfo << "firstGeo is " << firstGeo << "\n";

        std::pair< gsMatrix<T>, gsMatrix<T> > result1 = evalGeoAtSide(firstGeo, firstPS, samplingPoints);
        
        gsMatrix<T>& values1 = result1.first;
        gsInfo << "values1:\n" << values1 << "\n";
    
        gsMatrix<T>& derivatives1 = result1.second;
        gsInfo << "derivatives1:\n" << derivatives1 << "\n";

        
        
        // B. Second patch
        const gsGeometry<>& secondGeo = mp.patch(secondPS.patch);
        gsInfo << "secondGeo is " << secondGeo << "\n";
        
        std::pair< gsMatrix<T>, gsMatrix<T> > result2 = evalGeoAtSide(secondGeo, secondPS, samplingPoints);
        
        gsMatrix<T>& values2 = result2.first;
        gsInfo << "values2:\n" << values2 << "\n";
    
        gsMatrix<T>& derivatives2 = result2.second;
        gsInfo << "derivatives2:\n" << derivatives2 << "\n";
        

        // TODO: C. Check if vaules1 == values2 *or* if they are the same after flipping the columns
        //          In the second case, flip also the columns of derivatives2
        
        // TODO: D, Compute gluing data
        // ...
        
        
        // It is always the gluing data for the two corners
        gsMatrix<T> alpha1(1,2); // fillme
        gsMatrix<T> alpha2(1,2); // fillme
        gsMatrix<T> beta1(1,2);  // fillme
        gsMatrix<T> beta2(1,2);  // fillme
    }
    
    
    
    
    return 0;
    
}

