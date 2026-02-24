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
        
        const T tol = 1e-10;
        if ((values1 - values2).norm() < tol)
        {
            gsInfo << "Boundary values match directly.\n";
        }
        else if ((values1 - values2.rowwise().reverse()).norm() < tol)
        {
            gsInfo << "Boundary values match after flipping. Reversing values2 and derivatives2.\n";
            values2 = values2.rowwise().reverse().eval();
            derivatives2 = derivatives2.rowwise().reverse().eval();
        }
        else
        {
            gsInfo << "Boundary values do NOT match. Something is wrong!\n";
            return -1;
        }
        
     
        // TODO: C. Check if vaules1 == values2 *or* if they are the same after flipping the columns
        //          In the second case, flip also the columns of derivatives2
        //DONE ABOVE
        
        // TODO: D, Compute gluing data
        // ...
        
        
        // It is always the gluing data for the two corners
        /*gsMatrix<T> alpha1(1,2); // fillme
        gsMatrix<T> alpha2(1,2); // fillme
        gsMatrix<T> beta1(1,2);  // fillme
        gsMatrix<T> beta2(1,2);  // fillme*/

        
        // derivatives1 and derivatives2 are 4 x samplingPoints matrices.
        // For a 2D geometry, deriv_into gives:
        //   row 0: dF/du (x-component)
        //   row 1: dF/du (y-component)
        //   row 2: dF/dv (x-component)
        //   row 3: dF/dv (y-component)
        //
        // The normal direction for each patch side is the "u" direction here
        // (the direction crossing the interface), and "v" is along the interface.
        // We need to extract dF/du and dF/dv properly based on normDir and ambDir.
        
        const index_t normDir1 = firstPS.direction();
        const index_t ambDir1  = 1 - normDir1;
        const index_t normDir2 = secondPS.direction();
        const index_t ambDir2  = 1 - normDir2;
        
        // For deriv_into, the layout is:
        //   row 0: dF/d(xi_0), x-component
        //   row 1: dF/d(xi_0), y-component
        //   row 2: dF/d(xi_1), x-component
        //   row 3: dF/d(xi_1), y-component
        // So dF/d(xi_k) has x-comp at row 2*k and y-comp at row 2*k+1.
        
        // Helper lambda: extract dF/d(xi_dir) at corner idx from a derivatives matrix
        auto getPartial = [](const gsMatrix<T>& derivs, index_t dir, index_t idx) -> gsVector<T> {
            gsVector<T> v(2);
            v(0) = derivs(2*dir,     idx);
            v(1) = derivs(2*dir + 1, idx);
            return v;
        };
        
        // Helper lambda: 2D determinant det([a, b])
        auto det2D = [](const gsVector<T>& a, const gsVector<T>& b) -> T {
            return a(0)*b(1) - a(1)*b(0);
        };
        
        // Corner indices: first column (idx=0) and last column (idx=samplingPoints-1)
        const index_t c0 = 0;
        const index_t c1 = samplingPoints - 1;
        
        // --- Patch 1 (Left) ---
        // dFL/du = derivative in normal direction, dFL/dv = derivative in ambient direction
        gsVector<T> dFLdu_c0 = getPartial(derivatives1, normDir1, c0);
        gsVector<T> dFLdv_c0 = getPartial(derivatives1, ambDir1,  c0);
        gsVector<T> dFLdu_c1 = getPartial(derivatives1, normDir1, c1);
        gsVector<T> dFLdv_c1 = getPartial(derivatives1, ambDir1,  c1);
        
        // alphaL = det(dFL/du, dFL/dv) at corners, linearly interpolated
        // alpha1 stores [value_at_corner0, value_at_corner1]
        gsMatrix<T> alpha1(1, 2);
        alpha1(0, 0) = det2D(dFLdu_c0, dFLdv_c0);
        alpha1(0, 1) = det2D(dFLdu_c1, dFLdv_c1);
        
        // betaL = dot(dFL/du, dFL/dv) / dot(dFL/dv, dFL/dv) at corners
        gsMatrix<T> beta1(1, 2);
        beta1(0, 0) = dFLdu_c0.dot(dFLdv_c0) / dFLdv_c0.dot(dFLdv_c0);
        beta1(0, 1) = dFLdu_c1.dot(dFLdv_c1) / dFLdv_c1.dot(dFLdv_c1);
        
        // --- Patch 2 (Right) ---
        gsVector<T> dFRdu_c0 = getPartial(derivatives2, normDir2, c0);
        gsVector<T> dFRdv_c0 = getPartial(derivatives2, ambDir2,  c0);
        gsVector<T> dFRdu_c1 = getPartial(derivatives2, normDir2, c1);
        gsVector<T> dFRdv_c1 = getPartial(derivatives2, ambDir2,  c1);
        
        gsMatrix<T> alpha2(1, 2);
        alpha2(0, 0) = det2D(dFRdu_c0, dFRdv_c0);
        alpha2(0, 1) = det2D(dFRdu_c1, dFRdv_c1);
        
        gsMatrix<T> beta2(1, 2);
        beta2(0, 0) = dFRdu_c0.dot(dFRdv_c0) / dFRdv_c0.dot(dFRdv_c0);
        beta2(0, 1) = dFRdu_c1.dot(dFRdv_c1) / dFRdv_c1.dot(dFRdv_c1);
        
        gsInfo << "alpha1 (alphaL): " << alpha1 << "\n";
        gsInfo << "alpha2 (alphaR): " << alpha2 << "\n";
        gsInfo << "beta1  (betaL):  " << beta1  << "\n";
        gsInfo << "beta2  (betaR):  " << beta2  << "\n";
        
        // betanew = alphaL * betaR - alphaR * betaL (at each corner)
        gsMatrix<T> betanew(1, 2);
        betanew(0, 0) = alpha1(0, 0) * beta2(0, 0) - alpha2(0, 0) * beta1(0, 0);
        betanew(0, 1) = alpha1(0, 1) * beta2(0, 1) - alpha2(0, 1) * beta1(0, 1);
        gsInfo << "betanew: " << betanew << "\n";
        
        // --- Check AS-G1 condition at all sampling points ---
        // The condition: alphaR * dFL/du - alphaL * dFR/du + beta * dFL/dv == 0
        // where alpha and beta are linearly interpolated along the interface.
        gsInfo << "\nChecking AS-G1 condition at sampling points:\n";
        for (index_t i = 0; i < samplingPoints; ++i)
        {
            // Linear interpolation parameter t in [0,1]
            T t = static_cast<T>(i) / (samplingPoints - 1.);
            
            T aL = alpha1(0, 0) * (1. - t) + alpha1(0, 1) * t;
            T aR = alpha2(0, 0) * (1. - t) + alpha2(0, 1) * t;
            T bn = betanew(0, 0) * (1. - t) + betanew(0, 1) * t;
            
            gsVector<T> dFLdu_i = getPartial(derivatives1, normDir1, i);
            gsVector<T> dFLdv_i = getPartial(derivatives1, ambDir1,  i);
            gsVector<T> dFRdu_i = getPartial(derivatives2, normDir2, i);
            
            gsVector<T> residual = aR * dFLdu_i + aL * dFRdu_i + bn * dFLdv_i;
            
            gsInfo << "  Point " << i << " (t=" << t << "): residual = ["
                   << residual.transpose() << "], norm = " << residual.norm() << "\n";
        }
        
        // ...existing code...
    }
    
    
    
    
    return 0;
    
}

