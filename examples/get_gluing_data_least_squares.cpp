/** @file get_gluing_data_least_squares.cpp

    @brief Compute gluing data via least-squares fitting of linear
           alpha/beta functions to the AS-G1 condition.

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

    std::string geometry("domain2d/two_bilinear_patches.xml");
   
    gsCmdLine cmd("Compute gluing data via least-squares.");
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

    const index_t samplingPoints = 20;
    
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
        
        // TODO: D, Compute gluing data via least squares
        //
        // We seek linear gluing functions:
        //   alphaL(t) = aL0*(1-t) + aL1*t
        //   alphaR(t) = aR0*(1-t) + aR1*t
        //   beta(t)   = b0*(1-t)  + b1*t
        //
        // such that the AS-G1 condition
        //   alphaR(t) * dFL/dn - alphaL(t) * dFR/dn + beta(t) * dFL/dt = 0
        // is satisfied as well as possible (in the least-squares sense)
        // at all sampling points along the interface.
        //
        // Unknowns vector: x = [aL0, aL1, aR0, aR1, b0, b1]^T  (6 unknowns)
        // At each sampling point i with parameter t_i, the condition gives
        // 2 scalar equations (x- and y-components).
        //
        // For a given point i, let:
        //   n1 = dFL/dn,  n2 = dFR/dn,  tau = dFL/dt   (each a 2-vector)
        //   phi0 = (1 - t_i),  phi1 = t_i
        //
        // The condition becomes (per component k = 0,1):
        //   (aR0*phi0 + aR1*phi1) * n1[k]
        // - (aL0*phi0 + aL1*phi1) * n2[k]
        // + (b0*phi0  + b1*phi1)  * tau[k]  = 0
        //
        // Row for component k at point i in the matrix M:
        //   col 0 (aL0): -phi0 * n2[k]
        //   col 1 (aL1): -phi1 * n2[k]
        //   col 2 (aR0):  phi0 * n1[k]
        //   col 3 (aR1):  phi1 * n1[k]
        //   col 4 (b0):   phi0 * tau[k]
        //   col 5 (b1):   phi1 * tau[k]
        
        const index_t normDir1 = firstPS.direction();
        const index_t ambDir1  = 1 - normDir1;
        const index_t normDir2 = secondPS.direction();
        
        auto getPartial = [](const gsMatrix<T>& derivs, index_t dir, index_t idx) -> gsVector<T> {
            gsVector<T> v(2);
            v(0) = derivs(2*dir,     idx);
            v(1) = derivs(2*dir + 1, idx);
            return v;
        };
        
        // Build the (2*N) x 6 matrix M and the zero right-hand side
        const index_t N = samplingPoints;
        gsMatrix<T> M(2*N, 6);
        M.setZero();
        
        for (index_t i = 0; i < N; ++i)
        {
            T t_i  = static_cast<T>(i) / (N - 1.);
            T phi0 = 1. - t_i;
            T phi1 = t_i;
            
            gsVector<T> n1  = getPartial(derivatives1, normDir1, i);
            gsVector<T> n2  = getPartial(derivatives2, normDir2, i);
            gsVector<T> tau = getPartial(derivatives1, ambDir1,  i);
            
            for (index_t k = 0; k < 2; ++k)
            {
                index_t row = 2*i + k;
                M(row, 0) = -phi0 * n2(k);   // aL0
                M(row, 1) = -phi1 * n2(k);   // aL1
                M(row, 2) =  phi0 * n1(k);   // aR0
                M(row, 3) =  phi1 * n1(k);   // aR1
                M(row, 4) =  phi0 * tau(k);  // b0
                M(row, 5) =  phi1 * tau(k);  // b1
            }
        }
        
        // We want M * x = 0 with x != 0.
        // This is a homogeneous least-squares problem: find the right
        // singular vector of M corresponding to the smallest singular value.
        //
        // Equivalently, find the eigenvector of M^T M with smallest eigenvalue.
        
        gsMatrix<T> MtM = M.transpose() * M;  // 6x6
        typename gsMatrix<T>::SelfAdjEigenSolver eigSolver(MtM);
        
        gsVector<T> eigenvalues = eigSolver.eigenvalues();
        gsMatrix<T> eigenvectors = eigSolver.eigenvectors();
        
        gsInfo << "Eigenvalues of M^T M: " << eigenvalues.transpose() << "\n";
        
        // The solution is the eigenvector with the smallest eigenvalue
        gsVector<T> x = eigenvectors.col(0);  // sorted ascending by Eigen
        
        // Extract gluing data
        gsMatrix<T> alpha1(1, 2), alpha2(1, 2), beta(1, 2);
        alpha1(0, 0) = x(0);  // aL0
        alpha1(0, 1) = x(1);  // aL1
        alpha2(0, 0) = x(2);  // aR0
        alpha2(0, 1) = x(3);  // aR1
        beta(0, 0)   = x(4);  // b0
        beta(0, 1)   = x(5);  // b1
        
        // Normalize so that alpha1 has a nice scale (e.g., alpha1(0) = det of left Jacobian at corner 0)
        {
            gsVector<T> dFLdn_c0 = getPartial(derivatives1, normDir1, 0);
            gsVector<T> dFLdt_c0 = getPartial(derivatives1, ambDir1,  0);
            T detJL0 = dFLdn_c0(0)*dFLdt_c0(1) - dFLdn_c0(1)*dFLdt_c0(0);
            
            if (std::abs(alpha1(0, 0)) > 1e-14)
            {
                T scale = detJL0 / alpha1(0, 0);
                alpha1 *= scale;
                alpha2 *= scale;
                beta   *= scale;
            }
        }
        
        gsInfo << "\nalpha1 (alphaL): " << alpha1 << "\n";
        gsInfo << "alpha2 (alphaR): " << alpha2 << "\n";
        gsInfo << "beta:            " << beta   << "\n";
        
        // --- Check AS-G1 condition at all sampling points ---
        gsInfo << "\nChecking AS-G1 condition at sampling points:\n";
        T maxResNorm = 0;
        for (index_t i = 0; i < N; ++i)
        {
            T t = static_cast<T>(i) / (N - 1.);
            
            T aL = alpha1(0, 0) * (1. - t) + alpha1(0, 1) * t;
            T aR = alpha2(0, 0) * (1. - t) + alpha2(0, 1) * t;
            T bn = beta(0, 0)   * (1. - t) + beta(0, 1)   * t;
            
            gsVector<T> dFLdn_i = getPartial(derivatives1, normDir1, i);
            gsVector<T> dFLdt_i = getPartial(derivatives1, ambDir1,  i);
            gsVector<T> dFRdn_i = getPartial(derivatives2, normDir2, i);
            
            gsVector<T> residual = aR * dFLdn_i - aL * dFRdn_i + bn * dFLdt_i;
            
            maxResNorm = std::max(maxResNorm, residual.norm());
            
            if (N <= 20)
                gsInfo << "  Point " << i << " (t=" << t << "): residual = ["
                       << residual.transpose() << "], norm = " << residual.norm() << "\n";
        }
        gsInfo << "  Max residual norm = " << maxResNorm << "\n";
    
    }
    
    
    
    
    return 0;
    
}

