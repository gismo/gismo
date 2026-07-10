#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <gismo.h>
#include <gsModeling/gsAsG1GlobalBasis.hpp>

using namespace gismo;

template <typename T>
bool testFile(const std::string& filename, int refinements, int degree, T tolVal, T tolGrad)
{
    gsInfo << "Testing file: " << filename << "\n";

    // Read geometry
    typename gsMultiPatch<T>::uPtr mpPtr;
    try {
        mpPtr = gsReadFile<>(filename);
    } catch (...) {
        gsInfo << "  Error: Cannot read file.\n";
        return false;
    }
    if (!mpPtr) {
        gsInfo << "  Error: Null patch pointer.\n";
        return false;
    }
    gsMultiPatch<T>& mp = *mpPtr;
    mp.computeTopology();

    // Ensure degree >= targetDegree
    const short_t inputDeg = mp.patch(0).basis().degree(0);
    if (inputDeg < degree)
    {
        const short_t elev = degree - inputDeg;
        mp.degreeElevate(elev);
    }

    // Refine
    if (refinements > 0)
    {
        const short_t deg = mp.patch(0).basis().degree(0);
        const index_t mult = std::max<index_t>(deg - 1, 1);
        for (index_t i = 0; i < refinements; ++i)
            mp.uniformRefine(1, mult);
    }

    // Construct MultiBasis and GlobalBasis
    gsMultiBasis<T> mb(mp);
    gsAsG1GlobalBasis<T> globalBasis(mp, mb);
    gsSparseMatrix<T> m_M = globalBasis.getTransformationMatrix();

    // Export Matrix
    {
        std::string outName = filename + ".mtx";
        std::ofstream file(outName);
        if (file.is_open())
        {
            file << "%%MatrixMarket matrix coordinate real general\n";
            file << m_M.rows() << " " << m_M.cols() << " " << m_M.nonZeros() << "\n";
            for (int col = 0; col < m_M.cols(); ++col)
            {
                for (typename gsSparseMatrix<T>::InnerIterator it(m_M, col); it; ++it)
                {
                    file << (it.row() + 1) << " " << (it.col() + 1) << " " << it.value() << "\n";
                }
            }
            file.close();
            gsInfo << "  Exported matrix to " << outName << "\n";
        }
    }

    // Setup offsets
    int nPatches = mp.nPatches();
    std::vector<int> offsets(nPatches + 1, 0);
    for (int p = 0; p < nPatches; ++p)
        offsets[p+1] = offsets[p] + mb.basis(p).size();

    int nGlobal = m_M.cols();
    gsInfo << "  Global DOFs: " << nGlobal << "\n";
    T maxMM = m_M.coeffs().cwiseAbs().maxCoeff();
    gsInfo << "  Max m_M coeff: " << maxMM << "\n";
    if (maxMM > 1e10) {
        for (int i = 0; i < nGlobal; ++i) {
            gsVector<T> v = gsVector<T>::Zero(nGlobal);
            v(i) = 1.0;
            gsVector<T> c = m_M * v;
            T m = c.cwiseAbs().maxCoeff();
            if (m > 1e10) gsInfo << "    Huge coef on Global DOF " << i << ": " << m << "\n";
        }
    }

    // Numerical G1 smoothness verification along all interfaces
    const index_t nCheck = 21;
    T maxValErr = 0;
    T maxGradErr = 0;

    for (auto it = mp.iBegin(); it != mp.iEnd(); ++it)
    {
        const boundaryInterface& ifc = *it;
        const patchSide ps1 = ifc.first();
        const patchSide ps2 = ifc.second();

        const short_t normDir1 = ps1.direction(), tanDir1_ = 1 - normDir1;
        const short_t normDir2 = ps2.direction(), tanDir2_ = 1 - normDir2;
        const bool par1_ = ps1.parameter(), par2_ = ps2.parameter();
        const bool ifcFlipped = !ifc.dirOrientation(ps1, tanDir1_);

        gsMatrix<T> sup1 = mp.patch(ps1.patch).support();
        gsMatrix<T> sup2 = mp.patch(ps2.patch).support();
        const T ifcCoord1 = sup1(normDir1, par1_ ? 1 : 0);
        const T ifcCoord2 = sup2(normDir2, par2_ ? 1 : 0);
        const T t1a = sup1(tanDir1_, 0), t1b = sup1(tanDir1_, 1);
        const T t2a = sup2(tanDir2_, 0), t2b = sup2(tanDir2_, 1);

        for (index_t idx = 0; idx < nGlobal; ++idx)
        {
            gsVector<T> globalVec = gsVector<T>::Zero(nGlobal);
            globalVec(idx) = T(1);
            gsVector<T> coefs = m_M * globalVec;

            gsVector<T> c1 = coefs.segment(offsets[ps1.patch], mb.basis(ps1.patch).size());
            gsVector<T> c2 = coefs.segment(offsets[ps2.patch], mb.basis(ps2.patch).size());

            auto func1 = static_cast<const gsTensorBSplineBasis<2,T>&>(mb.basis(ps1.patch)).makeGeometry(c1);
            auto func2 = static_cast<const gsTensorBSplineBasis<2,T>&>(mb.basis(ps2.patch)).makeGeometry(c2);

            for (index_t i = 0; i < nCheck; ++i)
            {
                T s = T(i) / T(nCheck - 1);
                T t1 = t1a + s * (t1b - t1a);
                T s2 = ifcFlipped ? (1.0 - s) : s;
                T t2 = t2a + s2 * (t2b - t2a);

                gsMatrix<T> pt1(2, 1), pt2(2, 1);
                pt1(normDir1, 0) = ifcCoord1; pt1(tanDir1_, 0) = t1;
                pt2(normDir2, 0) = ifcCoord2; pt2(tanDir2_, 0) = t2;

                // Check C0
                gsMatrix<T> v1 = func1->eval(pt1);
                gsMatrix<T> v2 = func2->eval(pt2);
                T valErr = std::abs(v1(0, 0) - v2(0, 0));
                maxValErr = std::max(maxValErr, valErr);

                // Check G1
                gsMatrix<T> df1, df2, dG1, dG2;
                func1->deriv_into(pt1, df1);
                func2->deriv_into(pt2, df2);
                mp.patch(ps1.patch).deriv_into(pt1, dG1);
                mp.patch(ps2.patch).deriv_into(pt2, dG2);

                gsMatrix<T> J1(2,2), J2(2,2);
                J1(0,0) = dG1(0,0); J1(0,1) = dG1(1,0);
                J1(1,0) = dG1(2,0); J1(1,1) = dG1(3,0);
                J2(0,0) = dG2(0,0); J2(0,1) = dG2(1,0);
                J2(1,0) = dG2(2,0); J2(1,1) = dG2(3,0);

                gsVector<T> paramGrad1(2), paramGrad2(2);
                paramGrad1(0) = df1(0,0); paramGrad1(1) = df1(1,0);
                paramGrad2(0) = df2(0,0); paramGrad2(1) = df2(1,0);

                gsVector<T> physGrad1 = J1.inverse().transpose() * paramGrad1;
                gsVector<T> physGrad2 = J2.inverse().transpose() * paramGrad2;

                T gradErr = (physGrad1 - physGrad2).norm();
                maxGradErr = std::max(maxGradErr, gradErr);
            }
        }
    }

    gsInfo << "  Max C0 error: " << maxValErr << "\n";
    gsInfo << "  Max G1 error: " << maxGradErr << "\n";

    if (maxValErr < tolVal && maxGradErr < tolGrad)
    {
        gsInfo << "  Result: PASS\n\n";
        return true;
    }
    else
    {
        gsInfo << "  Result: FAIL\n\n";
        return false;
    }
}

#include <dirent.h>
#include <algorithm>

int main(int argc, char* argv[])
{
    using T = real_t;

    std::string folder("/home/fhasanova/gismo/filedata/domain2d/2patch");
    index_t refinements = 1;
    index_t degree = 3;
    T tolVal = 1e-8;
    T tolGrad = 1e-3;

    gsCmdLine cmd("Test gsAsG1GlobalBasis G1 continuity on a directory.");
    cmd.addString("d", "dir", "Directory containing XML files", folder);
    cmd.addInt("r", "refinements", "Number of refinements", refinements);
    cmd.addInt("p", "degree", "Target polynomial degree", degree);

    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    DIR *dir;
    struct dirent *ent;
    std::vector<std::string> xml_files;
    
    if ((dir = opendir(folder.c_str())) != NULL) {
        while ((ent = readdir(dir)) != NULL) {
            std::string fname = ent->d_name;
            if (fname.size() > 4 && fname.substr(fname.size()-4) == ".xml") {
                xml_files.push_back(fname);
            }
        }
        closedir(dir);
    } else {
        gsInfo << "Could not open directory " << folder << "\n";
        return 1;
    }
    
    std::sort(xml_files.begin(), xml_files.end());
    
    int passed = 0;
    int failed = 0;
    
    for (size_t i = 0; i < xml_files.size(); ++i) {
        std::string filepath = folder + "/" + xml_files[i];
        bool success = false;
        try {
            success = testFile<T>(filepath, refinements, degree, tolVal, tolGrad);
        } catch (const std::exception& e) {
            gsInfo << "  Exception: " << e.what() << "\n";
            success = false;
        } catch (...) {
            gsInfo << "  Unknown exception.\n";
            success = false;
        }
        
        if (success) passed++;
        else failed++;
    }
    
    gsInfo << "Summary: " << passed << " passed, " << failed << " failed.\n";
    return failed == 0 ? 0 : 1;
}
