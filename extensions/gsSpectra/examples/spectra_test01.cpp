/** @file

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
*/

#include <gismo.h>
#include <gsSpectra/gsSpectra.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    std::string matrixA, matrixB;

    gsCmdLine cmd("Spectra tester! If you want to provide matrices from file, use options -A and -B");
    cmd.addString( "A", "matA", "Input XML file for matrix A", matrixA );
    cmd.addString( "B", "matB", "Input XML file for matrix B", matrixB );

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsSparseMatrix<real_t> A, B;
    index_t sz;
    if (matrixA.empty() || matrixB.empty())
    {
        sz = 10;
        gsMatrix<real_t> Ad(sz,sz);
        Ad.setRandom();
        Ad = (Ad + Ad.transpose());

        B.resize(sz,sz);
        A = Ad.sparseView();

        // Define the B matrix, a tridiagonal matrix with 2 on the diagonal
        // and 1 on the subdiagonals
        for (index_t i = 0; i < sz; i++)
        {
            B.insert(i, i) = 2.0;
            if (i > 0)
                B.insert(i - 1, i) = 1.0;
            if (i < sz - 1)
                B.insert(i + 1, i) = 1.0;
        }

    }
    else
    {
        gsReadFile<>(matrixA,A);
        gsReadFile<>(matrixB,B);
        sz = A.cols();
    }

    gsInfo<<"Test for eigenvalue solvers.\n A is a symmetric matrix and B is positive definite\n";



    gsMatrix<> Ad = A.toDense();
    // Eigen::GeneralizedSelfAdjointEigenSolver< typename gsMatrix<>::Base >  eigSolverA;
    // eigSolverA.compute(Ad);
    // gsInfo << "Eigen dense solver:\n";
    // gsInfo << "Eigenvalues A*x=lambda*B*x:\n" << eigSolverA.eigenvalues().transpose() <<"\n\n";

    gsSpectraSymSolver<gsMatrix<real_t>> minev(Ad, math::floor(sz/2), sz);
    minev.compute(Spectra::SortRule::SmallestAlge,1000,1e-10,Spectra::SortRule::SmallestAlge);
    gsInfo << "Symmetric solver:\n";
    gsInfo << "Eigenvalues A*x=lambda*x:\n" << minev.eigenvalues().transpose() <<"\n\n";

    gsSpectraSymShiftSolver<gsMatrix<real_t>> minev2(Ad, math::floor(sz/2), sz,0.0);
    minev2.compute(Spectra::SortRule::SmallestAlge,1000,1e-10,Spectra::SortRule::SmallestAlge);
    gsInfo << "Symmetric solver:\n";
    gsInfo << "Eigenvalues A*x=lambda*x:\n" << minev2.eigenvalues().transpose() <<"\n\n";

    gsSpectraSymShiftSolver<gsSparseMatrix<real_t>> minev3(A, math::floor(sz/2), sz,0.0);
    minev3.compute(Spectra::SortRule::SmallestAlge,1000,1e-10,Spectra::SortRule::SmallestAlge);
    gsInfo << "Symmetric solver:\n";
    gsInfo << "Eigenvalues A*x=lambda*x:\n" << minev3.eigenvalues().transpose() <<"\n\n";

    // gsMatrix<> Bd = B.toDense();
    // Eigen::GeneralizedSelfAdjointEigenSolver< typename gsMatrix<>::Base >  eigSolverAB;
    // eigSolverAB.compute(Ad,Bd);
    // gsInfo << "Eigen dense solver:\n";
    // gsInfo << "Eigenvalues A*x=lambda*B*x:\n" << eigSolverAB.eigenvalues().transpose() <<"\n\n";

    gsSpectraGenSymSolver<gsSparseMatrix<real_t>,Spectra::GEigsMode::Cholesky> Chsolver(A,B,math::floor(sz/2),sz);
    Chsolver.compute(Spectra::SortRule::SmallestAlge,1000,1e-3,Spectra::SortRule::SmallestMagn);
    gsInfo << "General Symmetric solver, Cholesky:\n";
    gsInfo << "Eigenvalues A*x=lambda*B*x:\n" << Chsolver.eigenvalues().transpose() <<"\n\n";

    gsSpectraGenSymSolver<gsSparseMatrix<real_t>,Spectra::GEigsMode::RegularInverse> Rsolver(A,B,math::floor(sz/2),sz);
    Rsolver.compute(Spectra::SortRule::SmallestAlge,1000,1e-3,Spectra::SortRule::SmallestMagn);
    gsInfo << "General Symmetric solver, Regular Inverse:\n";
    gsInfo << "Eigenvalues A*x=lambda*B*x:\n" << Rsolver.eigenvalues().transpose() <<"\n\n";

    gsSpectraGenSymShiftSolver<gsSparseMatrix<real_t>,Spectra::GEigsMode::ShiftInvert> Ssolver(A,B,2,sz,1e-3);
    Ssolver.compute(Spectra::SortRule::LargestMagn,1000,1e-10,Spectra::SortRule::SmallestMagn);
    gsInfo << "General Symmetric Shift solver, Shift Invert (LM sorting):\n";
    gsInfo << "Eigenvalues A*x=lambda*B*x (shift=1):\n" << Ssolver.eigenvalues().transpose() <<"\n\n";
    Ssolver.compute(Spectra::SortRule::LargestAlge,1000,1e-10,Spectra::SortRule::SmallestMagn);
    gsInfo << "General Symmetric Shift solver, Shift Invert (LA sorting):\n";
    gsInfo << "Eigenvalues A*x=lambda*B*x (shift=1):\n" << Ssolver.eigenvalues().transpose() <<"\n\n";
    Ssolver.compute(Spectra::SortRule::SmallestMagn,1000,1e-10,Spectra::SortRule::SmallestMagn);
    gsInfo << "General Symmetric Shift solver, Shift Invert (SM sorting):\n";
    gsInfo << "Eigenvalues A*x=lambda*B*x (shift=1):\n" << Ssolver.eigenvalues().transpose() <<"\n\n";
    Ssolver.compute(Spectra::SortRule::SmallestAlge,1000,1e-10,Spectra::SortRule::SmallestMagn);
    gsInfo << "General Symmetric Shift solver, Shift Invert (SA sorting):\n";
    gsInfo << "Eigenvalues A*x=lambda*B*x (shift=1):\n" << Ssolver.eigenvalues().transpose() <<"\n\n";
    Ssolver.compute(Spectra::SortRule::BothEnds,1000,1e-10,Spectra::SortRule::SmallestMagn);
    gsInfo << "General Symmetric Shift solver, Shift Invert (BE sorting):\n";
    gsInfo << "Eigenvalues A*x=lambda*B*x (shift=1):\n" << Ssolver.eigenvalues().transpose() <<"\n\n";

    gsSpectraGenSymShiftSolver<gsSparseMatrix<real_t>,Spectra::GEigsMode::Cayley> Csolver(A,B,math::floor(sz/2),sz,1e-3);
    Csolver.compute(Spectra::SortRule::LargestAlge,1000,1e-3,Spectra::SortRule::SmallestMagn);
    gsInfo << "General Symmetric Shift solver, Cayley (LAsorting):\n";
    gsInfo << "Eigenvalues A*x=lambda*B*x (shift=1):\n" << Csolver.eigenvalues().transpose() <<"\n\n";
    Csolver.compute(Spectra::SortRule::LargestMagn,1000,1e-3,Spectra::SortRule::SmallestMagn);
    gsInfo << "General Symmetric Shift solver, Cayley (LMsorting):\n";
    gsInfo << "Eigenvalues A*x=lambda*B*x (shift=1):\n" << Csolver.eigenvalues().transpose() <<"\n\n";
    Csolver.compute(Spectra::SortRule::SmallestAlge,1000,1e-3,Spectra::SortRule::SmallestMagn);
    gsInfo << "General Symmetric Shift solver, Cayley (SAsorting):\n";
    gsInfo << "Eigenvalues A*x=lambda*B*x (shift=1):\n" << Csolver.eigenvalues().transpose() <<"\n\n";
    Csolver.compute(Spectra::SortRule::SmallestMagn,1000,1e-3,Spectra::SortRule::SmallestMagn);
    gsInfo << "General Symmetric Shift solver, Cayley (SMsorting):\n";
    gsInfo << "Eigenvalues A*x=lambda*B*x (shift=1):\n" << Csolver.eigenvalues().transpose() <<"\n\n";
    Csolver.compute(Spectra::SortRule::BothEnds,1000,1e-3,Spectra::SortRule::SmallestMagn);
    gsInfo << "General Symmetric Shift solver, Cayley (BEsorting):\n";
    gsInfo << "Eigenvalues A*x=lambda*B*x (shift=1):\n" << Csolver.eigenvalues().transpose() <<"\n\n";

    gsSpectraGenSymShiftSolver<gsSparseMatrix<real_t>,Spectra::GEigsMode::Buckling> Bsolver(B,A,math::floor(sz/2),sz,1);
    Bsolver.compute(Spectra::SortRule::SmallestAlge,1000,1e-3,Spectra::SortRule::SmallestMagn);
    gsInfo << "General Symmetric Shift solver, Buckling:\n";
    gsInfo << "Eigenvalues B*x=lambda*A*x (!) (shift=1):\n" << Bsolver.eigenvalues().transpose() <<"\n\n";

    return EXIT_SUCCESS;

}//main
