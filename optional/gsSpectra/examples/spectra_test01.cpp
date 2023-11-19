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

    index_t Nhalf = math::floor(sz/2);
    index_t N = std::min(10,Nhalf);

    gsMatrix<> Ad = A.toDense();
    gsSpectraSymSolver<gsMatrix<real_t>> minev(Ad, N, 2*N);
    minev.compute(Spectra::SortRule::SmallestAlge,1000,1e-10,Spectra::SortRule::SmallestAlge);
    gsInfo << "Symmetric solver ";
    if (minev.converged())
        gsInfo<<"(converged in "<<minev.num_iterations()<<"iterations)."<<"\n";
    else
        gsInfo<<"(did not converge)\n";
    gsInfo << "Eigenvalues A*x=lambda*x:\n" << minev.eigenvalues().transpose() <<"\n\n";

    gsSpectraSymShiftSolver<gsMatrix<real_t>> minev2(Ad, N, 2*N,0.0);
    minev2.compute(Spectra::SortRule::SmallestAlge,1000,1e-10,Spectra::SortRule::SmallestAlge);
    gsInfo << "Symmetric solver";
    if (minev2.converged())
        gsInfo<<"(converged in "<<minev2.num_iterations()<<"iterations)."<<"\n";
    else
        gsInfo<<"(did not converge)\n";

    gsInfo << "Eigenvalues A*x=lambda*x:\n" << minev2.eigenvalues().transpose() <<"\n\n";

    gsSpectraSymShiftSolver<gsSparseMatrix<real_t>> minev3(A, N, 2*N,0.0);
    minev3.compute(Spectra::SortRule::SmallestAlge,1000,1e-10,Spectra::SortRule::SmallestAlge);
    gsInfo << "Symmetric solver";
    if (minev3.converged())
        gsInfo<<"(converged in "<<minev3.num_iterations()<<"iterations)."<<"\n";
    else
        gsInfo<<"(did not converge)\n";

    gsInfo << "Eigenvalues A*x=lambda*x:\n" << minev3.eigenvalues().transpose() <<"\n\n";

    gsSpectraGenSymSolver<gsSparseMatrix<real_t>,Spectra::GEigsMode::Cholesky> Chsolver(A,B,N,2*N);
    Chsolver.compute(Spectra::SortRule::SmallestAlge,1000,1e-3,Spectra::SortRule::SmallestMagn);
    gsInfo << "General Symmetric solver, Cholesky: ";
    if (Chsolver.converged())
        gsInfo<<"(converged in "<<Chsolver.num_iterations()<<"iterations)."<<"\n";
    else
        gsInfo<<"(did not converge)\n";

    gsInfo << "Eigenvalues A*x=lambda*B*x:\n" << Chsolver.eigenvalues().transpose() <<"\n\n";

    gsSpectraGenSymSolver<gsSparseMatrix<real_t>,Spectra::GEigsMode::RegularInverse> Rsolver(A,B,N,2*N);
    Rsolver.compute(Spectra::SortRule::SmallestAlge,1000,1e-3,Spectra::SortRule::SmallestMagn);
    gsInfo << "General Symmetric solver, Regular Inverse: ";
    if (Rsolver.converged())
        gsInfo<<"(converged in "<<Rsolver.num_iterations()<<"iterations)."<<"\n";
    else
        gsInfo<<"(did not converge)\n";

    gsInfo << "Eigenvalues A*x=lambda*B*x:\n" << Rsolver.eigenvalues().transpose() <<"\n\n";

    gsSpectraGenSymShiftSolver<gsSparseMatrix<real_t>,Spectra::GEigsMode::ShiftInvert> Ssolver(A,B,N,2*N,1e-3);
    Ssolver.compute(Spectra::SortRule::LargestMagn,1000,1e-10,Spectra::SortRule::SmallestMagn);
    gsInfo << "General Symmetric Shift solver, Shift Invert (LM sorting): ";
    if (Ssolver.converged())
        gsInfo<<"(converged in "<<Ssolver.num_iterations()<<"iterations)."<<"\n";
    else
        gsInfo<<"(did not converge)\n";

    gsInfo << "Eigenvalues A*x=lambda*B*x (shift=1):\n" << Ssolver.eigenvalues().transpose() <<"\n\n";
    Ssolver.compute(Spectra::SortRule::LargestAlge,1000,1e-10,Spectra::SortRule::SmallestMagn);
    gsInfo << "General Symmetric Shift solver, Shift Invert (LA sorting): ";
    if (Ssolver.converged())
        gsInfo<<"(converged in "<<Ssolver.num_iterations()<<"iterations)."<<"\n";
    else
        gsInfo<<"(did not converge)\n";

    gsInfo << "Eigenvalues A*x=lambda*B*x (shift=1):\n" << Ssolver.eigenvalues().transpose() <<"\n\n";

    Ssolver.compute(Spectra::SortRule::SmallestMagn,1000,1e-10,Spectra::SortRule::SmallestMagn);
    gsInfo << "General Symmetric Shift solver, Shift Invert (SM sorting): ";
    if (Ssolver.converged())
        gsInfo<<"(converged in "<<Ssolver.num_iterations()<<"iterations)."<<"\n";
    else
        gsInfo<<"(did not converge)\n";

    gsInfo << "Eigenvalues A*x=lambda*B*x (shift=1):\n" << Ssolver.eigenvalues().transpose() <<"\n\n";
    Ssolver.compute(Spectra::SortRule::SmallestAlge,1000,1e-10,Spectra::SortRule::SmallestMagn);
    gsInfo << "General Symmetric Shift solver, Shift Invert (SA sorting): ";
    if (Ssolver.converged())
        gsInfo<<"(converged in "<<Ssolver.num_iterations()<<"iterations)."<<"\n";
    else
        gsInfo<<"(did not converge)\n";

    gsInfo << "Eigenvalues A*x=lambda*B*x (shift=1):\n" << Ssolver.eigenvalues().transpose() <<"\n\n";
    Ssolver.compute(Spectra::SortRule::BothEnds,1000,1e-10,Spectra::SortRule::SmallestMagn);
    gsInfo << "General Symmetric Shift solver, Shift Invert (BE sorting): ";
    if (Ssolver.converged())
        gsInfo<<"(converged in "<<Ssolver.num_iterations()<<"iterations)."<<"\n";
    else
        gsInfo<<"(did not converge)\n";

    gsInfo << "Eigenvalues A*x=lambda*B*x (shift=1):\n" << Ssolver.eigenvalues().transpose() <<"\n\n";

    gsSpectraGenSymShiftSolver<gsSparseMatrix<real_t>,Spectra::GEigsMode::Cayley> Csolver(A,B,N,2*N,1e-3);
    Csolver.compute(Spectra::SortRule::LargestAlge,1000,1e-3,Spectra::SortRule::SmallestMagn);
    gsInfo << "General Symmetric Shift solver, Cayley (LAsorting): ";
    if (Csolver.converged())
            gsInfo<<"(converged in "<<Csolver.num_iterations()<<"iterations)."<<"\n";
        else
            gsInfo<<"(did not converge)\n";

    gsInfo << "Eigenvalues A*x=lambda*B*x (shift=1):\n" << Csolver.eigenvalues().transpose() <<"\n\n";
    Csolver.compute(Spectra::SortRule::LargestMagn,1000,1e-3,Spectra::SortRule::SmallestMagn);
    gsInfo << "General Symmetric Shift solver, Cayley (LMsorting): ";
    if (Csolver.converged())
            gsInfo<<"(converged in "<<Csolver.num_iterations()<<"iterations)."<<"\n";
        else
            gsInfo<<"(did not converge)\n";

    gsInfo << "Eigenvalues A*x=lambda*B*x (shift=1):\n" << Csolver.eigenvalues().transpose() <<"\n\n";
    Csolver.compute(Spectra::SortRule::SmallestAlge,1000,1e-3,Spectra::SortRule::SmallestMagn);
    gsInfo << "General Symmetric Shift solver, Cayley (SAsorting): ";
    if (Csolver.converged())
            gsInfo<<"(converged in "<<Csolver.num_iterations()<<"iterations)."<<"\n";
        else
            gsInfo<<"(did not converge)\n";

    gsInfo << "Eigenvalues A*x=lambda*B*x (shift=1):\n" << Csolver.eigenvalues().transpose() <<"\n\n";
    Csolver.compute(Spectra::SortRule::SmallestMagn,1000,1e-3,Spectra::SortRule::SmallestMagn);
    gsInfo << "General Symmetric Shift solver, Cayley (SMsorting): ";
    if (Csolver.converged())
            gsInfo<<"(converged in "<<Csolver.num_iterations()<<"iterations)."<<"\n";
        else
            gsInfo<<"(did not converge)\n";

    gsInfo << "Eigenvalues A*x=lambda*B*x (shift=1):\n" << Csolver.eigenvalues().transpose() <<"\n\n";
    Csolver.compute(Spectra::SortRule::BothEnds,1000,1e-3,Spectra::SortRule::SmallestMagn);
    gsInfo << "General Symmetric Shift solver, Cayley (BEsorting): ";
    if (Csolver.converged())
            gsInfo<<"(converged in "<<Csolver.num_iterations()<<"iterations)."<<"\n";
        else
            gsInfo<<"(did not converge)\n";

    gsInfo << "Eigenvalues A*x=lambda*B*x (shift=1):\n" << Csolver.eigenvalues().transpose() <<"\n\n";

    gsSpectraGenSymShiftSolver<gsSparseMatrix<real_t>,Spectra::GEigsMode::Buckling> Bsolver(B,A,N,2*N,1);
    Bsolver.compute(Spectra::SortRule::SmallestAlge,1000,1e-3,Spectra::SortRule::SmallestMagn);
    gsInfo << "General Symmetric Shift solver, Buckling: ";
    gsInfo << "Eigenvalues B*x=lambda*A*x (!) (shift=1):\n" << Bsolver.eigenvalues().transpose() <<"\n\n";

    return EXIT_SUCCESS;

}//main
