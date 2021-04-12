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
    index_t sz = 10;
    gsMatrix<real_t> Ad(sz,sz);
    Ad.setRandom();
    Ad = (Ad + Ad.transpose());

    gsSparseMatrix<real_t> A, B(sz,sz);
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

    gsInfo<<"Test for eigenvalue solvers.\n A is a symmetric matrix and B is positive definite\n";

    gsSpectraSymSolver<gsMatrix<real_t>> minev(Ad, math::floor(sz/2), sz);
    minev.compute(Spectra::SortRule::SmallestAlge);
    gsInfo << "Symmetric solver:\n";
    gsInfo << "Eigenvalues A*x=lambda*x:\n" << minev.eigenvalues().transpose() <<"\n\n";

    gsSpectraGenSymSolver<gsSparseMatrix<real_t>,Spectra::GEigsMode::Cholesky> Chsolver(A,B,math::floor(sz/2),sz);
    Chsolver.compute(Spectra::SortRule::SmallestAlge,1000,1e-3);
    gsInfo << "General Symmetric solver, Cholesky:\n";
    gsInfo << "Eigenvalues A*x=lambda*B*x:\n" << Chsolver.eigenvalues().transpose() <<"\n\n";

    gsSpectraGenSymSolver<gsSparseMatrix<real_t>,Spectra::GEigsMode::RegularInverse> Rsolver(A,B,math::floor(sz/2),sz);
    Rsolver.compute(Spectra::SortRule::SmallestAlge,1000,1e-3);
    gsInfo << "General Symmetric solver, Regular Inverse:\n";
    gsInfo << "Eigenvalues A*x=lambda*B*x:\n" << Rsolver.eigenvalues().transpose() <<"\n\n";

    gsSpectraGenSymShiftSolver<gsSparseMatrix<real_t>,Spectra::GEigsMode::ShiftInvert> Ssolver(A,B,math::floor(sz/2),sz,1);
    Ssolver.compute(Spectra::SortRule::SmallestAlge,1000,1e-3);
    gsInfo << "General Symmetric Shift solver, Shift Invert:\n";
    gsInfo << "Eigenvalues A*x=lambda*B*x (shift=1):\n" << Ssolver.eigenvalues().transpose() <<"\n\n";

    gsSpectraGenSymShiftSolver<gsSparseMatrix<real_t>,Spectra::GEigsMode::Cayley> Csolver(A,B,math::floor(sz/2),sz,1);
    Csolver.compute(Spectra::SortRule::SmallestAlge,1000,1e-3);
    gsInfo << "General Symmetric Shift solver, Cayley:\n";
    gsInfo << "Eigenvalues A*x=lambda*B*x (shift=1):\n" << Csolver.eigenvalues().transpose() <<"\n\n";

    gsSpectraGenSymShiftSolver<gsSparseMatrix<real_t>,Spectra::GEigsMode::Buckling> Bsolver(B,A,math::floor(sz/2),sz,1);
    Bsolver.compute(Spectra::SortRule::SmallestAlge,1000,1e-3);
    gsInfo << "General Symmetric Shift solver, Buckling:\n";
    gsInfo << "Eigenvalues B*x=lambda*A*x (!) (shift=1):\n" << Bsolver.eigenvalues().transpose() <<"\n\n";

    return EXIT_SUCCESS;

}//main
