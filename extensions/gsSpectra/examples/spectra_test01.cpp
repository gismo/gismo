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
    gsMatrix<real_t> Ad(10,10), Bd(10,10);
    Ad.setRandom();
    Ad = (Ad + Ad.transpose()) * 0.5;
    Bd.setIdentity();
    Bd = (Bd + Bd.transpose()) * 0.5;

    gsSparseMatrix<real_t> A, B;
    A = Ad.sparseView();
    B = Bd.sparseView();


    gsSpectraSymSolver<gsMatrix<real_t>> minev(Ad, 5, 10);
    minev.compute(Spectra::SortRule::SmallestAlge);
    gsInfo << "Eigenvalues:" << minev.eigenvalues().transpose() <<"\n";

    gsSpectraGenSymSolver<gsSparseMatrix<real_t>,Spectra::GEigsMode::Cholesky> Chsolver(A,B,5,10);
    Chsolver.compute(Spectra::SortRule::SmallestAlge,1000,1e-3);
    gsInfo << "Eigenvalues:" << Chsolver.eigenvalues().transpose() <<"\n";

    gsSpectraGenSymSolver<gsSparseMatrix<real_t>,Spectra::GEigsMode::RegularInverse> Rsolver(A,B,5,10);
    Rsolver.compute(Spectra::SortRule::SmallestAlge,1000,1e-3);
    gsInfo << "Eigenvalues:" << Rsolver.eigenvalues().transpose() <<"\n";

    gsSpectraGenSymShiftSolver<gsSparseMatrix<real_t>,Spectra::GEigsMode::ShiftInvert> Ssolver(A,B,5,10,1);
    Ssolver.compute(Spectra::SortRule::SmallestAlge,1000,1e-3);
    gsInfo << "Eigenvalues:" << Ssolver.eigenvalues().transpose() <<"\n";

    gsSpectraGenSymShiftSolver<gsSparseMatrix<real_t>,Spectra::GEigsMode::Buckling> Bsolver(A,B,5,10,1);
    Bsolver.compute(Spectra::SortRule::SmallestAlge,1000,1e-3);
    gsInfo << "Eigenvalues:" << Bsolver.eigenvalues().transpose() <<"\n";

    gsSpectraGenSymShiftSolver<gsSparseMatrix<real_t>,Spectra::GEigsMode::Cayley> Csolver(A,B,5,10,1);
    Csolver.compute(Spectra::SortRule::SmallestAlge,1000,1e-3);
    gsInfo << "Eigenvalues:" << Csolver.eigenvalues().transpose() <<"\n";

    return EXIT_SUCCESS;

}//main
