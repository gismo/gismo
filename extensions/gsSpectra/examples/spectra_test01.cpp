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
    gsMatrix<real_t> A(10,10);
    A.setRandom();
    A = (A + A.transpose()) * 0.5;

    gsSpectraSymSolver<gsMatrix<>> minev(A, 5, 10);
    minev.compute(Spectra::SortRule::SmallestAlge);
    gsInfo << "Eigenvalues:" << minev.eigenvalues().transpose() <<"\n";
    return EXIT_SUCCESS;

}//main
