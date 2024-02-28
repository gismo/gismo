/** @file mex_uniformPointGrid.cpp

    @brief Mex adaptor for uniformPointGrid function

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Noertoft
*/

#include "mex.h"
#include <gismo.h>

using namespace std;
using namespace gismo;

typedef size_t gs_vec_sz;

// --------------------------------------------------------------------------
// "main" gateway function
void mexFunction ( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] ) 
{
    // Extract inputs
    real_t *plower = mxGetPr(prhs[0]);
    real_t *pupper = mxGetPr(prhs[1]);
    mwIndex numPoints = (mwIndex) mxGetScalar(prhs[2]);
    // Make copies of inputs

    GISMO_ENSURE(mxGetM(prhs[0])==1 && mxGetM(prhs[1])==1,"Input must have one row");
    GISMO_ENSURE(mxGetN(prhs[0])==mxGetN(prhs[1]),"Input vectors must have the same number of rows");
    size_t PARDIM = mxGetN(prhs[0]);

    gsVector<real_t> lower(PARDIM), upper(PARDIM);
    for (gs_vec_sz i=0; i<PARDIM; ++i) {
        lower(i) = *plower++;
        upper(i) = *pupper++;
    }
    // Call C++ code to generate result
    gsMatrix<real_t> pts = uniformPointGrid(lower, upper, numPoints);
    // Copy result for output
    plhs[0] = mxCreateDoubleMatrix(pts.rows(), pts.cols(), mxREAL);
    real_t *out = mxGetPr(plhs[0]);
    for (gs_vec_sz j=0; j<pts.cols(); ++j) 
    {
        for (gs_vec_sz i=0; i<pts.rows(); ++i) 
        {
            *out++ = pts(i,j);
        }
    }
}
