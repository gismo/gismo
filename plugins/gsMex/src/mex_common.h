/** @file mex_common.h

    @brief Mex common utility functions

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): O. Chanon, A. Mantzaflaris, P. Noertoft
*/

// to do : uint64_t --> uintptr_t

#pragma once

#include "mex.h"

#define CLASS_HANDLE_SIGNATURE 0xFF00F0A5

// This macro denotes the maximal length of strings passed to this mexFunction.
#define __MAXSTRLEN__ 128

// This macro denotes the parametric dimension
#define __DIM__ 2


// Copies the array pointed to by the input pointer and returns it
// as a gsMatrix.
template<class T>
gismo::gsAsConstMatrix<T> extractMatrixFromPointer(const mxArray *pnt) 
{
    return gismo::gsAsConstMatrix<real_t>(mxGetDoubles(pnt), mxGetM(pnt), mxGetN(pnt) );
}

// Copies the content of the input gsMatrix into a MATLAB real
// double matrix and returns a pointer to it. Note: memory is
// allocated (intended for output so it is not freed).
template<class T>
mxArray* createPointerFromMatrix(const gismo::gsMatrix<T> & mat) 
{
    const mwSize numRows = mat.rows(), numCols = mat.cols();
    mxArray *pnt = mxCreateDoubleMatrix(numRows,numCols,mxREAL);
    gismo::gsAsMatrix<real_t>(mxGetDoubles(pnt),numRows,numCols) =
        mat.template cast<real_t>();
    return pnt;
}

mxArray * createPointerFromStdVector(const std::vector<double>& v)
{
    mxArray * mx = mxCreateDoubleMatrix(1,v.size(), mxREAL);
    std::copy(v.begin(), v.end(), mxGetDoubles(mx));
    return mx;
}

template<class base> inline mxArray *convertPtr2Mat(base *ptr)
{
    mexLock();
    mxArray *out = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
    *(mxGetUint64s(out)) = reinterpret_cast<mxUint64>(ptr);
    return out;
}

template<class base> inline base *convertMat2HandlePtr(const mxArray *in)
{
    if (mxGetNumberOfElements(in) != 1 || mxGetClassID(in) != mxUINT64_CLASS || mxIsComplex(in))
        mexErrMsgTxt("mex_common.h (convertMat2HandlePtr): Input must be a real uint64 scalar.");
    base *ptr = reinterpret_cast<base *>(*(mxGetUint64s(in)));
    //if (!ptr->isValid())
    //    mexErrMsgTxt("Handle not valid.");
    return ptr;
}

template<class base> inline base *convertMat2Ptr(const mxArray *in)
{
    return convertMat2HandlePtr<base>(in);
}

template<class base> inline void destroyObject(const mxArray *in)
{
    delete convertMat2HandlePtr<base>(in);
    mexUnlock();
}

class mystream : public std::streambuf
{
protected:
    virtual std::streamsize xsputn(const char *s, std::streamsize n) { mexPrintf("%.*s", n, s); return n; }
    virtual int overflow(int c=EOF) { if (c != EOF) { mexPrintf("%.1s", &c); } return 1; }
};

class scoped_redirect_cout
{
public:
    scoped_redirect_cout() { old_buf = std::cout.rdbuf(); std::cout.rdbuf(&mout); }
    ~scoped_redirect_cout() { std::cout.rdbuf(old_buf); }
private:
    mystream mout;
    std::streambuf *old_buf;
};
static scoped_redirect_cout mycout_redirect;
