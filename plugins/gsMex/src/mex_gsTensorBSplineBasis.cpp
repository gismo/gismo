/** @file mex_gsTensorBSplineBasis.cpp

    @brief Mex adaptor for gsTensorBSplineBasis class

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): O. Chanon, A. Mantzaflaris, P. Noertoft
*/

#include <gismo.h>
#include "mex_common.h"

#define T real_t
#define Z index_t

using namespace gismo;

// --------------------------------------------------------------------------
// "main" gateway function
//
// plhs -> Matlab output
// prhs -> Matlab input
//
void mexFunction ( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    // --------------------------------------------------------------------------
    // Try command.
    char cmd[__MAXSTRLEN__];
    try  {

        // Fetch command string.
        if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
            throw("First input argument should be a command string"
                  "less than MAXSTRLEN characters long.");

        if (!strcmp(cmd,"constructor")) {

            // ----------------------------------------------------------------------
            // Constructors

            if (nrhs==3)
            {
                if (mxIsDouble(prhs[1]) && mxIsDouble(prhs[2])) // two knot vectors are directly passed from Matlab
                {
                    // argument 2: list of knots
                    if (mxIsScalar(prhs[1])) { // or not scalar
                        mexErrMsgIdAndTxt("MATLAB:mexcpp:typeargin", "Second argument has to be double vector.");
                    }
                    if (!mxGetM(prhs[1])==1) { // or not scalar
                        mexErrMsgIdAndTxt("MATLAB:mexcpp:typeargin", "Second argument must be a row vector.");
                    }

                    // argument 3: list of knots
                    if (mxIsScalar(prhs[2])) { // or not scalar
                        mexErrMsgIdAndTxt("MATLAB:mexcpp:typeargin", "Second argument has to be double vector.");
                    }
                    if (!mxGetM(prhs[2])==1) { // or not scalar
                        mexErrMsgIdAndTxt("MATLAB:mexcpp:typeargin", "Second argument must be a row vector.");
                    }

                    T *knots1 = mxGetDoubles(prhs[1]);
                    T *knots2 = mxGetDoubles(prhs[2]);

                    gsKnotVector<T> kv1(*knots1);
                    gsKnotVector<T> kv2(*knots2);
                    gsTensorBSplineBasis<2,T> * tbasis = new gsTensorBSplineBasis<2,T>(kv1,kv2);
                    plhs[0] = convertPtr2Mat<gsTensorBSplineBasis<2,T> >(tbasis);
                }
                else if (mxIsChar(prhs[1]) && mxIsChar(prhs[2])) // two knot vectors are passed as G+Smo types
                {
                    char constructSwitch1[__MAXSTRLEN__];
                    char constructSwitch2[__MAXSTRLEN__];
                    if (mxGetString(prhs[1], constructSwitch1, sizeof(constructSwitch1)))
                        throw("Second input argument should be a string"
                              "less than MAXSTRLEN characters long.");
                    if (mxGetString(prhs[2], constructSwitch2, sizeof(constructSwitch2)))
                        throw("Third input argument should be a string"
                              "less than MAXSTRLEN characters long.");
                    if (!strcmp(constructSwitch1,"gsKnotVector") && strcmp(constructSwitch2,"gsKnotVector"))
                    {
                        gsKnotVector<T> *kv1 = convertMat2Ptr<gsKnotVector<T> >(prhs[3]);
                        gsKnotVector<T> *kv2 = convertMat2Ptr<gsKnotVector<T> >(prhs[4]);
                        gsTensorBSplineBasis<2,T> * tbasis = new gsTensorBSplineBasis<2,T>(*kv1,*kv2);
                        plhs[0] = convertPtr2Mat<gsTensorBSplineBasis<2,T> >(tbasis);
                    }
                }

            }
            else
            {
                // INVALID constructor
                throw("Invalid construction.");
            }

        } else if (!strcmp(cmd,"destructor")) {

            // ----------------------------------------------------------------------
            // Destructor
            destroyObject<gsTensorBSplineBasis<__DIM__> >(prhs[1]);

        } else if (!strcmp(cmd,"degree")) {

            gsTensorBSplineBasis <__DIM__> *instance = convertMat2Ptr < gsTensorBSplineBasis < __DIM__ > > (prhs[1]);
            mwIndex deg = (mwIndex) mxGetScalar(prhs[2]);
            mxArray *out = mxCreateDoubleScalar((double)instance->degree(deg));
            plhs[0] = out;

        } else if (!strcmp(cmd,"eval")) {

            // ----------------------------------------------------------------------
            // eval(pts)

            gsTensorBSplineBasis<__DIM__> *instance = convertMat2Ptr<gsTensorBSplineBasis<__DIM__> >(prhs[1]);
            // Copy the input (FIXME: this should be avoided)
            const gsMatrix<real_t> pts = extractMatrixFromPointer<real_t>(prhs[2]);
            // Call the method
            gsMatrix<real_t> vals = instance->eval(pts);
            // Copy result to output (FIXME: this should be avoided)
            plhs[0] = createPointerFromMatrix<real_t>(vals);

        } else if (!strcmp(cmd,"evalSingle")) {

            // ----------------------------------------------------------------------
            // evalSingle(ind,pts)

            gsTensorBSplineBasis <__DIM__> *instance = convertMat2Ptr < gsTensorBSplineBasis < __DIM__ > > (prhs[1]);
            // Copy the input (FIXME: this should be avoided)
            const mwIndex ind = (mwIndex) * mxGetPr(prhs[2]);
            const gsMatrix <real_t> pts = extractMatrixFromPointer<real_t>(prhs[3]);
            // Call the method
            gsMatrix <real_t> vals = instance->evalSingle(ind-1, pts);
            // Copy result to output (FIXME: this should be avoided)
            plhs[0] = createPointerFromMatrix<real_t>(vals);

        } else if (!strcmp(cmd,"save")) {

            // ----------------------------------------------------------------------
            // save(file)

            gsTensorBSplineBasis <__DIM__> *instance = convertMat2Ptr < gsTensorBSplineBasis < __DIM__ > > (prhs[1]);
            char* input_buf = mxArrayToString(prhs[2]);
            // Save the THB-spline basis in the specified file
            std::string filename(input_buf); // Reading requires a std::string
            gsWrite(*instance, filename);

        } else if (!strcmp(cmd,"knots")) {

            // ----------------------------------------------------------------------
            // knots(level, direction)

            gsTensorBSplineBasis <__DIM__> *instance = convertMat2Ptr < gsTensorBSplineBasis < __DIM__ > > (prhs[1]);
            mwIndex direction = (mwIndex) mxGetScalar(prhs[2]);
            const gsKnotVector<>& kv = instance->knots(direction-1);
            plhs[0] = createPointerFromStdVector(kv);

        } else if (!strcmp(cmd,"active")) {

            // ----------------------------------------------------------------------
            // active(pts)

            gsTensorBSplineBasis<__DIM__> *instance = convertMat2Ptr<gsTensorBSplineBasis<__DIM__> >(prhs[1]);
            // Copy the input (FIXME: this should be avoided)
            gsMatrix<real_t> pts = extractMatrixFromPointer<real_t>(prhs[2]);
            // Call method
            const gsMatrix<unsigned> vals = instance->active(pts);
            // Copy the result for output (FIXME: this should be avoided)
            plhs[0] = createPointerFromMatrix<unsigned>(vals);

        } else {

            // ----------------------------------------------------------------------
            // Unknown command
            throw("unknown command.");

        }

        // ------------------------------------------------------------------------
        // That's it (command executed).
        return;

        // --------------------------------------------------------------------------
        // Catch (command failed).

    } catch (std::exception& e) { // Caught from library.
        std::string errMsg = std::string("\n  While executing the following command: ") + cmd + std::string("\n  The following exception/error ocurred: ") + e.what();
        mexErrMsgTxt(errMsg.c_str());
    } catch (const char* str) { // Caught from within this mexFunction.
        std::string errMsg = std::string("\n  While executing the following command: ") + cmd + std::string("\n  The following exception/error ocurred: ") + std::string(str);
        mexErrMsgTxt(errMsg.c_str());
    } catch (...) { // Something else went wrong
        std::string errMsg = std::string("\n  While executing the following command: ") + cmd + std:: string("\n  An error ocurred.");
        mexErrMsgTxt(errMsg.c_str());
    } // end try-catch

} // end mexFunction
