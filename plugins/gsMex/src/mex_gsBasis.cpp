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

        gsBasis<real_t> * instance = convertMat2Ptr<gsBasis<real_t> >(prhs[1]);
        if      (!strcmp(cmd,"dim"))
        {
            int val      = instance->dim();
            mxArray *out = mxCreateDoubleScalar((double)val);
            plhs[0]      = out;
        }
        else if (!strcmp(cmd,"numElements"))
        {
            int val      = instance->numElements();
            mxArray *out = mxCreateDoubleScalar((double)val);
            plhs[0]      = out;
        }
        else if (!strcmp(cmd,"size"))
        {
            int val      = instance->size();
            mxArray *out = mxCreateDoubleScalar((double)val);
            plhs[0]      = out;
        }
        else if (!strcmp(cmd,"support"))
        {
            gsMatrix<real_t> supp = instance->support();
            plhs[0] = createPointerFromMatrix(supp);
        }
        else if (!strcmp(cmd,"degree"))
        {
            mwIndex deg = (mwIndex) mxGetScalar(prhs[2]);
            mxArray *out = mxCreateDoubleScalar((double)instance->degree(deg));
            plhs[0] = out;
        }
        else if (!strcmp(cmd,"eval"))
        {
            // Copy the input (FIXME: this should be avoided)
            const gsMatrix<real_t> pts = extractMatrixFromPointer<real_t>(prhs[2]);
            // Call the method
            gsMatrix<real_t> vals = instance->eval(pts);
            // Copy result to output (FIXME: this should be avoided)
            plhs[0] = createPointerFromMatrix<real_t>(vals);
        }
        else if (!strcmp(cmd,"evalSingle"))
        {
            // Copy the input (FIXME: this should be avoided)
            const mwIndex ind = (mwIndex) * mxGetDoubles(prhs[2]);
            const gsMatrix <real_t> pts = extractMatrixFromPointer<real_t>(prhs[3]);
            // Call the method
            gsMatrix <real_t> vals = instance->evalSingle(ind-1, pts);
            // Copy result to output (FIXME: this should be avoided)
            plhs[0] = createPointerFromMatrix<real_t>(vals);
        }
        else if (!strcmp(cmd,"active"))
        {
            // Copy the input (FIXME: this should be avoided)
            gsMatrix<real_t> pts = extractMatrixFromPointer<real_t>(prhs[2]);
            // Call method
            const gsMatrix<unsigned> vals = instance->active(pts);
            // Copy the result for output (FIXME: this should be avoided)
            plhs[0] = createPointerFromMatrix<unsigned>(vals);
        } else if (!strcmp(cmd,"print")) {
            std::ostringstream a;
            instance->print(a);
            mexPrintf("%s\n", a.str().c_str());
        }
        else {

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
