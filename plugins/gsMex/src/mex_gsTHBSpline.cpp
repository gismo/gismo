/** @file mex_gsTHBSpline.cpp

    @brief Mex adaptor for gsTHBSpline class

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): O. Chanon, A. Mantzaflaris
*/

#include <gismo.h>
#include "mex_common.h"
#include <cstring>

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

            if (nrhs==3) {
                // constructor from 1 argument (+ 1 type switch)
                char constructSwitch[__MAXSTRLEN__];
                if (mxGetString(prhs[1], constructSwitch, sizeof(constructSwitch)))
                    throw("Second input argument should be a string"
                          "less than MAXSTRLEN characters long.");
                if (!strcmp(constructSwitch,"char")) {
                    // constructor ( char filename )
                    // Extract the input
                    char* input_buf = mxArrayToString(prhs[2]);
                    // Read the THB-spline basis from the specified file
                    std::string filename(input_buf); // Reading requires a std::string
                    gsFileData<real_t>  data( filename );
                    gsTHBSpline<__DIM__> * hbs = data.getFirst< gsTHBSpline<__DIM__> >().release();
                    plhs[0] = convertPtr2Mat<gsTHBSpline<__DIM__> >(hbs);
                    // Free the memory allocated by mxArrayToString
                    mxFree(input_buf);
                } else if (!strcmp(constructSwitch,"gsTensorBSpline")) {
                    // constructor ( gsTensorBSplineBasis )
                    gsTensorBSpline<__DIM__> *instance = convertMat2Ptr<gsTensorBSpline<__DIM__> >(prhs[2]);
                    plhs[0] = convertPtr2Mat<gsTHBSpline<__DIM__> >(new gsTHBSpline<__DIM__>(*instance));
                } else {
                    throw ("Invalid construction.");
                }
            } else {
                // INVALID constructor
                throw("Invalid construction.");
            }

        } else if (!strcmp(cmd,"destructor")) {

            // ----------------------------------------------------------------------
            // Destructor
            destroyObject<gsTHBSpline<__DIM__> >(prhs[1]);

        } else if (!strcmp(cmd,"accessor")) {

            // ----------------------------------------------------------------------
            // Accessor

            // Fetch instance and property to be accessed
            gsTHBSpline<__DIM__> *instance = convertMat2Ptr<gsTHBSpline<__DIM__> >(prhs[1]);
            char prop[__MAXSTRLEN__];
            if (mxGetString(prhs[2], prop, sizeof(prop)))
                throw("Third input argument should be a property string less than MAXSTRLEN characters long.");

            // Call method as specified by the input string
            if (!strcmp(prop,"dim")) {
                int val      = instance->parDim();
                mxArray *out = mxCreateDoubleScalar((double)val);
                plhs[0]      = out;
            } else if (!strcmp(prop,"size")) {
                int val      = instance->size();
                mxArray *out = mxCreateDoubleScalar((double)val);
                plhs[0]      = out;
            } else if (!strcmp(prop,"support")) {
                gsMatrix<real_t> supp = instance->support();
                plhs[0] = createPointerFromMatrix(supp);
            } else {
                // Unknown property
                throw("Third input argument contains an unknown property string.");
            }

        } else if (!strcmp(cmd,"eval")) {

            // ----------------------------------------------------------------------
            // eval(pts)

            gsTHBSpline<__DIM__> *instance = convertMat2Ptr<gsTHBSpline<__DIM__> >(prhs[1]);
            // Copy the input (FIXME: this should be avoided)
            const gsMatrix<real_t> pts = extractMatrixFromPointer<real_t>(prhs[2]);
            // Call the method
            gsMatrix<real_t> vals = instance->eval(pts);
            // Copy result to output (FIXME: this should be avoided)
            plhs[0] = createPointerFromMatrix<real_t>(vals);

        } else if (!strcmp(cmd,"deriv")) {

            // ----------------------------------------------------------------------
            // deriv(pts)

            gsTHBSpline<__DIM__> *instance = convertMat2Ptr<gsTHBSpline<__DIM__> >(prhs[1]);
            // Copy the input (FIXME: this should be avoided)
            const gsMatrix<real_t> pts = extractMatrixFromPointer<real_t>(prhs[2]);
            // Call the method
            gsMatrix<real_t> vals = instance->jacobian(pts);
            // Copy result to output (FIXME: this should be avoided)
            plhs[0] = createPointerFromMatrix<real_t>(vals);

        } else if (!strcmp(cmd,"hess")) {

            // ----------------------------------------------------------------------
            // hess(pts, dir)

            gsTHBSpline<__DIM__> *instance = convertMat2Ptr<gsTHBSpline<__DIM__> >(prhs[1]);
            // Copy the input (FIXME: this should be avoided)
            const gsMatrix<real_t> pts = extractMatrixFromPointer<real_t>(prhs[2]);
            mwIndex dir = (mwIndex) mxGetScalar(prhs[3]);
            // Call the method
            gsMatrix<real_t> vals = instance->hess(pts, dir-1);
            // Copy result to output (FIXME: this should be avoided)
            plhs[0] = createPointerFromMatrix<real_t>(vals);

        } else if (!strcmp(cmd,"active")) {

            // ----------------------------------------------------------------------
            // active(pts)

            gsTHBSpline<__DIM__> *instance = convertMat2Ptr<gsTHBSpline<__DIM__> >(prhs[1]);
            // Copy the input (FIXME: this should be avoided)
            gsMatrix<real_t> pts = extractMatrixFromPointer<real_t>(prhs[2]);
            // Call method
            const gsMatrix<unsigned> vals = instance->active(pts);
            // Copy the result for output (FIXME: this should be avoided)
            plhs[0] = createPointerFromMatrix<unsigned>(vals);

        } else if (!strcmp(cmd,"basis")) {

            // ----------------------------------------------------------------------
            // basis()

            gsTHBSpline<__DIM__> *instance = convertMat2Ptr<gsTHBSpline<__DIM__> >(prhs[1]);
            // Copy the result for output (FIXME: this should be avoided)
            gsTHBSplineBasis<__DIM__> * hbs = new gsTHBSplineBasis<__DIM__>(instance->basis());
            plhs[0] = convertPtr2Mat<gsTHBSplineBasis<__DIM__> >(hbs);

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
        std::string errMsg = std::string("\n  While executing the following command: ") + cmd + 
            std::string("\n  The following exception/error ocurred: ") + e.what();
        mexErrMsgTxt(errMsg.c_str());
    } catch (const char* str) { // Caught from within this mexFunction.
        std::string errMsg = std::string("\n  While executing the following command: ") + cmd + 
            std::string("\n  The following exception/error ocurred: ") + std::string(str);
        mexErrMsgTxt(errMsg.c_str());
    } catch (...) { // Something else went wrong
        std::string errMsg = std::string("\n  While executing the following command: ") + cmd + 
            std::string("\n  An error ocurred.");
        mexErrMsgTxt(errMsg.c_str());
    } // end try-catch

} // end mexFunction
