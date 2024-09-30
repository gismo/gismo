/** @file mex_gsTHBSplineBasis.cpp

    @brief Mex adaptor for gsTHBSplineBasis class

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): O. Chanon, A. Mantzaflaris, P. Noertoft
*/

#include <gismo.h>
#include "mex_common.h"

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
                    gsTHBSplineBasis<__DIM__> * hbs = data.getFirst< gsTHBSplineBasis<2> >().release();
                    plhs[0] = convertPtr2Mat<gsTHBSplineBasis<__DIM__> >(hbs);
                    // Free the memory allocated by mxArrayToString
                    mxFree(input_buf);
                } else if (!strcmp(constructSwitch,"gsTensorBSplineBasis")) {
                    // constructor ( gsTensorBSplineBasis )
                    gsTensorBSplineBasis<__DIM__> *instance = convertMat2Ptr<gsTensorBSplineBasis<__DIM__> >(prhs[2]);
                    plhs[0] = convertPtr2Mat<gsTHBSplineBasis<__DIM__> >(new gsTHBSplineBasis<__DIM__>(*instance));
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
            destroyObject<gsTHBSplineBasis<__DIM__> >(prhs[1]);

        } else if (!strcmp(cmd,"accessor")) {

            // ----------------------------------------------------------------------
            // Accessor

            // Fetch instance and property to be accessed
            gsTHBSplineBasis<__DIM__> *instance = convertMat2Ptr<gsTHBSplineBasis<__DIM__> >(prhs[1]);
            char prop[__MAXSTRLEN__];
            if (mxGetString(prhs[2], prop, sizeof(prop)))
                throw("Third input argument should be a property string less than MAXSTRLEN characters long.");

            // Call method as specified by the input string
            if (!strcmp(prop,"dim")) {
                int val      = instance->dim();
                mxArray *out = mxCreateDoubleScalar((double)val);
                plhs[0]      = out;
            } else if (!strcmp(prop,"numElements")) {
                int val      = instance->numElements();
                mxArray *out = mxCreateDoubleScalar((double)val);
                plhs[0]      = out;
            } else if (!strcmp(prop,"size")) {
                int val      = instance->size();
                mxArray *out = mxCreateDoubleScalar((double)val);
                plhs[0]      = out;
            } else if (!strcmp(prop,"treeSize")) {
                int val      = instance->treeSize();
                mxArray *out = mxCreateDoubleScalar((double)val);
                plhs[0]      = out;
            } else if (!strcmp(prop,"treeLeafSize")) {
                int val      = instance->tree().leafSize();
                mxArray *out = mxCreateDoubleScalar((double)val);
                plhs[0]      = out;
            } else if (!strcmp(prop,"support")) {
                gsMatrix<real_t> supp = instance->support();
                plhs[0] = createPointerFromMatrix(supp);
            } else if (!strcmp(prop,"maxLevel")) {
                int val      = instance->maxLevel()+1;
                mxArray *out = mxCreateDoubleScalar((double)val);
                plhs[0]      = out;
            }

        } else if (!strcmp(cmd,"treePrintLeaves")) {

            // ----------------------------------------------------------------------
            // tree().printLeaves()

            gsTHBSplineBasis <__DIM__> *instance = convertMat2Ptr < gsTHBSplineBasis < __DIM__ > > (prhs[1]);
            instance->tree().printLeaves();

        } else if (!strcmp(cmd,"degree")) {

            gsTHBSplineBasis <__DIM__> *instance = convertMat2Ptr < gsTHBSplineBasis < __DIM__ > > (prhs[1]);
            mwIndex deg = (mwIndex) mxGetScalar(prhs[2]);
            mxArray *out = mxCreateDoubleScalar((double)instance->degree(deg));
            plhs[0] = out;

        } else if (!strcmp(cmd,"eval")) {

            // ----------------------------------------------------------------------
            // eval(pts)

            gsTHBSplineBasis<__DIM__> *instance = convertMat2Ptr<gsTHBSplineBasis<__DIM__> >(prhs[1]);
            // Copy the input (FIXME: this should be avoided)
            const gsMatrix<real_t> pts = extractMatrixFromPointer<real_t>(prhs[2]);
            // Call the method
            gsMatrix<real_t> vals = instance->eval(pts);
            // Copy result to output (FIXME: this should be avoided)
            plhs[0] = createPointerFromMatrix<real_t>(vals);

        } else if (!strcmp(cmd,"evalSingle")) {

            // ----------------------------------------------------------------------
            // evalSingle(ind,pts)

            gsTHBSplineBasis <__DIM__> *instance = convertMat2Ptr < gsTHBSplineBasis < __DIM__ > > (prhs[1]);
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

            gsTHBSplineBasis <__DIM__> *instance = convertMat2Ptr < gsTHBSplineBasis < __DIM__ > > (prhs[1]);
            char* input_buf = mxArrayToString(prhs[2]);
            // Save the THB-spline basis in the specified file
            std::string filename(input_buf); // Reading requires a std::string
            gsWrite(*instance, filename);

        } else if (!strcmp(cmd,"knots")) {

            // ----------------------------------------------------------------------
            // knots(level, direction)

            gsTHBSplineBasis <__DIM__> *instance = convertMat2Ptr < gsTHBSplineBasis < __DIM__ > > (prhs[1]);
            mwIndex level = (mwIndex) mxGetScalar(prhs[2]);
            mwIndex direction = (mwIndex) mxGetScalar(prhs[3]);
            const gsKnotVector<>& kv = instance->tensorLevel(level-1).knots(direction-1);
            plhs[0] = createPointerFromStdVector(kv);

        } else if (!strcmp(cmd,"active")) {

            // ----------------------------------------------------------------------
            // active(pts)

            gsTHBSplineBasis<__DIM__> *instance = convertMat2Ptr<gsTHBSplineBasis<__DIM__> >(prhs[1]);
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
