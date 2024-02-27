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



template<int D> void
callMember(const gsTHBSplineBasis<D,double> & instance,  char * prop,
           mxArray* plhs[], const mxArray* prhs[])
{
            // Call method as specified by the input string
            if (!strcmp(prop,"dim")) {
                int val      = instance.dim();
                mxArray *out = mxCreateDoubleScalar((double)val);
                plhs[0]      = out;
            } else if (!strcmp(prop,"numElements")) {
                int val      = instance.numElements();
                mxArray *out = mxCreateDoubleScalar((double)val);
                plhs[0]      = out;
            } else if (!strcmp(prop,"size")) {
                int val      = instance.size();
                mxArray *out = mxCreateDoubleScalar((double)val);
                plhs[0]      = out;
            } else if (!strcmp(prop,"treeSize")) {
                int val      = instance.treeSize();
                mxArray *out = mxCreateDoubleScalar((double)val);
                plhs[0]      = out;
            } else if (!strcmp(prop,"treeLeafSize")) {
                int val      = instance.tree().leafSize();
                mxArray *out = mxCreateDoubleScalar((double)val);
                plhs[0]      = out;
            } else if (!strcmp(prop,"support")) {
                gsMatrix<real_t> supp = instance.support();
                plhs[0] = createPointerFromMatrix(supp);
            } else if (!strcmp(prop,"maxLevel")) {
                int val      = instance.maxLevel()+1;
                mxArray *out = mxCreateDoubleScalar((double)val);
                plhs[0]      = out;
            }
            else if (!strcmp(prop,"treePrintLeaves"))
            {
            // ----------------------------------------------------------------------
            // tree().printLeaves()
                instance.tree().printLeaves();
            }
            else if (!strcmp(prop,"degree"))
            {
                mwIndex deg = (mwIndex) mxGetScalar(prhs[2]);
                mxArray *out = mxCreateDoubleScalar((double)instance.degree(deg));
                plhs[0] = out;
            }
            else if (!strcmp(prop,"eval"))
            {
                // Copy the input (FIXME: this should be avoided)
                const gsMatrix<real_t> pts = extractMatrixFromPointer<real_t>(prhs[2]);
                // Call the method
                gsMatrix<real_t> vals = instance.eval(pts);
                // Copy result to output (FIXME: this should be avoided)
                plhs[0] = createPointerFromMatrix<real_t>(vals);
            }
            else if (!strcmp(prop,"evalSingle"))
            {
                // Copy the input (FIXME: this should be avoided)
                const mwIndex ind = (mwIndex) * mxGetDoubles(prhs[2]);
                const gsMatrix <real_t> pts = extractMatrixFromPointer<real_t>(prhs[3]);
                // Call the method
                gsMatrix <real_t> vals = instance.evalSingle(ind-1, pts);
                // Copy result to output (FIXME: this should be avoided)
                plhs[0] = createPointerFromMatrix<real_t>(vals);
            }
            else if (!strcmp(prop,"save"))
            {
                char* input_buf = mxArrayToString(prhs[2]);
                // Save the THB-spline basis in the specified file
                std::string filename(input_buf); // Reading requires a std::string
                gsWrite(instance, filename);
            }
            else if (!strcmp(prop,"knots"))
            {
                mwIndex level = (mwIndex) mxGetScalar(prhs[2]);
                mwIndex direction = (mwIndex) mxGetScalar(prhs[3]);
                const gsKnotVector<>& kv = instance.tensorLevel(level-1).knots(direction-1);
                plhs[0] = createPointerFromStdVector(kv);
            }
            else if (!strcmp(prop,"active"))
            {
                // Copy the input (FIXME: this should be avoided)
                gsMatrix<real_t> pts = extractMatrixFromPointer<real_t>(prhs[2]);
                // Call method
                const gsMatrix<unsigned> vals = instance.active(pts);
                // Copy the result for output (FIXME: this should be avoided)
                plhs[0] = createPointerFromMatrix<unsigned>(vals);
            }
            else
            {
                throw("call Member: unknown command.");
            }
}

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

        if (!strcmp(cmd,"constructor"))
        {
            // ----------------------------------------------------------------------
            // Constructors

            if (nrhs==3) {
                // constructor from 1 argument (+ 1 type switch)
                char constructSwitch[__MAXSTRLEN__];
                if (mxGetString(prhs[1], constructSwitch, sizeof(constructSwitch)))
                    throw("Second input argument should be a string"
                          "less than MAXSTRLEN characters long.");
                if (!strcmp(constructSwitch,"char"))
                {
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
                }
                else if (!strcmp(constructSwitch,"gsTensorBSplineBasis"))
                {
                    gsBasis<real_t> * instance = convertMat2Ptr<gsBasis<real_t> >(prhs[1]);
                    const int dim = instance->domainDim();
                    // Use preprocessor to shorten switch ??
                    switch (dim)
                    {
                    case 2:
                        plhs[0] = convertPtr2Mat<gsTHBSplineBasis<2> >( new gsTHBSplineBasis<2>(*instance) );
                        break;
                    case 3:
                        plhs[0] = convertPtr2Mat<gsTHBSplineBasis<3> >( new gsTHBSplineBasis<3>(*instance) );
                        break;
                    default:
                        GISMO_ERROR("Error in gsTHBSplineBasis constructor.");
                    }

                } else {
                    throw ("gsTHBSplineBasis: Invalid construction.");
                }
            } else {
                // INVALID constructor
                throw("Invalid construction.");
            }

        } else if (!strcmp(cmd,"destructor"))
        {

            // ----------------------------------------------------------------------
            // Destructor
            destroyObject<gsBasis<real_t> >(prhs[1]);

        }
        else
        {
            // ----------------------------------------------------------------------
            // Member functions
            // Fetch instance and property to be accessed
            gsBasis<real_t> * instance = convertMat2Ptr<gsBasis<real_t> >(prhs[1]);
            const int dim = instance->domainDim();
            switch (dim)
            {
            case 2:
                callMember(*convertMat2Ptr<gsTHBSplineBasis<2> >(prhs[1]), cmd, plhs, prhs);
                break;
            case 3:
                callMember(*convertMat2Ptr<gsTHBSplineBasis<3> >(prhs[1]), cmd, plhs, prhs);
                break;
            default:
                GISMO_ERROR("Error in gsTHBSplineBasis accessor.");
            }
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
