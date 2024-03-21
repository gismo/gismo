/** @file mex_gsKnotVector.cpp

    @brief Mex adaptor for gsKnotVector class

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
// prhs[0]: output
//
// prhs -> Matlab input
// prhs[0]: command name
// prhs[1]: pointer to matlab object
// prhs[2]: first argument
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

            if (nrhs==3)
            {
                // constructor from 1 argument (+ 1 type switch)
                char constructSwitch[__MAXSTRLEN__];
                if (mxGetString(prhs[1], constructSwitch, sizeof(constructSwitch)))
                    throw("Second input argument should be a string"
                          "less than MAXSTRLEN characters long.");
                if (!strcmp(constructSwitch,"double"))
                {
                    // argument 2: list of knots
                    if (!mxIsDouble(prhs[2]) || // not double
                         mxIsScalar(prhs[2]))
                    { // or not scalar
                        mexErrMsgIdAndTxt("MATLAB:mexcpp:typeargin", "Second argument has to be double vector.");
                    }
                    if (!mxGetM(prhs[2])==1)
                    { // or not scalar
                        mexErrMsgIdAndTxt("MATLAB:mexcpp:typeargin", "Second argument must be a row vector.");
                    }

                    T * knotptr = mxGetDoubles(prhs[2]);
                    gsKnotVector<T> * kv = new gsKnotVector<T>(-1,knotptr,knotptr+mxGetN(prhs[2]));
                    //std::vector<real_t> knots(knotptr,knotptr+mxGetN(prhs[2]));
                    //gsKnotVector<T> * kv = new gsKnotVector<T>(knots);
                    plhs[0] = convertPtr2Mat<gsKnotVector<T> >(kv);
                }
                else
                {
                    throw ("Invalid construction.");
                }
            } else
            {
                // INVALID constructor
                throw("Invalid construction.");
            }

        } else if (!strcmp(cmd,"destructor")) {

            // ----------------------------------------------------------------------
            // Destructor
            destroyObject<gsKnotVector<T>>(prhs[1]);

        // } else if (!strcmp(cmd,"accessor")) {

        //     // ----------------------------------------------------------------------
        //     // Accessor

        //     // Fetch instance and property to be accessed
        //     gsKnotVector<T> *instance = convertMat2Ptr< gsKnotVector<T> >(prhs[1]);
        //     char prop[__MAXSTRLEN__];
        //     if (mxGetString(prhs[2], prop, sizeof(prop)))
        //         throw("Third input argument should be a property string less than MAXSTRLEN characters long.");

        //     // Call method as specified by the input string
        //     if (!strcmp(prop,"dim")) {
        //         int val      = instance->dim();
        //         mxArray *out = mxCreateDoubleScalar((double)val);
        //         plhs[0]      = out;
        //     } else if (!strcmp(prop,"numElements")) {
        //         int val      = instance->numElements();
        //         mxArray *out = mxCreateDoubleScalar((double)val);
        //         plhs[0]      = out;
        //     } else if (!strcmp(prop,"size")) {
        //         int val      = instance->size();
        //         mxArray *out = mxCreateDoubleScalar((double)val);
        //         plhs[0]      = out;
        //     } else if (!strcmp(prop,"treeSize")) {
        //         int val      = instance->treeSize();
        //         mxArray *out = mxCreateDoubleScalar((double)val);
        //         plhs[0]      = out;
        //     } else if (!strcmp(prop,"treeLeafSize")) {
        //         int val      = instance->tree().leafSize();
        //         mxArray *out = mxCreateDoubleScalar((double)val);
        //         plhs[0]      = out;
        //     } else if (!strcmp(prop,"support")) {
        //         gsMatrix<T> supp = instance->support();
        //         plhs[0] = createPointerFromMatrix(supp);
        //     } else if (!strcmp(prop,"maxLevel")) {
        //         int val      = instance->maxLevel()+1;
        //         mxArray *out = mxCreateDoubleScalar((double)val);
        //         plhs[0]      = out;
        //     }

        // } else if (!strcmp(cmd,"treePrintLeaves")) {

        //     // ----------------------------------------------------------------------
        //     // tree().printLeaves()

        //     gsKnotVector<T> *instance = convertMat2Ptr <  gsKnotVector<T> > (prhs[1]);
        //     instance->tree().printLeaves();

        } else if (!strcmp(cmd,"degree")) {

            gsKnotVector<T> *instance = convertMat2Ptr< gsKnotVector<T> >(prhs[1]);
            mxArray *out = mxCreateDoubleScalar(instance->degree());
            plhs[0]      = out;

        // } else if (!strcmp(cmd,"eval")) {

        //     // ----------------------------------------------------------------------
        //     // eval(pts)

        //     gsKnotVector<T> *instance = convertMat2Ptr< gsKnotVector<T> >(prhs[1]);
        //     // Copy the input (FIXME: this should be avoided)
        //     const gsMatrix<T> pts = extractMatrixFromPointer<T>(prhs[2]);
        //     // Call the method
        //     gsMatrix<T> vals = instance->eval(pts);
        //     // Copy result to output (FIXME: this should be avoided)
        //     plhs[0] = createPointerFromMatrix<T>(vals);

        // } else if (!strcmp(cmd,"evalSingle")) {

        //     // ----------------------------------------------------------------------
        //     // evalSingle(ind,pts)

        //     gsKnotVector<T> *instance = convertMat2Ptr <  gsKnotVector<T> > (prhs[1]);
        //     // Copy the input (FIXME: this should be avoided)
        //     const mwIndex ind = (mwIndex) * mxGetPr(prhs[2]);
        //     const gsMatrix <T> pts = extractMatrixFromPointer<T>(prhs[3]);
        //     // Call the method
        //     gsMatrix <T> vals = instance->evalSingle(ind-1, pts);
        //     // Copy result to output (FIXME: this should be avoided)
        //     plhs[0] = createPointerFromMatrix<T>(vals);

        // } else if (!strcmp(cmd,"save")) {

        //     // ----------------------------------------------------------------------
        //     // save(file)

        //     gsKnotVector<T> *instance = convertMat2Ptr <  gsKnotVector<T> > (prhs[1]);
        //     char* input_buf = mxArrayToString(prhs[2]);
        //     // Save the THB-spline basis in the specified file
        //     std::string filename(input_buf); // Reading requires a std::string
        //     gsWrite(*instance, filename);

        // } else if (!strcmp(cmd,"knots")) {

        //     // ----------------------------------------------------------------------
        //     // knots(level, direction)

        //     gsKnotVector<T> *instance = convertMat2Ptr <  gsKnotVector<T> > (prhs[1]);
        //     mwIndex level = (mwIndex) mxGetScalar(prhs[2]);
        //     mwIndex direction = (mwIndex) mxGetScalar(prhs[3]);
        //     const gsKnotVector<>& kv = instance->tensorLevel(level-1).knots(direction-1);
        //     plhs[0] = createPointerFromStdVector(kv);

        // } else if (!strcmp(cmd,"active")) {

        //     // ----------------------------------------------------------------------
        //     // active(pts)

        //     gsKnotVector<T> *instance = convertMat2Ptr< gsKnotVector<T> >(prhs[1]);
        //     // Copy the input (FIXME: this should be avoided)
        //     gsMatrix<T> pts = extractMatrixFromPointer<T>(prhs[2]);
        //     // Call method
        //     const gsMatrix<unsigned> vals = instance->active(pts);
        //     // Copy the result for output (FIXME: this should be avoided)
        //     plhs[0] = createPointerFromMatrix<unsigned>(vals);
        } else if (!strcmp(cmd,"get")) {
            gsKnotVector<T> *instance = convertMat2Ptr< gsKnotVector<T> >(prhs[1]);
            plhs[0] = createPointerFromStdVector(*instance);
        } else if (!strcmp(cmd,"print")) {
            gsKnotVector<T> *instance = convertMat2Ptr< gsKnotVector<T> >(prhs[1]);
            std::ostringstream a;
            instance->print(a);
            mexPrintf("%s\n", a.str().c_str());
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
