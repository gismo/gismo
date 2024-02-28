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

template<int D> void
callMember(const gsTensorBSplineBasis<D,T> & instance,  char * prop,
           mxArray* plhs[], const mxArray* prhs[])
{
    if (!strcmp(prop,"save"))
    {
        char* input_buf = mxArrayToString(prhs[2]);
        // Save the THB-spline basis in the specified file
        std::string filename(input_buf); // Reading requires a std::string
        gsWrite(instance, filename);
    }
    else if (!strcmp(prop,"knots"))
    {
        mwIndex direction = (mwIndex) mxGetScalar(prhs[2]);
        const gsKnotVector<>& kv = instance.knots(direction-1);
        plhs[0] = createPointerFromStdVector(kv);
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

            if (nrhs>3)
            {
                const short_t DIM = nrhs-2;
                // GISMO_ENSURE(__DIM__==nrhs-2,"Number of input arguments should be "<<__DIM__-2);
                char constructSwitch[__MAXSTRLEN__];
                if (mxGetString(prhs[1], constructSwitch, sizeof(constructSwitch)))
                    throw("Second input argument should be a string"
                          "less than MAXSTRLEN characters long.");

                std::vector<gsKnotVector<T>> kvs(DIM);
                for (index_t d=0; d!=DIM; d++)
                {
                    if (!strcmp(constructSwitch,"gsKnotVector"))
                    {
                        // This does not work
                        kvs[d] = *convertMat2Ptr<gsKnotVector<T> >(prhs[2+d]);
                    }
                    else if (!strcmp(constructSwitch,"double"))
                    {
                        // argument 2: list of knots
                        if (mxIsScalar(prhs[2+d])) { // or not scalar
                            mexErrMsgIdAndTxt("MATLAB:mexcpp:typeargin", "Argument ",std::to_string(d)," has to be double vector.");
                        }
                        if (!mxGetM(prhs[2+d])==1) { // or not scalar
                            mexErrMsgIdAndTxt("MATLAB:mexcpp:typeargin", "Argument ",std::to_string(d)," must be a row vector.");
                        }
                        T *knotptr = mxGetDoubles(prhs[2+d]);
                        std::vector<real_t> knots(knotptr,knotptr+mxGetN(prhs[2+d]));
                        kvs[d] = gsKnotVector<T>(knots);
                    }
                    else
                    {
                        throw ("Invalid construction.");
                    }
                }

                switch (DIM)
                {
                case 2:
                {
                    gsTensorBSplineBasis<2> * tbasis = new gsTensorBSplineBasis<2>(kvs);
                    plhs[0] = convertPtr2Mat<gsTensorBSplineBasis<2> >(tbasis);
                    break;
                }
                case 3:
                {
                    gsTensorBSplineBasis<3> * tbasis = new gsTensorBSplineBasis<3>(kvs);
                    plhs[0] = convertPtr2Mat<gsTensorBSplineBasis<3> >(tbasis);
                    break;
                }
                default:
                    GISMO_ERROR("Error in gsTensorBSplineBasis constructor.");
                }
            }
            else
            {
                // INVALID constructor
                throw("Invalid construction.");
            }

        }
        else if (!strcmp(cmd,"destructor"))
        {

            // ----------------------------------------------------------------------
            // Destructor
            destroyObject<gsBasis<T> >(prhs[1]);

        }
        else
        {
            // ----------------------------------------------------------------------
            // Member functions
            // Fetch instance and property to be accessed
            gsBasis<T> * instance = convertMat2Ptr<gsBasis<T> >(prhs[1]);
            const int dim = instance->domainDim();
            switch (dim)
            {
            case 2:
                callMember(*convertMat2Ptr<gsTensorBSplineBasis<2> >(prhs[1]), cmd, plhs, prhs);
                break;
            case 3:
                callMember(*convertMat2Ptr<gsTensorBSplineBasis<3> >(prhs[1]), cmd, plhs, prhs);
                break;
            default:
                GISMO_ERROR("Error in gsTensorBSplineBasis accessor.");
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
