/** @file mex_gsTensorBSplineBasis.cpp

    @brief Mex adaptor for gsTensorBSplineBasis class

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): O. Chanon
*/

#include <gismo.h>
#include "mex_common.h"

using namespace gismo;

template<class T>
void mexFunctionTemplate ( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[], char* cmd);

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
        if (nrhs < 2 || mxGetString(prhs[0], cmd, sizeof(cmd)))
            throw("First input argument should be a command string"
                  "less than MAXSTRLEN characters long.");

        int8_t dimension = (int8_t) mxGetScalar(prhs[nrhs-1]);
       if (dimension == 2)
            mexFunctionTemplate < gsTensorBSplineBasis<2> > (nlhs, plhs, nrhs, prhs, cmd);
        else if (dimension == 3)
            mexFunctionTemplate < gsTensorBSplineBasis<3> > (nlhs, plhs, nrhs, prhs, cmd);
        else
            throw ("Dimension should be specified in the last argument"
                   "and must be equal to 2 or 3.");

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


template<class T>
void mexFunctionTemplate ( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[], char* cmd) {
        if (!strcmp(cmd,"constructor")) {

            // ----------------------------------------------------------------------
            // Constructors

            if (nrhs==4) {
                // constructor from 1 argument (+ 1 type switch)
                char constructSwitch[__MAXSTRLEN__];
                if (mxGetString(prhs[1], constructSwitch, sizeof(constructSwitch)))
                    throw("Second input argument should be a string"
                          "less than MAXSTRLEN characters long.");
                if (!strcmp(constructSwitch,"char")) {
                    // constructor ( char filename )
                    // Extract the input
                    char* input_buf = mxArrayToString(prhs[2]);
                    // Read the B-spline basis from the specified file
                    std::string filename(input_buf); // Reading requires a std::string
                    gsFileData<real_t>  data( filename );
                    T * bspb = data.getFirst<T>().release();
                    plhs[0] = convertPtr2Mat<T>(bspb);
                    // Free the memory allocated by mxArrayToString
                    mxFree(input_buf);
                }
                else if (!strcmp(constructSwitch,"gsTensorBSplineBasis")) {
                    // constructor ( gsTensorBSplineBasis )
                    T *instance = convertMat2Ptr<T>(prhs[2]);
                    plhs[0] = convertPtr2Mat<T>(new T(*instance));
                }
                else if (!strcmp(constructSwitch,"cell")) {
                    mwSize dim(mxGetNumberOfElements(prhs[2]));

                    // Set up and construct the knot vectors...
                    std::vector<gsKnotVector<real_t>> kts(dim);
                    for (mwIndex index=0; index<dim; index++) {
                        const mxArray *cell_element_ptr(mxGetCell(prhs[2], index));
                        if (cell_element_ptr == NULL)
                            throw ("Invalid construction: void cell element.");
                        else {
                            mwSize total_num_of_knots(mxGetNumberOfElements(cell_element_ptr));
                            double *pr(mxGetPr(cell_element_ptr));
                            std::vector<real_t> knots_dim(total_num_of_knots);
                            for (mwIndex index2=0; index2 < total_num_of_knots; index2++) {
                                knots_dim[index2] = *pr++;
                            }
                            kts[index] = gsKnotVector<>(knots_dim);
                        }
                    }
                    // ...a 2D-tensor-B-spline basis with this knot vector...
                    T tens(kts);

                    // ...and a 2D B-spline basis out of the tensor-B-spline basis.
                    T bspb( tens );
                    plhs[0] = convertPtr2Mat<T>(new T(bspb));
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
            destroyObject<T>(prhs[1]);

        } else if (!strcmp(cmd,"accessor")) {

            // ----------------------------------------------------------------------
            // Accessor

            // Fetch instance and property to be accessed
            T *instance = convertMat2Ptr<T>(prhs[1]);
            char prop[__MAXSTRLEN__];
            if (mxGetString(prhs[2], prop, sizeof(prop)))
                throw("Third input argument should be a property string less than MAXSTRLEN characters long.");

            // Call method as specified by the input string
            if (!strcmp(prop,"domainDim")) {
                int val      = instance->domainDim();
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
            } else if (!strcmp(prop,"support")) {
                gsMatrix<real_t> supp = instance->support();
                plhs[0] = createPointerFromMatrix(supp);
            }

        } else if (!strcmp(cmd,"degree")) {

            T *instance = convertMat2Ptr <T> (prhs[1]);
            mwIndex deg = (mwIndex) mxGetScalar(prhs[2]);
            mxArray *out = mxCreateDoubleScalar((double)instance->degree(deg-1));
            plhs[0] = out;

        } else if (!strcmp(cmd,"eval")) {

            // ----------------------------------------------------------------------
            // eval(pts)

            T *instance = convertMat2Ptr<T>(prhs[1]);
            // Copy the input (FIXME: this should be avoided)
            const gsMatrix<real_t> pts = extractMatrixFromPointer<real_t>(prhs[2]);
            // Call the method
            gsMatrix<real_t> vals = instance->eval(pts);
            // Copy result to output (FIXME: this should be avoided)
            plhs[0] = createPointerFromMatrix<real_t>(vals);

        } else if (!strcmp(cmd,"evalSingle")) {

            // ----------------------------------------------------------------------
            // evalSingle(ind,pts)

            T *instance = convertMat2Ptr <T> (prhs[1]);
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

            T *instance = convertMat2Ptr <T> (prhs[1]);
            char* input_buf = mxArrayToString(prhs[2]);
            // Save the B-spline basis in the specified file
            std::string filename(input_buf); // Reading requires a std::string
            gsWrite(*instance, filename);

        } else if (!strcmp(cmd,"knots")) {

            // ----------------------------------------------------------------------
            // knots(direction)

            T *instance = convertMat2Ptr <T> (prhs[1]);
            mwIndex direction = (mwIndex) mxGetScalar(prhs[2]);
            const gsKnotVector<>& kv = instance->knots(direction-1);
            plhs[0] = createPointerFromStdVector(kv);

        } else if (!strcmp(cmd,"active")) {

            // ----------------------------------------------------------------------
            // active(pts)

            T *instance = convertMat2Ptr<T>(prhs[1]);
            // Copy the input (FIXME: this should be avoided)
            gsMatrix<real_t> pts = extractMatrixFromPointer<real_t>(prhs[2]);
            // Call method
            gsMatrix<unsigned> vals = instance->active(pts);
            vals = vals + gsMatrix<unsigned>::Ones(vals.rows(),vals.cols());
            // Copy the result for output (FIXME: this should be avoided)
            plhs[0] = createPointerFromMatrix<unsigned>(vals);

        } else if (!strcmp(cmd,"uniformRefine")) {

            // ----------------------------------------------------------------------
            // uniformRefine(numKnots, mult)

            T *instance = convertMat2Ptr<T>(prhs[1]);
            // Copy the input (FIXME: this should be avoided)
            mwIndex numKnots = (mwIndex) mxGetScalar(prhs[2]);
            mwIndex mult = (mwIndex) mxGetScalar(prhs[3]);

            instance->uniformRefine(numKnots, mult);

        } else if (!strcmp(cmd,"uniformRefine_withCoefs")) {

            // ----------------------------------------------------------------------
            // uniformRefine_withCoefs(coefs, numKnots, mult)

            T *instance = convertMat2Ptr<T>(prhs[1]);
            // Copy the input (FIXME: this should be avoided)
            gsMatrix<real_t> coefs = extractMatrixFromPointer<real_t>(prhs[2]);
            mwIndex numKnots = (mwIndex) mxGetScalar(prhs[3]);
            mwIndex mult = (mwIndex) mxGetScalar(prhs[4]);

            instance->uniformRefine_withCoefs(coefs, numKnots, mult);
            plhs[0] = createPointerFromMatrix<real_t>(coefs);

        } else if (!strcmp(cmd,"refineElements")) {

            // ----------------------------------------------------------------------
            // refineElements(boxes)

            T *instance = convertMat2Ptr<T>(prhs[1]);
            // Copy the input (FIXME: this should be avoided)
            std::vector<unsigned int> boxes = extractStdVectorFromPointer<unsigned int>(prhs[2]);
            std::for_each(boxes.begin(), boxes.end(), [](unsigned int& d) { d-=1;});

            instance->refineElements(boxes);

        /* } else if (!strcmp(cmd,"refineElements_withCoefs")) {

            // ----------------------------------------------------------------------
            // refineElements_withCoefs(coefs, boxes)

            T *instance = convertMat2Ptr<T> (prhs[1]);
            // Copy the input (FIXME: this should be avoided)
            gsMatrix <real_t> coefs = extractMatrixFromPointer<real_t>(prhs[2]);
            std::vector<unsigned int> boxes = extractStdVectorFromPointer<unsigned int>(prhs[3]);

            std::for_each(boxes.begin(), boxes.end(), [](unsigned int &d) { d -= 1; });

            instance->refineElements_withCoefs(coefs, boxes); // VIRTUAL METHOD
            plhs[0] = createPointerFromMatrix<real_t>(coefs); */

        } else {

            // ----------------------------------------------------------------------
            // Unknown command
            throw("unknown command.");

        }

}
