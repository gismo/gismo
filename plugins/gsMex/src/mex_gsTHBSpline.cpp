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

template<class T, class Tbasis, class Ttensorbasis>
void mexFunctionTemplate ( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[], char* cmd);

template<class T>
void specializedSliceCoefs(mxArray* plhs[], const mxArray* prhs[]);

template<>
void specializedSliceCoefs<gsTHBSpline<1> >(mxArray* plhs[], const mxArray* prhs[]);

/*
template<class T>
void specializedComputeCoonsPatch(mxArray* plhs[], const mxArray* prhs[]);

template<>
void specializedComputeCoonsPatch<gsTHBSpline<1> >(mxArray* plhs[], const mxArray* prhs[]);
*/

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
    try {

        // Fetch command string.
        if (nrhs < 2 || mxGetString(prhs[0], cmd, sizeof(cmd)))
            throw ("First input argument should be a command string "
                   "less than MAXSTRLEN characters long, second input argument "
                   "should be an int8.");

        int8_t dimension = (int8_t) mxGetScalar(prhs[nrhs - 1]);
        if (dimension == 1)
            mexFunctionTemplate<gsTHBSpline<1>, gsTHBSplineBasis<1>, gsBSplineBasis<> >
                                                                          (nlhs, plhs, nrhs, prhs, cmd);
        else if (dimension == 2)
            mexFunctionTemplate<gsTHBSpline<2>, gsTHBSplineBasis<2>, gsTensorBSplineBasis <2> >
                                                                          (nlhs, plhs, nrhs, prhs, cmd);
        else if (dimension == 3)
            mexFunctionTemplate<gsTHBSpline<3>, gsTHBSplineBasis<3>, gsTensorBSplineBasis <3> >
                                                                          (nlhs, plhs, nrhs, prhs, cmd);
        else
            throw ("Dimension should be specified in the last argument "
                   "and must be equal to 1, 2 or 3.");

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

template<>
void specializedSliceCoefs<gsTHBSpline<1> >(mxArray* plhs[], const mxArray* prhs[]){
    throw("sliceCoefs cannot be called with a gsTHBSpline of parametric dimension 1.");
}

template<class T>
void specializedSliceCoefs(mxArray* plhs[], const mxArray* prhs[]){
    T *instance = convertMat2Ptr<T> (prhs[1]);
    typedef typename T::BoundaryGeometryType BoundaryGeometryType;

    // Copy the input (FIXME: this should be avoided)
    mwIndex dir_fixed = (mwIndex) mxGetScalar(prhs[2]);
    real_t par = (real_t) mxGetScalar(prhs[3]);
    BoundaryGeometryType* res = new BoundaryGeometryType();
    // Call method
    instance->slice(dir_fixed-1, par, *res);
    plhs[0] = createPointerFromMatrix(res->coefs());

}

/* template<>
void specializedComputeCoonsPatch<gsTHBSpline<1> >(mxArray* plhs[], const mxArray* prhs[]){
    throw("coonsPatch cannot be called with a gsTHBSpline of parametric dimension 1.");
}

template<class T>
void specializedComputeCoonsPatch(mxArray* plhs[], const mxArray* prhs[]){
    typedef typename T::BoundaryGeometryType BoundaryGeometryType;
    gsMultiPatch<> gs;
    BoundaryGeometryType *instance = convertMat2Ptr<BoundaryGeometryType>(prhs[1]);
    gs.addPatch(*instance);
    for (unsigned int i=2; i<=2*instance->parDim(); i++)
        gs.addPatch(*convertMat2Ptr<BoundaryGeometryType>(prhs[1]));
    gsCoonsPatch<double> coons(gs);
    gsGeometry<> res(coons.compute());
} */

template<class T, class Tbasis, class Ttensorbasis>
void mexFunctionTemplate ( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[], char* cmd) {
    if (!strcmp(cmd, "constructor")) {
        // ----------------------------------------------------------------------
        // Constructors

        if (nrhs == 4) {
            // constructor from 1 argument (+ 1 type switch)
            char constructSwitch[__MAXSTRLEN__];
            if (mxGetString(prhs[1], constructSwitch, sizeof(constructSwitch)))
                throw ("Second input argument should be a string"
                       "less than MAXSTRLEN characters long.");
            if (!strcmp(constructSwitch, "char")) {
                // constructor ( char filename )
                // Extract the input
                char *input_buf = mxArrayToString(prhs[2]);
                // Read the THB-spline basis from the specified file
                std::string filename(input_buf); // Reading requires a std::string
                gsFileData <real_t> data(filename);
                T *hbs = data.getFirst <T> ().release();
                plhs[0] = convertPtr2Mat<T> (hbs);
                // Free the memory allocated by mxArrayToString
                mxFree(input_buf);
            } else if (!strcmp(constructSwitch, "gsTHBSpline")) {
                // copy constructor ( T )
                T *instance = convertMat2Ptr<T> (prhs[2]);
                plhs[0] = convertPtr2Mat<T> (new T(*instance));

            /* } else if (!strcmp(constructSwitch, "gsTensorBSpline")) { // TODO
                // constructor ( gsTensorBSpline )
                Ttensorbasis *instance = convertMat2Ptr<Ttensorbasis> (prhs[2]);
                plhs[0] = convertPtr2Mat<T> (new T(*instance)); */
            } else {
                throw ("Invalid construction.");
            }
        } else if (nrhs == 6) {
            // constructor from 2 argument (+ 2 type switch)
            char constructSwitch1[__MAXSTRLEN__];
            char constructSwitch2[__MAXSTRLEN__];
            if (mxGetString(prhs[1], constructSwitch1, sizeof(constructSwitch1))) {
                throw ("Second input argument should be a string"
                       "less than MAXSTRLEN characters long.");
            }
            if (mxGetString(prhs[2], constructSwitch2, sizeof(constructSwitch2))) {
                throw ("Third input argument should be a string"
                       "less than MAXSTRLEN characters long.");
            }
            if (!(strcmp(constructSwitch1, "gsTHBSplineBasis") || strcmp(constructSwitch2, "double"))) {
                // constructor ( gsTHBSplineBasis, controlPts )
                const gsMatrix <real_t> coefs = extractMatrixFromPointer<real_t>(prhs[4]);
                Tbasis *instance = convertMat2Ptr<Tbasis> (prhs[3]);
                plhs[0] = convertPtr2Mat<T> (new T(*instance, coefs));
            } else {
                throw ("Invalid construction.");
            }
        } else {
            // INVALID constructor
            throw ("Invalid construction.");
        }

    } else if (!strcmp(cmd, "destructor")) {

        // ----------------------------------------------------------------------
        // Destructor
        destroyObject<T> (prhs[1]);

    } else if (!strcmp(cmd, "accessor")) {

        // ----------------------------------------------------------------------
        // Accessor

        // Fetch instance and property to be accessed
        T *instance = convertMat2Ptr<T> (prhs[1]);
        char prop[__MAXSTRLEN__];
        if (mxGetString(prhs[2], prop, sizeof(prop)))
            throw ("Third input argument should be a property string less than MAXSTRLEN characters long.");

        // Call method as specified by the input string
        if (!strcmp(prop, "parDim")) {
            int val = instance->parDim();
            mxArray *out = mxCreateDoubleScalar((double) val);
            plhs[0] = out;
        } else if (!strcmp(prop, "geoDim")) {
            int val = instance->geoDim();
            mxArray *out = mxCreateDoubleScalar((double) val);
            plhs[0] = out;
        } else if (!strcmp(prop, "size")) {
            int val = instance->size();
            mxArray *out = mxCreateDoubleScalar((double) val);
            plhs[0] = out;
        } else if (!strcmp(prop, "support")) {
            gsMatrix <real_t> supp = instance->support();
            plhs[0] = createPointerFromMatrix(supp);
        } else if (!strcmp(prop, "basis")) {
            // Copy the result for output (FIXME: this should be avoided)
            Tbasis *hbs = new Tbasis(instance->basis());
            plhs[0] = convertPtr2Mat<Tbasis> (hbs);
        } else if (!strcmp(prop, "coefs")) {
            // const gsMatrix<>& cc = instance->coefs();
            plhs[0] = createPointerFromMatrix(instance->coefs());
        } else {
            // Unknown property
            throw ("Third input argument contains an unknown property string.");
        }

    } else if (!strcmp(cmd, "eval")) {

        // ----------------------------------------------------------------------
        // eval(pts)
        T *instance = convertMat2Ptr<T> (prhs[1]);
        // Copy the input (FIXME: this should be avoided)
        const gsMatrix <real_t> pts = extractMatrixFromPointer<real_t>(prhs[2]);
        // Call the method
        gsMatrix <real_t> vals = instance->eval(pts);
        // Copy result to output (FIXME: this should be avoided)
        plhs[0] = createPointerFromMatrix<real_t>(vals);

    } else if (!strcmp(cmd, "jacobian")) {

        // ----------------------------------------------------------------------
        // deriv(pts)

        T *instance = convertMat2Ptr<T> (prhs[1]);
        // Copy the input (FIXME: this should be avoided)
        const gsMatrix <real_t> pts = extractMatrixFromPointer<real_t>(prhs[2]);
        // Call the method
        gsMatrix <real_t> vals = instance->jacobian(pts);
        // Copy result to output (FIXME: this should be avoided)
        plhs[0] = createPointerFromMatrix<real_t>(vals);

    } else if (!strcmp(cmd, "hess")) {

        // ----------------------------------------------------------------------
        // hess(pts, dir)

        T *instance = convertMat2Ptr<T> (prhs[1]);
        // Copy the input (FIXME: this should be avoided)
        const gsMatrix <real_t> pts = extractMatrixFromPointer<real_t>(prhs[2]);
        mwIndex dir = (mwIndex) mxGetScalar(prhs[3]);
        // Call the method
        gsMatrix <real_t> vals = instance->hess(pts, dir - 1);
        // Copy result to output (FIXME: this should be avoided)
        plhs[0] = createPointerFromMatrix<real_t>(vals);

    } else if (!strcmp(cmd, "sliceCoefs")) {

        // ----------------------------------------------------------------------
        // sliceCoefs(dir_fixed, par)
        specializedSliceCoefs<T>(plhs, prhs);

    } else if (!strcmp(cmd, "save")) {

        // ----------------------------------------------------------------------
        // save(file)

        T *instance = convertMat2Ptr<T>(prhs[1]);
        char *input_buf = mxArrayToString(prhs[2]);
        // Save the THB-spline in the specified file
        std::string filename(input_buf); // Reading requires a std::string
        gsWrite(*instance, filename);

 /*   } else if (!strcmp(cmd, "computeCoonsPatch")) {

        // ----------------------------------------------------------------------
        // computeCoonsPatch(boundaries)
        specializedComputeCoonsPatch<T>(plhs, prhs);
*/
    } else {

        // ----------------------------------------------------------------------
        // Unknown command
        throw ("unknown command.");

    }
}
