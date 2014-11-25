#include "mex.h"
#include <gismo.h>
#include "classHandle.hpp"
#include <cstring>

using namespace std;
using namespace gismo;

// This macro denotes the maximal length of strings passed to this mexFunction.
#define __MAXSTRLEN__ 128

// This macro denotes the parametric dimension
#define __DIM__ 2

// Copies the array pointed to by the input pointer and returns it
// as a gsMatrix.
template<class C>
gsMatrix<C> extractMatrixFromPointer(const mxArray *pnt) 
{
    return gsAsConstMatrix<real_t>(mxGetPr(pnt), mxGetM(pnt), mxGetN(pnt) );
}

// Copies the content of the input gsMatrix into a MATLAB real
// double matrix and returns a pointer to it. Note: memory is
// allocated (intended for output so it is not freed).
template<class C>
mxArray* createPointerFromMatrix( const gsMatrix<C> mat ) 
{
    //Note: there is the possibility to get the pointer to the matrix
    //data by: real_t * ptr = mat.data();
    const mwSize numRows = mat.rows(),
                 numCols = mat.cols();
    mxArray *pnt = mxCreateDoubleMatrix(numRows,numCols,mxREAL);
    real_t *out = mxGetPr(pnt);
    for (size_t j=0; j<numCols; ++j) 
    {
        for (size_t i=0; i<numRows; ++i) 
        {
            *out++ = mat(i,j);
        }
    }
  return pnt;
}


// --------------------------------------------------------------------------
// "main" gateway function
void mexFunction ( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

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
          gsTHBSplineBasis<2> * hbs = NULL;
          hbs = data.getFirst< gsTHBSplineBasis<2> >();
          gsTHBSplineBasis<__DIM__> * tmp = new gsTHBSplineBasis<__DIM__>(*hbs);
          plhs[0] = convertPtr2Mat<gsTHBSplineBasis<__DIM__> >(tmp);
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
        gsVector<real_t> supp = instance->support();
        plhs[0] = mxCreateDoubleMatrix(1, supp.size(), mxREAL);
        double *out = mxGetPr(plhs[0]);
        copy(supp.begin(), supp.end(), &out[0]);
      } else {
	// Unknown property
	throw("Third input argument contains an unknown property string.");
      }

    } else if (!strcmp(cmd,"treePrintLeaves")) {

      // ----------------------------------------------------------------------
      // tree().printLeaves()

      gsTHBSplineBasis<__DIM__> *instance = convertMat2Ptr<gsTHBSplineBasis<__DIM__> >(prhs[1]);
      instance->tree().printLeaves();

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

      gsTHBSplineBasis<__DIM__> *instance = convertMat2Ptr<gsTHBSplineBasis<__DIM__> >(prhs[1]);
      // Copy the input (FIXME: this should be avoided)
      const mwIndex     ind = (mwIndex) *mxGetPr(prhs[2]);
      const gsMatrix<real_t> pts = extractMatrixFromPointer<real_t>(prhs[3]);
      // Call the method
      gsMatrix<real_t> vals = instance->evalSingle(ind,pts);
      // Copy result to output (FIXME: this should be avoided)
      plhs[0] = createPointerFromMatrix<real_t>(vals);

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

  } catch (exception& e) { // Caught from library.
    string errMsg = string("\n  While executing the following command: ") + cmd + 
                    string("\n  The following exception/error ocurred: ") + e.what();
    mexErrMsgTxt(errMsg.c_str());
  } catch (const char* str) { // Caught from within this mexFunction.
    string errMsg = string("\n  While executing the following command: ") + cmd + 
                    string("\n  The following exception/error ocurred: ") + string(str);
    mexErrMsgTxt(errMsg.c_str());
  } catch (...) { // Something else went wrong
    string errMsg = string("\n  While executing the following command: ") + cmd + 
                    string("\n  An error ocurred.");
    mexErrMsgTxt(errMsg.c_str());
  } // end try-catch

} // end mexFunction
