
----------------------------------------------------------------------
 Introduction
----------------------------------------------------------------------

 * This is a MEX (Matlab EXecutable) interface for G+SMO.

----------------------------------------------------------------------
 Dependecies
----------------------------------------------------------------------

 To generate the .mex files, you need mex (comes with MATLAB).

 To use the .mex files in MATLAB, you need MATLAB.

----------------------------------------------------------------------
 Source tree
----------------------------------------------------------------------

 ./ : contains this file (README) and a license file for the
      classHandle.

 /include : contains relevant headers (currently only a header for
            safe memory treatment of classes, authored by Oliver
            Woodward, see license.txt).

 src/ : contains the mex source class definitions (gateways between
        C++ and MATLAB).

 m/   : constains mostly MATLAB class definitions file. Each of these
        call the relevant mex file. Also contains a test script 
        testTHBSplineBasis.m.

----------------------------------------------------------------------
 Using the Gismo class in MATLAB
----------------------------------------------------------------------

 * Make sure that the gs* methods are in the MATLAB path.

 * To construct an instance of gsTHBSplineBasis from an input file:
   > hbs = gsTHBSplineBasis(file);

 * To inspect the available methods for that class:
   > methods(hbs);

 * To get help on the method 'support':
   > help gsTHBSplineBasis/support

 * To test the gsTHBSplineBasis:
   > gsTestTHBSplineBasis

----------------------------------------------------------------------
 Notes
----------------------------------------------------------------------

 * All type-checking is done in the m-files, while practically no
   validations of the input and output arguments are performed in the
   mex-file 'work horse'.

----------------------------------------------------------------------
 Issues
----------------------------------------------------------------------

 * If you experience linking problems relating to libstdc++.so and
   GLIBCXX: 

   -> This could be because Matlab uses an old shared object file
      (e.g. /usr/local/matlab/R2012a/sys/os/glnx86/libstdc++.so.6). This
      could perhaps be resolved by replacing that link with a link to
      a newer shared object file in your installation
      (e.g. /usr/lib/i386-linux-gnu/libstdc++.so.6). Or even better:
      by enforcing the use of this at the time of linking.

 * If you experience errors like "Unknown symbol: _ZN" when running a
   mex function:

   -> In Linux, try updating the LD_LIBRARY_PATH to include the
      directory in which the gismo shared library (libgismo.so) is
      located. In bash:

      export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/shared/library/folder

      This pertains to the development phase. When gismo is installed,
      the shared library should be found without any difficulties.

----------------------------------------------------------------------
 Author: Peter Noertoft
----------------------------------------------------------------------
