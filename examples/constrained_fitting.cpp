/** @file constrained_fitting.cpp

    @brief Demonstration of a least-squares fitting (see
    fitting_example.cpp) with constraints. In this case, the values
    for v=0 and v=1 are prescribed by the user and only the remaining
    DOFs are computed by the LS. Cf. Prautzch, Boehm, Paluszny: Bezier
    and B-spline techniques, Section 4.7.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): D. Mokris
*/

#include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    // The following parameters are hard-wired but could easily be
    // made available through the command line analogously to
    // fitting_example.cpp.
    int iter = 5;
    real_t tolerance = 1e-5;
    real_t threshold = 1e-5;
    real_t v_min = -1;
    real_t v_max = 1;
    index_t num_int_knots = 0;
    index_t bound_mult = 4;
    real_t lambda = 1e-6;
    index_t extension = 2;
    real_t refPercent = 0.1;
    bool save = true;

    std::string uvxyz_file = "fitting/deepdrawingC.xml";

    // Read the fitting points and their parameters.
    gsFileData<> fd_in(uvxyz_file);
    gsMatrix<> uv, xyz;
    fd_in.getId<gsMatrix<> >(0, uv );
    fd_in.getId<gsMatrix<> >(1, xyz);

    // Design the curves prescribing the values along v=0 and v=1.
    gsMatrix<real_t> coefs0(5, 3);
    coefs0 << -50, -50, 25,
	-25, -50, 25,
	0, -50, 25,
	25, -50, 25,
	50, -50, 25;
    gsBSpline<> v0(-1, 1, 1, 3, coefs0);

    gsMatrix<real_t> coefs1(5, 3);
    coefs1 << -50, 50, 25,
	-25, 50, 25,
	0, 50, 25,
	25, 50, 25,
	50, 50, 25;
    gsBSpline<> v1(-1, 1, 1, 3, coefs1);

    // Find an appropriate basis. One could also modify the code by using
    // further parameters from fitting_example.cpp.
    // Note: u_knots are the knots in the u-direction, i.e., of the curve with v=0.
    gsKnotVector<> u_knots (v0.knots());
    GISMO_ASSERT(u_knots == v1.knots(), "The knot vectors in the v-direction differ.");
    gsKnotVector<real_t> v_knots(v_min, v_max, num_int_knots, bound_mult);
    gsTensorBSplineBasis<2> basis0(u_knots, v_knots);
    gsTHBSplineBasis<2> basis(basis0);

    // Extension, refin and lambda could also easily changed.
    std::vector<unsigned> ext;
    ext.push_back(extension);
    ext.push_back(extension);
    gsHFitting<2, real_t> fitting(uv, xyz, basis, refPercent, ext, lambda);

    // Prescribe the curves with v=0 and v=1.
    std::vector<gsBSpline<real_t> > prescribedCurves;
    prescribedCurves.push_back(v0);
    prescribedCurves.push_back(v1);
    std::vector<boxSide>     prescribedSides;
    prescribedSides.push_back(boxSide(3));
    prescribedSides.push_back(boxSide(4));
    fitting.setConstraints(prescribedSides, prescribedCurves);

    gsStopwatch time;
    for(int i = 0; i <= iter; i++)
    {
        gsInfo<<"----------------\n";
        gsInfo<<"Iteration "<<i<<".."<<"\n";

        time.restart();
        fitting.nextIteration(tolerance, threshold, prescribedSides);
        time.stop();
        gsInfo<<"Fitting time: "<< time <<"\n";
        gsInfo<<"Fitted with "<< fitting.result()->basis() <<"\n";
        gsInfo<<"Min distance : "<< fitting.minPointError() <<" / ";
        gsInfo<<"Max distance : "<< fitting.maxPointError() <<"\n";

	if(save)
	{
	    gsFileData<> fd_out;
	    fd_out << *fitting.result() ;
	    fd_out.dump("fitting_out" + util::to_string(i));
	}

        if ( fitting.maxPointError() < tolerance )
        {
            gsInfo<<"Error tolerance achieved after "<<i<<" iterations.\n";
            break;
        }
    }

    return 0;

}
