/** @file embeddedCurve_example.cpp

    @brief Example for embedded curve

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Chianese, H.M. Verhelst
*/

#include <iostream>
#include <gismo.h>

using namespace gismo;

template <class T>
gsMatrix<T> findRoot(const gsGeometry<T> &geometry, gsKnotVector<T> &kv, const gsMatrix<T> &value, const index_t i=0,
            const short_t dir=0, const index_t maxiter=50, const T epsilon=1e-4);

int main(int argc, char *argv[])
{
    // bool plot = false; // If set to true, paraview file is generated and launched on exit
    // bool trim = false; // If set to true, trim/merge operations are displayed
    // bool intersect = false; // If set to true, intersection example is displayed

    gsCmdLine cmd("TODO");
    // cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    // cmd.addSwitch("trim", "Basic trim/merge operations", trim);
    // cmd.addSwitch("intersect", "Intersection operations", intersect);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Superdomain
    gsKnotVector<> kv_ξ(0, 1, 3, 3); //along ξ dimension
    gsKnotVector<> kv_η(0, 1, 3, 3); //along η dimension

    gsMatrix<> coefs_s(36,2); //x,y;..
    coefs_s << 0,0,0.5,0,1.5,0,2.5,0,3.5,0,4,0,0,0.5,0.5,0.5,1.5,0.5,2.5,0.5,3.5,0.5,4,0.5,
               0,1.5,0.5,1.5,1.5,1.5,2.5,1.5,3.5,1.5,4,1.5,0,2.5,0.5,2.5,1.5,2.5,2.5,2.5,3.5,2.5,4,2.5,
               0,3.5,0.5,3.5,1.5,3.5,2.5,3.5,3.5,3.5,4,3.5,0,4,0.5,4,1.5,4,2.5,4,3.5,4,4,4;

    gsTensorBSplineBasis<2,real_t> basis(kv_ξ, kv_η);
    gsTensorBSpline<2, real_t>  surface(basis, coefs_s);

    // Print the superdomain
    gsInfo << "I am the superdomain" << "\n";
    gsInfo << "Knot vector along ξ dimension:" << surface.knots(0) << "\n";
    gsInfo << "Knot vector along η dimension:" << surface.knots(1) << "\n";
    gsInfo << "Control points:" << "\n";
    gsInfo << surface.coefs() << "\n";

    // Embedded entity in the superdomain parameter space
    gsKnotVector<> kv_c(0, 1, 2, 3); //start,end,interior knots, start/end multiplicity
    gsMatrix<> coefs_c(5, 2); //ξ,η;..
    coefs_c << 0, 0,
             0.2618, 0.053,
             0.5, 0.5,
             0.738, 0.9465,
             1,1;

    coefs_c.col(0) = coefs_c.col(0) * kv_ξ.last();
    coefs_c.col(1) = coefs_c.col(1) * kv_η.last();

    gsBSpline<> curve_param(kv_c, coefs_c);

    // Print the embedded entity
    gsInfo << "I am the embedded entity in the superdomain parameter space" << "\n";
    gsInfo << "Knot vector:" << curve_param.knots() << "\n";
    gsInfo << "Control points:" << "\n";
    gsInfo << curve_param.coefs() << "\n";
    
    // Search for θ values at the crossing between the embedded curve and the surface
    // isoparametric line at ξ = 0.50
    gsMatrix<> θ(1,1);
    θ << 0.25; //initial guess
    index_t maxiter = 50;
    index_t i = 2;
    real_t epsilon = 1e-8;
    short_t dir = 0;

    gsMatrix<> θ_crossing_ξ = findRoot(curve_param, kv_ξ, θ, i, dir, maxiter, epsilon);

    gsInfo << θ_crossing_ξ << "\n";

    //// gsInfo<<curveCoordinate(coefs);

    return 0;
}

template <class T>
gsMatrix<T> findRoot(const gsGeometry<T> &geometry, gsKnotVector<T> &kv, const gsMatrix<T> &value, const index_t i=0,
            const short_t dir=0, const index_t maxiter=50, const real_t epsilon=1e-4)
{
     GISMO_ASSERT(geometry.parDim() == 1, "The embedded geometry must be a curve (parameter dimension 1)");

     gsMatrix<T> xn = value; 
     for (index_t n=0; n < maxiter; n++ )
     {
         real_t Dfxn = geometry.deriv(xn)(dir,0);
         if (Dfxn == 0)
         {
             gsInfo <<"Zero derivative: no crossing found.\n";
             xn << math::sqrt(-1);
             return xn;
         }
         real_t fxn = geometry.eval(xn)(dir,0) - kv(i);
         if (abs(fxn) <= epsilon)
         {
             gsInfo <<"Crossing found after "<< n << " iterations.\n";
             return xn;
         }
         xn(0,0) = xn(0,0) - fxn/Dfxn;
     }
     gsInfo << "Max number of iterations exceeded: no crossing found.\n";
     xn << math::sqrt(-1);
     return xn;
};