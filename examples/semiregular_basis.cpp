
#include <iostream>
#include <gismo.h>

#include "gsSemiRegularBasis.h"

using namespace gismo;

int main(int argc, char* argv[])
{
    index_t p = 2, l = 2, s = 3;
    bool plot = false;
    gsCmdLine cmd("Tutorial on gsBasis class.");
    cmd.addInt   ("p", "degree", "Degree", p);
    cmd.addInt   ("l", "layers", "layers", l);
    cmd.addInt   ("s", "samples", "number of sample points", s);
    cmd.addSwitch("plot"   , "Plot the result", plot);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // ======================================================================
    // making the basis
    // ======================================================================

    gsSemiRegularBasis<2> BB(p, l);

    // ======================================================================
    // printing some information about the basis
    // ======================================================================

    // printing the basis
    gsInfo << "The basis: \n" << BB << "\n";


    // printing some properties of the basis
    gsInfo << "Dimension of the parameter space: " << BB.dim() << "\n"
           << "Number of basis functions: "        << BB.size() << "\n"
           << "\n";


    // support of the basis
    // (dim x 2 matrix, the parametric domain)
    //gsMatrix<> support = BB.support();


    // ======================================================================
    // basis evaluation
    // ======================================================================
    gsVector<> a(2); a.setZero();
    gsVector<> b(2); b.setOnes();
    gsVector<unsigned> np = uniformSampleCount(a,b, s);
    gsMatrix<> u = gsPointGrid(a, b, np);
    if (s<100)
        gsInfo << "Points (" << u.cols() << "): \n" << u << "\n\n";

    // indices of active (nonzero) functions at parameter u
    gsMatrix<index_t> active = BB.active(u);
    if (s<100)
        gsInfo << "#Active basis functions at u: \n"
               << active << "\n\n";

    // values of all active functions at u
    gsMatrix<> values = BB.eval(u);
    if (s<100)
        gsInfo << "Values at u : \n"
               << values << "\n\n";

    // values of single basis functions
    if (s<100)
    {
        for (index_t i = 0; i != active.rows(); i++)
        {
            gsMatrix<> val = BB.evalSingle(active(i), u);
            
            gsInfo << "basis fun. index:  " << active(i)
                   << "   value: " << val(0, 0) << "\n";
        }
        gsInfo << "\n";
    }

    // ----------------------------------------------------------------------
    // derivatives
    // ----------------------------------------------------------------------


    gsMatrix<> derivs = BB.deriv(u);
    if (s<100)
        gsInfo << "Derivatives at u : \n"
               << derivs << "\n\n";

    // ----------------------------------------------------------------------
    // second derivatives
    // ----------------------------------------------------------------------

    gsMatrix<> derivs2 = BB.deriv2(u);
    if (s<100)
        gsInfo << "Second derivatives at u : \n"
               << derivs2 << "\n\n";

    gsInfo << "\nFor more information about evaluation "
           << "(and order of derivatives) look at doxygen documentation."
           << "\n\n";

    // ----------------------------------------------------------------------
    // fitting
    // ----------------------------------------------------------------------

    gsFunctionExpr<real_t> ff("x", "x*(1-y)", "x*y^2", 2);
    gsMatrix<> fu = ff.eval(u);
    gsKnotVector<> knot_v(0,1,0,4); //multiplicity at the end
    gsTensorBSplineBasis<2> tp_bezier(knot_v,knot_v);
    auto control_points = tp_bezier.interpolateAtAnchors(ff.eval(tp_bezier.anchors()));

    gsInfo<<"Control points are = " <<control_points->coefs()<<"\n";


    
    gsFitting<> fit(u, fu, BB);
    fit.compute();
    gsInfo << *fit.result() <<"\n";

    auto res = fit.result()->eval(u);

    auto error_L2 = ( fu - res ).norm();

    gsInfo<<"Error is: "<<error_L2<< "\n";
    
    if (plot)
    {
        gsWriteParaview(BB, "srbasis", 10000);
        gsWriteParaviewPoints(fu, "srpoints");
        gsWriteParaview(*fit.result(), "srfitting", 10000);
        //gsFileManager::open("srbasis.pvd");
    }

    return EXIT_SUCCESS;
}
