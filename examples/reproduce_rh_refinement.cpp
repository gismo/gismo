/** @file reproduce_rh_refinement.cpp

    @brief rh-refinement study (replaces the shell example): combines the
           r-adaptivity of this paper (control-point relocation via the composed
           map G = S o sigma) with h-adaptivity (local THB-spline refinement),
           on the Section-4.3 spherical-patch Poisson problem with a sharp
           tanh boundary layer.

    Four strategies are compared on identical error-vs-DoF axes:
        (h)    uniform    : identity sigma, tensor-product basis, uniform refine
        (H)    adaptive-h : identity sigma, THB basis, adaptive (residual) refine
        (r)    relocate   : reparameterised sigma, tensor basis, uniform refine
        (rh)   combined   : reparameterised sigma, THB basis, adaptive refine

    The relocated map concentrates the parametric resolution at the layer a
    priori; THB then adds local degrees of freedom where the discrete error
    remains large. The comparison shows rh reaches a target accuracy with the
    fewest degrees of freedom.

    Usage:
        reproduce_rh_refinement -i <sphere_patch.xml> [-p 2] [-L 6]

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Ye Ji (y.ji-1@tudelft.nl)
*/

#include <gismo.h>
#include <gsNurbs/gsSquareDomain.h>
#include <gsOptim/gsOptim.h>
#include <gsAssembler/gsAdaptiveParametrization.h>
#include <gsAssembler/gsAdaptiveRefUtils.h>

using namespace gismo;

// Manufactured layer  u = tanh((1/2 - v)/t) + 1  at polar angle v = 1/2, width t
// (kappa = 1/t). Smaller t = sharper layer = stronger case for h-adaptivity.
static const std::string VSTR = "atan2(sqrt(x^2+y^2),z)";
static std::string g_t = "0.05";
static std::string exactStr() { return "tanh((1/2-("+VSTR+"))/"+g_t+") + 1"; }
static std::string sourceStr()
{
    return "u:=atan2(y,x); v:="+VSTR+"; w:=sqrt(x^2+y^2+z^2); t:="+g_t+";"
           "-(-cos(v)*(1 - tanh((1/2 - v)/t)^2)/t - 2*sin(v)*tanh((1/2 - v)/t)*(1 - tanh((1/2 - v)/t)^2)/t^2)/(w^2*sin(v))";
}

// Relocated (r-adaptive) or identity composition sigma of bi-degree (d,d), NxN.
static gsGeometry<real_t>::uPtr makeSigma(const gsGeometry<real_t>& S, index_t d, index_t N,
                                          bool reparam, real_t smoothing, real_t penalty)
{
    gsGeometry<>::uPtr comp = gsNurbsCreator<>::BSplineSquareDeg(d);
    const index_t nref=(index_t)std::round(std::log2((double)N));
    for (index_t r=0;r<nref;++r) comp->uniformRefine();
    gsSquareDomain<real_t> dom(*comp);
    dom.options().addSwitch("Slide","",false); dom.applyOptions();   // fixed boundary
    if (reparam && dom.nControls()>0)
    {
        gsFunctionExpr<> indicator(exactStr(), S.targetDim());
        gsOptim<real_t>::LBFGS opt; opt.options().setInt("MaxIterations",250);
        gsAdaptiveParametrization<real_t,MonitorMode::GradientBased> relo(dom,S,indicator,opt,false);
        relo.options().setReal("Smoothing",smoothing); relo.options().setReal("Penalty",penalty);
        relo.solve();
        const real_t mJ = relo.computeMinJacobian();
        GISMO_ENSURE(math::isfinite(mJ)&&mJ>0, "reparam folded (minJ="<<mJ<<")");
        gsInfo<<"    [relocated sigma: min det Jsigma = "<<mJ<<"]\n";
    }
    return dom.domain().clone();
}

// Assemble + solve the surface Poisson problem for composed geometry cmp and analysis mb.
// Returns L2, H1 error and (optionally) the element-wise squared L2 error for marking.
static real_t g_quA = 1.0;   // quadrature over-integration for the composed integrand
static void solveOnce(const gsMultiPatch<real_t>& cmp, gsMultiBasis<real_t>& mb,
                      real_t& l2, real_t& h1, index_t& dofs,
                      std::vector<real_t>* elErr)
{
    const index_t dim = cmp.patch(0).targetDim();
    gsFunctionExpr<> msFun(exactStr(),dim), fFun(sourceStr(),dim);
    gsComposedFunction<real_t> cms(cmp.patch(0), msFun), cf(cmp.patch(0), fFun);

    gsBoundaryConditions<> bc;
    for (short_t s=1;s<=4;++s) bc.addCondition(s, condition_type::dirichlet, &cms, 0, true);
    bc.setGeoMap(cmp);

    gsExprAssembler<> A(1,1);
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    A.setIntegrationElements(mb);                    // quadrature on the analysis elements
    A.options().setReal("quA",g_quA); A.options().setInt("quB",2);
    gsOptionList seOff; seOff.addSwitch("SameElement","",false);
    A.options().update(seOff, gsOptionList::addIfUnknown);   // required for composed geometry

    gsExprEvaluator<> ev(A);
    ev.options().setReal("quA",g_quA); ev.options().setInt("quB",2);
    ev.options().update(seOff, gsOptionList::addIfUnknown);

    geometryMap G = A.getMap(cmp);
    space u = A.getSpace(mb);
    auto fsrc = A.getCoeff(cf);
    auto uex  = ev.getVariable(cms);

    gsMatrix<real_t> solVector;
    solution usol = A.getSolution(u, solVector);

    u.setup(bc, dirichlet::l2Projection, 0);
    A.initSystem();
    A.assemble( igrad(u,G)*igrad(u,G).tr()*meas(G),  u*fsrc*meas(G) );
    // SPD stiffness -> conjugate gradients (far faster than QR as the DoF grow).
    gsSparseSolver<>::CGDiagonal solver(A.matrix());
    solver.setTolerance(1e-12); solver.setMaxIterations(10000);
    solVector = solver.solve(A.rhs());
    GISMO_ENSURE(solVector.allFinite(), "linear solve produced non-finite values");

    l2 = math::sqrt( ev.integral( (uex-usol).sqNorm()*meas(G) ) );
    h1 = l2 + math::sqrt( ev.integral( (igrad(uex,G)-igrad(usol,G)).sqNorm()*meas(G) ) );
    dofs = A.numDofs();

    if (elErr)
    {
        ev.integralElWise( (uex-usol).sqNorm()*meas(G) );
        *elErr = ev.elementwise();
    }
}

int main(int argc, char** argv)
{
    std::string input = "parametrization/monitor_example_sphere_patch_with_gradient_function.xml";
    index_t p = 2, nLoops = 6, compDeg = 2, compN = 4;
    real_t smoothing = 1e-2, penalty = 1e-4, refParam = 0.85;
    index_t refExt = 1, refCrit = 3;   // 3 = Doerfler/bulk marking
    index_t adaptMax = 40, targetDof = 5000;  // adaptive runs until DoF>=targetDof (or adaptMax steps)
    std::string tLayer = "0.05";
    gsCmdLine cmd("rh-refinement: uniform-h vs THB-h vs r vs rh on the Sec.4.3 Poisson layer.");
    cmd.addString("i","input","Sphere-patch xml",input);
    cmd.addString("t","layer","Layer width t=1/kappa (smaller=sharper)",tLayer);
    cmd.addReal("q","quA","Quadrature over-integration for composed integrand",g_quA);
    cmd.addInt ("p","degree","Analysis degree",p);
    cmd.addInt ("L","loops","Refinement steps per strategy",nLoops);
    cmd.addInt ("d","compDegree","Composition sigma degree",compDeg);
    cmd.addInt ("N","compN","Composition sigma elements per direction",compN);
    cmd.addReal("S","smoothing","Monitor smoothing",smoothing);
    cmd.addReal("P","penalty","Fold-barrier penalty",penalty);
    cmd.addReal("m","markParam","Marking parameter",refParam);
    cmd.addInt ("D","targetDof","Adaptive strategies run until this many DoF",targetDof);
    cmd.addInt ("","adaptMax","Max adaptive refinement steps",adaptMax);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    g_t = tLayer;
    GISMO_ENSURE(gsFileManager::fileExists(input), "Input not found: "<<input);
    gsInfo<<"Layer width t = "<<g_t<<" (kappa = 1/t)\n";

    gsFileData<> fd(input); gsMultiPatch<> geomMP;
    if (fd.hasLabel("geometry")) fd.getLabel("geometry",geomMP);
    else                         geomMP.addPatch(*fd.getFirst<gsGeometry<>>());
    const gsGeometry<>& S = geomMP.patch(0);

    // identity and relocated compositions
    gsGeometry<>::uPtr sigId = makeSigma(S,compDeg,compN,false,smoothing,penalty);
    gsInfo<<"Relocating sigma (r-adaptivity)...\n";
    gsGeometry<>::uPtr sigR  = makeSigma(S,compDeg,compN,true ,smoothing,penalty);

    // analysis-space factory: fresh tensor / THB basis on [0,1]^2
    gsKnotVector<> kv(0,1,3,p+1);    // start with a few elements
    gsTensorBSplineBasis<2> tp(kv,kv);

    struct Run { std::string tag; const gsGeometry<>* sigma; bool thb; };
    std::vector<Run> runs = {
        {"h  uniform ", sigId.get(), false},
        {"H  THB-adap", sigId.get(), true },
        {"r  reloc   ", sigR.get(),  false},
        {"rh combined", sigR.get(),  true }
    };

    gsInfo<<"\n"<<std::left<<std::setw(13)<<"strategy"<<std::setw(8)<<"step"
          <<std::setw(10)<<"DoF"<<std::setw(14)<<"L2"<<std::setw(14)<<"H1"<<"\n";
    gsInfo<<std::string(59,'-')<<"\n";

    for (const Run& run : runs)
    {
        gsMultiPatch<> cmp; cmp.addPatch(gsComposedGeometry<real_t>(*run.sigma, S));
        gsMultiBasis<> mb;
        if (run.thb) mb.addBasis(gsTHBSplineBasis<2,real_t>(tp).clone().release());
        else         mb.addBasis(tp.clone().release());

        // Uniform strategies quadruple the DoF per step, so nLoops steps suffice.
        // Adaptive strategies add few DoF per step, so we let them run until they reach
        // a comparable DoF range (targetDof) -- otherwise their curves stop far short of
        // the uniform ones and look truncated.
        const index_t maxSteps = run.thb ? adaptMax : nLoops;
        for (index_t step=0; step<maxSteps; ++step)
        {
            real_t l2,h1; index_t dofs; std::vector<real_t> elErr;
            solveOnce(cmp, mb, l2, h1, dofs, run.thb ? &elErr : nullptr);
            gsInfo<<std::left<<std::setw(13)<<run.tag<<std::setw(8)<<step
                  <<std::setw(10)<<dofs<<std::scientific<<std::setprecision(4)
                  <<std::setw(14)<<l2<<std::setw(14)<<h1<<std::fixed<<"\n";

            if (!run.thb && step+1==nLoops) break;
            if ( run.thb && (dofs>=targetDof || step+1==maxSteps)) break;
            if (run.thb)
            {
                std::vector<bool> marked;
                gsMarkElementsForRef(elErr, refCrit, refParam, marked);
                gsRefineMarkedElements(mb, marked, refExt);
            }
            else mb.uniformRefine();
        }
        gsInfo<<std::string(59,'-')<<"\n";
    }
    gsInfo<<"\nCompare L2/H1 vs DoF across strategies: rh should reach a given error with fewest DoF.\n";
    return 0;
}
