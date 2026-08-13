/** @file reproduce_rh_anisotropic.cpp

    @brief Anisotropic-feature companion to reproduce_rh_refinement.cpp.

    Isotropic (THB) h-refinement resolves an anisotropic feature inefficiently:
    it can only split an element into equal squares, never stretch it along the
    feature. Control-point relocation (r-adaptivity), by contrast, can align and
    stretch the mesh along the feature. This example makes that concrete: a sharp,
    CURVED ridge

        u = tanh( (y - m(x))/t ) + 1 ,   m(x) = 1/2 + 0.3*sin(pi*x),

    solved as a Poisson problem on the unit square. Because the ridge is curved,
    a low-DoF sigma can only align the mesh approximately; the residual is then
    resolved by THB. Four strategies are compared on error-vs-DoF axes:
        (h)  uniform    : identity sigma, tensor basis, uniform refine
        (H)  THB        : identity sigma, THB basis, adaptive refine
        (r)  relocate   : reparameterised sigma, tensor basis, uniform refine
        (rh) combined   : reparameterised sigma, THB basis, adaptive refine

    Usage: reproduce_rh_anisotropic [-t 0.03] [-p 2] [-L 6]
*/

#include <gismo.h>
#include <gsNurbs/gsSquareDomain.h>
#include <gsOptim/gsOptim.h>
#include <gsAssembler/gsAdaptiveParametrization.h>
#include <gsAssembler/gsAdaptiveRefUtils.h>

using namespace gismo;

static std::string g_t = "0.03";
// Curved anisotropic ridge: sharp across y = m(x), smooth along it.
static std::string exactStr()
{ return "tanh((y-(1/2+0.3*sin(pi*x)))/"+g_t+")+1"; }
static std::string sourceStr()
{
    // f = -Delta u = (2/t^2) sech^2(s) tanh(s) (1+m'^2) + (m''/t) sech^2(s),
    // s=(y-m)/t, m=1/2+0.3 sin(pi x), m'=0.3 pi cos(pi x), m''=-0.3 pi^2 sin(pi x).
    // Fully inlined (no local variables) for robust exprtk parsing.
    const std::string T  = g_t;
    const std::string s  = "((y-(1/2+0.3*sin(pi*x)))/"+T+")";
    const std::string th = "tanh"+s;
    const std::string se = "(1-"+th+"^2)";
    const std::string mp = "(0.3*pi*cos(pi*x))";
    const std::string mpp= "(-0.3*pi^2*sin(pi*x))";
    return "(2/("+T+")^2)*"+se+"*"+th+"*(1+"+mp+"^2) + ("+mpp+"/"+T+")*"+se;
}

static real_t g_quA = 2.0;

static gsGeometry<real_t>::uPtr makeSigma(const gsGeometry<real_t>& S, index_t d, index_t N,
                                          bool reparam, real_t smoothing, real_t penalty)
{
    gsGeometry<>::uPtr comp = gsNurbsCreator<>::BSplineSquareDeg(d);
    const index_t nref=(index_t)std::round(std::log2((double)N));
    for (index_t r=0;r<nref;++r) comp->uniformRefine();
    gsSquareDomain<real_t> dom(*comp);
    dom.options().addSwitch("Slide","",true); dom.applyOptions();   // sliding boundary
    if (reparam && dom.nControls()>0)
    {
        gsFunctionExpr<> indicator(exactStr(), S.targetDim());
        gsOptim<real_t>::LBFGS opt; opt.options().setInt("MaxIterations",250);
        gsAdaptiveParametrization<real_t,MonitorMode::GradientBased> relo(dom,S,indicator,opt,false);
        relo.options().setReal("Smoothing",smoothing); relo.options().setReal("Penalty",penalty);
        const gsVector<real_t> before = dom.getControls();
        relo.solve();
        const bool moved = (dom.getControls()-before).norm() > 1e-12;
        const real_t mJ = relo.computeMinJacobian();
        GISMO_ENSURE(math::isfinite(mJ)&&mJ>0, "reparam folded (minJ="<<mJ<<")");
        if (!moved) gsWarn<<"  relocation did not move; try a larger --penalty.\n";
        gsInfo<<"    [relocated sigma: min det Jsigma = "<<mJ<<"]\n";
    }
    return dom.domain().clone();
}

static void solveOnce(const gsMultiPatch<real_t>& cmp, gsMultiBasis<real_t>& mb,
                      real_t& l2, real_t& h1, index_t& dofs, std::vector<real_t>* elErr)
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

    A.setIntegrationElements(mb);
    A.options().setReal("quA",g_quA); A.options().setInt("quB",2);
    gsOptionList seOff; seOff.addSwitch("SameElement","",false);
    A.options().update(seOff, gsOptionList::addIfUnknown);
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
    gsSparseSolver<>::CGDiagonal solver(A.matrix());
    solver.setTolerance(1e-12); solver.setMaxIterations(10000);
    solVector = solver.solve(A.rhs());
    GISMO_ENSURE(solVector.allFinite(), "non-finite solve");

    l2 = math::sqrt( ev.integral( (uex-usol).sqNorm()*meas(G) ) );
    h1 = l2 + math::sqrt( ev.integral( (igrad(uex,G)-igrad(usol,G)).sqNorm()*meas(G) ) );
    dofs = A.numDofs();
    if (elErr) { ev.integralElWise( (uex-usol).sqNorm()*meas(G) ); *elErr = ev.elementwise(); }
}

int main(int argc, char** argv)
{
    index_t p = 2, nLoops = 6, compDeg = 2, compN = 8;
    real_t smoothing = 1e-2, penalty = 1e-2, refParam = 0.7;
    index_t refExt = 1, refCrit = 3, adaptMax = 40, targetDof = 5000;
    std::string tLayer = "0.03";
    gsCmdLine cmd("rh-refinement on an anisotropic curved ridge (planar Poisson).");
    cmd.addString("t","layer","Layer width t",tLayer);
    cmd.addReal ("q","quA","Quadrature over-integration",g_quA);
    cmd.addInt  ("p","degree","Analysis degree",p);
    cmd.addInt  ("L","loops","Uniform refinement steps",nLoops);
    cmd.addInt  ("d","compDegree","Composition degree",compDeg);
    cmd.addInt  ("N","compN","Composition elements/direction",compN);
    cmd.addReal ("S","smoothing","Monitor smoothing",smoothing);
    cmd.addReal ("P","penalty","Fold-barrier penalty",penalty);
    cmd.addReal ("m","markParam","Marking parameter",refParam);
    cmd.addInt  ("D","targetDof","Adaptive DoF target",targetDof);
    cmd.addInt  ("","adaptMax","Max adaptive steps",adaptMax);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    g_t = tLayer;
    gsInfo<<"Anisotropic curved ridge, layer width t = "<<g_t<<"\n";

    // Planar unit square as the "surface" S (identity 2->2 map).
    gsGeometry<>::uPtr S = gsNurbsCreator<>::BSplineSquare(1.0, 0.0, 0.0);

    gsGeometry<>::uPtr sigId = makeSigma(*S,compDeg,compN,false,smoothing,penalty);
    gsInfo<<"Relocating sigma (r-adaptivity)...\n";
    gsGeometry<>::uPtr sigR  = makeSigma(*S,compDeg,compN,true ,smoothing,penalty);

    gsKnotVector<> kv(0,1,3,p+1);
    gsTensorBSplineBasis<2> tp(kv,kv);

    struct Run { std::string tag; const gsGeometry<>* sigma; bool thb; };
    std::vector<Run> runs = {
        {"h  uniform ", sigId.get(), false},
        {"H  THB-adap", sigId.get(), true },
        {"r  reloc   ", sigR.get(),  false},
        {"rh combined", sigR.get(),  true }
    };

    gsInfo<<"\n"<<std::left<<std::setw(13)<<"strategy"<<std::setw(8)<<"step"
          <<std::setw(10)<<"DoF"<<std::setw(14)<<"L2"<<std::setw(14)<<"H1"<<"\n"
          <<std::string(59,'-')<<"\n";

    for (const Run& run : runs)
    {
        gsMultiPatch<> cmp; cmp.addPatch(gsComposedGeometry<real_t>(*run.sigma, *S));
        gsMultiBasis<> mb;
        if (run.thb) mb.addBasis(gsTHBSplineBasis<2,real_t>(tp).clone().release());
        else         mb.addBasis(tp.clone().release());

        const index_t maxSteps = run.thb ? adaptMax : nLoops;
        for (index_t step=0; step<maxSteps; ++step)
        {
            real_t l2,h1; index_t dofs; std::vector<real_t> elErr;
            solveOnce(cmp, mb, l2, h1, dofs, run.thb ? &elErr : nullptr);
            gsInfo<<std::left<<std::setw(13)<<run.tag<<std::setw(8)<<step<<std::setw(10)<<dofs
                  <<std::scientific<<std::setprecision(4)<<std::setw(14)<<l2<<std::setw(14)<<h1<<std::fixed<<"\n";
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
    return 0;
}
