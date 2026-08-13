/** @file export_rh_mesh.cpp

    @brief Exports Paraview (.vts + _mesh.vtp) files of the adapted THB meshes for
           the anisotropic curved-ridge problem of reproduce_rh_anisotropic.cpp,
           for visualizing the hierarchical structure:
             thb_mesh*  : identity sigma  (isotropic THB refinement)
             rh_mesh*   : relocated sigma (THB on the sigma-aligned mesh)
           Both are refined adaptively until ~targetDof degrees of freedom and
           written in physical coordinates (the mesh is the THB solution basis
           mapped through G = S o sigma).
*/

#include <gismo.h>
#include <gsNurbs/gsSquareDomain.h>
#include <gsOptim/gsOptim.h>
#include <gsAssembler/gsAdaptiveParametrization.h>
#include <gsAssembler/gsAdaptiveRefUtils.h>

using namespace gismo;

static std::string g_t = "0.03";
static std::string exactStr(){ return "tanh((y-(1/2+0.3*sin(pi*x)))/"+g_t+")+1"; }
static std::string sourceStr()
{
    const std::string T=g_t;
    const std::string s="((y-(1/2+0.3*sin(pi*x)))/"+T+")";
    const std::string th="tanh"+s, se="(1-"+th+"^2)";
    const std::string mp="(0.3*pi*cos(pi*x))", mpp="(-0.3*pi^2*sin(pi*x))";
    return "(2/("+T+")^2)*"+se+"*"+th+"*(1+"+mp+"^2) + ("+mpp+"/"+T+")*"+se;
}

static gsGeometry<real_t>::uPtr makeSigma(const gsGeometry<real_t>& S, index_t d, index_t N, bool reparam)
{
    gsGeometry<>::uPtr comp = gsNurbsCreator<>::BSplineSquareDeg(d);
    const index_t nref=(index_t)std::round(std::log2((double)N));
    for (index_t r=0;r<nref;++r) comp->uniformRefine();
    gsSquareDomain<real_t> dom(*comp);
    dom.options().addSwitch("Slide","",true); dom.applyOptions();
    if (reparam && dom.nControls()>0)
    {
        gsFunctionExpr<> indicator(exactStr(), S.targetDim());
        gsOptim<real_t>::LBFGS opt; opt.options().setInt("MaxIterations",250);
        gsAdaptiveParametrization<real_t,MonitorMode::GradientBased> relo(dom,S,indicator,opt,false);
        relo.options().setReal("Smoothing",1e-2); relo.options().setReal("Penalty",1e-2);
        relo.solve();
        gsInfo<<"  [relocated sigma: min det Jsigma = "<<relo.computeMinJacobian()<<"]\n";
    }
    return dom.domain().clone();
}

// Solve on composed geometry cmp with THB basis mb; returns solVector, dofs, elErr.
static gsMatrix<real_t> solveOnce(const gsMultiPatch<real_t>& cmp, gsMultiBasis<real_t>& mb,
                                  index_t& dofs, std::vector<real_t>& elErr, real_t quA)
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
    A.options().setReal("quA",quA); A.options().setInt("quB",2);
    gsOptionList seOff; seOff.addSwitch("SameElement","",false);
    A.options().update(seOff, gsOptionList::addIfUnknown);
    gsExprEvaluator<> ev(A);
    ev.options().setReal("quA",quA); ev.options().setInt("quB",2);
    ev.options().update(seOff, gsOptionList::addIfUnknown);

    geometryMap G = A.getMap(cmp);
    space u = A.getSpace(mb);
    auto fsrc = A.getCoeff(cf);
    auto uex  = ev.getVariable(cms);
    gsMatrix<real_t> solVector;
    solution usol = A.getSolution(u, solVector);

    u.setup(bc, dirichlet::l2Projection, 0);
    A.initSystem();
    A.assemble( igrad(u,G)*igrad(u,G).tr()*meas(G), u*fsrc*meas(G) );
    gsSparseSolver<>::CGDiagonal solver(A.matrix());
    solver.setTolerance(1e-12); solver.setMaxIterations(10000);
    solVector = solver.solve(A.rhs());
    dofs = A.numDofs();
    ev.integralElWise( (uex-usol).sqNorm()*meas(G) );
    elErr = ev.elementwise();
    return solVector;
}

int main(int argc, char** argv)
{
    index_t compDeg=2, compN=8, p=2, targetDof=1200, refExt=1, refCrit=3;
    real_t refParam=0.7, quA=2.0;
    std::string tLayer="0.03";
    gsCmdLine cmd("Export adapted THB meshes (THB vs rh) for the anisotropic ridge.");
    cmd.addString("t","layer","Layer width",tLayer);
    cmd.addInt("D","targetDof","Refine until this many DoF",targetDof);
    cmd.addReal("m","markParam","Marking parameter",refParam);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    g_t=tLayer;

    gsGeometry<>::uPtr S = gsNurbsCreator<>::BSplineSquare(1.0,0.0,0.0);
    gsGeometry<>::uPtr sigId = makeSigma(*S,compDeg,compN,false);
    gsInfo<<"Relocating sigma...\n";
    gsGeometry<>::uPtr sigR  = makeSigma(*S,compDeg,compN,true );

    gsKnotVector<> kv(0,1,3,p+1);
    gsTensorBSplineBasis<2> tp(kv,kv);

    struct Case { std::string name; const gsGeometry<>* sigma; };
    for (Case c : { Case{"thb_mesh", sigId.get()}, Case{"rh_mesh", sigR.get()} })
    {
        gsMultiPatch<> cmp; cmp.addPatch(gsComposedGeometry<real_t>(*c.sigma, *S));
        gsMultiBasis<> mb; mb.addBasis(gsTHBSplineBasis<2,real_t>(tp).clone().release());

        index_t dofs=0; gsMatrix<real_t> solVec; std::vector<real_t> elErr;
        for (index_t step=0; step<40; ++step)
        {
            solVec = solveOnce(cmp, mb, dofs, elErr, quA);
            if (dofs>=targetDof) break;
            std::vector<bool> marked;
            gsMarkElementsForRef(elErr, refCrit, refParam, marked);
            gsRefineMarkedElements(mb, marked, refExt);
        }

        // Reconstruct the solution as a multipatch on the THB basis and write it,
        // together with its mesh, in physical (G) coordinates.
        gsExprAssembler<> A(1,1);
        A.setIntegrationElements(mb);
        typedef gsExprAssembler<>::space space;
        space u = A.getSpace(mb);
        gsFunctionExpr<> msFun(exactStr(), cmp.patch(0).targetDim());
        gsComposedFunction<real_t> cms(cmp.patch(0), msFun);
        gsBoundaryConditions<> bc;
        for (short_t s=1;s<=4;++s) bc.addCondition(s, condition_type::dirichlet, &cms, 0, true);
        bc.setGeoMap(cmp);
        u.setup(bc, dirichlet::l2Projection, 0);
        gsExprAssembler<>::solution usol = A.getSolution(u, solVec);
        gsMultiPatch<> solMP; usol.extract(solMP);

        gsField<> field(cmp, solMP);
        gsWriteParaview<>(field, c.name, 60000, true);   // .vts (field) + _mesh.vtp
        gsInfo<<c.name<<": "<<dofs<<" DoF, "<<mb.basis(0).numElements()<<" elements -> "
              <<c.name<<".vts (+ _mesh.vtp)\n";
    }
    return 0;
}
