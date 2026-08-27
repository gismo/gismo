/** @file export_transient_snapshots.cpp

    @brief Exports Paraview (.vts + _mesh.vtp) snapshots of the r-adaptive (ALE)
           transient solution and the relocated mesh at selected time instants, for
           the moving-colatitude-band heat equation on the exact spherical patch
           (Sec.4.5). Used to build the time-series figure showing that the
           relocated mesh tracks the moving front.

    At each requested step k, the composed geometry G = S o sigma(t_k) (the
    relocated mesh) is written together with the solution field u^k, so a renderer
    can draw both the solution (coloured surface) and the mesh at that instant.

    Usage:
        export_transient_snapshots -i <sphere_patch.xml> [-t 0.05] [-K 16] [-N 24]
                                   [--compN 16] [--c0 0.35] [--c1 0.65] [--snaps 4,8,12,16]

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

static const std::string VSTR = "atan2(sqrt(x^2+y^2),z)";
static std::string g_t="0.05", g_c="0.35", g_cdot="0.30";
static std::string exactStr(){ return "tanh(("+g_c+"-("+VSTR+"))/"+g_t+")+1"; }
static std::string fStr()
{
    const std::string V ="(atan2(sqrt(x^2+y^2),z))";
    const std::string W2="(x^2+y^2+z^2)";
    const std::string Tc="tanh(("+g_c+"-"+V+")/"+g_t+")";
    const std::string g2="(1-"+Tc+"^2)";
    const std::string dudt="("+g2+"*"+g_cdot+"/"+g_t+")";
    const std::string mlap="( -( -cos"+V+"*"+g2+"/"+g_t
                          +" - 2*sin"+V+"*"+Tc+"*"+g2+"/"+g_t+"^2 )/("+W2+"*sin"+V+") )";
    return dudt+" + "+mlap;
}
static real_t g_quA = 2.0;

// Physical mesh velocity w(xi) = (G_new(xi) - G_prev(xi))/Dt in R^3.
class MeshVel : public gsFunction<real_t>
{
    const gsGeometry<real_t>& m_Gnew; const gsGeometry<real_t>& m_Gprev; real_t m_Dt;
public:
    MeshVel(const gsGeometry<real_t>& Gnew, const gsGeometry<real_t>& Gprev, real_t Dt)
        : m_Gnew(Gnew), m_Gprev(Gprev), m_Dt(Dt) {}
    GISMO_CLONE_FUNCTION(MeshVel)
    short_t domainDim() const { return 2; }
    short_t targetDim() const { return 3; }
    void eval_into(const gsMatrix<real_t>& u, gsMatrix<real_t>& res) const
    { gsMatrix<real_t> a,b; m_Gnew.eval_into(u,a); m_Gprev.eval_into(u,b); res=(a-b)/m_Dt; }
};

static gsGeometry<real_t>::uPtr makeIdentity(index_t d, index_t N)
{
    gsGeometry<>::uPtr comp = gsNurbsCreator<>::BSplineSquareDeg(d);
    const index_t nref=(index_t)std::round(std::log2((double)N));
    for (index_t r=0;r<nref;++r) comp->uniformRefine();
    return comp;
}
static gsGeometry<real_t>::uPtr relocateFrom(const gsGeometry<real_t>& init, const gsGeometry<real_t>& S,
                                             real_t smoothing, real_t penalty, index_t maxIt)
{
    gsGeometry<>::uPtr comp = init.clone();
    gsSquareDomain<real_t> dom(*comp);
    dom.options().addSwitch("Slide","",false); dom.applyOptions();
    if (dom.nControls()>0)
    {
        gsFunctionExpr<> indicator(exactStr(), S.targetDim());
        gsOptim<real_t>::LBFGS opt; opt.options().setInt("MaxIterations",maxIt);
        gsAdaptiveParametrization<real_t,MonitorMode::GradientBased> relo(dom,S,indicator,opt,false);
        relo.options().setReal("Smoothing",smoothing); relo.options().setReal("Penalty",penalty);
        relo.solve();
        const real_t mJ = relo.computeMinJacobian();
        GISMO_ENSURE(math::isfinite(mJ)&&mJ>0, "reloc folded (minJ="<<mJ<<")");
    }
    return dom.domain().clone();
}
static void aleStep(const gsMultiPatch<real_t>& cmp, gsMultiBasis<real_t>& mb,
                    const gsFunction<real_t>* wFun, const gsMultiPatch<real_t>& uPrevField,
                    real_t Dt, gsMultiPatch<real_t>& uNextField)
{
    const index_t dim = cmp.patch(0).targetDim();
    gsFunctionExpr<> uexFun(exactStr(),dim), fFun(fStr(),dim);
    gsComposedFunction<real_t> cuex(cmp.patch(0), uexFun), cf(cmp.patch(0), fFun);
    gsBoundaryConditions<> bc;
    for (short_t s=1;s<=4;++s) bc.addCondition(s, condition_type::dirichlet, &cuex, 0, true);
    bc.setGeoMap(cmp);
    gsExprAssembler<> A(1,1);
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;
    A.setIntegrationElements(mb);
    A.options().setReal("quA",g_quA); A.options().setInt("quB",2);
    gsOptionList seOff; seOff.addSwitch("SameElement","",false);
    A.options().update(seOff, gsOptionList::addIfUnknown);
    geometryMap G = A.getMap(cmp);
    space u = A.getSpace(mb);
    auto fsrc  = A.getCoeff(cf);
    auto uprev = A.getCoeff(uPrevField.patch(0));
    gsMatrix<real_t> solVector;
    solution usol = A.getSolution(u, solVector);
    u.setup(bc, dirichlet::l2Projection, 0);
    A.initSystem();
    if (wFun) { auto w = A.getCoeff(*wFun);
        A.assemble( u*u.tr()*(1.0/Dt)*meas(G) + igrad(u,G)*igrad(u,G).tr()*meas(G)
                  - (igrad(u,G)*w)*u.tr()*meas(G),  u*fsrc*meas(G) + u*uprev*(1.0/Dt)*meas(G) );
    } else {
        A.assemble( u*u.tr()*(1.0/Dt)*meas(G) + igrad(u,G)*igrad(u,G).tr()*meas(G),
                    u*fsrc*meas(G) + u*uprev*(1.0/Dt)*meas(G) );
    }
    gsSparseSolver<>::CGDiagonal solver(A.matrix());
    solver.setTolerance(1e-12); solver.setMaxIterations(20000);
    solVector = solver.solve(A.rhs());
    GISMO_ENSURE(solVector.allFinite(), "non-finite solve");
    usol.extract(uNextField);
}
// THB basis (identity geometry) refined toward the current band by ||grad_S u_ex||,
// mirroring the THB feature-tracking strategy of reproduce_transient_sphere.cpp.
static gsMultiBasis<real_t> buildTHBforFront(const gsMultiPatch<real_t>& cmp,
                                             const gsTensorBSplineBasis<2>& base,
                                             index_t nPass, real_t refParam, index_t refExt, index_t refCrit)
{
    gsMultiBasis<> mb; mb.addBasis(gsTHBSplineBasis<2,real_t>(base).clone().release());
    const index_t dim = cmp.patch(0).targetDim();
    for (index_t pass=0; pass<nPass; ++pass){
        gsFunctionExpr<> uexFun(exactStr(),dim);
        gsComposedFunction<real_t> cuex(cmp.patch(0), uexFun);
        gsExprAssembler<> A(1,1); A.setIntegrationElements(mb);
        A.options().setReal("quA",g_quA); A.options().setInt("quB",2);
        gsOptionList seOff; seOff.addSwitch("SameElement","",false);
        A.options().update(seOff, gsOptionList::addIfUnknown);
        gsExprEvaluator<> ev(A); ev.options().setReal("quA",g_quA); ev.options().setInt("quB",2);
        ev.options().update(seOff, gsOptionList::addIfUnknown);
        auto Gm = A.getMap(cmp);
        auto uex = ev.getVariable(cuex);
        ev.integralElWise( igrad(uex,Gm).sqNorm()*meas(Gm) );
        std::vector<real_t> elErr = ev.elementwise();
        std::vector<bool> marked; gsMarkElementsForRef(elErr, refCrit, refParam, marked);
        gsRefineMarkedElements(mb, marked, refExt);
    }
    return mb;
}

static void projectExact(const gsMultiPatch<real_t>& cmp, gsMultiBasis<real_t>& mb, gsMultiPatch<real_t>& field)
{
    const index_t dim = cmp.patch(0).targetDim();
    gsFunctionExpr<> uexFun(exactStr(),dim);
    gsComposedFunction<real_t> cuex(cmp.patch(0), uexFun);
    gsBoundaryConditions<> bc;
    for (short_t s=1;s<=4;++s) bc.addCondition(s, condition_type::dirichlet, &cuex, 0, true);
    bc.setGeoMap(cmp);
    gsExprAssembler<> A(1,1);
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;
    A.setIntegrationElements(mb);
    A.options().setReal("quA",g_quA); A.options().setInt("quB",2);
    gsOptionList seOff; seOff.addSwitch("SameElement","",false);
    A.options().update(seOff, gsOptionList::addIfUnknown);
    geometryMap G = A.getMap(cmp); space u = A.getSpace(mb); auto uc = A.getCoeff(cuex);
    gsMatrix<real_t> solVector; solution usol = A.getSolution(u, solVector);
    u.setup(bc, dirichlet::l2Projection, 0); A.initSystem();
    A.assemble( u*u.tr()*meas(G), u*uc*meas(G) );
    gsSparseSolver<>::CGDiagonal solver(A.matrix()); solver.setTolerance(1e-14); solver.setMaxIterations(20000);
    solVector = solver.solve(A.rhs()); usol.extract(field);
}

int main(int argc, char** argv)
{
    std::string input = "parametrization/monitor_example_sphere_patch_with_gradient_function.xml";
    index_t compDeg=2, compN=16, Na=24, K=16, relocIt=40, p=2, thbPass=4;
    real_t smoothing=1e-2, penalty=1e-4, T=1.0, c0=0.35, c1=0.65;
    std::string tLayer="0.05", snaps="4,8,12,16";
    gsCmdLine cmd("Export r-ALE transient solution+mesh snapshots on the sphere patch.");
    cmd.addString("i","input","Sphere-patch xml",input);
    cmd.addString("t","layer","Front width",tLayer);
    cmd.addString("","snaps","Comma-separated step indices to export",snaps);
    cmd.addInt("K","steps","Time steps",K);
    cmd.addInt("N","Na","Fixed analysis mesh per direction",Na);
    cmd.addInt("","compN","Sigma mesh per direction",compN);
    cmd.addInt("","relocIt","Warm-start relocation iters per step",relocIt);
    cmd.addInt("","thbPass","THB feature-refinement passes per step",thbPass);
    cmd.addReal("","c0","Front start colatitude",c0);
    cmd.addReal("","c1","Front end colatitude",c1);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    g_t=tLayer;
    const real_t Dt=T/K, cdot=(c1-c0)/T;
    { std::ostringstream o; o<<cdot; g_cdot=o.str(); }
    GISMO_ENSURE(gsFileManager::fileExists(input), "Input not found: "<<input);
    std::set<index_t> want;
    { std::stringstream ss(snaps); std::string tok; while(std::getline(ss,tok,',')) if(!tok.empty()) want.insert(std::stoi(tok)); }

    gsFileData<> fd(input); gsMultiPatch<> geomMP;
    if (fd.hasLabel("geometry")) fd.getLabel("geometry",geomMP);
    else                         geomMP.addPatch(*fd.getFirst<gsGeometry<>>());
    gsGeometry<>::uPtr S = geomMP.patch(0).clone();

    gsKnotVector<> kvA(0,1,Na-1,p+1); gsTensorBSplineBasis<2> tpA(kvA,kvA);
    gsGeometry<>::uPtr sigId = makeIdentity(compDeg, compN);

    gsMultiBasis<> mb; mb.addBasis(tpA.clone().release());
    g_c=std::to_string(c0);
    gsGeometry<>::uPtr sigPrev = relocateFrom(*sigId, *S, smoothing, penalty, 250);
    gsMultiPatch<> cmp0; cmp0.addPatch(gsComposedGeometry<real_t>(*sigPrev,*S));
    gsMultiPatch<> un; projectExact(cmp0, mb, un);
    gsInfo<<"Exporting r-ALE snapshots at steps "<<snaps<<" (K="<<K<<", t="<<g_t<<")\n";

    for (index_t k=1;k<=K;++k){
        const real_t ck=c0+(c1-c0)*real_t(k)/real_t(K); g_c=std::to_string(ck);
        gsGeometry<>::uPtr sigNew = relocateFrom(*sigPrev, *S, smoothing, penalty, relocIt);
        gsComposedGeometry<real_t> Gnew(*sigNew,*S), Gprev(*sigPrev,*S);
        MeshVel w(Gnew,Gprev,Dt);
        gsMultiPatch<> cmp; cmp.addPatch(gsComposedGeometry<real_t>(*sigNew,*S));
        gsMultiPatch<> unext; aleStep(cmp, mb, &w, un, Dt, unext);
        if (want.count(k)){
            gsField<> field(cmp, unext);
            const std::string name = "trans_snap_"+std::to_string(k);
            gsWriteParaview<>(field, name, 40000, true);   // .vts (u) + _mesh.vtp (relocated mesh)
            gsInfo<<"  step "<<k<<" (t="<<k*Dt<<", band centre "<<ck<<") -> "<<name<<".vts\n";
        }
        un = unext; sigPrev = sigNew->clone();
    }

    // ---------- THB (local h-refinement, re-adapt each step, identity geometry) ----------
    // Same feature-tracking strategy as reproduce_transient_sphere.cpp: each step a THB
    // mesh is re-adapted from a coarse base toward the current band, on the identity map.
    gsKnotVector<> kvB(0,1,3,p+1); gsTensorBSplineBasis<2> base(kvB,kvB);
    gsMultiPatch<> cmpT; cmpT.addPatch(gsComposedGeometry<real_t>(*sigId,*S));
    g_c=std::to_string(c0);
    gsMultiBasis<> mbT0 = buildTHBforFront(cmpT, base, thbPass, 0.85, 1, 3);
    gsMultiPatch<> unT; projectExact(cmpT, mbT0, unT);
    gsInfo<<"Exporting THB (local-refinement) snapshots at steps "<<snaps<<"\n";
    for (index_t k=1;k<=K;++k){
        const real_t ck=c0+(c1-c0)*real_t(k)/real_t(K); g_c=std::to_string(ck);
        gsMultiBasis<> mb = buildTHBforFront(cmpT, base, thbPass, 0.85, 1, 3);
        gsMultiPatch<> unext; aleStep(cmpT, mb, nullptr, unT, Dt, unext);
        if (want.count(k)){
            gsField<> field(cmpT, unext);
            const std::string name = "thb_snap_"+std::to_string(k);
            gsWriteParaview<>(field, name, 40000, true);   // .vts (u) + _mesh.vtp (THB mesh)
            gsInfo<<"  step "<<k<<" (band centre "<<ck<<") -> "<<name<<".vts ("
                  <<mb.basis(0).size()<<" DoF)\n";
        }
        unT = unext;
    }

    // ---------- rh (relocated sigma + per-step THB, ALE) ----------
    // Combined strategy of reproduce_transient_sphere.cpp: the composite map is
    // relocated to track the band AND a THB mesh is re-adapted on it each step.
    g_c=std::to_string(c0);
    gsGeometry<>::uPtr sigPrevRH = relocateFrom(*sigId, *S, smoothing, penalty, 250);
    gsMultiPatch<> cmpRH0; cmpRH0.addPatch(gsComposedGeometry<real_t>(*sigPrevRH,*S));
    gsMultiBasis<> mbRH0 = buildTHBforFront(cmpRH0, base, thbPass, 0.85, 1, 3);
    gsMultiPatch<> unRH; projectExact(cmpRH0, mbRH0, unRH);
    gsInfo<<"Exporting rh (relocated + THB) snapshots at steps "<<snaps<<"\n";
    for (index_t k=1;k<=K;++k){
        const real_t ck=c0+(c1-c0)*real_t(k)/real_t(K); g_c=std::to_string(ck);
        gsGeometry<>::uPtr sigNew = relocateFrom(*sigPrevRH, *S, smoothing, penalty, relocIt);
        gsComposedGeometry<real_t> Gnew(*sigNew,*S), Gprev(*sigPrevRH,*S);
        MeshVel w(Gnew,Gprev,Dt);
        gsMultiPatch<> cmp; cmp.addPatch(gsComposedGeometry<real_t>(*sigNew,*S));
        gsMultiBasis<> mb = buildTHBforFront(cmp, base, thbPass, 0.85, 1, 3);
        gsMultiPatch<> unext; aleStep(cmp, mb, &w, unRH, Dt, unext);
        if (want.count(k)){
            gsField<> field(cmp, unext);
            const std::string name = "rh_snap_"+std::to_string(k);
            gsWriteParaview<>(field, name, 40000, true);   // .vts (u) + _mesh.vtp (relocated+THB mesh)
            gsInfo<<"  step "<<k<<" (band centre "<<ck<<") -> "<<name<<".vts ("
                  <<mb.basis(0).size()<<" DoF)\n";
        }
        unRH = unext; sigPrevRH = sigNew->clone();
    }
    gsInfo<<"done.\n";
    return 0;
}
