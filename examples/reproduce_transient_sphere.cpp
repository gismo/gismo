/** @file reproduce_transient_sphere.cpp

    @brief Transient heat equation with a MOVING sharp front on the EXACT
           spherical patch (Sec.4.3 / Fig.10 geometry), backward Euler, comparing
           r-adaptivity (ALE on a fixed mesh) against uniform, THB, and rh.

    Surface heat equation  du/dt - Lap_S(u) = f  on the sphere patch, manufactured
    moving colatitude layer
        u(x,t) = tanh( (c(t) - vartheta)/tw ) + 1,   vartheta = colatitude,
        c(t)   = c0 + (c1-c0) t/T   (the latitude band sweeps in colatitude),
    with f = du/dt - Lap_S(u) (both analytic) and u imposed as a time-dependent
    Dirichlet datum. Lap_S is the Laplace-Beltrami of Sec.4.3 with 1/2 -> c(t).

    r-adaptivity: ALE backward Euler on a FIXED analysis basis whose composed
    geometry G = S o sigma is relocated to track the band each step (warm-started).
    Because the basis is fixed, the previous coefficient field carries over with no
    inter-mesh projection; the mesh motion enters through the PHYSICAL mesh velocity
        w(xi) = ( S(sigma_new(xi)) - S(sigma_prev(xi)) ) / Dt   in R^3,
    a tangential surface velocity (this is the surface analogue of the planar
    reproduce_transient_moving.cpp, where w reduced to Dsigma/Dt).

    Usage:
        reproduce_transient_sphere -i <sphere_patch.xml> [-t 0.05] [-K 20] [-N 24]
                                   [--c0 0.35] [--c1 0.7] [--Kreloc 1] [--thbPass 4]

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

using namespace gismo;

static const std::string VSTR = "atan2(sqrt(x^2+y^2),z)";   // colatitude vartheta
static std::string g_t="0.05", g_c="0.35", g_cdot="0.35";
static std::string exactStr(){ return "tanh(("+g_c+"-("+VSTR+"))/"+g_t+")+1"; }
// f = du/dt - Lap_S u, in a single exprtk declaration block (du/dt + (-Lap_S u)).
// -Lap_S u is the Sec.4.3 Laplace-Beltrami with the band centre 1/2 replaced by c(t).
static std::string fStr()
{
    const std::string V ="(atan2(sqrt(x^2+y^2),z))";   // colatitude vv
    const std::string W2="(x^2+y^2+z^2)";              // radius^2 (=ww^2)
    const std::string Tc="tanh(("+g_c+"-"+V+")/"+g_t+")";
    const std::string g2="(1-"+Tc+"^2)";               // sech^2
    const std::string dudt="("+g2+"*"+g_cdot+"/"+g_t+")";
    const std::string mlap="( -( -cos"+V+"*"+g2+"/"+g_t
                          +" - 2*sin"+V+"*"+Tc+"*"+g2+"/"+g_t+"^2 )/("+W2+"*sin"+V+") )";
    return dudt+" + "+mlap;                            // du/dt + (-Lap_S u)
}
static real_t g_quA = 2.0;

// Physical mesh velocity w(xi) = (G_new(xi) - G_prev(xi))/Dt in R^3, G = S o sigma.
class MeshVel : public gsFunction<real_t>
{
    const gsGeometry<real_t>& m_Gnew;
    const gsGeometry<real_t>& m_Gprev;
    real_t m_Dt;
public:
    MeshVel(const gsGeometry<real_t>& Gnew, const gsGeometry<real_t>& Gprev, real_t Dt)
        : m_Gnew(Gnew), m_Gprev(Gprev), m_Dt(Dt) {}
    GISMO_CLONE_FUNCTION(MeshVel)
    short_t domainDim() const { return 2; }
    short_t targetDim() const { return 3; }
    void eval_into(const gsMatrix<real_t>& u, gsMatrix<real_t>& res) const
    {
        gsMatrix<real_t> a,b; m_Gnew.eval_into(u,a); m_Gprev.eval_into(u,b);
        res = (a-b)/m_Dt;
    }
};

static gsGeometry<real_t>::uPtr makeIdentity(index_t d, index_t N)
{
    gsGeometry<>::uPtr comp = gsNurbsCreator<>::BSplineSquareDeg(d);
    const index_t nref=(index_t)std::round(std::log2((double)N));
    for (index_t r=0;r<nref;++r) comp->uniformRefine();
    return comp;
}

// Relocate sigma to track the current front (gradient monitor, fixed boundary).
static gsGeometry<real_t>::uPtr relocateFrom(const gsGeometry<real_t>& init, const gsGeometry<real_t>& S,
                                             real_t smoothing, real_t penalty, index_t maxIt)
{
    gsGeometry<>::uPtr comp = init.clone();
    gsSquareDomain<real_t> dom(*comp);
    dom.options().addSwitch("Slide","",false); dom.applyOptions();   // fixed boundary (Sec.4.3)
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

// One backward-Euler (ALE) step on the surface. wFun != null -> moving mesh.
static void aleStep(const gsMultiPatch<real_t>& cmp, gsMultiBasis<real_t>& mb,
                    const gsFunction<real_t>* wFun, const gsMultiPatch<real_t>& uPrevField,
                    real_t Dt, gsMultiPatch<real_t>& uNextField, real_t& l2err, real_t& h1err)
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
    gsExprEvaluator<> ev(A);
    ev.options().setReal("quA",g_quA); ev.options().setInt("quB",2);
    ev.options().update(seOff, gsOptionList::addIfUnknown);

    geometryMap G = A.getMap(cmp);
    space u = A.getSpace(mb);
    auto fsrc  = A.getCoeff(cf);
    auto uprev = A.getCoeff(uPrevField.patch(0));
    auto uex   = ev.getVariable(cuex);
    gsMatrix<real_t> solVector;
    solution usol = A.getSolution(u, solVector);

    u.setup(bc, dirichlet::l2Projection, 0);
    A.initSystem();
    if (wFun) {
        auto w = A.getCoeff(*wFun);
        A.assemble( u*u.tr()*(1.0/Dt)*meas(G)
                  + igrad(u,G)*igrad(u,G).tr()*meas(G)
                  - (igrad(u,G)*w)*u.tr()*meas(G),
                    u*fsrc*meas(G) + u*uprev*(1.0/Dt)*meas(G) );
    } else {
        A.assemble( u*u.tr()*(1.0/Dt)*meas(G)
                  + igrad(u,G)*igrad(u,G).tr()*meas(G),
                    u*fsrc*meas(G) + u*uprev*(1.0/Dt)*meas(G) );
    }
    gsSparseSolver<>::CGDiagonal solver(A.matrix());
    solver.setTolerance(1e-12); solver.setMaxIterations(20000);
    solVector = solver.solve(A.rhs());
    GISMO_ENSURE(solVector.allFinite(), "non-finite solve");
    l2err = math::sqrt( ev.integral( (uex-usol).sqNorm()*meas(G) ) );
    h1err = l2err + math::sqrt( ev.integral( (igrad(uex,G)-igrad(usol,G)).sqNorm()*meas(G) ) );
    usol.extract(uNextField);
}

static void projectExact(const gsMultiPatch<real_t>& cmp, gsMultiBasis<real_t>& mb,
                         gsMultiPatch<real_t>& field)
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
    geometryMap G = A.getMap(cmp);
    space u = A.getSpace(mb);
    auto uc = A.getCoeff(cuex);
    gsMatrix<real_t> solVector;
    solution usol = A.getSolution(u, solVector);
    u.setup(bc, dirichlet::l2Projection, 0);
    A.initSystem();
    A.assemble( u*u.tr()*meas(G), u*uc*meas(G) );
    gsSparseSolver<>::CGDiagonal solver(A.matrix());
    solver.setTolerance(1e-14); solver.setMaxIterations(20000);
    solVector = solver.solve(A.rhs());
    usol.extract(field);
}

// THB basis (identity geometry) refined toward the current band by ||grad_S u_ex||.
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

int main(int argc, char** argv)
{
    std::string input = "parametrization/monitor_example_sphere_patch_with_gradient_function.xml";
    index_t p=2, compDeg=2, compN=16, Na=24, K=20, relocIt=40, thbPass=4, Kreloc=1;
    real_t smoothing=1e-2, penalty=1e-4, T=1.0, c0=0.35, c1=0.70;
    std::string tLayer="0.05";
    gsCmdLine cmd("Transient moving-front heat equation on the exact sphere patch.");
    cmd.addString("i","input","Sphere-patch xml",input);
    cmd.addString("t","layer","Front width tw",tLayer);
    cmd.addInt("K","steps","Time steps",K);
    cmd.addInt("N","Na","Fixed analysis mesh per direction",Na);
    cmd.addInt("","compN","Deformation sigma mesh per direction",compN);
    cmd.addInt("","relocIt","Warm-start relocation iters per step",relocIt);
    cmd.addInt("","Kreloc","Relocate sigma every Kreloc steps (static in between)",Kreloc);
    cmd.addReal("","c0","Front start colatitude",c0);
    cmd.addReal("","c1","Front end colatitude",c1);
    cmd.addInt("","thbPass","THB feature-refinement passes per step",thbPass);
    cmd.addReal("S","smoothing","Monitor smoothing theta",smoothing);
    cmd.addReal("P","penalty","Fold-barrier penalty",penalty);
    cmd.addReal("T","final","Final time",T);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    g_t=tLayer;
    const real_t Dt=T/K, cdot=(c1-c0)/T;
    { std::ostringstream o; o<<cdot; g_cdot=o.str(); }
    GISMO_ENSURE(gsFileManager::fileExists(input), "Input not found: "<<input);

    gsFileData<> fd(input); gsMultiPatch<> geomMP;
    if (fd.hasLabel("geometry")) fd.getLabel("geometry",geomMP);
    else                         geomMP.addPatch(*fd.getFirst<gsGeometry<>>());
    gsGeometry<>::uPtr S = geomMP.patch(0).clone();

    gsInfo<<"Transient heat on sphere patch: band width tw="<<g_t<<", colatitude "<<c0<<"->"<<c1
          <<", "<<K<<" steps, Dt="<<Dt<<", analysis "<<Na<<"x"<<Na<<", sigma "<<compN<<"x"<<compN<<"\n";

    gsKnotVector<> kvA(0,1,Na-1,p+1);
    gsTensorBSplineBasis<2> tpA(kvA,kvA);
    gsGeometry<>::uPtr sigId = makeIdentity(compDeg, compN);

    std::vector<real_t> errU(K,0), errR(K,0), errT(K,0), errRH(K,0);
    std::vector<real_t> errU1(K,0), errR1(K,0), errT1(K,0), errRH1(K,0);
    std::vector<index_t> dofT(K,0), dofRH(K,0);
    real_t wallU=0, wallR=0, wallT=0, wallRH=0;
    const index_t dofFixed=(Na+p)*(Na+p);

    // ---------- static uniform (identity sigma, w=0) ----------
    {
        gsStopwatch clk;
        gsMultiPatch<> cmp; cmp.addPatch(gsComposedGeometry<real_t>(*sigId,*S));
        gsMultiBasis<> mb; mb.addBasis(tpA.clone().release());
        g_c=std::to_string(c0);
        gsMultiPatch<> un; projectExact(cmp, mb, un);
        for (index_t k=1;k<=K;++k){
            const real_t ck=c0+(c1-c0)*real_t(k)/real_t(K); g_c=std::to_string(ck);
            gsMultiPatch<> unext; real_t e,eh; aleStep(cmp, mb, nullptr, un, Dt, unext, e, eh);
            errU[k-1]=e; errU1[k-1]=eh; un=unext;
        }
        wallU=clk.stop();
    }
    // ---------- r-adaptivity (ALE, fixed mesh, relocated sigma) ----------
    {
        gsStopwatch clk;
        gsMultiBasis<> mb; mb.addBasis(tpA.clone().release());
        g_c=std::to_string(c0);
        gsGeometry<>::uPtr sigPrev = relocateFrom(*sigId, *S, smoothing, penalty, 250);
        gsMultiPatch<> cmp0; cmp0.addPatch(gsComposedGeometry<real_t>(*sigPrev,*S));
        gsMultiPatch<> un; projectExact(cmp0, mb, un);
        for (index_t k=1;k<=K;++k){
            const real_t ck=c0+(c1-c0)*real_t(k)/real_t(K); g_c=std::to_string(ck);
            gsMultiPatch<> unext; real_t e,eh;
            if ((k-1)%Kreloc==0){
                gsGeometry<>::uPtr sigNew = relocateFrom(*sigPrev, *S, smoothing, penalty, relocIt);
                gsComposedGeometry<real_t> Gnew(*sigNew,*S), Gprev(*sigPrev,*S);
                MeshVel w(Gnew,Gprev,Dt);
                gsMultiPatch<> cmp; cmp.addPatch(gsComposedGeometry<real_t>(*sigNew,*S));
                aleStep(cmp, mb, &w, un, Dt, unext, e, eh);
                sigPrev = sigNew->clone();
            } else {
                gsMultiPatch<> cmp; cmp.addPatch(gsComposedGeometry<real_t>(*sigPrev,*S));
                aleStep(cmp, mb, nullptr, un, Dt, unext, e, eh);
            }
            errR[k-1]=e; errR1[k-1]=eh; un=unext;
        }
        wallR=clk.stop();
    }
    // ---------- THB (feature-tracking, re-adapt each step, identity geometry) ----------
    {
        gsStopwatch clk;
        gsKnotVector<> kvB(0,1,3,p+1); gsTensorBSplineBasis<2> base(kvB,kvB);
        gsMultiPatch<> cmp; cmp.addPatch(gsComposedGeometry<real_t>(*sigId,*S));
        g_c=std::to_string(c0);
        gsMultiBasis<> mb0 = buildTHBforFront(cmp, base, thbPass, 0.85, 1, 3);
        gsMultiPatch<> un; projectExact(cmp, mb0, un);
        for (index_t k=1;k<=K;++k){
            const real_t ck=c0+(c1-c0)*real_t(k)/real_t(K); g_c=std::to_string(ck);
            gsMultiBasis<> mb = buildTHBforFront(cmp, base, thbPass, 0.85, 1, 3);
            gsMultiPatch<> unext; real_t e,eh; aleStep(cmp, mb, nullptr, un, Dt, unext, e, eh);
            errT[k-1]=e; errT1[k-1]=eh; dofT[k-1]=mb.basis(0).size(); un=unext;
        }
        wallT=clk.stop();
    }
    // ---------- rh (relocated sigma + per-step THB, ALE) ----------
    {
        gsStopwatch clk;
        gsKnotVector<> kvB(0,1,3,p+1); gsTensorBSplineBasis<2> base(kvB,kvB);
        g_c=std::to_string(c0);
        gsGeometry<>::uPtr sigPrev = relocateFrom(*sigId, *S, smoothing, penalty, 250);
        gsMultiPatch<> cmp0; cmp0.addPatch(gsComposedGeometry<real_t>(*sigPrev,*S));
        gsMultiBasis<> mb0 = buildTHBforFront(cmp0, base, thbPass, 0.85, 1, 3);
        gsMultiPatch<> un; projectExact(cmp0, mb0, un);
        for (index_t k=1;k<=K;++k){
            const real_t ck=c0+(c1-c0)*real_t(k)/real_t(K); g_c=std::to_string(ck);
            gsGeometry<>::uPtr sigCur; bool moved=false;
            if ((k-1)%Kreloc==0){
                gsGeometry<>::uPtr sigNew = relocateFrom(*sigPrev, *S, smoothing, penalty, relocIt);
                sigCur = sigNew->clone(); moved=true;
            } else sigCur = sigPrev->clone();
            gsComposedGeometry<real_t> Gcur(*sigCur,*S), Gprev(*sigPrev,*S);
            MeshVel w(Gcur,Gprev,Dt);
            gsMultiPatch<> cmp; cmp.addPatch(gsComposedGeometry<real_t>(*sigCur,*S));
            gsMultiBasis<> mb = buildTHBforFront(cmp, base, thbPass, 0.85, 1, 3);
            gsMultiPatch<> unext; real_t e,eh;
            aleStep(cmp, mb, moved? &w : nullptr, un, Dt, unext, e, eh);
            if (moved) sigPrev = sigCur->clone();
            errRH[k-1]=e; errRH1[k-1]=eh; dofRH[k-1]=mb.basis(0).size(); un=unext;
        }
        wallRH=clk.stop();
    }

    gsInfo<<"\n"<<std::left<<std::setw(6)<<"step"<<std::setw(8)<<"t"
          <<std::setw(14)<<"uniform L2"<<std::setw(14)<<"r(ALE) L2"<<std::setw(14)<<"THB L2"
          <<std::setw(14)<<"rh L2"<<std::setw(9)<<"THB DoF"<<"\n"<<std::string(74,'-')<<"\n";
    real_t mu=0,mr=0,mt=0,mrh=0,mu1=0,mr1=0,mt1=0,mrh1=0; real_t adofT=0,adofRH=0;
    for (index_t k=0;k<K;++k){
        gsInfo<<std::left<<std::setw(6)<<(k+1)<<std::fixed<<std::setprecision(3)<<std::setw(8)<<((k+1)*Dt)
              <<std::scientific<<std::setprecision(3)<<std::setw(14)<<errU[k]<<std::setw(14)<<errR[k]
              <<std::setw(14)<<errT[k]<<std::setw(14)<<errRH[k]<<std::fixed<<std::setw(9)<<dofT[k]<<"\n";
        mu=math::max(mu,errU[k]); mr=math::max(mr,errR[k]); mt=math::max(mt,errT[k]); mrh=math::max(mrh,errRH[k]);
        mu1=math::max(mu1,errU1[k]); mr1=math::max(mr1,errR1[k]); mt1=math::max(mt1,errT1[k]); mrh1=math::max(mrh1,errRH1[k]);
        adofT+=dofT[k]; adofRH+=dofRH[k];
    }
    adofT/=K; adofRH/=K;
    gsInfo<<std::string(74,'-')<<"\n";
    gsInfo<<"DoF:  uniform/r (fixed) = "<<dofFixed<<",  THB (avg) = "<<std::fixed<<std::setprecision(0)<<adofT
          <<",  rh (avg) = "<<adofRH<<"\n";
    gsInfo<<"max H1 over time:  uniform = "<<std::scientific<<std::setprecision(3)<<mu1
          <<",  r(ALE) = "<<mr1<<",  THB = "<<mt1<<",  rh = "<<mrh1<<"\n";
    gsInfo<<"max L2 over time:  uniform = "<<std::scientific<<std::setprecision(3)<<mu
          <<",  r(ALE) = "<<mr<<",  THB = "<<mt<<",  rh = "<<mrh<<"\n";
    gsInfo<<std::fixed<<std::setprecision(3)
          <<"wall-clock:  uniform = "<<wallU<<" s,  r(ALE) = "<<wallR<<" s,  THB = "<<wallT<<" s,  rh = "<<wallRH<<" s\n";
    return 0;
}
