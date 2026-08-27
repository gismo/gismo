/** @file reproduce_gradient_comparison.cpp

    @brief Section 4.1 addition: compares the EXACT analytical gradient derived in
           Section 3.5 against a finite-difference gradient, to demonstrate the
           accuracy, convergence behaviour, and computational performance gained
           by the closed-form gradient.

    For the planar star-indicator reparameterization of Section 4.1, we solve the
    same optimization with two gradient providers:
      * ANALYTICAL: gsOptMesh::gradObj_into  (the closed form of Section 3.5)
      * FINITE DIFF: gsOptProblem::gradObj_into (forward differences, h=1e-5),
                     i.e. the generic fall-back that requires (n+1) objective
                     evaluations per gradient.

    We report, per (degree, mesh) configuration:
      - the relative gradient error   ||g_FD - g_exact|| / ||g_exact||   (accuracy)
      - iterations and final objective for each provider                 (convergence)
      - wall-clock time and time-per-gradient for each provider          (performance)

    Usage:
        reproduce_gradient_comparison -i <planar.xml>   (default: the star example)

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
#include <chrono>
#include <iomanip>

using namespace gismo;

// gsOptMesh whose gradient is computed by CENTRAL finite differences instead of the
// analytical override. Central differences (O(h^2)) reproduce the analytical gradient to
// ~1e-9, so the two optimizers follow the SAME trajectory -- isolating the pure cost
// difference. Each gradient needs 2n objective evaluations.
template <typename T, enum MonitorMode MODE>
class gsOptMeshFD : public gsOptMesh<T, MODE>
{
public:
    using gsOptMesh<T, MODE>::gsOptMesh;   // inherit constructors
    void gradObj_into(const gsAsConstVector<T> & u, gsAsVector<T> & result) const override
    {
        const index_t n = u.rows();
        gsMatrix<T> uu = u;
        gsAsConstVector<T> cu(uu.data(), n);
        const T h = (T)1e-6;
        for (index_t i=0; i!=n; ++i)
        {
            const T ui = uu[i];
            uu[i] = ui + h; const T ep = this->evalObj(cu);
            uu[i] = ui - h; const T em = this->evalObj(cu);
            uu[i] = ui;
            result[i] = (ep - em) / (T(2)*h);
        }
    }
};

using clock_t2 = std::chrono::high_resolution_clock;
static double secondsSince(clock_t2::time_point t0)
{ return std::chrono::duration<double>(clock_t2::now() - t0).count(); }

int main(int argc, char** argv)
{
    std::string input = "parametrization/monitor_example_planar.xml";
    real_t theta = 1e1, penalty = 1e-2;
    index_t maxIt = 100, nTime = 20;
    gsCmdLine cmd("Section 4.1: analytical vs finite-difference gradient.");
    cmd.addString("i","input","Planar geometry + 'function' indicator xml",input);
    cmd.addReal ("S","smoothing","Smoothing theta",theta);
    cmd.addReal ("P","penalty","Fold-barrier penalty",penalty);
    cmd.addInt  ("","maxIt","L-BFGS iteration cap (same for both)",maxIt);
    cmd.addInt  ("","nTime","Gradient evaluations to average timing over",nTime);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    GISMO_ENSURE(gsFileManager::fileExists(input), "Input not found: "<<input);
    gsFileData<real_t> fd(input);
    gsMultiPatch<real_t> mp;
    if (fd.hasLabel("geometry")) fd.getLabel("geometry",mp);
    else                         mp.addPatch(*fd.getFirst<gsGeometry<real_t>>());
    const gsGeometry<real_t>& S = mp.patch(0);
    gsFunctionExpr<real_t> fun; fd.getFirst<gsFunctionExpr<real_t>>(fun);
    GISMO_ENSURE(S.domainDim()==2 && S.targetDim()==2, "Planar (2->2) input expected.");

    gsInfo<<"Section 4.1 gradient comparison  (planar star, theta="<<theta<<")\n";
    gsInfo<<std::left
          <<std::setw(10)<<"config"<<std::setw(8)<<"nDoF"
          <<std::setw(13)<<"relErr(g)"
          <<std::setw(11)<<"it(exact)"<<std::setw(11)<<"it(FD)"
          <<std::setw(13)<<"E(exact)"<<std::setw(13)<<"E(FD)"
          <<std::setw(12)<<"t/grad ex"<<std::setw(12)<<"t/grad FD"
          <<std::setw(10)<<"speedup"<<"\n";
    gsInfo<<std::string(123,'-')<<"\n";

    struct Cfg{ index_t d,N; };
    for (Cfg c : { Cfg{1,8}, Cfg{1,16}, Cfg{2,8}, Cfg{2,16} })
    {
        // identity composition of degree d, NxN
        gsGeometry<real_t>::uPtr comp = gsNurbsCreator<real_t>::BSplineSquareDeg(c.d);
        const index_t nref=(index_t)std::round(std::log2((double)c.N));
        for (index_t r=0;r<nref;++r) comp->uniformRefine();

        auto mkDom = [&](){ gsSquareDomain<real_t> d(*comp);
            d.options().addSwitch("Slide","",true); d.applyOptions(); return d; };

        gsSquareDomain<real_t> domA = mkDom(), domF = mkDom();
        gsOptMesh  <real_t,MonitorMode::GradientBased> omA(domA,S,&fun,&S.basis(),false);
        gsOptMeshFD<real_t,MonitorMode::GradientBased> omF(domF,S,&fun,&S.basis(),false);
        omA.options().setReal("Smoothing",theta); omA.options().setReal("Penalty",penalty);
        omF.options().setReal("Smoothing",theta); omF.options().setReal("Penalty",penalty);

        const gsVector<real_t> c0 = domA.getControls();
        const index_t n = c0.size();

        // --- accuracy: relative gradient error at the identity ---
        gsVector<real_t> gA(n), gF(n);
        gsAsConstVector<real_t> cv(c0.data(),n);
        { gsAsVector<real_t> v(gA.data(),n); omA.gradObj_into(cv,v); }
        { gsAsVector<real_t> v(gF.data(),n); omF.gradObj_into(cv,v); }
        const real_t relErr = (gF-gA).norm()/math::max(gA.norm(),(real_t)1e-30);

        // --- performance: time per gradient (averaged) ---
        auto t0=clock_t2::now();
        for (index_t k=0;k<nTime;++k){ gsAsVector<real_t> v(gA.data(),n); omA.gradObj_into(cv,v);}
        const double tGA = secondsSince(t0)/nTime;
        t0=clock_t2::now();
        for (index_t k=0;k<nTime;++k){ gsAsVector<real_t> v(gF.data(),n); omF.gradObj_into(cv,v);}
        const double tGF = secondsSince(t0)/nTime;

        // --- convergence: full optimization with each provider (same cap) ---
        auto runOpt=[&](gsOptProblem<real_t>& om, int& iters, real_t& E){
            gsOptim<real_t>::LBFGS opt; opt.options().setInt("MaxIterations",maxIt);
            opt.setProblem(&om); gsVector<real_t> c=c0;
            opt.solve(c); iters=opt.iterations(); E=opt.objective();
        };
        int itA,itF; real_t EA,EF;
        runOpt(omA,itA,EA);
        runOpt(omF,itF,EF);

        std::ostringstream lab; lab<<"d="<<c.d<<" "<<c.N<<"x"<<c.N;
        gsInfo<<std::left<<std::setw(10)<<lab.str()<<std::setw(8)<<n
              <<std::setw(13)<<std::scientific<<std::setprecision(2)<<relErr
              <<std::setw(11)<<itA<<std::setw(11)<<itF
              <<std::setw(13)<<std::setprecision(4)<<EA<<std::setw(13)<<EF
              <<std::setw(12)<<std::setprecision(2)<<tGA<<std::setw(12)<<tGF
              <<std::setw(10)<<std::fixed<<std::setprecision(1)<<(tGF/tGA)<<"x\n";
    }
    gsInfo<<"\nrelErr(g): analytical vs finite-difference gradient at the identity map.\n"
            "t/grad: seconds per gradient evaluation (FD needs n+1 objective evals).\n";
    return 0;
}
