/** @file export_sphere_rh_mesh.cpp

    @brief Exports Paraview (.vts + _mesh.vtp) of the adapted THB meshes on the
           EXACT spherical patch (Sec.4.3 / Fig.10 geometry), for the colatitude
           tanh layer u_ex = tanh(kappa*(1/2 - vartheta)) + 1, kappa = 1/t:
             sphere_thb_mesh* : identity sigma  (isotropic THB refinement)
             sphere_rh_mesh*  : relocated sigma (THB on the sigma-aligned mesh)
           Both are refined adaptively to ~targetDof degrees of freedom and written
           in physical (surface) coordinates: the mesh is the THB analysis basis
           mapped through the composite map G = S o sigma onto the sphere, coloured
           by the projected solution. Surface analogue of export_rh_mesh.cpp.

    Usage:
        export_sphere_rh_mesh -i <sphere_patch.xml> [-t 0.025] [-N 8] [-D 1300]

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

// Colatitude tanh layer at vartheta = 1/2, width t = 1/kappa (as in reproduce_rh_refinement.cpp).
static const std::string VSTR = "atan2(sqrt(x^2+y^2),z)";
static std::string g_t = "0.025";
static std::string exactStr() { return "tanh((1/2-("+VSTR+"))/"+g_t+") + 1"; }
static std::string sourceStr()
{
    return "u:=atan2(y,x); v:="+VSTR+"; w:=sqrt(x^2+y^2+z^2); t:="+g_t+";"
           "-(-cos(v)*(1 - tanh((1/2 - v)/t)^2)/t - 2*sin(v)*tanh((1/2 - v)/t)*(1 - tanh((1/2 - v)/t)^2)/t^2)/(w^2*sin(v))";
}

// Identity or gradient-relocated composition sigma of bi-degree (d,d), NxN, fixed boundary.
static gsGeometry<real_t>::uPtr makeSigma(const gsGeometry<real_t>& S, index_t d, index_t N, bool reparam)
{
    gsGeometry<>::uPtr comp = gsNurbsCreator<>::BSplineSquareDeg(d);
    const index_t nref=(index_t)std::round(std::log2((double)N));
    for (index_t r=0;r<nref;++r) comp->uniformRefine();
    gsSquareDomain<real_t> dom(*comp);
    dom.options().addSwitch("Slide","",false); dom.applyOptions();   // fixed boundary (Sec.4.3)
    if (reparam && dom.nControls()>0)
    {
        gsFunctionExpr<> indicator(exactStr(), S.targetDim());
        gsOptim<real_t>::LBFGS opt; opt.options().setInt("MaxIterations",250);
        gsAdaptiveParametrization<real_t,MonitorMode::GradientBased> relo(dom,S,indicator,opt,false);
        relo.options().setReal("Smoothing",1e-2); relo.options().setReal("Penalty",1e-4);
        relo.solve();
        const real_t mJ = relo.computeMinJacobian();
        GISMO_ENSURE(math::isfinite(mJ)&&mJ>0, "reparam folded (minJ="<<mJ<<")");
        gsInfo<<"  [relocated sigma: min det Jsigma = "<<mJ<<"]\n";
    }
    return dom.domain().clone();
}

// Solve the surface Poisson problem on composed geometry cmp with THB basis mb.
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
    std::string input = "parametrization/monitor_example_sphere_patch_with_gradient_function.xml";
    index_t compDeg=2, compN=8, p=2, targetDof=1300, refExt=1, refCrit=3;
    real_t refParam=0.85, quA=2.0;
    std::string tLayer="0.025";
    gsCmdLine cmd("Export adapted THB meshes (THB vs rh) on the exact sphere patch.");
    cmd.addString("i","input","Sphere-patch xml",input);
    cmd.addString("t","layer","Layer width t=1/kappa",tLayer);
    cmd.addInt("N","compN","Composition sigma elements per direction",compN);
    cmd.addInt("D","targetDof","Refine until this many DoF",targetDof);
    cmd.addReal("m","markParam","Marking parameter (Doerfler bulk fraction)",refParam);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    g_t=tLayer;
    GISMO_ENSURE(gsFileManager::fileExists(input), "Input not found: "<<input);

    gsFileData<> fd(input); gsMultiPatch<> geomMP;
    if (fd.hasLabel("geometry")) fd.getLabel("geometry",geomMP);
    else                         geomMP.addPatch(*fd.getFirst<gsGeometry<>>());
    const gsGeometry<>& S = geomMP.patch(0);

    gsGeometry<>::uPtr sigId = makeSigma(S,compDeg,compN,false);
    gsInfo<<"Relocating sigma (r-adaptivity)...\n";
    gsGeometry<>::uPtr sigR  = makeSigma(S,compDeg,compN,true );

    gsKnotVector<> kv(0,1,3,p+1);
    gsTensorBSplineBasis<2> tp(kv,kv);

    struct Case { std::string name; const gsGeometry<>* sigma; };
    for (Case c : { Case{"sphere_thb_mesh", sigId.get()}, Case{"sphere_rh_mesh", sigR.get()} })
    {
        gsMultiPatch<> cmp; cmp.addPatch(gsComposedGeometry<real_t>(*c.sigma, S));
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

        // Reconstruct the solution on the THB basis and write it with its mesh,
        // in physical (surface) coordinates through G = S o sigma.
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
        gsWriteParaview<>(field, c.name, 60000, true);   // physical: .vts + _mesh.vtp

        // Parametric-domain view: same THB analysis mesh mapped through sigma into
        // S's parametric domain [0,1]^2 (flat). geometry = sigma (2->2), so the mesh
        // shows the isotropic THB squares (identity sigma) or the sigma-aligned,
        // concentrated cells (relocated sigma) that the composition then wraps onto S.
        gsMultiPatch<> cmpP; cmpP.addPatch(*c.sigma);
        gsField<> fieldP(cmpP, solMP);
        gsWriteParaview<>(fieldP, c.name+"_param", 60000, true);

        gsInfo<<c.name<<": "<<dofs<<" DoF, "<<mb.basis(0).numElements()<<" elements -> "
              <<c.name<<".vts and "<<c.name<<"_param.vts (+ _mesh.vtp)\n";
    }
    gsInfo<<"wrote sphere_{thb,rh}_mesh.vts (+ _mesh.vtp)\n";
    return 0;
}
