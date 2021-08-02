/** @file poisson2_example.cpp

    @brief Tutorial on how to use expression assembler to solve the Poisson equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozillexL.org/MPL/2.0/.

    Author(s): exL. Mantzaflaris
*/

//! [Include namespace]
#include <gismo.h>
#include <gsAssembler/gsAdaptiveRefUtils.h>

using namespace gismo;
//! [Include namespace]

template<typename T>
class gsElementErrorPlotter : public gsFunction<T>
{
public:
    gsElementErrorPlotter(const gsBasis<T>& mp, const std::vector<T>& errors ) : m_mp(mp),m_errors(errors)
    {

    }

    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& res) const
    {
        // Initialize domain element iterator -- using unknown 0
        res.setZero(1,u.cols());
        for(index_t i=0; i<u.cols();++i)
        {
            int iter =0;
            // Start iteration over elements

            typename gsBasis<T>::domainIter domIt = m_mp.makeDomainIterator();
            for (; domIt->good(); domIt->next() )
            {
                 bool flag = true;
                const gsVector<T>& low = domIt->lowerCorner();
                const gsVector<T>& upp = domIt->upperCorner();


                for(int d=0; d<domainDim();++d )
                {
                    if(low(d)> u(d,i) || u(d,i) > upp(d))
                    {
                        flag = false;
                        break;
                    }
                }
                if(flag)
                {
                     res(0,i) = m_errors.at(iter);
                     break;
                }
                iter++;
            }
        }
    }

    short_t domainDim() const { return m_mp.dim();}

private:
    const gsBasis<T>& m_mp;
    const std::vector<T>& m_errors;
};


int main(int argc, char *argv[])
{
    //! [Parse command line]
    index_t est = 0;
    bool plot = false;
    index_t numRefine  = 5;
    index_t numElevate = 0;
    bool last = false;
    std::string fn("pde/poisson2d_bvp2.xml");

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addInt( "E", "errorEstimation",
                "Error estimation method; 0 = strong, 1 = weak", est );
    cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement", last);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read input file]

    gsVector<> pt(2);
    pt.setConstant(0.5);

    gsFileData<> fd(fn);
    gsInfo << "Loaded file "<< fd.lastPath() <<"\n";

    // gsMultiPatch<> mp;
    // mp = *gsNurbsCreator<>::BSplineRectangle(-0.5,-0.5,0.5,0.5);

    // gsFunctionExpr<> f("0",2);

    // gsBoundaryConditions<> bc;
    // bc.addCondition(boundary::north, condition_type::dirichlet, 0 );
    // bc.addCondition(boundary::east, condition_type::dirichlet, 0 );
    // bc.addCondition(boundary::south, condition_type::dirichlet, 0 );
    // bc.addCondition(boundary::west, condition_type::dirichlet, 0 );

    gsMultiPatch<> mp;
    fd.getId(0, mp); // id=0: Multipatch domain

    // // Cast all patches of the mp object to THB splines
    // gsMultiPatch<> mpBspline;
    // // gsTensorBSpline<2,real_t> *geo;
    // gsTHBSpline<2,real_t> thb;
    // for (index_t k=0; k!=mpBspline.nPatches(); ++k)
    // {
    //     gsTensorBSpline<2,real_t> *geo = dynamic_cast< gsTensorBSpline<2,real_t> * > (&mpBspline.patch(k));
    //     thb = gsTHBSpline<2,real_t>(*geo);
    //     mp.addPatch(thb);
    // }

    gsFunctionExpr<> f;
    fd.getId(1, f); // id=1: source function
    gsInfo<<"Source function "<< f << "\n";

    gsBoundaryConditions<> bc;
    fd.getId(2, bc); // id=2: boundary conditions
    bc.setGeoMap(mp);
    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";

    gsOptionList Aopt;
    fd.getId(4, Aopt); // id=4: assembler options

    //! [Read input file]

    //! [Refinement]

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    mp.degreeElevate(numElevate);
    // h-refine each basis
    for (int r =0; r < numRefine-1; ++r)
        mp.uniformRefine();

    // gsTensorBSpline<2,real_t> *geo;
    gsTHBSpline<2,real_t> thb;
    for (index_t k=0; k!=mp.nPatches(); ++k)
    {
        gsTensorBSpline<2,real_t> *geo = dynamic_cast< gsTensorBSpline<2,real_t> * > (&mp.patch(k));
        thb = gsTHBSpline<2,real_t>(*geo);
        mp.patch(k) = thb;
    }

    numRefine = 0;

    gsMultiBasis<> basisL(mp);
    gsMultiBasis<> basisH = basisL;
    basisH.degreeElevate();
    gsInfo<<"Basis Primal: "<<basisL.basis(0)<<"\n";
    gsInfo<<"Basis Dual:   "<<basisH.basis(0)<<"\n";


    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< basisL.minCwiseDegree() <<"\n";
    //! [Refinement]

    //! [Problem setup]
    gsExprAssembler<> exL(1,1);
    gsExprAssembler<> exH(1,1);
    exL.setOptions(Aopt);
    exH.setOptions(Aopt);

    //gsInfo<<"Active options:\n"<< exL.options() <<"\n";
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // Elements used for numerical integration
    exL.setIntegrationElements(basisL);
    gsExprEvaluator<> evL(exL);

    exH.setIntegrationElements(basisH);
    gsExprEvaluator<> evH(exH);

    // Set the geometry map
    geometryMap G = exL.getMap(mp);
    geometryMap H = exH.getMap(mp);

    // Set the discretization space
    space u = exL.getSpace(basisL);
    space v = exH.getSpace(basisH);

    // Set the source term
    auto ff = exL.getCoeff(f, G);
    auto gg = exH.getCoeff(f, H);

    // Recover manufactured solution
    // gsFunctionExpr<> ms("((x-0.5)*(x+0.5)*(y-0.5)*(y+0.5)) / ((x+2/3)*y+3/4)",2);
    //gsInfo<<"Exact solution: "<< ms << "\n";
    gsFunctionExpr<> msP;
    fd.getId(3, msP); // id=3: reference solution primal
    gsFunctionExpr<> msD;
    fd.getId(9, msD); // id=9: reference solution dual

    auto primal_exL = evL.getVariable(msP, G);
    auto dual_exL   = evL.getVariable(msD, G);
    auto primal_exH = evH.getVariable(msP, H);
    auto dual_exH   = evH.getVariable(msD, H);

    // Solution vector and solution variable
    gsMatrix<> primalL, primalH, primalLp, dualL, dualLp, dualH;

    // Solutions and their projections
        // PDE solution on low-order mesh
    solution uL = exL.getSolution(u, primalL);
        // Dual solution on low-order mesh
    solution zL = exL.getSolution(u, dualL);
        // Dual solution projection on high-order mesh
    solution zH = exH.getSolution(v, dualH);

    // zH2
    gsMultiPatch<> zH2_mp(mp);//just initialize for not being empty
    auto zH2 = exL.getCoeff(zH2_mp);
    // zL2
    gsMultiPatch<> zL2_mp(mp);//just initialize for not being empty
    auto zL2 = exH.getCoeff(zL2_mp);
    // uL2
    gsMultiPatch<> uL2_mp(mp);//just initialize for not being empty
    auto uL2 = exH.getCoeff(uL2_mp);


    gsSparseSolver<>::CGDiagonal solver;

    u.setup(bc, dirichlet::interpolation, 0);
    v.setup(bc, dirichlet::interpolation, 0);

    exL.initSystem();
    exH.initSystem();

    //! [Problem setup]
    gsInfo<< "NumDofs Primal: "<<exL.numDofs() <<std::flush;

    // [ORIGINAL PROBLEM]
    // Compute the system matrix and right-hand side
    exL.assemble( igrad(u, G) * igrad(u, G).tr() * meas(G), u * ff * meas(G) );

    // Enforce Neumann conditions to right-hand side
    // auto g_N = exL.getBdrFunction();
    // exL.assembleRhsBc(u * g_N * nv(G).norm(), bc.neumannSides() );


    //gsInfo<<"Sparse Matrix:\n"<< exL.matrix().toDense() <<"\n";
    //gsInfo<<"Rhs vector:\n"<< exL.rhs().transpose() <<"\n";

    // gsDebugVar(exL.matrix().toDense());
    // gsDebugVar(exL.rhs().transpose());

    gsInfo<< ".\n" <<std::flush;// Assemblying done

    solver.compute( exL.matrix() );
    primalL = solver.solve(exL.rhs());

    uL.extract(uL2_mp);// this updates uL2 variable

    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        // ev.options().setSwitch("plot.elements", true);
        gsWriteParaview( mp, "mp",1000,true,true);
        gsWriteParaview( basisL.basis(0), "basis");
        evL.writeParaview( uL   , G, "solution");
        evH.writeParaview( uL2   ,H, "solution2");
        evL.writeParaview( primal_exL   , G, "solution_exact");
    }

    gsInfo<<"Objective function errors J(u)-J(u_h)\n";
    // gsInfo<<"\t exact: "<<ev.integral(u_ex*meas(G))-ev.integral(u_sol*meas(G))<<"\n";
    gsInfo<<"\t exact: "<<evL.integral((primal_exL-uL)*meas(G))<<"\n";
    // [!ORIGINAL PROBLEM]

    // [Low-order Dual PROBLEM]
    // Compute the system matrix and right-hand side
    exL.initSystem();
    exL.assemble( igrad(u, G) * igrad(u, G).tr() * meas(G), u * uL * meas(G) );

    // Enforce Neumann conditions to right-hand side
    // auto g_N = exL.getBdrFunction();
    // exL.assembleRhsBc(u * g_N * nv(G).norm(), bc.neumannSides() );

    solver.compute( exL.matrix() );
    dualL = solver.solve(exL.rhs());

    zL.extract(zL2_mp);// this updates zL2 variable

    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        // ev.options().setSwitch("plot.elements", true);
        evL.writeParaview( zL   , G, "dualL");
        evH.writeParaview( zL2   ,H, "dualL2");
        evL.writeParaview( dual_exL   , G, "dual_exact");
    }

    // [!Low-order Dual PROBLEM]
    // Initialize the system
    exH.initSystem();

    gsInfo<< "NumDofs Dual: "<<exH.numDofs() <<std::flush;

    // [DUAL PROBLEM]
    // Compute the system matrix and right-hand side
    exH.assemble( igrad(v, H) * igrad(v, H).tr() * meas(H), v * uL2 * meas(H) );

    gsInfo<< ".\n" <<std::flush;// Assemblying done

    solver.compute( exH.matrix() );
    dualH = solver.solve(exH.rhs());

    gsMatrix<> dualH0;
    solution zH0 = exH.getSolution(v,dualH0);

    zH.extractFull(dualH0);
    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        // ev.options().setSwitch("plot.elements", true);
        evH.writeParaview( zH   , H, "dualH");
    }
    // [!DUAL PROBLEM]

    zH.extract(zH2_mp);// this updates zH2 variable

    gsDebug<<"zL "<<evL.eval(zL,pt)<<"\n";
    gsDebug<<"zL "<<evH.eval(zL,pt)<<"\n";
    gsDebug<<"zH "<<evL.eval(zH2,pt)<<"\n";
    gsDebug<<"zH "<<evH.eval(zH,pt)<<"\n"; // Different from the above
    gsDebug<<"\n";
    gsDebug<<"grad zL "<<evL.eval(grad(zL),pt)<<"\n";
    gsDebug<<"grad zL "<<evH.eval(grad(zL2),pt)<<"\n";
    gsDebug<<"grad zH "<<evL.eval(grad(zH2),pt)<<"\n";
    gsDebug<<"grad zH "<<evH.eval(grad(zH),pt)<<"\n";
    gsDebug<<"\n";
    gsDebug<<"grad uL "<<evL.eval(grad(zH2),pt)<<"\n";
    gsDebug<<"\n";
    // gsDebug<<"lapl zL "<<evL.eval(lapl(zL),pt)<<"\n";
    // gsDebug<<"lapl zL "<<evH.eval(lapl(zL2),pt)<<"\n";
    // gsDebug<<"lapl zH "<<evL.eval(lapl(zH2),pt)<<"\n";
    // gsDebug<<"lapl zH "<<evH.eval(lapl(zH),pt)<<"\n";
    // gsDebug<<"\n";
    // Integrals via the assembler and partition of unity.
    space v0 = exH.getSpace(basisH); // full basis
    space u0 = exL.getSpace(basisL); // full basis

    gsBoundaryConditions<> nobc;

    u0.setup(nobc,dirichlet::interpolation,0);
    v0.setup(nobc,dirichlet::interpolation,0);

    exL.initSystem();

    gsDebugVar(evL.eval(zL,pt));
    gsDebugVar(evL.eval(u0,pt));
    gsDebugVar(evL.eval(zL.val() * u0,pt));

    exL.assemble(zL.val() * u0);

    gsDebug<<"int zL "<<evL.integral(zL)<<"; "<<exL.rhs().sum()<<"\n";

    exH.initSystem();
    exH.assemble(zL2.val() * v0);
    gsDebug<<"int zH "<<evH.integral(zL2)<<"; "<<exH.rhs().sum()<<"\n";

    exL.initSystem();
    exL.assemble(zH2.val() * u0);
    gsDebug<<"int zH "<<evL.integral(zH2)<<"; "<<exL.rhs().sum()<<"\n";

    exH.initSystem();
    exH.assemble(zH.val() * v0);
    gsDebug<<"int zH "<<evH.integral(zH)<<"; "<<exH.rhs().sum()<<"\n";

    exL.initSystem();
    exL.assemble(u0 * grad(zH2)*grad(uL).tr());
    gsDebug<<"int grad-norm "<<evL.integral(grad(zH2)*grad(uL).tr())<<"; "<<exL.rhs().sum()<<"\n";

    // exL.initSystem(false);
    // exL.assemble(u * grad(zH2)*grad(zH2).tr());
    // gsDebug<<"int grad-norm "<<evL.integral(grad(zH2)*grad(zH2).tr())<<"; "<<exL.rhs().sum()<<"\n";



    gsInfo<<"Objective function errors J(u)-J(u_h)\n";
    real_t error = evL.integral((0.5*primal_exL*primal_exL)*meas(G)) - evL.integral((0.5*uL*uL)*meas(G));
    real_t errest = evL.integral( (zH2-zL) * ff * meas(G)-(((igrad(zH2) - igrad(zL))*igrad(uL).tr()) ) * meas(G));

    // real_t errest = evL.integral( (zH2-zL) * ff * meas(G)-(((igrad(zH2) - igrad(zL))*igrad(uL).tr()) ) * meas(G));
    gsInfo<<"\texact:\t"   <<error<<"\n"
            "\testimate:\t"<<errest<<"\n"
            "\teff:\t"<<errest/error<<"\n";

    // ELEMENT WISE ERROR ESTIMATION
    if (est==0)
    {
        // gsInfo<<evL.eval((ff - lapl(uL)) * (zH2 - zL)*meas(G),pt)<<"\n";
        gsInfo<<evL.eval((ff - lapl(uL)) * (zH2 - zL)*meas(G),pt)<<"\n";



        evL.integralElWise((ff - lapl(uL)) * (zH2 - zL)*meas(G));
        // evH.integralElWise((gg - lapl(uL)) * (zH2 - zL)*meas(G));

        gsVector<> elementNorms = evL.allValues().transpose();
        std::vector<real_t> errors;
        errors.resize(elementNorms.size());
        gsVector<>::Map(&errors[0],elementNorms.size() ) = elementNorms;

        gsInfo<< "  Result (global)    : "<< elementNorms.sum() <<"\n";

        gsElementErrorPlotter<real_t> err_eh(basisH.basis(0),errors);

        const gsField<> elemError_eh( mp, err_eh, false );
        gsWriteParaview<>( elemError_eh, "error_elem_eh", 1000);


        MarkingStrategy adaptRefCrit = PUCA;
        const real_t adaptRefParam = 0.8;
        std::vector<bool> elMarked( errors.size() );

        gsMarkElementsForRef( errors, adaptRefCrit, adaptRefParam, elMarked);
        for (std::vector<bool>::const_iterator i = elMarked.begin(); i != elMarked.end(); ++i)
            gsInfo << *i << ' ';
        gsInfo<<"\n";

        // Refine the marked elements with a 1-ring of cells around marked elements
        gsRefineMarkedElements( mp, elMarked, 1 );
        gsWriteParaview(basisL.basis(0),"basisL",1000);
        gsWriteParaview(mp,"mp",1000,true);


    }
    // FUNCTION WISE ERROR ESTIMATION
    else if (est==1)
    {
        exL.initSystem();
        // ( grad(uLp) * ( grad(v) * (zH - zLp) + v * ( grad(zH) - grad(zLp) ) ) - gg * v * (zH - zLp)  );
        auto lhs = ((zH2 - zL) * grad(uL) * grad(u0).tr()).tr() + u0 * ( jac(zH2) * grad(uL).tr() - grad(zL) * grad(uL).tr() ); // + v * ( grad(zH) - grad(zLp) ) * grad(uLp).tr();
        // auto lhs2 = u0 * ( jac(zH2) * grad(uL).tr() - grad(zL) * grad(uL).tr() ) ;
        auto rhs = ff.val() * (zH2 - zL).val() * u0;

        evL.integralElWise((zH2-zL) * ff * meas(G)-(((igrad(zH2) - igrad(zL))*igrad(uL).tr()) ) * meas(G));
        std::vector<real_t> errors = evL.elementwise();


        gsInfo<< "  Result (global)    : "<< gsAsVector<real_t>(errors).sum()<<"\n";

        gsElementErrorPlotter<real_t> err_eh(basisH.basis(0),errors);

        const gsField<> elemError_eh( mp, err_eh, false );
        gsWriteParaview<>( elemError_eh, "error_elem_eh", 10000);

        MarkingStrategy adaptRefCrit = PUCA;
        const real_t adaptRefParam = 0.8;
        std::vector<bool> elMarked( errors.size() );

        for(index_t k = 0; k < errors.size(); k++)
            if(errors[k] < 0)
                errors[k] *= -1;


        gsMarkElementsForRef( errors, adaptRefCrit, adaptRefParam, elMarked);
        gsInfo<<"index\terror\trefined?\n";
        for (index_t k=0; k!= errors.size(); k++)
        {
            gsInfo << k<<"\t"<<errors[k] <<"\t"<<elMarked[k] << "\n";

        }
        // Refine the marked elements with a 1-ring of cells around marked supports
        gsWriteParaview(mp,"mp_old",1000,true);
        gsWriteParaview(mp.basis(0),"basisL_old",1000);
        gsRefineMarkedElements( mp, elMarked, 0 );
        basisL = gsMultiBasis<>(mp);
        gsWriteParaview(basisL.basis(0),"basisL",1000);
        gsWriteParaview(mp,"mp",1000,true);

    }
    // FUNCTION WISE ERROR ESTIMATION
    else if (est==2)
    {
        exL.initSystem();
        // ( grad(uLp) * ( grad(v) * (zH - zLp) + v * ( grad(zH) - grad(zLp) ) ) - gg * v * (zH - zLp)  );
        auto lhs = ((zH2 - zL) * grad(uL) * grad(u0).tr()).tr() + u0 * ( jac(zH2) * grad(uL).tr() - grad(zL) * grad(uL).tr() ); // + v * ( grad(zH) - grad(zLp) ) * grad(uLp).tr();
        // auto lhs2 = u0 * ( jac(zH2) * grad(uL).tr() - grad(zL) * grad(uL).tr() ) ;
        auto rhs = ff.val() * (zH2 - zL).val() * u0;

        gsMatrix<> res;
        exL.assemble((lhs - rhs)*meas(G));
        res = exL.rhs();

        gsInfo<< "  Result (global)    : "<< res.sum()<<"\n";

        MarkingStrategy adaptRefCrit = PUCA;
        const real_t adaptRefParam = 0.8;
        std::vector<bool> funMarked( res.size() );
        std::vector<real_t> errors( res.size() );
        gsVector<>::Map(&errors[0],res.size() ) = res;

        for(index_t k = 0; k < errors.size(); k++)
            if(errors[k] < 0)
                errors[k] *= -1;


        gsMarkElementsForRef( errors, adaptRefCrit, adaptRefParam, funMarked);
        gsInfo<<"index\terror\trefined?\n";
        for (index_t k=0; k!= errors.size(); k++)
        {
            gsInfo << k<<"\t"<<errors[k] <<"\t"<<funMarked[k] << "\n";

        }
        // Refine the marked elements with a 1-ring of cells around marked supports
        gsWriteParaview(mp,"mp_old",1000,true);
        gsWriteParaview(mp.basis(0),"basisL_old",1000);
        gsRefineMarkedFunctions( mp, funMarked, 0 );
        basisL = gsMultiBasis<>(mp);
        gsWriteParaview(basisL.basis(0),"basisL",1000);
        gsWriteParaview(mp,"mp",1000,true);

    }
    else
        GISMO_ERROR("Estimation method unknown");

    return EXIT_SUCCESS;

}// end main
