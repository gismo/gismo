/** @file composed_domain_poisson.cpp

    @brief Tutorial on how to use expression assembler to solve the Poisson equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

//! [Include namespace]
#include <gismo.h>
#include <gsNurbs/gsSquareDomain.h>

using namespace gismo;
//! [Include namespace]

template <class T>
class gsSquareDomain1D : public gsFunction<T>
{
    using Base = gsFunction<T> ;

public:
    gsSquareDomain1D()
    {
        index_t nInterior = 2;
        index_t degree = 2;
        gsKnotVector<T> kv(0,1,nInterior,degree+1);

        m_basis = gsBSplineBasis<T>(kv);
        m_domain = gsBSpline<T>(m_basis,m_basis.anchors().transpose());

        m_mapper = gsDofMapper(m_domain.basis(),m_domain.targetDim());

        gsMatrix<index_t> boundary = m_domain.basis().allBoundary();
        for (index_t a = 0; a!=boundary.rows(); a++)
            for (index_t d = 0; d!=m_domain.targetDim(); d++)
                m_mapper.eliminateDof(boundary(a,0),0,d);
        m_mapper.finalize();

        m_parameters.resize(m_mapper.freeSize());
        // std::vector<index_t> i(m_mapper.freeSize());
        // std::vector<index_t> j(m_mapper.freeSize());
        for (index_t k = 0; k!=m_domain.coefs().rows(); k++)
            for (index_t d = 0; d!=m_domain.targetDim(); d++)
                if (m_mapper.is_free(k,0,d))
                {
                    m_parameters[m_mapper.index(k,0,d)] = m_domain.coefs()(k,d);
                    // i[m_mapper.index(k,0,d)] = k; // i index of free entries
                    // j[m_mapper.index(k,0,d)] = d; // j index of free entries
                }

        // This is a way to cast only the free coefficients to a vector, and change an entry of that vector.
        // However, it cannot be used in ''gsVector<T> & controls() override { return m_parameters; };''
        //
        // gsDebugVar(m_domain.coefs()(i,j).diagonal()(0));
        // m_domain.coefs()(i,j).diagonal()(0) = 0.5;
        // gsDebugVar(m_domain.coefs()(i,j).diagonal()(0));
    }

    const gsBSpline<T> & domain() const
    {
        return m_domain;
    }

    gsMatrix<T> support() const override
    {
        return m_domain.support();
    }

    short_t domainDim() const override
    {
        return m_domain.domainDim();
    }

    short_t targetDim() const override
    {
        return m_domain.domainDim();
    }

    void eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const override
    {
        m_domain.eval_into(u,result);
    }

    void updateGeom()
    {
        for (index_t k = 0; k!=m_domain.coefs().rows(); k++)
            for (index_t d = 0; d!=m_domain.targetDim(); d++)
            {
                if (m_mapper.is_free(k,0,d))
                    m_domain.coefs()(k,d) = m_parameters[m_mapper.index(k,0,d)];
            }
    }

    /// Returns the controls of the function
    const gsVector<T> & controls() const override { return m_parameters; };
          gsVector<T> & controls()       override { return m_parameters; };

    /// Returns the number of controls of the function
    size_t nControls() const override
    {
        return m_mapper.freeSize();
    }

    /// Returns the control derivative
    virtual void control_deriv_into(const gsMatrix<T> & points, gsMatrix<T> & result) const override
    {
        gsMatrix<> tmp;

        result.resize(targetDim()*nControls(),points.cols());
        result.setZero();
        for (index_t p = 0; p!=points.cols(); p++)
        {
            gsAsMatrix<> res = result.reshapeCol(p,nControls(),targetDim());
            for (index_t k = 0; k!=m_domain.coefs().rows(); k++)
                for (index_t d = 0; d!=m_domain.targetDim(); d++)
                    if (m_mapper.is_free(k,0,d))
                    {
                        m_domain.basis().evalSingle_into(k,points.col(p),tmp); // evaluate basis function k
                        res(m_mapper.index(k,0,d),d) = tmp(0,0); // tmp is a single value (1 point, 1 basis function)
                    }
        }
    }

protected:
    gsBSplineBasis<T> m_basis;
    gsBSpline<T> m_domain;
    gsDofMapper m_mapper;
    gsVector<T> m_parameters;
};

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 1;
    index_t numElevate = 0;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    // Define a square domain: this is the intermediate map
    gsSquareDomain1D<real_t> domain;

    gsKnotVector<> kv(0,1,0,2);
    gsBSplineBasis<> basis(kv);
    gsMatrix<> coefs = basis.anchors().transpose();
    gsBSpline<> bspline(basis,coefs);
    bspline.embed(2);

    if (numElevate!=0)
        bspline.degreeElevate(numElevate);

    // h-refine
    for (int r =0; r < numRefine; ++r)
        bspline.uniformRefine();


    // Make composed geometry and basis
    const gsBasis<> & tbasis = bspline.basis(); // basis(u,v) -> deriv will give dphi/du ,dphi/dv
    const gsGeometry<> & tgeom = bspline; //G(u,v) -> deriv will give dG/du, dG/dv

    gsMatrix<> pars = domain.controls();
//    gsDebugVar(pars);
    pars *= 0.5;
    // pars(0,0) -= 0.1;
    domain.controls() = pars.col(0);
    domain.updateGeom();

    // Define a composite basis and composite geometry
    // The basis is composed by the square domain
    gsComposedBasis<real_t> cbasis(domain,tbasis); // basis(u,v) = basis(sigma(xi,eta)) -> deriv will give dphi/dxi, dphi/deta
    // The geometry is defined using the composite basis and some coefficients
    gsComposedGeometry<real_t> cgeom(cbasis, tgeom.coefs()); // G(u,v) = G(sigma(xi,eta))  -> deriv will give dG/dxi, dG/deta

    if (plot)
    {
        gsWriteParaview(cgeom,"cgeom",100);
        gsWriteParaview(cbasis,"cbasis",100);
    }

    gsMultiPatch<> mp;
    mp.addPatch(cgeom);

    gsMultiBasis<> dbasis(mp, true);

    //! [Refinement]

     // Source function:
     // gsFunctionExpr<> f("pi*pi*sin(pi*x)",2);

     // Exact solution
     gsFunctionExpr<> ms("sin(pi*x)",2);


    // Source function:
//    gsFunctionExpr<> f("2*pi^2*cos(pi*x)*cos(pi*y)",2);

    // Exact solution
//    gsFunctionExpr<> ms("cos(pi*x)*cos(pi*y)",2);

    gsBoundaryConditions<> bc;
    // bc.addCondition(boundary::side::west ,condition_type::dirichlet,&ms);
    // bc.addCondition(boundary::side::east ,condition_type::dirichlet,&ms);
    // bc.addCondition(boundary::side::south,condition_type::dirichlet,&ms);
    // bc.addCondition(boundary::side::north,condition_type::dirichlet,&ms);
    bc.setGeoMap(mp);

    //! [Problem setup]
    gsExprAssembler<> A(1,1);

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);

    // Set the geometry map
    geometryMap G = A.getMap(mp);

    // Set the discretization space
    space u = A.getSpace(dbasis);

    // // Set the source term
    // auto ff = A.getCoeff(f, G);

    // Recover manufactured solution
    auto u_ex = ev.getVariable(ms, G);

    // Solution vector and solution variable
    gsMatrix<> solVector;
    solution u_sol = A.getSolution(u, solVector);

    //! [Problem setup]

    gsSparseSolver<>::CGDiagonal solver;

    u.setup(bc, dirichlet::homogeneous, 0);

    // Initialize the system
    A.initSystem();

    gsInfo<< "A.numDofs() = " << A.numDofs() <<std::flush;

    // Compute the system matrix and right-hand side
   A.assemble(
       u * u.tr() * meas(G) //matrix
       ,
       u * u_ex * meas(G) //rhs vector
       );

// //  grad(u)*jac(G).ginv()
//   A.assemble(
//       (grad(u)*(jac(G).ginv().tr())) * (grad(u)*(jac(G).ginv().tr())).tr() * meas(G) //matrix
//       ,
//       u * ff * meas(G) //rhs vector
//   );

    solver.compute( A.matrix() );
    solVector = solver.solve(A.rhs());

    // Compute the error
    real_t L2err = math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) );
    gsInfo<<"\nL2 error = "<<L2err<<"\n";

    //! [Export visualization in ParaView]
    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";

        gsParaviewCollection collection("ParaviewOutput/solution", &ev);
        collection.options().setInt("numPoints", 50);
        collection.options().setInt("precision", 12);
//        collection.options().setInt("plotElements.resolution", 16);
        collection.newTimeStep(&mp);
        collection.addField(u_sol,"numerical solution");
        collection.addField(u_ex, "exact solution");
        collection.addField((u_ex-u_sol).sqNorm(), "error");
        collection.saveTimeStep();
        collection.save();

        // gsFileManager::open("ParaviewOutput/solution.pvd");
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]




//    gsDebug<<"GEOMETRY==================================================================\n";
//
//    gsMatrix<> point(2,1);
//    point.col(0) << 0.5, 0.25;
//    gsDebugVar(ev.eval(jac(G), point));
//
////    gsDebugVar( ev.integral(jac(G).det()) );
////    gsDebugVar( ev.integral( meas(G) ) );
//
//    gsVector<> pp = point.col(0);
//    gsMatrix<> ev1, ev2;
//    gsMatrix<> der;
//
//    real_t delta = 1e-6;
//    gsVector<> pt1 = pp, pt2 = pp;
//    pt1.at(1) += delta;
//    pt2.at(1) -= delta;
//    cgeom.eval_into(pt1,ev1);
//    cgeom.eval_into(pt2,ev2);
//    gsVector<> der_eta = (ev1-ev2)/(2*delta);
//    gsDebugVar(der_eta);
//
//    pt1 = pp, pt2 = pp;
//    pt1.at(0) += delta;
//    pt2.at(0) -= delta;
//    cgeom.eval_into(pt1, ev1);
//    cgeom.eval_into(pt2, ev2);
//    gsVector<> der_xi = (ev1-ev2)/(2*delta);
//    gsDebugVar(der_xi);
//
//    cgeom.deriv_into(pp,der);
////  gsDebugVar(der);
//    gsDebugVar(der.reshape(2,2));

//    gsDebug<<"BASIS==================================================================\n";
//    gsMatrix<index_t> act;
//    cbasis.active_into(pp,act);
//    cbasis.evalSingle_into(act(0,0),pt1,ev1);
//    cbasis.evalSingle_into(act(0,0),pt2,ev2);
//
//    gsVector<> der2 = (ev1-ev2)/(2*delta);
//    gsDebugVar(der2);
//
//    cbasis.derivSingle_into(act(0,0),pp,der);
//    gsDebugVar(der);

    return EXIT_SUCCESS;

}// end main
