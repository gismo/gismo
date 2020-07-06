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

template<class T>
class gsDWRHelper
{
    public:
        /// Shared pointer for solutionFunction
        typedef memory::shared_ptr< gsDWRHelper > Ptr;

        /// Unique pointer for solutionFunction
        typedef memory::unique_ptr< gsDWRHelper > uPtr;

    gsDWRHelper(const gsMultiBasis<>& basisL,
                const gsMultiBasis<>& basisH,
                const gsMultiPatch<> & mp
                ) : m_basisL(basisL), m_basisH(basisH), m_patches(mp) { this->initialize(); }

    protected:
        void initialize()
        {
            gsExprAssembler<> assembler(1,1);
            assembler.setIntegrationElements(m_basisH); //  is this correct?
            assembler.getSpace(m_basisH);
            assembler.initSystem();
            m_nH = assembler.numDofs();

            assembler.setIntegrationElements(m_basisL); //  is this correct?
            assembler.getSpace(m_basisL);
            assembler.initSystem();
            m_nL = assembler.numDofs();
        }

        void computeFromTo(const gsMultiBasis<>& basis1, const gsMultiBasis<>& basis2, gsSparseMatrix<T> & result)
        {
            gsExprAssembler<> assembler(1,1);
            geometryMap G = assembler.getMap(m_patches);

            assembler.setIntegrationElements(basis1); //  is this correct?
            space u1 = assembler.getSpace(basis1);
            space u2 = assembler.getTestSpace(u1 , basis2);
            assembler.initSystem();
            assembler.assemble(u2 * u1.tr() * meas(G));
            result = assembler.matrix();
        }

    public:
        void computeLL()
        {
            if ((m_matrixLL.rows()==0) || (m_matrixLL.cols()==0))
                computeFromTo(m_basisL, m_basisL, m_matrixLL);
            // else
            //     gsInfo<<"LL matrix already computed\n";
        }
        void computeHH()
        {
            if ((m_matrixHH.rows()==0) || (m_matrixHH.cols()==0))
                computeFromTo(m_basisH, m_basisH, m_matrixHH);
            // else
            //     gsInfo<<"HH matrix already computed\n";
        }
        void computeHL()
        {
            if ((m_matrixLH.rows()==0) || (m_matrixLH.cols()==0))
                computeFromTo(m_basisH, m_basisL, m_matrixHL);
            else
                m_matrixHL = m_matrixLH.transpose();
        }
        void computeLH()
        {
            if ((m_matrixHL.rows()==0) || (m_matrixHL.cols()==0))
                computeFromTo(m_basisL, m_basisH, m_matrixLH);
            else
                m_matrixLH = m_matrixHL.transpose();
        }

        const gsSparseMatrix<T> & matrixLL() {this->computeLL(); return m_matrixLL; }
        const gsSparseMatrix<T> & matrixHH() {this->computeHH(); return m_matrixHH; }
        const gsSparseMatrix<T> & matrixHL() {this->computeHL(); return m_matrixHL; }
        const gsSparseMatrix<T> & matrixLH() {this->computeLH(); return m_matrixLH; }

        void projectVector(const gsMatrix<T> & vectorIn, gsMatrix<T> & vectorOut)
        {
            size_t rows = vectorIn.rows();
            if ( rows == m_nL ) // from low to high
            {
                this->computeHH();
                this->computeLH();
                m_solver.compute(m_matrixHH);
                vectorOut = m_solver.solve(m_matrixLH*vectorIn);
            }
            else if (rows == m_nH) // from high to low
            {
                this->computeLL();
                this->computeHL();
                m_solver.compute(m_matrixLL);
                vectorOut = m_solver.solve(m_matrixHL*vectorIn);
            }
            else
                gsInfo<<"WARNING: cannot project vector, size mismatch; rows = "<<rows<<"; numDofs L = "<<m_nL<<"; numDofs H = "<<m_nL<<"\n";
        }


    protected:
        const gsMultiBasis<T> & m_basisL;
        const gsMultiBasis<T> & m_basisH;
        const gsMultiPatch<T> & m_patches;
        gsSparseMatrix<T> m_matrixHH, m_matrixLL, m_matrixHL, m_matrixLH;
        gsSparseSolver<>::CGDiagonal m_solver;

        size_t m_nL, m_nH;

        typedef gsExprAssembler<>::geometryMap geometryMap;
        typedef gsExprAssembler<>::variable    variable;
        typedef gsExprAssembler<>::space       space;
        typedef gsExprAssembler<>::solution    solution;


}; // class definition ends

template<class T>
class solutionFunction : public gismo::gsFunction<T>
{
  // Computes pressure and optionally (to do) traction on any point on a geometry
  //
protected:
    const gsMultiPatch<> & _mp;
    const gsExprAssembler<> & _A; // Expression Assembler
    const gsExprAssembler<>::solution & _soln;
    mutable index_t m_pIndex;

public:
    /// Shared pointer for solutionFunction
    typedef memory::shared_ptr< solutionFunction > Ptr;

    /// Unique pointer for solutionFunction
    typedef memory::unique_ptr< solutionFunction > uPtr;

    solutionFunction(   const gsMultiPatch<> & mp,
                        const gsExprAssembler<> & A,
                     gsExprAssembler<>::solution & soln
                     ) : _mp(mp), _A(A), _soln(soln), _mm_piece(nullptr){ }

    GISMO_CLONE_FUNCTION(solutionFunction)

    short_t domainDim() const {return 2;}

    short_t targetDim() const {return 1;}

    mutable solutionFunction<T> * _mm_piece; // todo: improve the way pieces are accessed

    const gsFunction<T> & piece(const index_t p) const
    {
        // delete m_piece;
        // m_piece = new gsMaterialMatrix(m_patches->piece(k), *m_thickness, *m_YoungsModulus, *m_PoissonRatio);
        _mm_piece = new solutionFunction(*this);
        _mm_piece->setPatch(p);
        return *_mm_piece;
    }

    ~solutionFunction() { delete _mm_piece; }

    void setPatch(index_t p) {m_pIndex = p; }

    void eval_into(const gsMatrix<>& u, gsMatrix<>& res) const
    {
        // u is an d x n matrix of points with dimension d and n points
        // res2 is temporary vector for pressure components
        // res3 is vector with pressure and traction.
        // Data is stored as: [ tractionX1,tractionX2,tractionX3,...
        //                      tractionY1,tractionY2,tractionY3,...
        //                      pressure1, pressure2, pressure3 ,... ]

        gsExprEvaluator<T> _ev(_A);

        // Pre-allocate memory for pressure.
        res.resize(this->targetDim(), u.cols());
        res.setZero();
        // Create objects for parametric points (v) and for pressure-ONLY result (res2)
        gsVector <T> v;
        gsVector <index_t> patches;
        gsMatrix<T> points;
        _mp.locatePoints(u, patches, points);

        for( index_t k = 0; k != u.cols(); k++)
        {
            v = points.col(k);
            res.col(k) = _ev.eval(_soln, v, patches.at(k));
        }
    }
};

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
    // bc.setMap(mp);
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
    // Set the degree of the higher-order basis one higher.
    gsMultiPatch<> mpH0 = mp;
    mpH0    .degreeElevate(1);

    // Cast all patches of the mp object to THB splines
    gsMultiPatch<> mpL, mpH;
    // gsTensorBSpline<2,real_t> *geo;
    gsTHBSpline<2,real_t> thb;
    for (index_t k=0; k!=mp.nPatches(); ++k)
    {
        gsTensorBSpline<2,real_t> *geoL = dynamic_cast< gsTensorBSpline<2,real_t> * > (&mp.patch(k));
        thb = gsTHBSpline<2,real_t>(*geoL);
        mpL.addPatch(thb);

        gsTensorBSpline<2,real_t> *geoH = dynamic_cast< gsTensorBSpline<2,real_t> * > (&mpH0.patch(k));
        thb = gsTHBSpline<2,real_t>(*geoH);
        mpH.addPatch(thb);
    }

    numRefine = 0;

    gsMultiBasis<> basisL(mpL);
    gsMultiBasis<> basisH(mpH);
    gsInfo<<"Basis Primal: "<<basisL.basis(0)<<"\n";
    gsInfo<<"Basis Dual:   "<<basisH.basis(0)<<"\n";


    gsInfo << "Patches: "<< mpH.nPatches() <<", degree: "<< basisL.minCwiseDegree() <<"\n";
    //! [Refinement]

    //! [Problem setup]
    gsExprAssembler<> exL(1,1);
    gsExprAssembler<> exH(1,1);
    exL.setOptions(Aopt);
    exH.setOptions(Aopt);

    //gsInfo<<"Active options:\n"<< exL.options() <<"\n";
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // Elements used for numerical integration
    exL.setIntegrationElements(basisH);
    gsExprEvaluator<> evL(exL);

    exH.setIntegrationElements(basisH);
    gsExprEvaluator<> evH(exH);

    // Set the geometry map
    geometryMap G = exL.getMap(mp);
    geometryMap H = exH.getMap(mp);

    // Set the discretization space
    space u = exL.getSpace(basisL);
    space v = exH.getSpace(basisH);

    u.setInterfaceCont(0);
    u.addBc( bc.get("Dirichlet") );

    v.setInterfaceCont(0);
    v.addBc( bc.get("Dirichlet") );

    // Set the source term
    variable ff = exL.getCoeff(f, G);
    variable gg = exH.getCoeff(f, H);

    // Recover manufactured solution
    // gsFunctionExpr<> ms("((x-0.5)*(x+0.5)*(y-0.5)*(y+0.5)) / ((x+2/3)*y+3/4)",2);
    //gsInfo<<"Exact solution: "<< ms << "\n";
    gsFunctionExpr<> msP;
    fd.getId(3, msP); // id=3: reference solution primal
    gsFunctionExpr<> msD;
    fd.getId(9, msD); // id=9: reference solution dual

    variable primal_exL = evL.getVariable(msP, G);
    variable dual_exL   = evL.getVariable(msD, G);
    variable primal_exH = evH.getVariable(msP, H);
    variable dual_exH   = evH.getVariable(msD, H);

    // Solution vector and solution variable
    gsMatrix<> primalL, primalH, primalLp, dualL, dualLp, dualH;

    // Solutions and their projections
        // PDE solution on low-order mesh
    solution uL = exL.getSolution(u, primalL);
        // PDE solution projection on high-order mesh
    solution uH = exH.getSolution(v, primalH);

        // Dual solution on low-order mesh
    solution zL = exL.getSolution(u, dualL);
        // Dual solution projection on high-order mesh
    solution zH = exH.getSolution(v, dualH);

    // zH2
    gsMultiPatch<> zH2_mp(mpL);//just initialize for not being empty
    variable zH2 = exL.getCoeff(zH2_mp);
    // zH2
    gsMultiPatch<> zL2_mp(mpL);//just initialize for not being empty
    variable zL2 = exH.getCoeff(zL2_mp);

    gsSparseSolver<>::CGDiagonal solver;

    exL.initSystem();
    exH.initSystem();

    //! [Problem setup]
    gsInfo<< "NumDofs Primal: "<<exL.numDofs() <<std::flush;

    // [ORIGINAL PROBLEM]
    // Compute the system matrix and right-hand side
    exL.assemble( igrad(u, G) * igrad(u, G).tr() * meas(G), u * ff * meas(G) );

    // Enforce Neumann conditions to right-hand side
    variable g_N = exL.getBdrFunction();
    exL.assembleRhsBc(u * g_N.val() * nv(G).norm(), bc.neumannSides() );
    //gsInfo<<"Sparse Matrix:\n"<< exL.matrix().toDense() <<"\n";
    //gsInfo<<"Rhs vector:\n"<< exL.rhs().transpose() <<"\n";

    // gsDebugVar(exL.matrix().toDense());
    // gsDebugVar(exL.rhs().transpose());

    gsInfo<< ".\n" <<std::flush;// Assemblying done

    solver.compute( exL.matrix() );
    primalL = solver.solve(exL.rhs());

    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        // ev.options().setSwitch("plot.elements", true);
        evL.writeParaview( uL   , G, "solution");
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
    // variable g_N = exL.getBdrFunction();
    // exL.assembleRhsBc(u * g_N.val() * nv(G).norm(), bc.neumannSides() );

    solver.compute( exL.matrix() );
    dualL = solver.solve(exL.rhs());

    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        // ev.options().setSwitch("plot.elements", true);
        evL.writeParaview( zL   , G, "dualL");
        evL.writeParaview( dual_exL   , G, "dual_exact");

    }

    // [!Low-order Dual PROBLEM]

    // [ Project solutions ]
    space v0 = exH.getSpace(basisH);
    space u0 = exL.getSpace(basisL);
    gsMatrix<> dualL0, primalL0;
    zL.extractFull(dualL0);
    uL.extractFull(primalL0);

    gsDWRHelper<real_t> L2projector(basisL,basisH,mp);
    L2projector.projectVector(dualL0,dualLp);
    L2projector.projectVector(primalL0,primalLp);

    solution zLp = exH.getSolution(v0, dualLp);
    solution uLp = exH.getSolution(v0, primalLp);
    exH.initSystem();

    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        // ev.options().setSwitch("plot.elements", true);
        evH.writeParaview( zLp   , H, "dualLp");
        evH.writeParaview( uLp   , H, "primalLp");
    }
    // [! Project solutions ]

    // Initialize the system
    exH.initSystem();

    gsInfo<< "NumDofs Dual: "<<exH.numDofs() <<std::flush;

    // [DUAL PROBLEM]
    // Compute the system matrix and right-hand side
    exH.assemble( igrad(v, H) * igrad(v, H).tr() * meas(H), v * uLp * meas(H) );

    gsInfo<< ".\n" <<std::flush;// Assemblying done

    solver.compute( exH.matrix() );
    dualH = solver.solve(exH.rhs());

    gsMatrix<> dualH0;
    solution zH0 = exH.getSolution(v0,dualH0);

    zH.extractFull(dualH0);
    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        // ev.options().setSwitch("plot.elements", true);
        evH.writeParaview( zH0   , H, "dualH");
    }
    // [!DUAL PROBLEM]

    zH.extract(zH2_mp);// this updates zH2 variable
    zL.extract(zL2_mp);// this updates zH2 variable

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
    gsDebug<<"lapl zL "<<evL.eval(slapl(zL),pt)<<"\n";
    gsDebug<<"lapl zL "<<evH.eval(lapl(zL2),pt)<<"\n";
    gsDebug<<"lapl zH "<<evL.eval(lapl(zH2),pt)<<"\n";
    gsDebug<<"lapl zH "<<evH.eval(slapl(zH),pt)<<"\n";
    gsDebug<<"\n";
    // Integrals via the assembler and partition of unity.

    exL.initSystem();
    exL.assemble(zL.val() * u0);
    gsDebug<<"int zL "<<evL.integral(zL.val())<<"; "<<exL.rhs().sum()<<"\n";

    exH.initSystem();
    exH.assemble(zL2.val() * v0);
    gsDebug<<"int zH "<<evH.integral(zL2)<<"; "<<exH.rhs().sum()<<"\n";

    exL.initSystem();
    exL.assemble(zH2.val() * u0);
    gsDebug<<"int zH "<<evL.integral(zH2)<<"; "<<exL.rhs().sum()<<"\n";

    exH.initSystem();
    exH.assemble(zH.val() * v0);
    gsDebug<<"int zH "<<evH.integral(zH.val())<<"; "<<exH.rhs().sum()<<"\n";

    gsInfo<<"Objective function errors J(u)-J(u_h)\n";
    real_t error = evL.integral((0.5*primal_exL*primal_exL)*meas(G)) - evL.integral((0.5*uL*uL)*meas(G));
    real_t errest = evH.integral( (zH-zLp) * gg * meas(H)-(((igrad(zH) - igrad(zLp))*igrad(uLp).tr()) ) * meas(H));
    // real_t errest = evL.integral( (zH2-zL) * ff * meas(G)-(((igrad(zH2) - igrad(zL))*igrad(uL).tr()) ) * meas(G));
    gsInfo<<"\texact:\t"   <<error<<"\n"
            "\testimate:\t"<<errest<<"\n"
            "\teff:\t"<<errest/error<<"\n";

    gsInfo<<"Computation errors (L2-norm)\n";
    real_t err1, err2;
    err1 = math::sqrt(evL.integral((primal_exL-uL).sqNorm()*meas(G)));
    err2 = math::sqrt(evH.integral((primal_exH-uLp).sqNorm()*meas(H)));
    gsInfo<<"\t (1) primal exact: "<<err1<<"\n"
            "\t (2) primal proje: "<<err2<<"\n"
            "\t ((1)-(2)/(1)    : "<<(err1-err2)/err1<<"\n";

    err1 = math::sqrt(evL.integral((dual_exL-zL).sqNorm()*meas(G)));
    err2 = math::sqrt(evH.integral((dual_exH-zLp).sqNorm()*meas(H)));
    gsInfo<<"\t (1) dual (L) exact: "<<err1<<"\n"
            "\t (2) dual (L) proje: "<<err2<<"\n"
            "\t ((1)-(2)/(1)      : "<<(err1-err2)/err1<<"\n"
            "\t dual (H) exact: "<<math::sqrt(evH.integral((dual_exH-zH).sqNorm()*meas(H)))<<"\n";

    // ELEMENT WISE ERROR ESTIMATION
    if (est==0)
    {
        gsInfo<<evL.eval((ff - slapl(uL)) * (zH2.val() - zL.val())*meas(G),pt)<<"\n";
        gsInfo<<evH.eval((gg - slapl(uLp)) * (zH.val() - zLp.val())*meas(H),pt)<<"\n";



        evL.integralElWise((ff - slapl(uL)) * (zH2.val() - zL.val())*meas(G));
        // evH.integralElWise((gg - slapl(uL)) * (zH2.val() - zL.val())*meas(G));

        // evH.integralElWise((gg - slapl(uLp)) * (zH - zLp)*meas(H));
        // evH.integralElWise( (zH-zLp) * gg * meas(H)-(((igrad(zH) - igrad(zLp))*igrad(uLp).tr()) ) * meas(H));
        gsVector<> elementNorms = evL.allValues().transpose();
        std::vector<real_t> errors;
        errors.resize(elementNorms.size());
        gsVector<>::Map(&errors[0],elementNorms.size() ) = elementNorms;

        gsInfo<< "  Result (global)    : "<< elementNorms.sum() <<"\n";
    return 0;

        gsElementErrorPlotter<real_t> err_eh(basisH.basis(0),errors);

        const gsField<> elemError_eh( mp, err_eh, false );
        gsWriteParaview<>( elemError_eh, "error_elem_eh", 1000);


        MarkingStrategy adaptRefCrit = PUCA;
        const real_t adaptRefParam = 0.9;
        std::vector<bool> elMarked( errors.size() );

        gsMarkElementsForRef( errors, adaptRefCrit, adaptRefParam, elMarked);
        for (std::vector<bool>::const_iterator i = elMarked.begin(); i != elMarked.end(); ++i)
            gsInfo << *i << ' ';
        gsInfo<<"\n";

        // Refine the marked elements with a 1-ring of cells around marked elements
        gsRefineMarkedElements( mpH, elMarked, 1 );
        gsRefineMarkedElements( mpL, elMarked, 1 );

        gsDebugVar(mpH.basis(0));
        gsDebugVar(mpL.basis(0));


    }
    // FUNCTION WISE ERROR ESTIMATION
    else if (est==1)
    {
        exH.initSystem(true);
        // exH.assemble( grad(uLp) * ( grad(v0) * (zH - zLp) + v0 * ( grad(zH) - grad(zLp) ) ) - gg * v0 * (zH - zLp)  );

        auto lhs = ((zH - zLp) * grad(uLp) * grad(v0).tr()).tr(); // + v0 * ( grad(zH) - grad(zLp) ) * grad(uLp).tr();
        auto lhs2 = v0 * ( grad(zH) - grad(zLp) ) * grad(uLp).tr();
        auto rhs = gg.val() * (zH - zLp).val() * v0;

        gsDebug<<evH.eval(zH,pt);


        gsMatrix<> res;
        exH.assemble(lhs*meas(H));
        res = exH.rhs();

        exH.assemble(lhs2*meas(H));
        res += exH.rhs();

        exH.assemble(rhs*meas(H));
        res += exH.rhs();

        gsDebugVar(res);

        gsInfo<< "  Result (global)    : "<< res.sum()<<"\n";

        MarkingStrategy adaptRefCrit = PUCA;
        const real_t adaptRefParam = 0.9;
        std::vector<bool> funMarked( res.size() );

        std::vector<real_t> errors( res.size() );
        gsVector<>::Map(&errors[0],res.size() ) = res;


        for (std::vector<real_t>::const_iterator i = errors.begin(); i != errors.end(); ++i)
            gsInfo << *i << "\n";
        gsInfo<<"\n";

        gsMarkElementsForRef( errors, adaptRefCrit, adaptRefParam, funMarked);
        for (std::vector<bool>::const_iterator i = funMarked.begin(); i != funMarked.end(); ++i)
            gsInfo << *i << ' ';
        gsInfo<<"\n";

        // Refine the marked elements with a 1-ring of cells around marked elements
        gsRefineMarkedFunctions( mpH, funMarked, 1 );
        gsDebugVar(mpH.basis(0));
    gsWriteParaview(mpH,"mpH",1000,true);


        gsMultiPatch<> tmp = mpH;
        tmp.patch(0).degreeReduce(1);
        mpL = tmp;
        // gsRefineMarkedFunctions( mpL, funMarked, 1 );
        gsDebugVar(mpL.basis(0));

    }
    // FUNCTION WISE ERROR ESTIMATION
    else if (est==2)
    {
        exL.initSystem(true);
        exH.initSystem(true);


        // exH.assemble( grad(uLp) * ( grad(v0) * (zH - zLp) + v0 * ( grad(zH) - grad(zLp) ) ) - gg * v0 * (zH - zLp)  );

        auto lhs = ((zH - zL) * grad(uL) * grad(u0).tr()).tr(); // + v0 * ( grad(zH) - grad(zLp) ) * grad(uLp).tr();
        auto lhs2 = u0 * ( grad(zH) - grad(zL) ) * grad(uL).tr();
        auto rhs = ff.val() * (zH - zL).val() * u0;

        gsDebug<<evL.eval(zL,pt)<<"\n";

        gsDebug<<evL.eval(zH,pt)<<"\n";

        gsDebug<<evH.integral(zH)<<"\n";

        exL.assemble(zH.val() * u0 * meas(G));
        gsDebugVar(exL.rhs().sum());
        gsDebug<<evH.integral(zH)<<"\n";
        exH.assemble(zH.val() * v0 * meas(H));
        gsDebugVar(exH.rhs().sum());

        exL.assemble(zL.val() * u0 * meas(G));
        gsDebugVar(exL.rhs().sum());


    //     gsDebug<<evL.eval(zL,pt)<<"\n";
    //     gsDebug<<evL.eval(zH.temp(),pt)<<"\n";
    //     gsDebug<<evH.eval(zH,pt)<<"\n";

    //     gsMatrix<> res;
    //     exL.assemble(lhs*meas(G));
    //     res = exL.rhs();

    //     exL.assemble(lhs2*meas(G));
    //     res += exL.rhs();

    //     exL.assemble(rhs*meas(G));
    //     res += exL.rhs();

    //     gsDebugVar(res);

    //     gsInfo<< "  Result (global)    : "<< res.sum()<<"\n";

    //     MarkingStrategy adaptRefCrit = PUCA;
    //     const real_t adaptRefParam = 0.9;
    //     std::vector<bool> funMarked( res.size() );

    //     std::vector<real_t> errors( res.size() );
    //     gsVector<>::Map(&errors[0],res.size() ) = res;


    //     for (std::vector<real_t>::const_iterator i = errors.begin(); i != errors.end(); ++i)
    //         gsInfo << *i << "\n";
    //     gsInfo<<"\n";

    //     gsMarkElementsForRef( errors, adaptRefCrit, adaptRefParam, funMarked);
    //     for (std::vector<bool>::const_iterator i = funMarked.begin(); i != funMarked.end(); ++i)
    //         gsInfo << *i << ' ';
    //     gsInfo<<"\n";

    //     // Refine the marked elements with a 1-ring of cells around marked elements
    //     gsRefineMarkedFunctions( mpL, funMarked, 1 );
    //     gsDebugVar(mpL.basis(0));
    // gsWriteParaview(mpL,"mpH",1000,true);


    //     gsMultiPatch<> tmp = mpH;
    //     tmp.patch(0).degreeReduce(1);
    //     mpL = tmp;
    //     // gsRefineMarkedFunctions( mpL, funMarked, 1 );
    //     gsDebugVar(mpL.basis(0));

    }
    else
        GISMO_ERROR("Estimation method unknown");




    gsWriteParaview(mpH,"mpH",1000,true);
    gsWriteParaview(mpL,"mpL",1000,true);

    return EXIT_SUCCESS;

}// end main
