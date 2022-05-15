/** @file gsNewton_test.cpp

    @brief Tests the newton iteration

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
 **/

//#define TEST_INFO

#include "gismo_unittest.h"

#ifdef GISMO_KLSHELL
#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/getMaterialMatrix.h>
#endif

SUITE(gsKLShell_test)
{
#ifdef GISMO_KLSHELL

    gsVector<real_t> modal_numerical(bool composite)
    {
        // Input options
        int numElevate  = 2;
        int numHref     = 4;

        real_t thickness = 0.01;
        real_t E_modulus = 1e5;
        real_t Density = 1e0;
        real_t PoissonRatio = 0.3;

        gsMultiPatch<> mp;

        mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
        mp.addAutoBoundaries();
        mp.embed(3);

        for(index_t i = 0; i< numElevate; ++i)
            mp.patch(0).degreeElevate();    // Elevate the degree

        // h-refine
        for(index_t i = 0; i< numHref; ++i)
            mp.patch(0).uniformRefine();

        gsMultiBasis<> dbasis(mp);

        // Boundary conditions
        gsBoundaryConditions<> BCs;

        // Plate
        // Pinned-Pinned-Pinned-Pinned
            // Left
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
            // Right
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
            // Top
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
            // Bottom
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 0 - x
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        BCs.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        BCs.setGeoMap(mp);

        // Initialise solution object
        gsMultiPatch<> mp_def = mp;

        // Linear isotropic material model
        gsVector<> tmp(3);
        tmp << 0, 0, 0;
        gsConstantFunction<> force(tmp,3);
        gsFunctionExpr<> t(std::to_string(thickness), 3);
        gsFunctionExpr<> E(std::to_string(E_modulus),3);
        gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
        gsConstantFunction<> rho(Density,3);

        // Linear anisotropic material model
        real_t pi = math::atan(1)*4;
        index_t kmax = 5;

        std::vector<gsFunctionSet<> * > Gs(kmax);
        std::vector<gsFunctionSet<> * > Ts(kmax);
        std::vector<gsFunctionSet<> * > Rs(kmax);
        std::vector<gsFunctionSet<> * > Phis(kmax);

        Rs[0] = Rs[1] = Rs[2] = Rs[3] = Rs[4] = &rho;


        gsMatrix<> Gmat = gsCompositeMatrix(E_modulus,E_modulus,0.5 * E_modulus / (1+PoissonRatio),PoissonRatio,PoissonRatio);
        Gmat.resize(Gmat.rows()*Gmat.cols(),1);
        gsConstantFunction<> Gfun(Gmat,3);
        Gs[0] = Gs[1] = Gs[2] = Gs[3] = Gs[4] = &Gfun;

        gsConstantFunction<> phi1, phi2, phi3, phi4, phi5;
        phi1.setValue(0/kmax * pi / 2.0,3);
        phi2.setValue(1/kmax * pi / 2.0,3);
        phi3.setValue(2/kmax * pi / 2.0,3);
        phi4.setValue(3/kmax * pi / 2.0,3);
        phi5.setValue(4/kmax * pi / 2.0,3);

        Phis[0] = &phi1;
        Phis[1] = &phi2;
        Phis[2] = &phi3;
        Phis[3] = &phi4;
        Phis[4] = &phi5;

        gsConstantFunction<> thicks(thickness/kmax,3);
        Ts[0] = Ts[1] = Ts[2] = Ts[3] = Ts[4] = &thicks;

        std::vector<gsFunction<>*> parameters;
        gsMaterialMatrixBase<real_t>* materialMatrix;

        gsOptionList options;

        if (composite)
        {
            materialMatrix = new gsMaterialMatrixComposite<3,real_t>(mp,Ts,Gs,Phis,Rs);
        }
        else
        {
            parameters.resize(2);
            parameters[0] = &E;
            parameters[1] = &nu;
            options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
            options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
            materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);
        }

        // options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
        // options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",static_cast<int>(!composite));
        // materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);

        gsThinShellAssemblerBase<real_t>* assembler;
        assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,BCs,force,materialMatrix);

        assembler->assemble();
        gsSparseMatrix<> K =  assembler->matrix();
        assembler->assembleMass();
        gsSparseMatrix<> M =  assembler->matrix();

        Eigen::GeneralizedSelfAdjointEigenSolver< gsMatrix<real_t>::Base >  eigSolver;
        eigSolver.compute(K,M);
        gsMatrix<> values  = eigSolver.eigenvalues();
        gsMatrix<> vectors = eigSolver.eigenvectors();

        values = values.cwiseSqrt();
        values = values.col(0).head(10);

        delete materialMatrix;
        delete assembler;

        return values;
    }

    gsVector<real_t> modal_analytical()
    {
        real_t thickness = 0.01;
        real_t E_modulus = 1e5;
        real_t Density = 1e0;
        real_t PoissonRatio = 0.3;

        real_t D = E_modulus*math::pow(thickness,3)/(12*(1-math::pow(PoissonRatio,2)));

        std::vector<real_t> omegas;
        for (index_t m=1; m!=10; m++)
          for (index_t n=1; n!=10; n++)
            omegas.push_back((math::pow(m/1.0,2)+math::pow(n/1.0,2))*math::pow(3.1415926535,2)*math::sqrt(D / (Density * thickness)));

        std::sort(omegas.begin(),omegas.end());
        omegas.resize(10);
        gsAsVector<> analytical(omegas);

        return analytical;
    }

    TEST(modal) // Declares test
    {
        // UNITTEST_TIME_CONSTRAINT(1000);// this will produce failure if test takes more than 1 sec

        gsVector<> an = modal_analytical();
        gsVector<> num, relError;
        real_t error;

        std::vector<bool> composite { true, false };

        for (std::vector<bool>::iterator comp = composite.begin(); comp!=composite.end(); comp++)
        {
            num = modal_numerical(*comp);
            relError = (num - an).array()/an.array();
            error = relError.norm();

            CHECK(error < 1e-3);
        }
    }

    std::pair<real_t,real_t> balloon_numerical(index_t material, index_t impl)
    {
        //! [Parse command line]
        index_t numRefine  = 1;
        index_t numElevate = 1;

        real_t E_modulus = 1.0;
        real_t Density = 1.0;
        real_t Ratio = 7.0;

        real_t thickness = 0.1;
        real_t mu = 4.225e5;

        real_t alpha1,alpha2,alpha3,mu1,mu2,mu3;
        alpha1 = 1.3;
        mu1    = 6.3e5/4.225e5*mu;
        alpha2 = 5.0;
        mu2    = 0.012e5/4.225e5*mu;
        alpha3 = -2.0;
        mu3    = -0.1e5/4.225e5*mu;

        real_t PoissonRatio = 0.5;
        E_modulus = 2*mu*(1+PoissonRatio);

        //! [Read input file]
        gsMultiPatch<> mp, mp_def;

        mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
        mp.embed(3);

        gsReadFile<>("surface/eighth_sphere.xml", mp);

        for(index_t i = 0; i< numElevate; ++i)
          mp.patch(0).degreeElevate();    // Elevate the degree

        // h-refine
        for(index_t i = 0; i< numRefine; ++i)
          mp.patch(0).uniformRefine();

        mp_def = mp;

        //! [Refinement]
        gsMultiBasis<> dbasis(mp);

        gsBoundaryConditions<> bc;
        bc.setGeoMap(mp);

        GISMO_ASSERT(mp.targetDim()==3,"Geometry must be surface (targetDim=3)!");
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 2 - z
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 2 - z

        bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        // Symmetry in x-direction:
        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 );
        bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 1 );
        bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 );

        // Symmetry in y-direction:
        bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 0 );
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 );
        bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );

        // Pressure
        real_t pressure = 10e3;

        //! [Refinement]

        // Linear isotropic material model
        gsVector<> tmp(3);
        tmp.setZero();
        gsConstantFunction<> force(tmp,3);
        gsConstantFunction<> pressFun(pressure,3);
        gsFunctionExpr<> t(std::to_string(thickness), 3);
        gsFunctionExpr<> E(std::to_string(E_modulus),3);
        gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
        gsFunctionExpr<> rho(std::to_string(Density),3);
        gsConstantFunction<> ratio(Ratio,3);

        gsConstantFunction<> alpha1fun(alpha1,3);
        gsConstantFunction<> mu1fun(mu1,3);
        gsConstantFunction<> alpha2fun(alpha2,3);
        gsConstantFunction<> mu2fun(mu2,3);
        gsConstantFunction<> alpha3fun(alpha3,3);
        gsConstantFunction<> mu3fun(mu3,3);

        std::vector<gsFunction<>*> parameters(3);
        parameters[0] = &E;
        parameters[1] = &nu;
        parameters[2] = &ratio;
        gsMaterialMatrixBase<real_t>* materialMatrix;

        if (material==4)
        {
            parameters.resize(8);
            parameters[0] = &E;
            parameters[1] = &nu;
            parameters[2] = &mu1fun;
            parameters[3] = &alpha1fun;
            parameters[4] = &mu2fun;
            parameters[5] = &alpha2fun;
            parameters[6] = &mu3fun;
            parameters[7] = &alpha3fun;
        }

        gsOptionList options;
        if      (material==0)
        {
            GISMO_ERROR("This test is not available for SvK models");
        }
        else
        {
            options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",material);
            options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",false);
            options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",impl);
            materialMatrix = getMaterialMatrix<3,real_t>(mp,t,parameters,rho,options);
        }

        gsThinShellAssemblerBase<real_t>* assembler;
        assembler = new gsThinShellAssembler<3, real_t, true >(mp,dbasis,bc,force,materialMatrix);

        assembler->setPressure(pressFun);

        // Function for the Jacobian
        typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>    Jacobian_t;
        typedef std::function<gsVector<real_t> (gsVector<real_t> const &) >         Residual_t;
        Jacobian_t Jacobian = [&assembler,&mp_def](gsVector<real_t> const &x)
        {
          assembler->constructSolution(x,mp_def);
          assembler->assembleMatrix(mp_def);
          gsSparseMatrix<real_t> m = assembler->matrix();
          return m;
        };
        // Function for the Residual
        Residual_t Residual = [&assembler,&mp_def](gsVector<real_t> const &x)
        {
          assembler->constructSolution(x,mp_def);
          assembler->assembleVector(mp_def);
          return assembler->rhs();
        };

        // Define Matrices
        assembler->assemble();

        gsSparseMatrix<> matrix = assembler->matrix();
        gsVector<> vector = assembler->rhs();

        // Solve linear problem
        gsVector<> solVector;
        gsSparseSolver<>::CGDiagonal solver;
        solver.compute( matrix );
        solVector = solver.solve(vector);

        gsVector<real_t> updateVector = solVector;
        gsVector<real_t> resVec = Residual(solVector);
        gsSparseMatrix<real_t> jacMat;
        for (index_t it = 0; it != 100; ++it)
        {
            jacMat = Jacobian(solVector);
            solver.compute(jacMat);
            updateVector = solver.solve(resVec); // this is the UPDATE
            solVector += updateVector;

            resVec = Residual(solVector);

            if (updateVector.norm() < 1e-6)
                break;
            else if (it+1 == it)
                gsWarn<<"Maximum iterations reached!\n";
        }

        mp_def = assembler->constructSolution(solVector);

        gsMultiPatch<> deformation = mp_def;
        for (size_t k = 0; k != mp_def.nPatches(); ++k)
            deformation.patch(k).coefs() -= mp.patch(k).coefs();

        gsMatrix<> pt(2,2);
        pt.col(0)<<0.0,1.0;
        pt.col(1)<<1.0,1.0;


        gsMatrix<> lambdas = assembler->computePrincipalStretches(pt,mp_def,0);

        real_t tol = 10e-10;
        GISMO_ASSERT((lambdas.col(0)-lambdas.col(1)).norm() < tol, "Stretches must be equal over the balloon");

        real_t r = (mp_def.patch(0).eval(pt).col(0)).norm();

        // Get the total force on the tension boundary
        real_t P = pressure * assembler->getArea(mp) / assembler->getArea(mp_def);

        std::pair<real_t,real_t> result;
        result.first = P;
        result.second = r;

        delete materialMatrix;
        delete assembler;

        return result;
    }

    real_t balloon_analytical(index_t material, index_t impl, real_t r)
    {
        real_t Pan;
        real_t R = 10.;

        real_t Ratio = 7.0;

        real_t thickness = 0.1;
        real_t mu = 4.225e5;

        real_t alpha1,alpha2,alpha3,mu1,mu2,mu3;
        alpha1 = 1.3;
        mu1    = 6.3e5/4.225e5*mu;
        alpha2 = 5.0;
        mu2    = 0.012e5/4.225e5*mu;
        alpha3 = -2.0;
        mu3    = -0.1e5/4.225e5*mu;

        real_t lambda = r/R;

        if      (material==1)
        {
            // Pan = 2*(thickness/R)*(mu*(1.0/lambdas(0)-lambdas(0)));
            Pan = 2*(thickness/R)*(mu*(math::pow(lambda,2-3)-math::pow(lambda,-2*2-3)));
        }
        else if (material==3)
        {
            real_t c2 = 1.0 / (Ratio+1);
            real_t c1 = 1.0 - c2;
            real_t m1 = c1*mu;
            real_t m2 = -c2*mu;
            real_t a1 = 2;
            real_t a2 = -2;
            Pan = 2*(thickness/R)*(m1*(math::pow(lambda,a1-3)-math::pow(lambda,-2*a1-3))+m2*(math::pow(lambda,a2-3)-math::pow(lambda,-2*a2-3)));
        }
        else if (material==4)
        {
            Pan=2*(thickness/R)*(
                mu1*(math::pow(lambda,alpha1-3)-math::pow(lambda,-2*alpha1-3))
                +mu2*(math::pow(lambda,alpha2-3)-math::pow(lambda,-2*alpha2-3))
                +mu3*(math::pow(lambda,alpha3-3)-math::pow(lambda,-2*alpha3-3)) );
        }
        else
            GISMO_ERROR("Material not treated");

        return Pan;
    }

    TEST(balloon) // Declares test
    {
        // UNITTEST_TIME_CONSTRAINT(1000);// this will produce failure if test takes more than 1 sec

        real_t P, Pan, rnum;

        std::vector<index_t> materials{ 1,3,4 };
        std::vector<index_t> implementations{ 1,2,3 };

        std::pair<real_t,real_t> num;
        real_t tol = 1e-3;
        for (std::vector<index_t>::iterator mat = materials.begin(); mat!=materials.end(); mat++)
        {
            for (std::vector<index_t>::iterator impl = implementations.begin(); impl!=implementations.end(); impl++)
            {
                if (*mat==4 && *impl!=3) continue;

                num = balloon_numerical(*mat,*impl);
                P = num.first;
                rnum = num.second;

                Pan = balloon_analytical(*mat,*impl,rnum);

                gsTestInfo<<"L error = "<<std::abs(P-Pan)/Pan<<"\n";
                CHECK((std::abs(P-Pan)/Pan < tol));
            }
        }
    }

    std::pair<real_t,real_t> UAT_numerical(index_t material, index_t impl, bool Compressibility)
    {
        //! [Parse command line]
        index_t numRefine  = 1;
        index_t numElevate = 1;

        real_t E_modulus = 1.0;
        real_t PoissonRatio;
        real_t Density = 1.0;
        real_t Ratio = 7.0;

        real_t mu = 1.5e6;
        real_t thickness = 0.001;

        real_t alpha1,alpha2,alpha3,mu1,mu2,mu3;
        alpha1 = 1.3;
        mu1    = 6.3e5/4.225e5*mu;
        alpha2 = 5.0;
        mu2    = 0.012e5/4.225e5*mu;
        alpha3 = -2.0;
        mu3    = -0.1e5/4.225e5*mu;

        if (!Compressibility)
          PoissonRatio = 0.5;
        else
          PoissonRatio = 0.45;

        E_modulus = 2*mu*(1+PoissonRatio);

        //! [Parse command line]

        //! [Read input file]
        gsMultiPatch<> mp, mp_def;

        mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree

        if (numElevate!=0)
            mp.degreeElevate(numElevate);

        // h-refine
        for (int r =0; r < numRefine; ++r)
            mp.uniformRefine();

        mp_def = mp;

        //! [Refinement]
        gsMultiBasis<> dbasis(mp);

        gsBoundaryConditions<> bc;
        bc.setGeoMap(mp);

        gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();

        real_t lambda = 2.0;
        gsConstantFunction<> displx(lambda-1.0,2);

        GISMO_ASSERT(mp.targetDim()==2,"Geometry must be planar (targetDim=2)!");
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0 );

        bc.addCondition(boundary::east, condition_type::dirichlet, &displx, 0, false, 0 );

        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 1 );

        //! [Refinement]

        // Linear isotropic material model
        gsVector<> tmp(2);
        tmp.setZero();
        gsConstantFunction<> force(tmp,2);
        gsFunctionExpr<> t(std::to_string(thickness),2);
        gsFunctionExpr<> E(std::to_string(E_modulus),2);
        gsFunctionExpr<> nu(std::to_string(PoissonRatio),2);
        gsFunctionExpr<> rho(std::to_string(Density),2);
        gsConstantFunction<> ratio(Ratio,2);

        gsConstantFunction<> alpha1fun(alpha1,2);
        gsConstantFunction<> mu1fun(mu1,2);
        gsConstantFunction<> alpha2fun(alpha2,2);
        gsConstantFunction<> mu2fun(mu2,2);
        gsConstantFunction<> alpha3fun(alpha3,2);
        gsConstantFunction<> mu3fun(mu3,2);

        std::vector<gsFunction<>*> parameters(3);
        parameters[0] = &E;
        parameters[1] = &nu;
        parameters[2] = &ratio;
        gsMaterialMatrixBase<real_t>* materialMatrix;

        if (material==4)
        {
            parameters.resize(8);
            parameters[0] = &E;
            parameters[1] = &nu;
            parameters[2] = &mu1fun;
            parameters[3] = &alpha1fun;
            parameters[4] = &mu2fun;
            parameters[5] = &alpha2fun;
            parameters[6] = &mu3fun;
            parameters[7] = &alpha3fun;
        }

        gsOptionList options;
        if      (material==0)
        {
            GISMO_ERROR("This test is not available for SvK models");
        }
        else
        {
            options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",material);
            options.addSwitch("Compressibility","Compressibility: (false): Imcompressible | (true): Compressible",Compressibility);
            options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",impl);
            materialMatrix = getMaterialMatrix<2,real_t>(mp,t,parameters,rho,options);
        }

        gsThinShellAssemblerBase<real_t>* assembler;
        assembler = new gsThinShellAssembler<2, real_t, false >(mp,dbasis,bc,force,materialMatrix);

        assembler->setPointLoads(pLoads);

        // Function for the Jacobian
        typedef std::function<gsSparseMatrix<real_t> (gsVector<real_t> const &)>    Jacobian_t;
        typedef std::function<gsVector<real_t> (gsVector<real_t> const &) >         Residual_t;
        Jacobian_t Jacobian = [&assembler,&mp_def](gsVector<real_t> const &x)
        {
          assembler->constructSolution(x,mp_def);
          assembler->assembleMatrix(mp_def);
          gsSparseMatrix<real_t> m = assembler->matrix();
          return m;
        };
        // Function for the Residual
        Residual_t Residual = [&assembler,&mp_def](gsVector<real_t> const &x)
        {
          assembler->constructSolution(x,mp_def);
          assembler->assembleVector(mp_def);
          return assembler->rhs();
        };

        // Define Matrices
        assembler->assemble();

        gsSparseMatrix<> matrix = assembler->matrix();
        gsVector<> vector = assembler->rhs();

        // Solve linear problem
        gsVector<> solVector;
        gsSparseSolver<>::CGDiagonal solver;
        solver.compute( matrix );
        solVector = solver.solve(vector);

        gsVector<real_t> updateVector = solVector;
        gsVector<real_t> resVec = Residual(solVector);
        gsSparseMatrix<real_t> jacMat;
        for (index_t it = 0; it != 100; ++it)
        {
            jacMat = Jacobian(solVector);
            solver.compute(jacMat);
            updateVector = solver.solve(resVec); // this is the UPDATE
            solVector += updateVector;

            resVec = Residual(solVector);

            if (updateVector.norm() < 1e-6)
                break;
            else if (it+1 == it)
                gsWarn<<"Maximum iterations reached!\n";
        }

        mp_def = assembler->constructSolution(solVector);

        gsMultiPatch<> deformation = mp_def;
        for (size_t k = 0; k != mp_def.nPatches(); ++k)
            deformation.patch(k).coefs() -= mp.patch(k).coefs();

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Check solutions
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // NOTE: all analytical solutions for compressible materials are fixed for displ=1; (lambda=2)

        // Compute stretches (should be the same everywhere)
        // Ordering: lambda(0) < lambda(1); lambda(2) is ALWAYS the through-thickness stretch
        gsVector<> pt(2);
        pt<<1,0;
        gsMatrix<> lambdas = assembler->computePrincipalStretches(pt,mp_def,0);

        // Get the total force on the tension boundary
        patchSide ps(0,boundary::east);
        gsMatrix<> forceVector = assembler->boundaryForceVector(mp_def,ps,0);
        real_t sideForce = forceVector.sum();
        real_t S   = -sideForce / (thickness*lambdas(0)*lambdas(2));
        real_t L   = lambdas(0);

        std::pair<real_t,real_t> result;
        result.first = L;
        result.second = S;

        delete materialMatrix;
        delete assembler;

        return result;
    }

    std::pair<real_t,real_t> UAT_analytical(index_t material, index_t impl, bool Compressibility)
    {
        real_t PoissonRatio;
        real_t Ratio = 7.0;

        real_t mu = 1.5e6;

        real_t alpha1,alpha2,alpha3,mu1,mu2,mu3;
        alpha1 = 1.3;
        mu1    = 6.3e5/4.225e5*mu;
        alpha2 = 5.0;
        mu2    = 0.012e5/4.225e5*mu;
        alpha3 = -2.0;
        mu3    = -0.1e5/4.225e5*mu;

        if (!Compressibility)
          PoissonRatio = 0.5;
        else
          PoissonRatio = 0.45;

        real_t lambda = 2.0;

        real_t San,J,K,Lan;
        if      (material==1 && Compressibility)
        {
            K = 2*mu*(1+PoissonRatio)/(3-6*PoissonRatio);
            J = 1.105598565;// specific for lambda==2!!
            San = lambda*(0.5*mu*(-(2*(math::pow(lambda,2)+2*J/lambda))/(3*math::pow(J,2./3.)*lambda)+2*lambda/math::pow(J,2./3.))+0.25*K*(2*math::pow(J,2)/lambda-2./lambda))/J;
            Lan = math::pow(J/lambda,0.5);
        }
        else if (material==1 && !Compressibility)
        {
            San = mu * (lambda*lambda - 1/lambda);
            Lan = math::pow(1./lambda,0.5);
        }
        else if (material==3 && Compressibility)
        {
            real_t c2 = 1.0 / (Ratio+1);
            real_t c1 = 1.0 - c2;
            K = 2*mu*(1+PoissonRatio)/(3-6*PoissonRatio);
            J = 1.099905842;// specific for lambda==2!!
            San = lambda*(0.5*c1*mu*(-(2*(math::pow(lambda,2)+2*J/lambda))/(3*math::pow(J,2./3.)*lambda)+2*lambda/math::pow(J,2./3.))+0.5*c2*mu*(-(4*(2*lambda*J+math::pow(J,2)/math::pow(lambda,2)))/(3*math::pow(J,4./3.)*lambda)+4/math::pow(J,1./3.))+0.25*K*(2*math::pow(J,2)/lambda-2/lambda))/J;
            Lan = math::pow(J/lambda,0.5);
        }
        else if (material==3 && !Compressibility)
        {
            real_t c2 = 1.0 / (Ratio+1);
            real_t c1 = 1.0 - c2;
            San =-mu*(c2*lambda*lambda+c2/lambda+c1)/lambda+lambda*(c1*lambda*mu+2*c2*mu);
            Lan = math::pow(1./lambda,0.5);
        }
        else if (material==4 && Compressibility)
        {
            K = 2*mu*(1+PoissonRatio)/(3-6*PoissonRatio);
            J = 1.088778638;// specific for lambda==2!!
            San = 1./J* (lambda *( mu1*(2*math::pow(lambda/math::pow(J,1./3.),alpha1)*alpha1/(3*lambda)-2*math::pow(math::pow(J/lambda,0.5)/math::pow(J,1./3.),alpha1)*alpha1/(3*lambda))/alpha1+mu2*(2*math::pow(lambda/math::pow(J,1./3.),alpha2)*alpha2/(3*lambda)-2*math::pow(math::pow(J/lambda,0.5)/math::pow(J,1./3.),alpha2)*alpha2/(3*lambda))/alpha2+mu3*(2*math::pow(lambda/math::pow(J,1./3.),alpha3)*alpha3/(3*lambda)-2*math::pow(math::pow(J/lambda,0.5)/math::pow(J,1./3.),alpha3)*alpha3/(3*lambda))/alpha3+0.25*K*(2*math::pow(J,2)/lambda-2/lambda) ) );
            Lan = math::pow(J/lambda,0.5);
        }
        else if (material==4 && !Compressibility)
        {
            San =-mu1*math::pow((1./lambda),0.5*alpha1)-mu2*math::pow((1./lambda),0.5*alpha2)-mu3*math::pow((1./lambda),0.5*alpha3)+mu1*math::pow(lambda,alpha1)+mu2*math::pow(lambda,alpha2)+mu3*math::pow(lambda,alpha3);
            Lan = math::pow(1./lambda,0.5);
        }
        else
            GISMO_ERROR("Material not treated");

        std::pair<real_t,real_t> result;
        result.first = Lan;
        result.second = San;
        return result;
    }

    // Choose among various shell examples, default = Thin Plate
    int main(int argc, char *argv[])
    {

        real_t S, San,L,Lan;

        std::vector<index_t> materials{ 1,3,4 };
        std::vector<index_t> implementations{ 1,2,3 };
        std::vector<bool> compressibility{ true,false };

        std::pair<real_t,real_t> num, an;
        real_t tol = 1e-9;
        for (std::vector<index_t>::iterator mat = materials.begin(); mat!=materials.end(); mat++)
        {
            for (std::vector<index_t>::iterator impl = implementations.begin(); impl!=implementations.end(); impl++)
            {
                for (std::vector<bool>::iterator comp = compressibility.begin(); comp!=compressibility.end(); comp++)
                {
                    if (*mat==4 && *impl!=3) continue;
                    num = numerical(*mat,*impl,*comp);
                    L = num.first;
                    S = num.second;

                    an = analytical(*mat,*impl,*comp);
                    Lan = an.first;
                    San = an.second;


                    if ( (std::abs(L-Lan)/Lan < tol) && (std::abs(S-San)/San < tol) )
                        gsInfo<<"Passed\n";
                    else
                        gsInfo<<"Failed; L error = "<<std::abs(L-Lan)/Lan<<"\t S error = "<<std::abs(S-San)/San<<"\n";
                }
            }
        }

        return EXIT_SUCCESS;

    }// end main


    TEST(UAT) // Declares test
    {
        // UNITTEST_TIME_CONSTRAINT(1000);// this will produce failure if test takes more than 1 sec

        real_t S, San,L,Lan;

        std::vector<index_t> materials{ 1,3,4 };
        std::vector<index_t> implementations{ 1,2,3 };
        std::vector<bool> compressibility{ true,false };

        std::pair<real_t,real_t> num, an;
        real_t tol = 1e-9;
        for (std::vector<index_t>::iterator mat = materials.begin(); mat!=materials.end(); mat++)
        {
            for (std::vector<index_t>::iterator impl = implementations.begin(); impl!=implementations.end(); impl++)
            {
                for (std::vector<bool>::iterator comp = compressibility.begin(); comp!=compressibility.end(); comp++)
                {
                    if (*mat==4 && *impl!=3) continue;
                    num = UAT_numerical(*mat,*impl,*comp);
                    L = num.first;
                    S = num.second;

                    an = UAT_analytical(*mat,*impl,*comp);
                    Lan = an.first;
                    San = an.second;

                    gsTestInfo<<"L error = "<<std::abs(P-Pan)/Pan<<"\n";
                    gsTestInfo<<"S error = "<<std::abs(S-San)/San<<"\n";
                    CHECK((std::abs(L-Lan)/Lan < tol));
                    CHECK((std::abs(S-San)/San < tol));
                }
            }
        }

        return EXIT_SUCCESS;
    }

#endif
}
