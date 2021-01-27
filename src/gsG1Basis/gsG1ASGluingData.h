//
// Created by afarahat on 3/23/20.
//


#pragma once

//#include <gsG1Basis/gsG1ASGluingDataAssembler.h>
#include <gsG1Basis/gsGluingData.h>
#include <gsG1Basis/gsG1ASGluingDataVisitorGlobal.h>
# include <gsG1Basis/gsG1OptionList.h>



namespace gismo
{

template<class T, class Visitor = gsG1ASGluingDataVisitorGlobal<T>>
class gsG1ASGluingData : public gsGluingData<T>
{
public:

    //Empty constructor will set the gluing data for the boundary edges
    gsG1ASGluingData()
    {   setGDEdge();
//        gsInfo << "Solution: " << sol << "\n";
//        gsInfo << "Solution Beta: " << solBeta << "\n";
    }


    gsG1ASGluingData(gsMultiPatch<T> const & mp,
                 gsMultiBasis<T> & mb)
        : gsGluingData<T>(mp, mb)
    {
    // Solve the system for alpha_L and alpha_R (also for beta, which will be splitted)
        refresh();
        assemble();
        solve();
    // Solve the system for beta_L and beta_R
        refreshBeta();
        assembleBeta();
        solveBeta();

//        gsInfo << "Solution: " << sol << "\n";
//        gsInfo << "Solution Beta: " << solBeta << "\n";

        AScondition(mp);
    }


    gsMatrix<> evalAlpha_R(gsMatrix<> points)
    {
        gsMatrix<> ones(1, points.cols());
        ones.setOnes();
        return sol.row(0) * ( ones - points ) + sol.row(1) * points;
    }

    gsMatrix<> evalAlpha_L(gsMatrix<> points)
    {
        gsMatrix<> ones(1, points.cols());
        ones.setOnes();
        return sol.row(2) * ( ones - points ) + sol.row(3) * points;
    }


    gsMatrix<> evalBeta_R(gsMatrix<> points)
    {
        gsMatrix<> ones(1, points.cols());
        ones.setOnes();
        return solBeta.row(0) * ( ones - points ) + solBeta.row(1) * points;
    }

    gsMatrix<> evalBeta_L(gsMatrix<> points)
    {
        gsMatrix<> ones(1, points.cols());
        ones.setOnes();
        return solBeta.row(2) * ( ones - points ) + solBeta.row(3) * points;
    }


    gsMatrix<> getSol()
    {
        return sol;
    }

    gsMatrix<> getSolBeta()
    {
        return solBeta;
    }

protected:

    gsSparseSystem<> mSys;
    gsSparseSystem<> mSysBeta;
    gsMatrix<> dirichletDofs;
    gsMatrix<> dirichletDofsBeta;



    gsMatrix<> sol; // In order, it contains: alpha_0L, alpha_1L, alpha_0R, alpha_1R, beta_0, beta_1, beta_2
                    // to construct the linear combination of the GD:
                    // alpha_L = ( 1 - t ) * alpha_0L + alpha_1L * t
                    // alpha_R = ( 1 - t ) * alpha_0R + alpha_1R * t
                    //beta = ( 1 - t )^2 * beta_0 + 2 * t * ( 1 - t ) * beta_1 + t^2 * beta_2

    gsMatrix<> solBeta;

    void refresh()
    {
        gsVector<index_t> size(1);
        size << 7;

        gsDofMapper map(size);
        map.finalize();

        gsSparseSystem<> sys(map);
        mSys = sys;
    }

    void refreshBeta()
    {
        gsVector<index_t> size(1);
        size << 4;

        gsDofMapper mapBeta(size);
        mapBeta.finalize();

        gsSparseSystem<> sysBeta(mapBeta);
        mSysBeta = sysBeta;

    }


    void assemble()
    {
        mSys.reserve(49, 1); // Reserve for the matrix 7x7 values

        dirichletDofs.setZero(mSys.colMapper(0).boundarySize(),1);

        // Assemble volume integrals
        Visitor visitor;
        apply(visitor);

        mSys.matrix().makeCompressed();

    }

    void assembleBeta()
    {
        mSysBeta.reserve(16, 1); // Reserve for the matrix 4x4 values

        dirichletDofsBeta.setZero(mSysBeta.colMapper(0).boundarySize(), 1);


        // Assemble volume integrals
        Visitor visitorBeta;
        applyBeta(visitorBeta);

        mSysBeta.matrix().makeCompressed();



    }

    void apply(Visitor visitor)
    {
    #pragma omp parallel
        {
            Visitor
    #ifdef _OPENMP
            // Create thread-private visitor
            visitor_(visitor);
            const int tid = omp_get_thread_num();
            const int nt  = omp_get_num_threads();
    #else
                &visitor_ = visitor;
    #endif

            gsQuadRule<T> quRule ; // Quadrature rule
            gsMatrix<T> quNodes  ; // Temp variable for mapped nodes
            gsVector<T> quWeights; // Temp variable for mapped weights

            const gsBasis<T> & basis = this->m_mb[0].basis(0).component(1); // = 0

            // Initialize reference quadrature rule and visitor data
            visitor_.initialize(basis,quRule);

            // Initialize domain element iterator -- using unknown 0
            typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator(boundary::none);

    #ifdef _OPENMP
            for ( domIt->next(tid); domIt->good(); domIt->next(nt) )
    #else
            for (; domIt->good(); domIt->next() )
    #endif
            {
                // Map the Quadrature rule to the element
                quRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights );

                // Perform required evaluations on the quadrature nodes
                visitor_.evaluate(quNodes, this->m_mp);

                // Assemble on element
                visitor_.assemble(*domIt, quWeights);

                // Push to global matrix and right-hand side vector
    #pragma omp critical(localToGlobal)
                visitor_.localToGlobal( dirichletDofs, mSys); // omp_locks inside

            }

        }//omp parallel

    }

    void applyBeta(Visitor visitorBeta)
    {
#pragma omp parallel
        {
            Visitor
#ifdef _OPENMP
            // Create thread-private visitor
            visitor_Beta(visitorBeta);
            const int tid = omp_get_thread_num();
            const int nt  = omp_get_num_threads();
#else
                &visitor_Beta = visitorBeta;
#endif

            gsQuadRule<T> quRule ; // Quadrature rule
            gsMatrix<T> quNodes  ; // Temp variable for mapped nodes
            gsVector<T> quWeights; // Temp variable for mapped weights

            const gsBasis<T> & basis = this->m_mb[0].basis(0).component(1); // = 0

            // Initialize reference quadrature rule and visitor data
            visitor_Beta.initialize(basis,quRule);

            // Initialize domain element iterator -- using unknown 0
            typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator(boundary::none);



#ifdef _OPENMP
            for ( domIt->next(tid); domIt->good(); domIt->next(nt) )
#else

            for (; domIt->good(); domIt->next() )
#endif
            {

                // Map the Quadrature rule to the element
                quRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights );

                // Perform required evaluations on the quadrature nodes
                visitor_Beta.evaluateBeta(quNodes, this->m_mp, sol);

                // Assemble on element
                visitor_Beta.assembleBeta(*domIt, quWeights);

                // Push to global matrix and right-hand side vector
#pragma omp critical(localToGlobal)
                visitor_Beta.localToGlobalBeta( dirichletDofsBeta, mSysBeta); // omp_locks inside

            }

        }//omp parallel

    }

    void solve()
    {
        gsSparseSolver<>::CGDiagonal solver;

//        gsInfo << "Matrix: " << mSys.matrix() << "\n";

        solver.compute(mSys.matrix());
        sol = solver.solve(mSys.rhs()); // My solution

//        gsInfo << "Rhs: " << mSys.rhs() << "\n";

    }


    void solveBeta()
    {
        gsSparseSolver<>::CGDiagonal solver;

//        gsInfo << "Matrix Beta: " << mSysBeta.matrix() << "\n";

        solver.compute(mSysBeta.matrix());
        solBeta = solver.solve(mSysBeta.rhs()); // My solution

    }

    void setGDEdge()
    {
        gsMatrix<> solTMP(7, 1);
        gsMatrix<> solBetaTMP(4, 1);

        solTMP(0, 0) = 1;
        solTMP(1, 0) = 1;
        solTMP(2, 0) = 1;
        solTMP(3, 0) = 1;
        solTMP(4, 0) = 0;
        solTMP(5, 0) = 0;
        solTMP(6, 0) = 0;

        solBetaTMP(0, 0) = 0;
        solBetaTMP(1, 0) = 0;
        solBetaTMP(2, 0) = 0;
        solBetaTMP(3, 0) = 0;

        sol = solTMP;
        solBeta = solBetaTMP;
    }



    void AScondition(gsMultiPatch<T> const & mp)
    {
        index_t p_size = 10;
        gsMatrix<> points(1, p_size);

        gsVector<> vec;
        vec.setLinSpaced(p_size,0,1);
        points = vec.transpose();

        gsGeometry<> & FR = mp.patch(0);
        gsGeometry<> & FL = mp.patch(1);

        gsMatrix<> DuFR, DvFR, DvFL;

        gsMatrix<> ones(1, points.cols());
        ones.setOnes();

        gsMatrix<> pointV(FR.parDim(), points.cols());
        pointV.setZero();
        pointV.row(1) = points;

        gsMatrix<> pointU(FL.parDim(), points.cols());
        pointU.setZero();
        pointU.row(0) = points;

        gsMatrix<> cond(1, p_size);

        gsMatrix<> alpha_R = sol.row(0) * ( ones - points ) + sol.row(1) * points;
        gsMatrix<> alpha_L = sol.row(2) * ( ones - points ) + sol.row(3) * points;
        gsMatrix<> beta = sol.row(4) * ( points.cwiseProduct(points) - 2 * points + ones ) + 2 * sol.row(5) * points.cwiseProduct( ones - points ) + sol.row(6) * points.cwiseProduct(points);


        refreshBeta();
        assembleBeta();
        solveBeta();

        gsMatrix<> beta_R = solBeta.row(0) * ( ones - points ) + solBeta.row(1) * points;
        gsMatrix<> beta_L = solBeta.row(2) * ( ones - points ) + solBeta.row(3) * points;



        for(index_t i = 0; i < points.cols(); i++)
        {
            DuFR = FR.jacobian(pointV.col(i)).col(0);
            DvFR = FR.jacobian(pointV.col(i)).col(1); // Same as DuFL

            DvFL = FL.jacobian(pointU.col(i)).col(1);

            cond.col(i) = alpha_R.col(i);
//                .cwiseProduct(DvFL) + alpha_L.col(i).cwiseProduct(DuFR) + beta.col(i).cwiseProduct(DvFR) ;

//            cond.col(i) = beta.col(i) - (alpha_R.col(i).cwiseProduct(beta_L.col(i)) + alpha_L.col(i).cwiseProduct(beta_R.col(i)));

//            gsInfo << "cond dim: " << cond.dim() << "\n";

            gsInfo << "Condition col " << i << ": " << cond.col(i) << "\n";


        }
    }


    }; // class gsGluingData

} // namespace gismo


