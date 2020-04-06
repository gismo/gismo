//
// Created by afarahat on 3/23/20.
//


#pragma once

#include <gsG1Basis/gsG1ASGluingDataAssembler.h>
#include <gsG1Basis/gsGluingData.h>
#include <gsG1Basis/gsG1ASGluingDataVisitorGlobal.h>
# include <gsG1Basis/gsG1OptionList.h>



namespace gismo
{

template<class T, class Visitor = gsG1ASGluingDataVisitorGlobal<T>>
class gsG1ASGluingData : public gsGluingData<T>
{
public:
    gsG1ASGluingData()
    { }

    gsG1ASGluingData(gsMultiPatch<T> const & mp,
                 gsMultiBasis<T> & mb)
        : gsGluingData<T>(mp, mb)
    {
        refresh();
        assemble();
        solve();
//        AScondition(mp);
    }





protected:

gsSparseSystem<> mSys;
gsMatrix<> dirichletDofs;

gsMapData<T> t;


gsMatrix<> sol; // In order, it contains: alpha_1L, alpha_0R, alpha_1R, beta_0, beta_1, beta_2 (alpha_0L already setted to zero in the system)
                // to construct the linear combination of the GD:
                // alpha_L = ( 1 - t ) * 1 + alpha_1L * t
                // alpha_R = ( 1 - t ) * alpha_0R + alpha_1R * t
                //beta = ( 1 - t )^2 * beta_0 + 2 * t * ( 1 - t ) * beta_1 + t^2 * beta_2



void refresh()
{
    gsVector<> size(1);
    size << 6;

    gsDofMapper map(size);
    map.finalize();

    gsSparseSystem<> sys(map);
    mSys = sys;
}


void assemble()
{
    mSys.reserve(36, 1); // Reserve for the matrix 6x6 values

    dirichletDofs.setZero(mSys.colMapper(0).boundarySize());

    // Assemble volume integrals
    Visitor visitor;
    apply(visitor);

    mSys.matrix().makeCompressed();

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

        //const gsGeometry<T> & patch = m_geo.patch(patchIndex); // 0 = patchindex

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

void solve()
{
    gsSparseSolver<>::CGDiagonal solver;


    solver.compute(mSys.matrix());
    sol = solver.solve(mSys.rhs()); // My solution

    gsInfo << "Solution: " << sol << "\n";
}


void AScondition(gsMultiPatch<T> const & mp)
{
    index_t p_size = 10000;
    gsMatrix<> points(1, p_size);
    points.setRandom();

    gsVector<> vec;
    vec.setLinSpaced(p_size,0,1);
    points = vec.transpose();

    gsGeometry<> & FR = mp.patch(0);
    gsGeometry<> & FL = mp.patch(1);

    gsMatrix<> pointV(FR.parDim(), points.cols());
    pointV.setZero();
    pointV.row(1) = points;

    gsMatrix<> pointU(FL.parDim(), points.cols());
    pointU.setZero();
    pointU.row(0) = points;

    gsMatrix<> DuFR = FR.jacobian(pointV).col(0);
    gsMatrix<> DvFR = FR.jacobian(pointV).col(1); // Same as DuFL

    gsMatrix<> DvFL = FL.jacobian(pointU).col(1);

    gsMatrix<> ones(1, points.cols());
    ones.setOnes();

    gsMatrix<> alpha_L = ( ones - points ) + sol.row(0) * points;
    gsMatrix<> alpha_R = sol.row(1) * ( ones - points ) + sol.row(2) * points;
    gsMatrix<> beta = sol.row(3) * ( points.cwiseProduct(points) - 2 * points + ones ) + 2 * sol.row(4) * points.cwiseProduct( ones - points ) + sol.row(5) * points.cwiseProduct(points);

    gsMatrix<> cond = alpha_R * DvFL + alpha_L * DuFR + beta * DvFR;

    gsInfo << "Condition 1: " << cond << "\n";
}


}; // class gsGluingData

} // namespace gismo


