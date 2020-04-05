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
    }





protected:

gsSparseSystem<> mSys;
gsMatrix<> dirichletDofs;


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
    gsMatrix<> sol;

    solver.compute(mSys.matrix());
    sol = solver.solve(mSys.rhs()); // My solution

    gsInfo << "Solution: " << sol << "\n";
}



}; // class gsGluingData

} // namespace gismo


