/** @file gsPoissonAssembler.hpp

    @brief Provides assembler implementation for the Poisson equation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss, A. Mantzaflaris, J. Sogn
*/


#include <gsAssembler/gsVisitorPoisson.h> // Stiffness volume integrals
#include <gsAssembler/gsVisitorNeumann.h> // Neumann boundary integrals
#include <gsAssembler/gsVisitorNitsche.h> // Nitsche boundary integrals
#include <gsAssembler/gsVisitorDg.h>      // DG interface integrals

namespace gismo
{

template <class T>
gsOptionList gsPoissonAssembler<T>::defaultOptions()
{
    gsOptionList options = gsAssembler<T>::defaultOptions();
    options.update( gsVisitorDg<T>::defaultOptions(), gsOptionList::addIfUnknown );
    options.update( gsVisitorNitsche<T>::defaultOptions(), gsOptionList::addIfUnknown );
    return options;
}

template<class T>
void gsPoissonAssembler<T>::refresh()
{
    // We use predefined helper which initializes the system matrix
    // rows and columns using the same test and trial space
    Base::scalarProblemGalerkinRefresh();
}

template<class T>
void gsPoissonAssembler<T>::assemble()
{
    GISMO_ASSERT(m_system.initialized(), 
                 "Sparse system is not initialized, call initialize() or refresh()");

    // Reserve sparse system
    m_system.reserve(m_bases[0], m_options, this->pde().numRhs());

    //Get sparsity pattern
    initMatrix();

    // Compute the Dirichlet Degrees of freedom (if needed by m_options)
    Base::computeDirichletDofs();

    // Clean the sparse system
   // m_system.setZero(); //<< this call leads to a quite significant performance degrade!

    // Assemble volume integrals
    Base::template push<gsVisitorPoisson<T> >();

    // Enforce Neumann boundary conditions
    Base::template push<gsVisitorNeumann<T> >(m_pde_ptr->bc().neumannSides() );

    switch (m_options.getInt("DirichletStrategy"))
    {
        case dirichlet::penalize:
            Base::penalizeDirichletDofs();
            break;
        case dirichlet::nitsche:
            Base::template push<gsVisitorNitsche<T> >(m_pde_ptr->bc().dirichletSides());
            break;
        default:
            break;
    }

    if ( m_options.getInt("InterfaceStrategy") == iFace::dg )
        Base::template pushInterface<gsVisitorDg<T> >();
    
    // Assembly is done, compress the matrix
    Base::finalize();
}

template <class T>
void gsPoissonAssembler<T>::initMatrix()
{
    const gsDofMapper & rowMap = m_system.rowMapper(0);
    GISMO_ASSERT( &rowMap == &m_system.colMapper(0), "Error");
    std::vector< std::map<index_t,bool> > colv(rowMap.freeSize());//m_system.matrix().cols()

#pragma omp parallel
{
    gsMatrix<index_t> actives;
    gsMatrix<T> cpt;
    //gsSparseEntries<T> entries;

#ifdef _OPENMP
    const int tid = omp_get_thread_num();
    const int nt  = omp_get_num_threads();
#   pragma omp for
#endif
    for (size_t np = 0; np < m_pde_ptr->domain().nPatches(); ++np)
    {
        // Initialize domain element iterator -- using unknown 0
        const gsBasis<T> & basis = m_bases[0][np];

        typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator();

        // Start iteration over elements
#ifdef _OPENMP
        for ( domIt->next(tid); domIt->good(); domIt->next(nt) )
//#       pragma omp critical(localMat)
#else
        for (; domIt->good(); domIt->next() )
#endif
        {
            cpt = domIt->center;
            actives = basis.active(cpt);
            m_system.mapColIndices(actives, np, actives);

            const index_t numActive = actives.rows();


            GISMO_ASSERT( m_system.matrix().cols() == m_system.rhs().rows(), "gsSparseSystem is not allocated");

            for (index_t i = 0; i < numActive; ++i)
                //for (index_t i = 0; i != sz; ++i)
            {
                const int ii =  actives(i);
                if ( rowMap.is_free_index(actives.at(i)) )
                {
#                   pragma omp parallel for
                    for (index_t j = 0; j < numActive; ++j)
                    {
                        const int jj = actives(j);
                        if ( rowMap.is_free_index(actives.at(j)) )
                        {
                            // Matrix is symmetric, we store only lower
                           // triangular part
                            auto it = colv[ii].find(jj);
                            if (colv[ii].end()==it)
                            {
#                               pragma omp critical (colv_indices)
//                            m_system.matrix().coeffRef(ii, jj) += 0.0;
                            //entries.add(ii,jj, 0.0 );
                                colv[ii][jj] = true;
                            }
                        }
                    }
                }
            }
        }
    }//for patches

#pragma omp barrier

#   pragma omp parallel for
    for (index_t c = 0; c!=(index_t)colv.size();++c)
        for ( const auto & t : colv[c] )
            m_system.matrix().coeffRef(c, t.first) = (T)(0);
    
}//omp parallel

    Base::finalize();
}


}// namespace gismo
