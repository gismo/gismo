/** @file gsNormL2.h

    @brief Computes the L2 norm.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/


#pragma once

#include <gsG1Basis/gsVisitorGlobalGDNorm.h>


namespace gismo
{


/** @brief The gsNormL2 class provides the functionality
 * to calculate the L2 - norm between a field and a function.
 *
 * \ingroup Assembler
*/
template <class T = real_t, class Visitor = gsVisitorGlobalGDNorm<T> >
class gsGlobalGDNorm
{

public:

    gsGlobalGDNorm(const gsMultiPatch<> & multiPatch,
            gsBSplineBasis<> & basis,
            const gsMatrix<> & sol)
        :  m_mp(multiPatch), m_basis(basis), m_sol(sol)
    {
        p = 2;
    }


public:

    /// @brief Returns the computed norm value
    T value() const { return m_value; }

    void compute()
    {
        m_value = T(0.0);

#pragma omp parallel
        {
#ifdef _OPENMP
            // Create thread-private visitor
        Visitor visitor;
        const int tid = omp_get_thread_num();
        const int nt  = omp_get_num_threads();
#else
            Visitor visitor;
#endif
            gsMatrix<T> quNodes; // Temp variable for mapped nodes
            gsVector<T> quWeights; // Temp variable for mapped weights
            gsQuadRule<T> QuRule; // Reference Quadrature rule

            // Evaluation flags for the Geometry map
            unsigned evFlags(0);


            // Obtain an integration domain
            const gsBasis<T> & dom = m_mp.basis(0).component(1); // v

            // Initialize visitor
            visitor.initialize(dom, QuRule, evFlags);

            typename gsBasis<T>::domainIter domIt = dom.makeDomainIterator(boundary::none);

            // TODO: optimization of the assembling routine, it's too slow for now
            // Start iteration over elements
#ifdef _OPENMP
            for ( domIt->next(tid); domIt->good(); domIt->next(nt) )
#else
            for (; domIt->good(); domIt->next() )
#endif
            {

                // Map the Quadrature rule to the element
                QuRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);

                // Evaluate on quadrature points
                visitor.evaluate(m_mp, m_basis, m_sol, quNodes);

#pragma omp critical
                {
                    visitor.compute(*domIt, quWeights, m_value);
                    //const T result = visitor.compute(*domIt, *geoEval, quWeights, m_value);
                };
            }


        }//omp parallel
        m_value = takeRoot(m_value);

    }

protected:

    inline T takeRoot(const T v)
    {
        switch (p)
        {
            case 0: // infinity norm
            case 1:
                return v;
            case 2:
                return math::sqrt(v);
            default:
                return math::pow(v, static_cast<T>(1)/p );
        }
    }



protected:

    const gsMultiPatch<T> m_mp;

    gsBSplineBasis<> m_basis;

    const gsMatrix<T> m_sol;

protected:

    index_t p;
    T              m_value;     // the total value of the norm




};



} // namespace gismo

