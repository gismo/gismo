/** @file gsApproxNormL2.h

    @brief Computes the L2 norm.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/


#pragma once

#include <gsG1Basis/gsNorm.h>
#include <gsG1Basis/gsVisitorApproxNormL2.h>

namespace gismo
{


/** @brief The gsNormL2 class provides the functionality
 * to calculate the L2 - norm between a field and a function.
 *
 * \ingroup Assembler
*/
template <class T = real_t, class Visitor = gsVisitorApproxNormL2<T> >
class gsApproxNormL2
{

public:

    gsApproxNormL2(const gsMultiPatch<> & multiPatch,
            const gsMultiPatch<T> & approx)
        :  patchesPtr( &multiPatch), m_approx(&approx)
    {
        p = 2;
    }


public:

    /// @brief Returns the computed norm value
    T value() const { return m_value; }

    void compute(gsOptionList optionList)
    {
        boxSide side = boundary::none;

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
            const gsBasis<T> & dom = patchesPtr->patch(0).basis();

            // Initialize visitor
            visitor.initialize(dom, QuRule, evFlags);

            // Initialize geometry evaluator
            typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, patchesPtr->patch(0)));

            typename gsBasis<T>::domainIter domIt = dom.makeDomainIterator(side);

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
                visitor.evaluate(patchesPtr, *geoEval, dom, m_approx, quNodes, optionList);
#pragma omp critical
                {
                    //visitor.compute(*domIt, *geoEval, quWeights, m_value);
                    visitor.compute(*domIt, *geoEval, quWeights, m_value);
                }
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

    const gsMultiPatch<T> * patchesPtr;

    const gsMultiPatch<T> * m_approx;

protected:

    T              m_value;     // the total value of the norm
    index_t p;



};


} // namespace gismo

