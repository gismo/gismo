/** @file gsVisitorNormL2BasisFunction.h

    @brief Computes the L2 norm, modified for G1 Basis functions.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris & P. Weinmueller
*/

#include<gsG1Basis/Norm/gsVisitorNormL2BasisFunction.h>
#include <gsCore/gsGeometryEvaluator.h>

#pragma once

namespace gismo
{

/** @brief The gsSeminormH1 class provides the functionality
 * to calculate the H1 - seminorm between a field and a function.
 *
 * \ingroup Assembler
*/
template <class T, class Visitor = gsVisitorNormL2BasisFunction<T> >
class gsNormL2BasisFunction
{

public:

    gsNormL2BasisFunction(const gsMultiPatch<T> & multiPatch,
                 gsMultiBasis<T> & multiBasis,
                 const gsBSplineBasis<> basis_plus,
                 const gsBSplineBasis<> basis_minus,
                 const gsApproxGluingData<T> gD,
                 const index_t uv)
        : patchesPtr( &multiPatch ), basisPtr( &multiBasis ),
          m_basis_plus(basis_plus), m_basis_minus(basis_minus), m_gD(gD), m_uv(uv)
    {
    }

public:
    /// @brief Returns the computed norm value
    T value() const { return m_value; }

    void compute(bool storeElWise = false)
    {
        boxSide side = boundary::none;

        if ( storeElWise )
            m_elWise.clear();

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

            for (size_t pn = 0; pn < patchesPtr->nPatches(); ++pn)// for all patches
            {
                // Obtain an integration domain
                const gsBasis<T> & dom = basisPtr->basis(0);

                // Initialize visitor
                visitor.initialize(dom, QuRule, evFlags);

                // Initialize geometry evaluator
                typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, patchesPtr->patch(pn)));

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
                    visitor.evaluate(*geoEval, quNodes, patchesPtr, basisPtr, m_basis_plus, m_basis_minus, m_gD, m_uv, pn);

                    // Accumulate value from the current element (squared)
                    T temp = 0.0;
                    //visitor.compute(*domIt, *geoEval, quWeights, m_value);
                    const T result = visitor.compute(*domIt, *geoEval, quWeights, temp);
#pragma omp critical
                    {
                        m_value += result;
                    }
                    //if (storeElWise)
                    //   m_elWise.push_back(takeRoot(result));
                }
            }
        }//omp parallel
        m_value = takeRoot(m_value);

    }



    inline T takeRoot(const T v) { return math::sqrt(v);}




protected:

    const gsMultiPatch<T> * patchesPtr;

    gsMultiBasis<T> * basisPtr;

    const gsBSplineBasis<> m_basis_plus;
    const gsBSplineBasis<> m_basis_minus;
    const gsApproxGluingData<T> m_gD;
    const index_t m_uv;


protected:

    std::vector<T> m_elWise;    // vector of the element-wise values of the norm
    T              m_value;     // the total value of the norm


};


} // namespace gismo

