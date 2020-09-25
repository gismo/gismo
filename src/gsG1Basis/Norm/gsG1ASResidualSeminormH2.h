/** @file gsSeminormH2.h

    @brief Computes the H2 seminorm, modified for g1 Basis functions.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Farahat
*/


#pragma once

#include <gsG1Basis/Norm/gsNorm.h>
# include <gsG1Basis/Norm/gsG1ASVisitorResidualSeminormH2.h>


namespace gismo
{

/** @brief The gsSeminormH2 class provides the functionality
 * to calculate the H2 - seminorm between a field and a function.
 *
 * \ingroup Assembler
*/
template <class T, class Visitor = gsG1ASVisitorResidualSeminormH2<T> >
class gsG1ASResidualSeminormH2
{

public:

    gsG1ASResidualSeminormH2(const gsMultiPatch<T> & multiPatch,
                             std::vector<gsSparseMatrix<>> & _field1,
                             std::vector<gsMultiBasis<>> & _func2,
                             bool _f2param = false)
        : patchesPtr( &multiPatch ),
          sparseMatrix(&_field1), basisVec(&_func2), f2param(_f2param)
    {
    }


public:

    /// @brief Returns the computed norm value
    T value() const { return m_value; }

    void compute(std::vector<gsG1System<real_t>> & g1SysVec, bool storeElWise = false)
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
                visitor.initialize(basisVec, QuRule, evFlags);

                // Initialize geometry evaluator
                typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, patchesPtr->patch(pn)));

                const gsBasis<T> & dom = basisVec->at(1).basis(pn);

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
                    visitor.evaluate(*geoEval, sparseMatrix, basisVec, g1SysVec, quNodes);

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

    const std::vector<gsSparseMatrix<>> * sparseMatrix;

    const std::vector<gsMultiBasis<>> * basisVec;

private:
    bool f2param;// not used yet

protected:

    std::vector<T> m_elWise;    // vector of the element-wise values of the norm
    T              m_value;     // the total value of the norm


};


} // namespace gismo

