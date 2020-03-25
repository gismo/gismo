/** @file gsSeminormH2.h

    @brief Computes the H2 seminorm, modified for g1 Basis functions.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger, P. Weinmueller
*/


#pragma once

#include <gsG1Basis/gsNorm.h>
# include <gsG1Basis/gsVisitorSeminormH2.h>


namespace gismo
{

/** @brief The gsSeminormH2 class provides the functionality
 * to calculate the H2 - seminorm between a field and a function.
 *
 * \ingroup Assembler
*/
template <class T, class Visitor = gsVisitorSeminormH2<T> >
class gsSeminormH2
{

public:

    gsSeminormH2(const gsField<T> & _field1,
                 const gsFunction<T> & _func2,
                 const std::vector< gsMultiPatch<>> & _field2,
                 bool _f2param = false)
        :patchesPtr( &_field1.patches() ),
         field1(&_field1), func2(&_func2), m_G1Basis(_field2)
    {
        g1basis = true;
    }

    gsSeminormH2(const gsField<T> & _field1,
                 const gsFunction<T> & _func2,
                 bool _f2param = false)
        : patchesPtr( &_field1.patches() ),
          field1(&_field1), func2(&_func2), f2param(_f2param)
    {
        g1basis = false;
    }

    gsSeminormH2(const gsField<T> & _field1)
        : patchesPtr( &_field1.patches() ),
          field1(&_field1), f2param(false)
    { }

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
            Visitor visitor(m_G1Basis);
            const int tid = omp_get_thread_num();
            const int nt  = omp_get_num_threads();
#else
            Visitor visitor(m_G1Basis);
#endif

            gsMatrix<T> quNodes; // Temp variable for mapped nodes
            gsVector<T> quWeights; // Temp variable for mapped weights
            gsQuadRule<T> QuRule; // Reference Quadrature rule

            // Evaluation flags for the Geometry map
            unsigned evFlags(0);

            for (size_t pn = 0; pn < patchesPtr->nPatches(); ++pn)// for all patches
            {

                const gsFunction<T> & func1 = field1->function(pn);
                const gsFunction<T> & func2p = func2->function(pn);

                // Obtain an integration domain
                const gsBasis<T> & dom = field1->isParametrized() ?
                                         field1->igaFunction(pn).basis() : field1->patch(pn).basis();

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
                    visitor.evaluate(*geoEval, func1, func2p, quNodes);

                    // Accumulate value from the current element (squared)
#pragma omp critical(compute)
                    visitor.compute(*domIt, *geoEval, quWeights, m_value);
                    //const T result = visitor.compute(*domIt, *geoEval, quWeights, m_value);

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

    const gsField<T>    * field1;

    const gsFunctionSet<T> * func2;

private:

    gsMatrix<T> f1ders, f1ders2, f2ders2;
    gsMatrix<T> f1pders2;

    bool f2param;// not used yet
    bool g1basis;

protected:
    std::vector< gsMultiPatch<>> m_G1Basis;

protected:

    std::vector<T> m_elWise;    // vector of the element-wise values of the norm
    T              m_value;     // the total value of the norm


};


} // namespace gismo

