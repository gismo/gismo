/** @file gsNormL2.h

    @brief Computes the L2 norm.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/


#pragma once

#include <gsG1Basis/gsNorm.h>
#include <gsG1Basis/gsVisitorNormL2.h>

# include<chrono>

namespace gismo
{


/** @brief The gsNormL2 class provides the functionality
 * to calculate the L2 - norm between a field and a function.
 *
 * \ingroup Assembler
*/
template <int p, class T = real_t, class Visitor = gsVisitorNormL2<T> >
class gsNormL
{

public:

    gsNormL(const gsMultiPatch<> & multiPatch,
            const gsSparseMatrix<T> & _field1,
            const gsFunction<T> & _func2,
            bool _f2param = false)
        :  patchesPtr( &multiPatch),
           sparseMatrix(&_field1), func2(&_func2), f2param(_f2param)
    {

    }


public:

    /// @brief Returns the computed norm value
    T value() const { return m_value; }

    void compute(gsVector<> numBasisFunctions, bool storeElWise = false)
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
                const gsFunction<T> & func2p = func2->function(pn);

                // Obtain an integration domain
                const gsBasis<T> & dom = patchesPtr->patch(pn).basis();

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
                    visitor.evaluate(*geoEval, func2p, dom, sparseMatrix, numBasisFunctions, quNodes);

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

    const gsSparseMatrix<T>    * sparseMatrix;

    const gsFunctionSet<T> * func2;

private:

    bool f2param;

protected:
    std::vector< gsMultiPatch<>> m_G1Basis;


protected:

    std::vector<T> m_elWise;    // vector of the element-wise values of the norm
    T              m_value;     // the total value of the norm


};


template <class T>
class gsNormL2 : public gsNormL<2,T,gsVisitorNormL2<T>>
{
public:
gsNormL2(const gsMultiPatch<> & multiPatch,
         const gsSparseMatrix<T> & sparseMatrix,
         const gsFunction<T> & _func2,
         bool _f2param = false)
    : gsNormL<2,T,gsVisitorNormL2<T>>(multiPatch, sparseMatrix, _func2, _f2param)
{ }

};



} // namespace gismo

