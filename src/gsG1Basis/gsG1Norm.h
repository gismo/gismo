/** @file gsSeminormH1.h
    @brief Computes the H1 norm, modified for G1 Basis functions.
    This file is part of the G+Smo library.
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    Author(s): A. Mantzaflaris & P. Weinmueller
*/

#include <gsG1Basis/gsVisitorG1Norm.h>


#pragma once

namespace gismo
{

/** @brief The gsSeminormH1 class provides the functionality
 * to calculate the H1 - seminorm between a field and a function.
 *
 * \ingroup Assembler
*/
    template <class T, class Visitor = gsVisitorG1Norm<T> >
    class gsG1Norm
    {

    public:

        gsG1Norm(const gsMultiPatch<T> & multiPatch,
                     const gsMultiBasis<> & multiBasis,
                     const gsMultiPatch<T> & mpSol,
                     const gsMatrix<T> & g1Solution,
                     const gsFunctionWithDerivatives<T> &exactSolution)
                : patchesPtr( &multiPatch ), basisPtr( &multiBasis ),
                  mpSolPtr(&mpSol), g1Sol(&g1Solution), exactSol(exactSolution)
        {
        }

    public:
        /// @brief Returns the computed norm value
        T valueL2() const { return m_valueL2; }
        T valueH1() const { return m_valueH1; }
        T valueH2() const { return m_valueH2; }

        void compute(bool storeElWise = false)
        {
            boxSide side = boundary::none;

            if ( storeElWise )
                m_elWise.clear();

            m_valueL2 = T(0.0);
            m_valueH1 = T(0.0);
            m_valueH2 = T(0.0);
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

                // G1 MultiBasis
                gsG1MultiBasis<T> g1MultiBasis(*patchesPtr, *basisPtr); // Maybe earlier? Only need for #patches>2

                for (size_t pn = 0; pn < patchesPtr->nPatches(); ++pn)// for all patches
                {
                    //const gsFunction<T> & func2p = exactSol->function(pn);

                    // Obtain an integration domain
                    const gsBasis<T> & dom = basisPtr->basis(pn);

                    // Initialize visitor
                    visitor.initialize(dom, QuRule, evFlags);

                    // Initialize geometry evaluator
                    //typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, patchesPtr->patch(pn)));
                    const gsGeometry<T> & patch = patchesPtr->patch(pn);

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
                        visitor.evaluate(patch, basisPtr, mpSolPtr, g1Sol, exactSol, g1MultiBasis, quNodes);

                        // Accumulate value from the current element (squared)
                        T temp = 0.0;
                        T temp2 = 0.0;
                        T temp3 = 0.0;
                        //visitor.compute(*domIt, *geoEval, quWeights, m_value);
                        visitor.compute(*domIt, quWeights, temp, temp2, temp3);
#pragma omp critical
                        {
                            m_valueL2 += temp;
                            m_valueH1 += temp2;
                            m_valueH2 += temp3;
                        }
                        //if (storeElWise)
                        //   m_elWise.push_back(takeRoot(result));
                    }
                }
            }//omp parallel

            m_valueL2 = takeRoot(m_valueL2);
            m_valueH1 = takeRoot(m_valueH1);
            m_valueH2 = takeRoot(m_valueH2);

        }



        inline T takeRoot(const T v) { return math::sqrt(v);}




    protected:

        const gsMultiPatch<T> * patchesPtr;

        const gsMultiBasis<T> * basisPtr;

        const gsMultiPatch<T> * mpSolPtr;

        const gsMatrix<T> * g1Sol;

        const gsFunctionWithDerivatives<T> & exactSol;

    private:

        bool f2param;

    protected:

        std::vector<T> m_elWise;    // vector of the element-wise values of the norm
        T              m_valueL2, m_valueH1, m_valueH2;     // the total value of the norm


    };


} // namespace gismo