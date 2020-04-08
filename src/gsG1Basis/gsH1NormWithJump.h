/** @file gsSeminormH1.h

    @brief Computes the H1 norm with the jump; for approx. g1 Basis.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
*/


#pragma once

namespace gismo
{

/** @brief The gsSeminormH1 class provides the functionality
 * to calculate the H1 - seminorm between a field and a function.
 *
 * \ingroup Assembler
*/
template <class T>
class gsH1NormWithJump
{

public:

    gsH1NormWithJump(const gsMultiPatch<> & multiPatch,
                     const gsSparseMatrix<T> & _field1)
        : m_mp(multiPatch) , m_sparse(&_field1)
    {
        m_value.setZero(m_mp.interfaces().size(),1);
    }


public:

    void compute(gsVector<> numBasisFunctions, gsVector<> numInterfaceFunctions)
    {

        for (size_t numInt = 0; numInt < m_mp.interfaces().size(); numInt++ )
        {
            T value = T(0.0);

            gsMatrix<T> quNodes_L; // Temp variable for mapped nodes
            gsVector<T> quWeights; // Temp variable for mapped weights
            gsQuadRule<T> QuRule_L; // Reference Quadrature rule

            gsMatrix<T> quNodes_R; // Temp variable for mapped nodes
            gsQuadRule<T> QuRule_R; // Reference Quadrature rule

            // Evaluation flags for the Geometry map
            unsigned evFlags(0);

            boundaryInterface & iFace = m_mp.interfaces()[numInt];

            index_t L = iFace.second().patch;
            index_t R = iFace.first().patch;

            // Obtain an integration domain
            const gsBasis<T> & dom_L = m_mp.patch(L).basis();
            const gsBasis<T> & dom_R = m_mp.patch(R).basis();

            boxSide side_L = iFace.second();
            boxSide side_R = iFace.first();

            // Initialize visitor
            initializeb(dom_L, QuRule_L, evFlags, side_L);
            initializeb(dom_R, QuRule_R, evFlags, side_R);

            // Initialize geometry evaluator
            typename gsGeometryEvaluator<T>::uPtr geoEval_L(getEvaluator(evFlags, m_mp.patch(L)));
            // Initialize geometry evaluator
            typename gsGeometryEvaluator<T>::uPtr geoEval_R(getEvaluator(evFlags, m_mp.patch(R)));

            const gsAffineFunction<> ifaceMap(m_mp.getMapForInterface(iFace));

            typename gsBasis<T>::domainIter domIt = dom_L.makeDomainIterator(side_R);
            for (; domIt->good(); domIt->next())
            {
                // Map the Quadrature rule to the element
                gsMatrix<T> domItCorner(2,2);
                domItCorner.col(0) = domIt->lowerCorner();
                domItCorner.col(1) = domIt->upperCorner();

                QuRule_R.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes_R, quWeights);
                ifaceMap.eval_into(domItCorner,domItCorner);
                if (domItCorner(1-side_L.direction(),0) > domItCorner(1-side_L.direction(),1)) // integral border switched
                {
                    gsMatrix<T> temp_domItCorner = domItCorner;
                    domItCorner.col(1) = temp_domItCorner.col(0);
                    domItCorner.col(0) = temp_domItCorner.col(1);
                }
                QuRule_L.mapTo(domItCorner.col(0), domItCorner.col(1), quNodes_L, quWeights);


                // Evaluate on quadrature points
                evaluateb(*geoEval_L, quNodes_L, dom_L, *geoEval_R, quNodes_R, dom_R, numBasisFunctions, numInterfaceFunctions, numInt, m_sparse);

                // Accumulate value from the current element (squared)
                computeb(*geoEval_L, side_L, *geoEval_R, quWeights, value);

            }


            m_value[numInt] = takeRoot(value);
        }

    }

    gsVector<T> value() const { return m_value; }

protected:

    void initializeb(const gsBasis<T> & basis,
                     gsQuadRule<T> & rule,
                     unsigned      & evFlags,
                     boxSide side) // replace with geoEval ?
    {
        // Setup Quadrature
        const unsigned d = basis.dim();
        const int dir = side.direction();
        gsVector<index_t> numQuadNodes( d );
        for (unsigned i = 0; i < d; ++i)
            numQuadNodes[i] = basis.degree(i) + 1;
        numQuadNodes[dir] = 1;

        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        // Set Geometry evaluation flags
        evFlags = NEED_VALUE|NEED_MEASURE|NEED_JACOBIAN|NEED_GRAD_TRANSFORM;
    }

    // Evaluate on element.
    void evaluateb(gsGeometryEvaluator<T> & geoEval,
                   gsMatrix<T>            & quNodes,
                   const gsBasis<T> & basis,
                   gsGeometryEvaluator<T> & geoEval_R,
                   gsMatrix<T>            & quNodes_R,
                   const gsBasis<T> & basis_R,
                   const gsVector<T> & numBasisFunctions,
                   const gsVector<T> & numInterfaceFunctions,
                   size_t numInt,
                   const gsSparseMatrix<T> * sol_sparse)
    {

        // Evaluate first function
        gsMatrix<unsigned> actives;
        gsMatrix<T> bGrads;

        basis.active_into(quNodes.col(0), actives);

        // Evaluate basis functions on element
        basis.deriv_into(quNodes,bGrads);

        f1ders.setZero(2,actives.rows());
        for (index_t i = numInterfaceFunctions[numInt]; i < numInterfaceFunctions[numInt+1]; i++)
            for (index_t j = 0; j < actives.rows(); j++)
                f1ders += sol_sparse->at(i,numBasisFunctions[geoEval.id()] + actives.at(j)) * bGrads.block(2*j,0,2,f1ders.dim().second);

        geoEval.evaluateAt(quNodes);


        // Evaluate second function
        gsMatrix<unsigned> actives2;
        gsMatrix<T> bGrads2;

        basis_R.active_into(quNodes_R.col(0), actives2);

        // Evaluate basis functions on element
        basis_R.deriv_into(quNodes_R,bGrads2);

        f2ders.setZero(2,actives2.rows());
        for (index_t i = numInterfaceFunctions[numInt]; i < numInterfaceFunctions[numInt+1]; i++)
            for (index_t j = 0; j < actives2.rows(); j++)
                f2ders += sol_sparse->at(i,numBasisFunctions[geoEval_R.id()] + actives2.at(j)) * bGrads2.block(2*j,0,2,f2ders.dim().second);

        geoEval_R.evaluateAt(quNodes_R);
    }

    // assemble on element
    inline T computeb(gsGeometryEvaluator<T> & geoEval,
                      boxSide side,
                      gsGeometryEvaluator<T> & geoEval_R,
                      gsVector<T> const      & quWeights,
                      T & accumulated)
    {

        T sum(0.0);
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            //const T d = element.dim();
            // Transform the gradients
            geoEval.transformGradients(k, f1ders, f1pders);
            geoEval_R.transformGradients(k, f2ders, f2pders);

            // Compute the unit normal
            gsVector<T> unormal;

            geoEval.outerNormal(k, side, unormal);

            const T weight = quWeights[k] * unormal.norm() ;

            // f2ders : N X 1
            sum += weight * ( f1pders - f2pders ).squaredNorm() ;
        }
        accumulated += sum;

        return sum;
    }

    inline T takeRoot(const T v) { return math::sqrt(v);}



private:
    gsMatrix<T> f1ders, f2ders;
    gsMatrix<T> f1pders, f2pders; // f2pders only needed if f2param = true

    gsMultiPatch<T> m_mp;
    const gsSparseMatrix<T> * m_sparse;

protected:
    gsVector<T> m_value;     // the total value of the norm


};


} // namespace gismo