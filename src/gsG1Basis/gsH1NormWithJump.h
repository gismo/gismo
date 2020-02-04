/** @file gsSeminormH1.h

    @brief Computes the H1 norm.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include<gsG1Basis/gsNorm.h>

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
    typedef std::multimap<index_t, std::map<index_t, std::map<index_t, gsMultiPatch<real_t>>>> typedef_g1;

public:

    gsH1NormWithJump(const gsField<T> & _field1,
                     const  std::vector< gsMultiPatch<>>& _field2,
                     boundaryInterface f,
                     bool _f2param = false)
        : m_field1(_field1) , m_G1Basis(_field2), iFace(f)
    {
        g1basis = true;
        g1basis_mp = false;
    }

    gsH1NormWithJump(const gsField<T> & _field1,
                     const typedef_g1 & _field2,
                     const boundaryInterface f,
                     bool _f2param = false)
        : m_field1(_field1) , m_G1Basis_mp(_field2), iFace(f)
    {
        g1basis = false;
        g1basis_mp = true;
    }


public:

    void compute(bool storeElWise = false)
    {
        if ( storeElWise )
        {
            // m_elWise.reserve( ..
            m_elWise.clear();
        }

        m_value = T(0.0);

        gsMatrix<T> quNodes_L  ; // Temp variable for mapped nodes
        gsVector<T> quWeights; // Temp variable for mapped weights
        gsQuadRule<T> QuRule_L; // Reference Quadrature rule

        gsMatrix<T> quNodes_R  ; // Temp variable for mapped nodes
        gsQuadRule<T> QuRule_R; // Reference Quadrature rule

        // Evaluation flags for the Geometry map
        unsigned evFlags(0);

        index_t L = iFace.second().patch, R = iFace.first().patch;

        const gsFunction<T> & func1  = m_field1.function(L);
        const gsFunction<T> & func2  = m_field1.function(R);

        // Obtain an integration domain
        const gsBasis<T> & dom = m_field1.isParametrized() ?
                                 m_field1.igaFunction(L).basis() : m_field1.patch(L).basis();

        boxSide side_L = iFace.second();
        boxSide side_R = iFace.first();

        // Initialize visitor
        initializeb(dom, QuRule_L, evFlags, side_L);
        initializeb(dom, QuRule_R, evFlags, side_R);

        // Initialize geometry evaluator
        typename gsGeometryEvaluator<T>::uPtr geoEval_L(getEvaluator(evFlags, m_field1.patches().patch(L)));
        // Initialize geometry evaluator
        typename gsGeometryEvaluator<T>::uPtr geoEval_R(getEvaluator(evFlags, m_field1.patches().patch(R)));

        typename gsBasis<T>::domainIter domIt_L = dom.makeDomainIterator(side_L);
        typename gsBasis<T>::domainIter domIt_R = dom.makeDomainIterator(side_R);
        domIt_R->good();
        for (; domIt_L->good(); domIt_L->next())
        {

            // Map the Quadrature rule to the element
            QuRule_L.mapTo( domIt_L->lowerCorner(), domIt_L->upperCorner(), quNodes_L, quWeights );
            QuRule_R.mapTo( domIt_R->lowerCorner(), domIt_R->upperCorner(), quNodes_R, quWeights );

            // Evaluate on quadrature points
            evaluateb(*geoEval_L, func1, quNodes_L, *geoEval_R, func2, quNodes_R);

            // Accumulate value from the current element (squared)
            const T result = computeb(*domIt_L, *geoEval_L, side_L ,*domIt_R, *geoEval_R, quWeights, m_value);
            if ( storeElWise )
                m_elWise.push_back( takeRoot(result) );
            domIt_R->next();
        }



        m_value = takeRoot(m_value);

    }

    T value() const { return m_value; }

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
                   const gsFunction<T>    & _func1,
                   gsMatrix<T>            & quNodes,
                   gsGeometryEvaluator<T> & geoEval_R,
                   const gsFunction<T>    & _func1_R,
                   gsMatrix<T>            & quNodes_R)
    {
        // Evaluate first function
        _func1.deriv_into(quNodes, f1ders);

        if (g1basis)
        {
            index_t n = m_G1Basis.at(geoEval.id()).nPatches();

            for (index_t i = 0; i < n; i++)
            {
                f1ders += m_G1Basis.at(geoEval.id()).patch(i).deriv(quNodes);
            }
        }

        if (g1basis_mp)
        {
            std::multimap<index_t, std::map<index_t, std::map<index_t, gsMultiPatch<real_t>>>>::iterator i_face;
            for (i_face=m_G1Basis_mp.equal_range(geoEval.id()).first; i_face!=m_G1Basis_mp.equal_range(geoEval.id()).second; ++i_face)
            {
                std::map<index_t, std::map<index_t, gsMultiPatch<real_t>>>::iterator i_side;
                for (i_side = i_face->second.begin(); i_side != i_face->second.end(); ++i_side)
                {
                    if (i_side->first == iFace.second().m_index) //left
                    {
                        for (std::map<index_t, gsMultiPatch<real_t>>::iterator i_mp = i_side->second.begin();
                             i_mp != i_side->second.end(); ++i_mp)
                        {
                            gsMultiPatch<real_t> mp_side = i_mp->second;
                            for (unsigned j = 0; j < mp_side.nPatches(); j++)
                            {
                                f1ders += mp_side.patch(j).deriv(quNodes);
                            }
                        }
                    }
                }
            }
        }

        // Evaluate second function
        geoEval.evaluateAt(quNodes);

        // Evaluate second function
        _func1_R.deriv_into(quNodes_R, f2ders);

        if (g1basis)
        {
            index_t n = m_G1Basis.at(geoEval_R.id()).nPatches();

            for (index_t i = 0; i < n; i++)
            {
                f2ders += m_G1Basis.at(geoEval_R.id()).patch(i).deriv(quNodes_R);
            }
        }
        if (g1basis_mp)
        {
            std::multimap<index_t, std::map<index_t, std::map<index_t, gsMultiPatch<real_t>>>>::iterator i_face;
            for (i_face=m_G1Basis_mp.equal_range(geoEval_R.id()).first; i_face!=m_G1Basis_mp.equal_range(geoEval_R.id()).second; ++i_face)
            {
                std::map<index_t, std::map<index_t, gsMultiPatch<real_t>>>::iterator i_side;
                for (i_side = i_face->second.begin(); i_side != i_face->second.end(); ++i_side)
                {

                    if (i_side->first == iFace.first().m_index) //right
                    {
                        for (std::map<index_t, gsMultiPatch<real_t>>::iterator i_mp = i_side->second.begin();
                             i_mp != i_side->second.end(); ++i_mp)
                        {
                            gsMultiPatch<real_t> mp_side = i_mp->second;
                            for (unsigned j = 0; j < mp_side.nPatches(); j++)
                            {
                                f2ders += mp_side.patch(j).deriv(quNodes_R);
                            }
                        }
                    }
                }
            }
        }

        // Evaluate second function
        geoEval_R.evaluateAt(quNodes_R);

    }

    // assemble on element
    inline T computeb(gsDomainIterator<T>    & ,
                      gsGeometryEvaluator<T> & geoEval,
                      boxSide side,
                      gsDomainIterator<T>    & ,
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

    bool f2param;
    bool g1basis;
    bool g1basis_mp;

private:
    gsVector<T> unormal;

protected:
    const gsField<T>    m_field1;
    std::vector<T> m_elWise;    // vector of the element-wise values of the norm
    T              m_value;     // the total value of the norm

    std::vector< gsMultiPatch<>> m_G1Basis;
    typedef_g1 m_G1Basis_mp;

    boundaryInterface iFace;
};


} // namespace gismo