/** @file gsNormL2Boundary.h

    @brief Computes the L2 norm only on the boudary (surface integral).

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, J. Sogn
*/


#pragma once

#include <gsAssembler/gsNorm.h>

namespace gismo
{

/** @brief The gsNormL2 class provides the functionality
 * to calculate the L2 - norm on a given boudary (surface integral)
 * between a field and a function.
 *
 * \ingroup Assembler
*/
template <class T>
class gsNormL2Boundary : public gsNorm<T>
{
    friend class gsNorm<T>;

    typedef gsNorm<T> Base;
    
public:

    gsNormL2Boundary(const gsField<T> & _field1,
                     const gsFunction<T> & _func2,
                     bool _f2param = false)
    : gsNorm<T>(_field1,_func2), f2param(_f2param)
    { 
        
    }

    explicit gsNormL2Boundary(const gsField<T> & _field1)
    : gsNorm<T>(_field1), f2param(false)
    {

    }

public:
    
    T compute(bool storeElWise = false)
    {
        if ( storeElWise )
        {
            // m_elWise.reserve( ..
            m_elWise.clear();
        }

        m_value = T(0.0);

        for (typename gsMultiPatch<T>::const_biterator bit =
                 patchesPtr->bBegin(); bit != patchesPtr->bEnd(); ++bit)
        {
            side = bit->side();
            this->apply1(*this, storeElWise, bit->patch, bit->side() );
        }

        m_value = takeRoot(m_value);
        return this->m_value;
    }

    inline T takeRoot(const T v) { return math::sqrt(v);}

protected:

    void initialize(const gsBasis<T> & basis,
                    gsQuadRule<T> & rule,
                    unsigned      & evFlags) // replace with geoEval ?
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
        evFlags = NEED_MEASURE| NEED_VALUE;
    }
    
    // Evaluate on element.
    inline void evaluate(gsGeometryEvaluator<T> & geoEval,
                         const gsFunction<T>    & _func1,
                         const gsFunction<T>    & _func2,
                         gsMatrix<T>            & quNodes)
    {
        // Evaluate first function
        _func1.eval_into(quNodes, f1vals);
        
        // Compute geometry related values
        geoEval.evaluateAt(quNodes);

        // Evaluate second function (defined of physical domain)
        _func2.eval_into(geoEval.values(), f2vals);
        
        // ** Evaluate function v
        //gsMatrix<T> f2val = func2Param ? _func2.eval(quNodes)
        //: _func2.eval( geoEval->values() );
    }
    
    // assemble on element
    inline T compute(gsDomainIterator<T>    & element, 
                     gsGeometryEvaluator<T> & geoEval,
                     gsVector<T> const      & quWeights,
                     T & accumulated)
    {
        T sum(0.0);
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            geoEval.outerNormal(k, side, unormal);
            const T weight = quWeights[k] * unormal.norm();
            sum +=weight * ( f1vals.col(k) - f2vals.col(k) ).squaredNorm();
        }
        accumulated += sum;
        return sum;
    }
    
private:

    using Base::m_elWise;
    using Base::m_value;
    using Base::patchesPtr;
    
    gsMatrix<T> f1vals, f2vals;
    gsVector<T> unormal;
    bool f2param;

protected:

    boxSide side;

};


} // namespace gismo

