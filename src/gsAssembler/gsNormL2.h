/** @file gsNormL2.h

    @brief Computes the L2 norm.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/


#pragma once

#include <gsAssembler/gsNorm.h>

namespace gismo
{


/** @brief The gsNormL2 class provides the functionality
 * to calculate the L2 - norm between a field and a function.
 *
 * \ingroup Assembler
*/
template <int p, class T = real_t>
class gsNormL : public gsNorm<T>
{
    typedef gsNorm<T> Base;
    friend class gsNorm<T>;

public:

    gsNormL(const gsField<T> & _field1,
             const gsFunction<T> & _func2,
             bool _f2param = false) 
    : Base(_field1,_func2), f2param(_f2param)
    { 
        
    }

    gsNormL(const gsField<T> & _field1)
    : Base(_field1), f2param(false)
    {

    }

public:
    
    T compute(bool storeElWise = false)
    {
        this->apply(*this,storeElWise);
        return this->m_value;
    }


protected:

    void initialize(const gsBasis<T> & basis,
                    gsQuadRule<T> & rule,
                    unsigned      & evFlags) // replace with geoEval ?
    {
        // Setup Quadrature
        const unsigned d = basis.dim();
        gsVector<index_t> numQuadNodes( d );
        for (unsigned i = 0; i < d; ++i)
            numQuadNodes[i] = basis.degree(i) + 2;
        
        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        // Set Geometry evaluation flags
        evFlags = NEED_MEASURE| NEED_VALUE;
    }
    
    // Evaluate on element.
    void evaluate(gsGeometryEvaluator<T> & geoEval,
                  const gsFunction<T>    & _func1,
                  const gsFunction<T>    & _func2,
                  gsMatrix<T>            & quNodes)
    {
        // Evaluate first function
        _func1.eval_into(quNodes, f1vals);
        
        // Compute geometry related values
        geoEval.evaluateAt(quNodes);
        
        // Evaluate second function
        _func2.eval_into( f2param ? quNodes : geoEval.values() , f2vals);        
    }
    
    // assemble on element
    T compute(gsDomainIterator<T>    & ,
              gsGeometryEvaluator<T> & geoEval,
              gsVector<T> const      & quWeights,
              T & accumulated)
    {
        T sum(0.0);
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            const T weight = quWeights[k] * geoEval.measure(k);
            switch (p)
            {
            case 0: // infinity norm
                // .template lpNorm<Eigen::Infinity>();
                sum = (f1vals - f2vals).array().abs().maxCoeff();
                accumulated = math::max(accumulated, sum);
                return sum;
                break;
            case 1:
                sum += weight * ( f1vals.col(k) - f2vals.col(k) ).template lpNorm<1>();
                break;
            case 2:
                sum += weight * ( f1vals.col(k) - f2vals.col(k) ).squaredNorm();
                break;
            default:
                sum += weight * ( f1vals.col(k) - f2vals.col(k) ).array().abs().pow(p).sum();
                //sum += weight * ( f1vals.col(k) - f2vals.col(k) ).template lpNorm<p>().squared();
            }
        }

        accumulated += sum;
        return sum;
    }

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

private:

    gsMatrix<T> f1vals, f2vals;

    bool f2param;
};


template <class T>
class gsNormL2 : public gsNormL<2,T>
{
public:
    gsNormL2(const gsField<T> & _field1,
             const gsFunction<T> & _func2,
             bool _f2param = false) 
    : gsNormL<2,T>(_field1, _func2, _f2param)
    { }

    gsNormL2(const gsField<T> & _field1)
    : gsNormL<2,T>(_field1)
    { }
};



} // namespace gismo

