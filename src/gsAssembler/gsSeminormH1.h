/** @file gsSeminormH1.h

    @brief Computes the H1 norm.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#include<gsAssembler/gsNorm.h>

#pragma once

namespace gismo
{

/** @brief The gsSeminormH1 class provides the functionality
 * to calculate the H1 - seminorm between a field and a function.
 *
 * \ingroup Assembler
*/
template <class T>
class gsSeminormH1 : public gsNorm<T>
{
    friend  class gsNorm<T>;

public:

    gsSeminormH1(const gsField<T> & _field1,
             const gsFunction<T> & _func2,
             bool _f2param = false) 
    : gsNorm<T>(_field1,_func2), f2param(_f2param)
    { }

    gsSeminormH1(const gsField<T> & _field1) 
    : gsNorm<T>(_field1), f2param(false)
    { }

public:
    
    T compute(bool storeElWise = false)
    {
        this->apply(*this,storeElWise);
        return m_value;
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
            numQuadNodes[i] = basis.degree(i) + 1;

        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here
        
        // Set Geometry evaluation flags
        evFlags = NEED_MEASURE|NEED_VALUE|NEED_GRAD_TRANSFORM;
    }
    
    // Evaluate on element.
    void evaluate(gsGeometryEvaluator<T> & geoEval,
                  const gsFunction<T>    & _func1,
                  const gsFunction<T>    & _func2,
                  gsMatrix<T>            & quNodes)
    {
        // Evaluate first function
        _func1.deriv_into(quNodes, f1ders);

        // Evaluate second function
        geoEval.evaluateAt(quNodes);
        _func2.deriv_into( f2param ? quNodes : geoEval.values() , f2ders);
    }

    // assemble on element
    inline T compute(gsDomainIterator<T>    & ,
                     gsGeometryEvaluator<T> & geoEval,
                     gsVector<T> const      & quWeights,
                     T & accumulated)
    {
        T sum(0.0);
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Transform the gradients
            geoEval.transformGradients(k, f1ders, f1pders);

            // Transform the gradients, if func2 is defined on the parameter space (f2param = true)
            if(f2param)
                geoEval.transformGradients(k, f2ders, f2pders);

            // old
            //if ( f2param )
            //f2ders.col(k)=geoEval.gradTransforms().block(0, k*d,d,d) * f2ders.col(k);// to do: generalize

            const T weight = quWeights[k] *  geoEval.measure(k);

            if(!f2param) // standard case: func2 defined on physical space
            {
                // for each k: put the gradients into the columns (as in f1pders)
                gsMatrix<T> f2dersk = f2ders.col(k);
                f2dersk.resize(gsNorm<T>::func2->domainDim(), gsNorm<T>::func2->targetDim());

                sum += weight * (f1pders - f2dersk).squaredNorm();
            }
            else // case: func2 defined on parameter space
                sum += weight * (f1pders - f2pders).squaredNorm();
        }
        accumulated += sum;
        return sum;
    }
    
    inline T takeRoot(const T v) { return math::sqrt(v);}

private:
    using gsNorm<T>::m_value;
    using gsNorm<T>::m_elWise;

    gsMatrix<T> f1ders, f2ders;
    gsMatrix<T> f1pders, f2pders; // f2pders only needed if f2param = true

    bool f2param;
};


} // namespace gismo

