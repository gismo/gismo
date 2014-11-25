

#pragma once

#include <gsAssembler/gsNorm.h>

namespace gismo
{

template <class T>
class gsNormL2 : public gsNorm<T>
{
    friend gsNorm<T>;

public:

    gsNormL2(const gsField<T> & _field1,
             const gsFunction<T> & _func2,
             bool _f2param = false) 
    : gsNorm<T>(_field1,_func2), f2param(_f2param)
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
            numQuadNodes[i] = basis.degree(i) + 1;
        
        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        // Set Geometry evaluation flags
        evFlags = NEED_MEASURE| NEED_VALUE;
    }
    
    // Evaluate on element.
    inline void evaluate(gsGeometryEvaluator<T> & geoEval,
                         const gsGeometry<T>    & func1,
                         const gsFunction<T>    & func2,
                         gsMatrix<T>            & quNodes)
    {
        // Evaluate first function
        func1.eval_into(quNodes, f1vals);
        
        // Compute geometry related values
        geoEval.evaluateAt(quNodes);
        
        // Evaluate second function (defined of physical domain)
        func2.eval_into(geoEval.values(), f2vals);
    }
    
    // assemble on element
    inline T compute(gsDomainIterator<T>    & element, 
                     gsGeometryEvaluator<T> & geoEval,
                     gsVector<T> const      & quWeights)
    {
        T sum(0.0);
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            const T weight = quWeights[k] * geoEval.measure(k);
            sum +=weight * ( f1vals.col(k) - f2vals.col(k) ).squaredNorm();
        }

        return sum;
    }
    
private:

    gsMatrix<T> f1vals, f2vals;

    bool f2param;
};


} // namespace gismo

