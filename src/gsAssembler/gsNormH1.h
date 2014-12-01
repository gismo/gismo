

#pragma once

namespace gismo
{


template <class T>
class gsNormH1 : public gsNorm<T>
{
    friend  class gsNorm<T>;

public:

    gsNormH1(const gsField<T> & _field1,
             const gsFunction<T> & _func2,
             bool _f2param = false) 
    : gsNorm<T>(_field1,_func2), f2param(_f2param)
    { 
        
    }

    gsNormH1(const gsField<T> & _field1) 
    : gsNorm<T>(_field1,gsConstantFunction<T>(T(0.0), 1)), f2param(false)
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
                  const gsGeometry<T>    & func1,
                  const gsFunction<T>    & func2,
                  gsMatrix<T>            & quNodes)
    {
        // Evaluate first function
        func1.deriv_into(quNodes, f1ders);
        // get the gradients to columns
        f1ders.transposeInPlace();
        f1ders.resize(quNodes.rows(), quNodes.cols() );

        // Evaluate second function (defined of physical domain)
        geoEval.evaluateAt(quNodes);
        func2.deriv_into(geoEval.values(), f2ders);
        // get the gradients to columns
        f2ders.transposeInPlace();
        f2ders.resize(quNodes.rows(), quNodes.cols() );
    }

    // assemble on element
    inline T compute(gsDomainIterator<T>    & element, 
                     gsGeometryEvaluator<T> & geoEval,
                     gsVector<T> const      & quWeights)
    {
        T sum(0.0);
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Transform the gradients
            geoEval.transformGradients(k, f1ders, f1pders);
            //if ( f2Param )
            //f2ders.col(k)=geoEval.gradTransforms().block(0, k*d,d,d) * f2ders.col(k);// to do: generalize
            
            const T weight = quWeights[k] *  geoEval.measure(k);
            sum += weight * (f1pders - f2ders.col(k)).squaredNorm();
        }
        return sum;
    }
    
private:

    using gsNorm<T>::m_value;
    using gsNorm<T>::m_elWise;

    gsMatrix<T> f1ders, f2ders;
    gsMatrix<T> f1pders;

    bool f2param;// not used yet
};


} // namespace gismo

