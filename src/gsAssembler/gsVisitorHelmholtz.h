
#pragma once


#include <gsAssembler/gsVisitorPoisson.h>// Stiffness volume integrals

namespace gismo
{

template <class T>
class gsVisitorHelmholtz : gsVisitorPoisson<T> // we inherit in order to reuse code
{
public:

    /// Constructor with the right hand side function of the poisson equation
    gsVisitorHelmholtz(const gsFunction<T> & rhs) : gsVisitorPoisson<T> (rhs)
    { }

    inline void assemble(gsDomainIterator<T>    & element,// UPDATE 
                         gsGeometryEvaluator<T> & geoEval,
                         const index_t k, T weight       )
    {
        const typename gsMatrix<T>::Block bVals  = basisData.topRows(numActive);
        const typename gsMatrix<T>::Block bGrads =
            basisData.bottomRows( numActive * element.dim() );

        for (index_t k = 0; k < quNodes.cols(); ++k) // loop over quadrature nodes
        {
            // Multiply weight by the geometry measure
            const T weight = quWeights[k] * geoEval.measure(k);
            
            // Compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, bGrads, physGrad);
            
            localRhs.noalias() += weight * ( bVals.col(k) * rhsVals.col(k).transpose() ) ;

            localMat.noalias() += weight * (physGrad.transpose() * physGrad
                               + bVals.col(k) * bVals.col(k).transpose() ); //Plus a mass term
        }

protected:

//using..

};


} // namespace gismo

