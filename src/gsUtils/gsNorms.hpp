
#pragma once

#include <gsUtils/gsNorms.h>

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsFunction.h>
#include <gsCore/gsConstantFunction.h>
#include <gsCore/gsField.h>
#include <gsCore/gsGeometry.h>
#include <gsUtils/gsPointGrid.h>

#include <gsCore/gsDomainIterator.h>
#include <gsTensor/gsTensorDomainIterator.h>
#include <gsTensor/gsTensorDomainBoundaryIterator.h>

namespace gismo
{


// L2 norm ///////////////////////////////////////////////////////////////////////////


template <typename T>
T computeL2Norm(const gsGeometry<T>& geo, const gsFunction<T>& u, bool isParametrized_u, int numEvals)
{
    // note: by specifying the constant function as parametric, we can avoid the geometry evaluations
    // if u is parametric too
    return computeL2Distance(geo, u, isParametrized_u, gsConstantFunction<T>(T(0), 1), true, numEvals);
}

template <typename T>
T computeL2Norm(const gsField<T>& u, int numSamples)
{
    return computeL2Distance(u, gsConstantFunction<T>(T(0), 1), true, numSamples);
}

template <typename T>
T computeL2Distance(const gsGeometry<T>& geo, const gsFunction<T>& u, bool isParametrized_u, const gsFunction<T>& v, bool isParametrized_v, int numEvals)
{
    GISMO_ASSERT( u.targetDim() == v.targetDim(), "Functions need to have same target dimension");

    const int d = geo.parDim();
    assert( d == geo.geoDim() );

    // compute the tensor Gauss rule
    gsMatrix<T> nodes;
    gsVector<T> weights;
    gsMatrix<T> range = geo.basis().support();

    // Number of nodes of the underlying Gauss rule to use
    const index_t nodesPerInterval = 1;
    const int nodesPerElement  = math::ipow(nodesPerInterval, d);
    const int numElements      = (numEvals + nodesPerElement - 1) / nodesPerElement;
    std::vector< std::vector<T> > intervals;
    uniformIntervals<T>(range.col(0), range.col(1), intervals, numElements);

    // perform the quadrature
    gsGaussRule<T> QuRule( gsVector<index_t>::Constant(d,nodesPerInterval) );
    const int numPts = QuRule.numNodes();

    gsTensorDomainIterator<T> domIt(intervals);
    gsMatrix<T> geo_pts, geo_jac, u_val, v_val;
    T sum = 0.0;

    for (; domIt.good(); domIt.next() )
    {
        // Map the Quadrature rule to the element
        QuRule.mapTo( domIt.lowerCorner(), domIt.upperCorner(), nodes, weights );

        // only compute the geometry points if either function is not parametrized
        geo_pts =  (!isParametrized_u || !isParametrized_v) ?
                    geo.eval(nodes) : gsMatrix<T>();
        geo_jac = geo.jacobian(nodes);
        
        // evaluate u and v
        u_val = isParametrized_u ? u.eval(nodes) : u.eval(geo_pts);
        v_val = isParametrized_v ? v.eval(nodes) : v.eval(geo_pts);
        
        for (index_t k = 0; k < numPts; ++k)
        {
            const T funcDet = fabs( geo_jac.block(0, k*d, d,d).determinant() ) ;
            const gsVector<T> diff = u_val.col(k) - v_val.col(k);
            sum += weights[k] * funcDet * diff.dot(diff);
        }
    }

    return math::sqrt(sum);
}


template <typename T>
T computeL2Distance(const gsField<T>& u, const gsFunction<T>& v, bool isParametrized_v, int numEvals)
{
    T dist = T();

    for (int i = 0; i < u.nPatches(); ++i)
    {
        T curDist = computeL2Distance( u.patch(i), u.function(i), u.isParametrized(), v, isParametrized_v, numEvals);
        dist += curDist * curDist;
    }

    return math::sqrt(dist);
}


template <typename T>
T computeL2Distance(const gsField<T>& u, const gsField<T>& v, int numEvals)
{
    T dist = T();
    assert( u.nPatches() == v.nPatches() );

    for (int i = 0; i < u.nPatches(); ++i)
    {
        assert( &u.patch(i) == &v.patch(i) );
        T curDist = computeL2Distance( u.patch(i), u.function(i), u.isParametrized(), v.function(i), v.isParametrized(), numEvals);
        dist += curDist * curDist;
    }

    return std::sqrt(dist);
}



// iga L2 norm ///////////////////////////////////////////////////////////////////////////

template <typename T>
T igaL2DistanceOnElt( const std::auto_ptr< gsGeometryEvaluator<T> > & geoEval ,
                      const std::auto_ptr< gsGeometryEvaluator<T> > & funcEval,
                      const gsFunction<T>& v,
                      const bool & v_isParam,
                      const typename gsBasis<T>::domainIter & domIt)
{

    // compute image of Gauss nodes under geometry mapping as well as Jacobians
    geoEval->evaluateAt(domIt->quNodes);
    funcEval->evaluateAt(domIt->quNodes);
    const gsMatrix<T> & func_vals = funcEval->values();

    // Evaluate function v
    gsMatrix<T> v_val = v_isParam ? v.eval(domIt->quNodes)
                                  : v.eval( geoEval->values() );

    // perform the quadrature
    T sum(0.0);

    for (index_t k = 0; k < domIt->numQuNodes(); ++k) // loop over quadrature nodes
    {
        const T weight = domIt->quWeights[k] * fabs( geoEval->measure(k) );
        const gsVector<T> diff = func_vals.col(k) - v_val.col(k);
        sum += weight * diff.dot(diff);
    }
    return sum;

}

template <typename T>
T igaL2DistanceOnElt( const std::auto_ptr< gsGeometryEvaluator<T> > & geoEval ,
                      const gsFunction<T> & func,
                      const gsFunction<T>& v,
                      const bool & v_isParam,
                      const typename gsBasis<T>::domainIter & domIt)
{

    // compute image of Gauss nodes under geometry mapping as well as Jacobians
    geoEval->evaluateAt(domIt->quNodes);
    gsMatrix<T> func_vals;
    func.eval_into(domIt->quNodes, func_vals);

    // Evaluate function v
    gsMatrix<T> v_val = v_isParam ? v.eval(domIt->quNodes)
                                  : v.eval( geoEval->values() );

    // perform the quadrature
    T sum(0.0);

    for (index_t k = 0; k < domIt->numQuNodes(); ++k) // loop over quadrature nodes
    {
        const T weight = domIt->quWeights[k] * fabs( geoEval->measure(k) );
        const gsVector<T> diff = func_vals.col(k) - v_val.col(k);
        sum += weight * diff.dot(diff);
    }
    return sum;

}


template <typename T>
T igaL2Distance(const gsGeometry<T>& patch, 
                const gsGeometry<T>& func,
                const gsFunction<T>& v,
                bool v_isParam)
{
    std::auto_ptr< gsGeometryEvaluator<T> > geoEval ( patch.evaluator(NEED_VALUE   |
                                                                      NEED_MEASURE ));
    // assuming real-valued function
    std::auto_ptr< gsGeometryEvaluator<T> > funcEval ( func.evaluator(NEED_VALUE) );
    
    // degree of the underlying Gauss rule to use
    gsVector<int> numNodes( func.parDim() );
    for (index_t i = 0; i < numNodes.size(); ++i)
        numNodes[i] = func.basis().degree(i) + 1;
    //gsVector<int> numNodes = gsGaussAssembler<T>::getNumIntNodesFor( func.basis() );

    T sum(0);
    typename gsBasis<T>::domainIter domIt = func.basis().makeDomainIterator();
    domIt->computeQuadratureRule( numNodes );
    for (; domIt->good(); domIt->next())
    {
        sum += igaL2DistanceOnElt( geoEval, funcEval, v, v_isParam, domIt);
    }
    return math::sqrt(sum);
}

template <typename T>
T igaL2Distance(const gsGeometry<T>& patch,
                const gsFunction<T>& func,
                const gsFunction<T>& v,
                const gsBasis<T>& B,
                bool v_isParam)
{
    std::auto_ptr< gsGeometryEvaluator<T> > geoEval ( patch.evaluator(NEED_VALUE   |
                                                                      NEED_MEASURE ));
    // degree of the underlying Gauss rule to use
    gsVector<int> numNodes( B.dim() );
    for (index_t i = 0; i < numNodes.size(); ++i)
        numNodes[i] = B.degree(i) + 1;
    //gsVector<int> numNodes = gsGaussAssembler<T>::getNumIntNodesFor( B );

    T sum(0);
    typename gsBasis<T>::domainIter domIt = B.makeDomainIterator();
    domIt->computeQuadratureRule( numNodes );
    for (; domIt->good(); domIt->next())
    {
        sum += igaL2DistanceOnElt( geoEval, func, v, v_isParam, domIt);
    }
    return math::sqrt(sum);
}

template <typename T>
gsMatrix<T> igaL2DistanceEltWiseSq(const gsGeometry<T>& patch,
                                   const gsGeometry<T>& func,
                                   const gsFunction<T>& v,
                                   bool v_isParam)
{
    const int d = func.parDim();

    std::auto_ptr< gsGeometryEvaluator<T> > geoEval ( patch.evaluator(NEED_VALUE   |
                                                                      NEED_MEASURE ));
    // assuming real-valued function
    std::auto_ptr< gsGeometryEvaluator<T> > funcEval ( func.evaluator(NEED_VALUE) );

    // degree of the underlying Gauss rule to use
    gsVector<int> numNodes( func.parDim() );
    for (index_t i = 0; i < numNodes.size(); ++i)
        numNodes[i] = func.basis().degree(i) + 1;
    //gsVector<int> numNodes = gsGaussAssembler<T>::getNumIntNodesFor( func.basis() );

    typename gsBasis<T>::domainIter domIt = func.basis().makeDomainIterator();
    domIt->computeQuadratureRule( numNodes );

    // Matrix to store the element-wise errors
    // in each line:
    // first entry: error
    // second to last entry: lower and upper
    gsMatrix<T> Errs( func.basis().numElements(), 1 + 2*d );

    gsVector<T> lowC;
    gsVector<T> uppC;

    index_t EltIdx = 0;
    for (; domIt->good(); domIt->next())
    {
        // compute the L2-distance
        Errs(EltIdx,0) = igaL2DistanceOnElt( geoEval, funcEval, v, v_isParam, domIt);

        // store coordinates of the element
        lowC = domIt->lowerCorner();
        uppC = domIt->upperCorner();
        for( int i=0; i < d; i++)
        {
            Errs( EltIdx , i+1) = lowC[i];
            Errs( EltIdx , i+1+d) = uppC[i];
        }

        EltIdx += 1;
    }
    return Errs;
}

// if u is given as a field, the information is split into
// "function"- and "geometry"-information and passed to igaL2Distance.
template <typename T>
T igaFieldL2Distance(const gsField<T>& u, const gsFunction<T>& v, bool v_isParam)
{
    T dist(0);

    for (int i = 0; i < u.nPatches(); ++i)
    {
        // extract the "function"-part of the gsField
        const gsGeometry<T> & func  = static_cast<const gsGeometry<T> &>( u.function(i) );
        // call igaL2Distance( patch, func, v, v_isParam)
        const T curDist = igaL2Distance( u.patch(i), func, v, v_isParam);
        dist += curDist * curDist;
    }

    return math::sqrt(dist);
}

template <typename T>
T igaFieldL2Distance(const gsField<T>& u, const gsFunction<T>& v, const gsMultiBasis<T>& B, bool v_isParam)
{
    T dist(0);

    for (int i = 0; i < u.nPatches(); ++i)
    {
        // extract the "function"-part of the gsField
        const gsFunction<T> & func  = u.function(i);
        // call igaL2Distance( patch, func, v, v_isParam)
        const T curDist = igaL2Distance( u.patch(i), func, v, B[i], v_isParam);
        dist += curDist * curDist;
    }

    return math::sqrt(dist);
}

template <typename T>
gsVector< gsMatrix<T> > igaFieldL2DistanceEltWiseSq(const gsField<T>& u, const gsFunction<T>& v, bool v_isParam)
{
    gsVector< gsMatrix<T> > Errs( u.nPatches() );

    for (int i = 0; i < u.nPatches(); ++i)
    {
        // extract the "function"-part of the gsField
        const gsGeometry<T> & func  = static_cast<gsGeometry<T> &>( u.function(i) );
        // call igaL2DistanceEltWiseSq( patch, func, v, v_isParam)
        // to get the element-wise squared L2-error on patch i
        Errs[i] = igaL2DistanceEltWiseSq( u.patch(i), func, v, v_isParam);
    }

    return Errs;
}


// H1 norm ///////////////////////////////////////////////////////////////////////////



template <typename T>
T igaH1DistanceOnElt( const std::auto_ptr< gsGeometryEvaluator<T> > & geoEval ,
                      const std::auto_ptr< gsGeometryEvaluator<T> > & funcEval,
                      const gsFunction<T>& v,
                      const bool & v_isParam,
                      const typename gsBasis<T>::domainIter & domIt,
                      const int d)
{
    // compute image of Gauss nodes under geometry mapping as well as Jacobians
    geoEval->evaluateAt(domIt->quNodes);
    funcEval->evaluateAt(domIt->quNodes);
    gsMatrix<T> physGrad_f, func_ders = funcEval->jacobians();// (!) coping

    // get the gradients to columns
    func_ders.transposeInPlace();
    func_ders.resize(d, domIt->numQuNodes() );

    // Evaluate function v
    gsMatrix<T> v_ders = v_isParam ? v.deriv(domIt->quNodes)
                                   : v.deriv( geoEval->values() );

    // get the gradients to columns
    v_ders.transposeInPlace();
    v_ders.resize(d,domIt->numQuNodes() );

    T sum(0.0);
    // perform the quadrature
    for (index_t k = 0; k < domIt->numQuNodes(); ++k) // loop over quadrature nodes
    {
        // Transform the gradients
        geoEval->transformGradients(k, func_ders, physGrad_f);
        if ( v_isParam )
            v_ders.col(k)=geoEval->gradTransforms().block(0, k*d,d,d) * v_ders.col(k);// to do: generalize

        const T weight = domIt->quWeights[k] * fabs( geoEval->measure(k) );
        sum += weight * (physGrad_f - v_ders.col(k)).squaredNorm();
    }
    return sum;
}



template <typename T>
T igaH1Distance(const gsGeometry<T>& patch, 
                const gsGeometry<T>& func,
                const gsFunction<T>& v,
                bool v_isParam)
{
    const int d = func.parDim();

    std::auto_ptr< gsGeometryEvaluator<T> > geoEval ( patch.evaluator(NEED_VALUE    |
                                                                      NEED_GRAD_TRANSFORM  |
                                                                      NEED_MEASURE ));
    // assuming real-valued function
    std::auto_ptr< gsGeometryEvaluator<T> > funcEval ( func.evaluator(NEED_JACOBIAN) );
    
    // degree of the underlying Gauss rule to use
    gsVector<int> numNodes( func.parDim() );
    for (index_t i = 0; i < numNodes.size(); ++i)
        numNodes[i] = func.basis().degree(i) + 1;
    //gsVector<int> numNodes = gsGaussAssembler<T>::getNumIntNodesFor( func.basis() );

    T sum(0);
    typename gsBasis<T>::domainIter domIt = func.basis().makeDomainIterator();
    domIt->computeQuadratureRule( numNodes );
    for (; domIt->good(); domIt->next())
    {
        sum += igaH1DistanceOnElt( geoEval, funcEval, v, v_isParam, domIt, d);
    }
    return math::sqrt(sum);
}

////////////

template <typename T>
T igaH1Distance(const gsGeometry<T>& patch,
                const gsFunction<T>& func,
                const gsFunction<T>& v,
                const gsBasis<T>& B,
                bool v_isParam)
{
    std::auto_ptr< gsGeometryEvaluator<T> > geoEval ( patch.evaluator(NEED_VALUE   |
                                                                      NEED_GRAD_TRANSFORM |
                                                                      NEED_MEASURE ));
    // degree of the underlying Gauss rule to use
    gsVector<int> numNodes( B.dim() );
    for (index_t i = 0; i < numNodes.size(); ++i)
        numNodes[i] = B.degree(i) + 1;
    //gsVector<int> numNodes = gsGaussAssembler<T>::getNumIntNodesFor( B );

    T sum(0);
    typename gsBasis<T>::domainIter domIt = B.makeDomainIterator();
    domIt->computeQuadratureRule( numNodes );
    for (; domIt->good(); domIt->next())
    {
        sum += igaH1DistanceOnElt( geoEval, func, v, v_isParam, domIt);
    }
    return math::sqrt(sum);
}

template <typename T>
T igaFieldH1Distance(const gsField<T>& u, const gsFunction<T>& v, const gsMultiBasis<T>& B, bool v_isParam)
{
    T dist(0);


    for (int i = 0; i < u.nPatches(); ++i)
    {
        // extract the "function"-part of the gsField
        const gsFunction<T> & func  = u.function(i);
        // call igaL2Distance( patch, func, v, v_isParam)
        const T curDist = igaH1Distance( u.patch(i), func, v, B[i], v_isParam);
        dist += curDist * curDist;
    }

    return math::sqrt(dist);
}
template <typename T>
T igaH1DistanceOnElt( const std::auto_ptr< gsGeometryEvaluator<T> > & geoEval ,
                      const gsFunction<T> & func,
                      const gsFunction<T>& v,
                      const bool & v_isParam,
                      const typename gsBasis<T>::domainIter & domIt)
{
    const int d = func.domainDim();

    // compute image of Gauss nodes under geometry mapping as well as Jacobians
    geoEval->evaluateAt(domIt->quNodes);
    gsMatrix<T> func_ders;
    func.deriv_into(domIt->quNodes, func_ders);

    // get the gradients to columns
    func_ders.transposeInPlace();
    func_ders.resize(d, domIt->numQuNodes() );

    // Evaluate function v
    gsMatrix<T> v_ders = v_isParam ? v.deriv(domIt->quNodes)
                                   : v.deriv( geoEval->values() );

    // get the gradients to columns
    v_ders.transposeInPlace();
    v_ders.resize(d,domIt->numQuNodes() );

    // perform the quadrature
    gsMatrix<T> physGrad_f;
    T sum(0.0);
    for (index_t k = 0; k < domIt->numQuNodes(); ++k) // loop over quadrature nodes
    {
        // Transform the gradients
        geoEval->transformGradients(k, func_ders, physGrad_f);
        if ( v_isParam )
            v_ders.col(k)=geoEval->gradTransforms().block(0, k*d,d,d) * v_ders.col(k);// to do: generalize

        const T weight = domIt->quWeights[k] * fabs( geoEval->measure(k) );
        sum += weight * (physGrad_f - v_ders.col(k)).squaredNorm();
    }
    return sum;

}
////////////

template <typename T>
gsMatrix<T> igaH1DistanceEltWiseSq(const gsGeometry<T>& patch,
                                   const gsGeometry<T>& func,
                                   const gsFunction<T>& v,
                                   bool v_isParam)
{
    const int d = func.parDim();

    std::auto_ptr< gsGeometryEvaluator<T> > geoEval ( patch.evaluator(NEED_VALUE    |
                                                                      NEED_GRAD_TRANSFORM  |
                                                                      NEED_MEASURE ));
    // assuming real-valued function
    std::auto_ptr< gsGeometryEvaluator<T> > funcEval ( func.evaluator(NEED_JACOBIAN) );

    // degree of the underlying Gauss rule to use
    gsVector<int> numNodes( func.parDim() );
    for (index_t i = 0; i < numNodes.size(); ++i)
        numNodes[i] = func.basis().degree(i) + 1;
    //gsVector<int> numNodes = gsGaussAssembler<T>::getNumIntNodesFor( func.basis() );

    typename gsBasis<T>::domainIter domIt = func.basis().makeDomainIterator();
    domIt->computeQuadratureRule( numNodes );

    // Matrix to store the element-wise errors
    // in each line:
    // first entry: error
    // second to last entry: lower and upper
    gsMatrix<T> Errs( func.basis().numElements(), 1 + 2*d );

    gsVector<T> lowC;
    gsVector<T> uppC;

    index_t EltIdx = 0;
    for (; domIt->good(); domIt->next())
    {
        // compute the H1-distance
        Errs(EltIdx,0) = igaH1DistanceOnElt( geoEval, funcEval, v, v_isParam, domIt, d);

        // store coordinates of the element
        lowC = domIt->lowerCorner();
        uppC = domIt->upperCorner();
        for( int i=0; i < d; i++)
        {
            Errs( EltIdx , i+1) = lowC[i];
            Errs( EltIdx , i+1+d) = uppC[i];
        }

        EltIdx += 1;
    }
    return Errs;
}

template <typename T>
T igaFieldH1Distance(const gsField<T>& u, const gsFunction<T>& v, bool v_isParam)
{
    T dist(0);

    for (int i = 0; i < u.nPatches(); ++i)
    {
        const gsGeometry<T> & func  = static_cast<gsGeometry<T> &>( u.function(i) );
        const T curDist = igaH1Distance( u.patch(i), func, v, v_isParam);
        dist += curDist * curDist;
    }

    return math::sqrt(dist);
}






template <typename T>
gsVector< gsMatrix<T> > igaFieldH1DistanceEltWiseSq(const gsField<T>& u, const gsFunction<T>& v, bool v_isParam)
{
    index_t N = u.nPatches();
    gsVector< gsMatrix<T> > Errs(N);

    for (int i = 0; i < u.nPatches(); ++i)
    {
        // extract the "function"-part of the gsField
        const gsGeometry<T> & func  = static_cast<gsGeometry<T> &>( u.function(i) );
        // call igaL2DistanceEltWiseSq( patch, func, v, v_isParam)
        // to get the element-wise squared L2-error on patch i
        Errs[i] = igaH1DistanceEltWiseSq( u.patch(i), func, v, v_isParam);
    }

    return Errs;
}



/////// DG norm ///////////////////////////////////////////////////////////

template <typename T>
T igaDGDistanceJump(const gsGeometry<T>& patch1, const gsGeometry<T>& patch2,
                    const gsGeometry<T>& func1,  const gsGeometry<T>& func2, // approximati solution
                    const gsFunction<T>& v1, const gsFunction<T>& v2,	// exact solution
                    const boundaryInterface & bi, // interface
                    const T mu,
                    bool v_isParam)
{
    const int d = func1.parDim();
    GISMO_ASSERT ( d == patch1.geoDim(), "Dimension mismatch" );

    std::auto_ptr< gsGeometryEvaluator<T> > geoEval1 ( patch1.evaluator(NEED_VALUE   |
                                                                        NEED_MEASURE ));

    std::auto_ptr< gsGeometryEvaluator<T> > geoEval2 ( patch2.evaluator(NEED_VALUE   |
                                                                        NEED_MEASURE ));
    // assuming real-valued function
    std::auto_ptr< gsGeometryEvaluator<T> > funcEval1 ( func1.evaluator(NEED_VALUE) );
    
    std::auto_ptr< gsGeometryEvaluator<T> > funcEval2 ( func2.evaluator(NEED_VALUE) );
    
    const boxSide side1 = bi.first();
    const boxSide side2 = bi.second();

    // "DG method not implemented yet for non-matching interfaces");
    
    // assumes matching orientation
    // degree of the underlying Gauss rule to use
    gsVector<int> intNodes1 ( func1.basis().dim() );
    const int dir1 = bi.first().direction();
    for (int i = 0; i < dir1; ++i)
        intNodes1[i] = ( func1.basis().degree(i) + func2.basis().degree(i) + 2 )/ 2 ;
    intNodes1[dir1] = 1;
    for (int i = dir1+1; i < func1.basis().dim(); ++i)
        intNodes1[i] = ( func1.basis().degree(i) + func2.basis().degree(i) + 2 )/ 2 ;

    gsVector<int> intNodes2 ( func2.basis().dim() );
    const int dir2 = bi.second().direction();
    for (int i = 0; i < dir2; ++i)
        intNodes2[i] = ( func1.basis().degree(i) + func2.basis().degree(i) + 2 )/ 2 ;
    intNodes2[dir2] = 1;
    for (int i = dir2+1; i < func1.basis().dim(); ++i)
        intNodes2[i] = ( func1.basis().degree(i) + func2.basis().degree(i) + 2 )/ 2 ;

    
    // Temporaries
    gsVector<T> unormal(d);

    T sum(0);
    // iterator on grid cells on the "right"
    typename gsDomainIterator<T>::uPtr domIter2= func2.basis().makeDomainIterator(side2);

    const int bSize1      = func1.basis().numElements( bi.first() .side() );
    const int bSize2      = func2.basis().numElements( bi.second().side() );
    const int ratio = bSize1 / bSize2;

    int count = 0;
    // iterate over all boundary grid cells on the "left"
    for (typename gsDomainIterator<T>::uPtr domIter1 = func2.basis().makeDomainIterator(side1);
         domIter1->good(); domIter1->next())
    {
        count++;
        // Compute the quadrature rule on both sides
        domIter1->computeQuadratureRule(intNodes1);
        domIter2->computeQuadratureRule(intNodes2);
        
        // compute image of Gauss nodes under geometry mapping as well
        // as Jacobians
        geoEval1->evaluateAt(domIter1->quNodes);
        geoEval2->evaluateAt(domIter2->quNodes);
        
        funcEval1->evaluateAt(domIter1->quNodes);
        funcEval2->evaluateAt(domIter2->quNodes);
        
        gsMatrix<T> func1_vals = funcEval1->values();// (!) coping
        gsMatrix<T> func2_vals = funcEval2->values();// (!) coping

        // exact solution
        gsMatrix<T> v1_vals = v_isParam ? v1.eval(domIter1->quNodes)
                                        : v1.eval( geoEval1->values() );
        
        gsMatrix<T> v2_vals = v_isParam ? v2.eval(domIter2->quNodes)
                                        : v2.eval( geoEval2->values() );

        for (index_t k=0; k!= domIter1->numQuNodes(); ++k)
        {
            // Compute the outer normal vector from patch1
            geoEval1->outerNormal(k, side1, unormal);
            
            // Integral transformation and quadarature weight (patch1)
            // assumed the same on both sides
            const T fff = mu * domIter1->quWeights[k] *  unormal.norm()*(1./domIter1->getCellSize()+1./domIter2->getCellSize());
            
            const T diff = (func1_vals(0,k) - v1_vals(0,k)) - (func2_vals(0,k) - v2_vals(0,k)) ;
            sum += fff * diff*diff;
        }
        if ( count % ratio == 0 ) // next master element ?
        {
            domIter2->next();
        }
    }
    return math::sqrt(sum);
}

/////////
template <typename T>
T igaFieldDGDistance(const gsField<T>& u, const gsFunction<T>& v, bool v_isParam)
{
    T curDist=0;
    T dist = igaFieldH1Distance(u, v, v_isParam);
    dist *= dist; // square

    gsMultiPatch<T> mp = u.patches();
    // iterate over all interfaces
    for ( typename gsMultiPatch<T>::const_iiterator it = mp.iBegin();
          it != mp.iEnd(); ++it ) // *it ---> interface
    {
        const gsGeometry<T> & func1  = static_cast<gsGeometry<T> &>( u.function(it->first().patch) );
        const gsGeometry<T> & func2  = static_cast<gsGeometry<T> &>( u.function(it->second().patch) );

        // get penalty parametr
        //const T h1 = math::pow( (T) func1.basis().size(), -1.0 / func1.basis().dim() );
        //const T h2 = math::pow( (T) func2.basis().size(), -1.0 / func2.basis().dim() );
        const T bdeg = (T)func1.basis().degree(0);
        T mu = ( (bdeg+func1.basis().dim())* (bdeg+1) * 2.0 );

        const bool reverse = func1.basis().numElements(it->first() .side() ) <
                func2.basis().numElements(it->second().side() ) ;

        const boundaryInterface & iFace =( reverse ? it->getInverse() : *it );

        if(!reverse)
        {
            curDist= igaDGDistanceJump( mp.patch(iFace.first().patch),
                                        mp.patch(iFace.second().patch),
                                        func1,
                                        func2,
                                        v, v,
                                        iFace,
                                        mu, // mu
                                        v_isParam);
        }
        else
        {
            curDist= igaDGDistanceJump( mp.patch(iFace.first().patch),
                                        mp.patch(iFace.second().patch),
                                        func2,
                                        func1,
                                        v, v,
                                        iFace,
                                        mu, // mu
                                        v_isParam);
        }

       dist += curDist * curDist;
    }
    return math::sqrt(dist);
}


// Maximum norm //////////////////////////////////////////////////////////////////////


template <typename T>
T computeMaximumNorm(const gsFunction<T>& f, const gsVector<T>& lower, const gsVector<T>& upper, int numSamples)
{
    gsMatrix<T> points = uniformPointGrid(lower, upper, numSamples);

    gsMatrix<T> values;
    f.eval_into( points, values );

    return values.array().abs().maxCoeff();
}

template <typename T>
T computeMaximumNorm(const gsGeometry<T>& geo, const gsFunction<T>& u, bool isParametrized_u, int numSamples)
{
    // note: by specifying the constant function as parametric, we can avoid the geometry evaluations
    // if u is parametric too
    return computeMaximumDistance(geo, u, isParametrized_u, gsConstantFunction<T>(T(0), 1), true, numSamples);
}

template <typename T>
T computeMaximumNorm(const gsField<T>& u, int numSamples)
{
    return computeMaximumDistance(u, gsConstantFunction<T>(T(0), 1), true, numSamples);
}

template <typename T>
T computeMaximumDistance(const gsFunction<T>& f1, const gsFunction<T>& f2, const gsVector<T>& lower, const gsVector<T>& upper, int numSamples)
{
    GISMO_ASSERT( f1.domainDim() == f2.domainDim(), "Functions need to have same domain dimension");
    GISMO_ASSERT( f1.targetDim() == f2.targetDim(), "Functions need to have same target dimension");

    gsMatrix<T> points = uniformPointGrid(lower, upper, numSamples);

    gsMatrix<T> values1, values2;
    f1.eval_into( points, values1 );
    f2.eval_into( points, values2 );

    return (values1 - values2).array().abs().maxCoeff();
}

template <typename T>
T computeMaximumDistance(const gsGeometry<T>& geo, const gsFunction<T>& u, bool isParametrized_u, const gsFunction<T>& v, bool isParametrized_v, int numSamples)
{
    GISMO_ASSERT( u.targetDim() == v.targetDim(), "Functions need to have same target dimension( got "<<u.targetDim() <<" and "<<v.targetDim() );

    // compute the point grid
    gsMatrix<T> range = geo.basis().support();
    gsVector<T> lower = range.col(0), upper = range.col(1);
    gsMatrix<T> points = uniformPointGrid(lower, upper, numSamples);

    // only compute the geometry points if either function is not parametrized
    gsMatrix<T> geo_pts =  (!isParametrized_u || !isParametrized_v) ?
                geo.eval(points) : gsMatrix<T>();

    // evaluate u and v
    gsMatrix<T> u_val = isParametrized_u ? u.eval(points) : u.eval(geo_pts),
            v_val = isParametrized_v ? v.eval(points) : v.eval(geo_pts);

    return (u_val - v_val).array().abs().maxCoeff();
}


template <typename T>
T computeMaximumDistance(const gsField<T>& u, const gsFunction<T>& v, bool isParametrized_v, int numSamples)
{
    T dist = T();

    for (int i = 0; i < u.nPatches(); ++i)
    {
        T curDist = computeMaximumDistance( u.patch(i), u.function(i), u.isParametrized(), v, isParametrized_v, numSamples );
        dist = std::max(dist, curDist);
    }

    return dist;
}


} // namespace gismo
