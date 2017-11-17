/** @file gsBasis.hpp

    @brief Provides implementation of Basis default operatiions.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsBasisFun.h>
#include <gsCore/gsDomainIterator.h>
#include <gsCore/gsBoundary.h>
#include <gsCore/gsGeometry.h>

namespace gismo
{

template<class T>
gsBasis<T>::gsBasis()
{ }

template<class T>
gsBasis<T>::gsBasis(const gsBasis& other)
{ }
        
template<class T>
gsBasis<T>::~gsBasis()
{ }

template<class T>
gsBasisFun<T> gsBasis<T>::function(unsigned i) const
{
    return gsBasisFun<T>(*this,i);
}


// Evaluates a linear combination of basis functions (default implementation)
template<class T>
void gsBasis<T>::evalFunc_into(const gsMatrix<T> &u, 
                           const gsMatrix<T> & coefs, 
                           gsMatrix<T>& result) const 
{
    gsMatrix<T> B ;
    gsMatrix<unsigned> actives;

    // compute function values
    this->eval_into(u,B);
    // compute active functions
    this->active_into(u,actives);

    // compute result as linear combination of
    // "coefs(actives)" and B
    linearCombination_into( coefs, actives, B, result );
}


// Evaluates the Jacobian of the function given by coefs (default implementation)
// For each point, result contains a geomDim x parDim matrix block containing the Jacobian matrix
template<class T>
void gsBasis<T>::jacobianFunc_into(const gsMatrix<T> &u, const gsMatrix<T> & coefs, gsMatrix<T>& result) const 
{  
    const index_t n = coefs.cols();
    const index_t numPts = u.cols();       // at how many points to evaluate the gradients
    const index_t pardim = this->dim();

    result.setZero( n, pardim * numPts );
    gsMatrix<T> B;
    gsMatrix<unsigned> ind;

    this->deriv_into(u,B);
    this->active_into(u,ind);
    const index_t numAct=ind.rows();

    for (index_t p = 0; p < numPts; ++p) // p = point
        for (index_t c=0; c<n; ++c )     // c = component
            for ( index_t a=0; a< numAct ; ++a ) // a = active function
            {
                result.block(c,p*pardim,  1, pardim).noalias() +=   
                    coefs(ind(a,p), c) * B.block(a*pardim, p, pardim, 1).transpose();  
            }
}

// Evaluates the partial derivatives of the function given by coefs (default implementation)
// For each point, result contains one column with stacked gradients of the components
template<class T>
void gsBasis<T>::derivFunc_into(const gsMatrix<T> &u, const gsMatrix<T> & coefs, gsMatrix<T>& result) const 
{
    gsMatrix<T>        B;
    gsMatrix<unsigned> actives;

    // compute first derivatives
    this->deriv_into(u,B);
    // compute active functions
    this->active_into(u,actives);

    // compute result as linear combination of
    // "coefs(actives)" and B
    linearCombination_into( coefs, actives, B, result );
}

// Evaluates the second derivatives of the function given by coefs (default implementation)
template<class T>
void gsBasis<T>::deriv2Func_into(const gsMatrix<T> & u,
                                 const gsMatrix<T> & coefs, 
                                 gsMatrix<T>& result) const 
{  
    gsMatrix<T> B;
    gsMatrix<unsigned> actives;

    // compute second derivatives
    this->deriv2_into(u,B);
    // compute active functions
    this->active_into(u,actives);

    // compute result as linear combination of
    // "coefs(actives)" and B
    linearCombination_into( coefs, actives, B, result );
}

template<class T>
void gsBasis<T>::evalAllDersFunc_into(const gsMatrix<T> &u,
                                      const gsMatrix<T> & coefs,
                                      const unsigned n,
                                      std::vector< gsMatrix<T> >& result) const
{
    // resize result so that it will hold
    // function values and up to the n-th derivatives
    result.resize(n+1);

    // B will contain the derivatives up to order n
    std::vector< gsMatrix<T> >B;
    // actives will contain the indices of the basis functions
    // which are active at the evaluation points
    gsMatrix<unsigned> actives;

    this->evalAllDers_into(u,n,B);
    this->active_into(u,actives);

    // for derivatives 0 to n, evaluate the function by linear combination
    // of coefficients with the respective function values/derivatives
    for( unsigned i = 0; i <= n; i++)
        linearCombination_into( coefs, actives, B[i], result[i] );
}


template<class T>
void gsBasis<T>::linearCombination_into(const gsMatrix<T> & coefs,
                                        const gsMatrix<unsigned> & actives,
                                        const gsMatrix<T> & values,
                                        gsMatrix<T> & result)
{
    const index_t numPts = values.cols() ;
    const index_t tarDim = coefs.cols()  ;
    const index_t stride = values.rows() / actives.rows();

    GISMO_ASSERT( actives.rows() * stride == values.rows(), 
                  "Number of values and actives does not fit together");

    result.resize( tarDim * stride, numPts );
    result.setZero();

    for ( index_t pt = 0; pt < numPts; ++pt ) // For pt, i.e., for every column of u
        for ( index_t i = 0; i < index_t( actives.rows() ); ++i )  // for all nonzero basis functions
            for ( index_t c = 0; c < tarDim; ++c )      // for all components of the geometry
            {
                result.block( stride * c, pt, stride, 1).noalias() +=
                    coefs( actives(i,pt), c) * values.block( stride * i, pt, stride, 1);
            }
}


template<class T>
inline gsMatrix<T> gsBasis<T>::laplacian(const gsMatrix<T> & u ) const
{
    gsMatrix<T> tmp;
    this->deriv2_into(u,tmp);
    return tmp.colwise().sum();
}

template<class T> inline
void gsBasis<T>::collocationMatrix(const gsMatrix<T> & u, gsSparseMatrix<T> & result) const 
{
    result.resize( u.cols(), this->size() );
    gsMatrix<T> ev;
    gsMatrix<unsigned> act;

    eval_into  (u.col(0), ev);
    active_into(u.col(0), act);
    result.reservePerColumn( act.rows() );
    for (index_t i=0; i!=act.rows(); ++i)
        result.insert(0, act.at(i) ) = ev.at(i);

    for (index_t k=1; k!=u.cols(); ++k)
    {
        eval_into  (u.col(k), ev );
        active_into(u.col(k), act);
        for (index_t i=0; i!=act.rows(); ++i)
            result.insert(k, act.at(i) ) = ev.at(i);
    }

    result.makeCompressed();
}

template<class T> inline
memory::unique_ptr<gsGeometry<T> > gsBasis<T>::interpolateData( gsMatrix<T> const& vals,
                                         gsMatrix<T> const& pts) const
{
    GISMO_ASSERT (dim()  == pts.rows() , "Wrong dimension of the points("<<
                  pts.rows()<<", expected "<<dim() <<").");
    GISMO_ASSERT (this->size() == pts.cols() , "Expecting as many points as the basis functions." );
    GISMO_ASSERT (this->size() == vals.cols(), "Expecting as many values as the number of points." );

    gsSparseMatrix<T>  Cmat;
    collocationMatrix(pts, Cmat);
    gsMatrix<T> x ( this->size(), vals.rows());

    // typename gsSparseSolver<T>::BiCGSTABIdentity solver( Cmat );
    // typename gsSparseSolver<T>::BiCGSTABDiagonal solver( Cmat );
    // typename gsSparseSolver<T>::QR solver( Cmat );
    typename gsSparseSolver<T>::BiCGSTABILUT solver( Cmat );

    // Solves for many right hand side  columns
    x =  solver.solve( vals.transpose() );
  
    // gsInfo <<"gs Interpolate error : " << solver.error() << std::"\n";
    // gsInfo <<"gs Interpolate iters : " << solver.iterations() << std::"\n";
    // gsInfo <<"intpl sol : " << x.transpose() << std::"\n";
  
    return makeGeometry( give(x) );
}

template<class T> inline
memory::unique_ptr<gsGeometry<T> > gsBasis<T>::interpolateAtAnchors(gsMatrix<T> const & vals) const
{
    GISMO_ASSERT (this->size() == vals.cols(), 
                  "Expecting as many values as the number of basis functions." );
    gsMatrix<T> pts;
    anchors_into(pts);
    return interpolateData(vals, pts);
}

/*
template<class T> inline
gsGeometry<T> * gsBasis<T>::interpolateAtAnchors(gsFunction<T> const & func) const
{
    gsMatrix<T> pts, vals;
    anchors_into(pts);
    func.eval_into(pts, vals);
    return interpolateData(vals, pts);
}

template<class T>
gsGeometry<T> * gsBasis<T>::projectL2(gsFunction<T> const & func) const
{

    return NULL;
}
*/

template<class T> inline
void gsBasis<T>::anchors_into(gsMatrix<T>& result) const
{ GISMO_NO_IMPLEMENTATION }

template<class T> inline
void gsBasis<T>::anchor_into(unsigned i, gsMatrix<T>& result) const
{ GISMO_NO_IMPLEMENTATION }


template<class T> inline
void gsBasis<T>::connectivityAtAnchors(gsMesh<T> & mesh) const
{ 
    gsMatrix<T> nodes = anchors();
    nodes.transposeInPlace();// coefficient vectors have ctrl points at rows
    connectivity(nodes, mesh); 
}

template<class T>
void gsBasis<T>::connectivity(const gsMatrix<T> & nodes, gsMesh<T> & mesh) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::active_into(const gsMatrix<T> & u, gsMatrix<unsigned>& result) const
{ GISMO_NO_IMPLEMENTATION }

template <class T>
bool gsBasis<T>::isActive(const unsigned i, const gsVector<T> & u) const 
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::numActive_into(const gsMatrix<T> & u, gsVector<unsigned>& result) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::activeCoefs_into(const gsVector<T> & u, const gsMatrix<T> & coefs, 
                                  gsMatrix<T>& result) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
gsMatrix<unsigned>
gsBasis<T>::allBoundary( ) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
gsMatrix<unsigned>
gsBasis<T>::boundaryOffset(boxSide const & s,unsigned offset) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
unsigned
gsBasis<T>::functionAtCorner(boxCorner const & c) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
gsBasis<T> * gsBasis<T>::boundaryBasis(boxSide const & s) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
gsMatrix<T> gsBasis<T>::support() const 
{ GISMO_NO_IMPLEMENTATION }

template<class T>
gsMatrix<T> gsBasis<T>::support(const unsigned & i) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
gsMatrix<T> gsBasis<T>::supportInterval(unsigned dir) const 
{ return support().row(dir); }

template<class T>
void gsBasis<T>::eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const
    { GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::evalSingle_into(unsigned i, const gsMatrix<T> & u, gsMatrix<T>& result) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::derivSingle_into(unsigned i, 
                                  const gsMatrix<T> & u, 
                                  gsMatrix<T>& result ) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::deriv2_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::deriv2Single_into(unsigned i, 
                                   const gsMatrix<T> & u, 
                                   gsMatrix<T>& result ) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::evalAllDers_into(const gsMatrix<T> & u, int n, 
                                  std::vector<gsMatrix<T> >& result) const
{
    result.resize(n+1);
    
    switch(n)
    {
    case 0:
        eval_into(u, result[0]);
        break;
    case 1:
        eval_into (u, result[0]);
        deriv_into(u, result[1]);
        break;
    case 2:
        eval_into  (u, result[0]);
        deriv_into (u, result[1]);
        deriv2_into(u, result[2]);
        break;
    default:
        GISMO_ERROR("evalAllDers implemented for order upto 2<"<<n<< " for "<<*this);
        break;
    }
}

template<class T>
void gsBasis<T>::evalAllDersSingle_into(unsigned i, const gsMatrix<T> & u, 
                                        int n, gsMatrix<T>& result) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::evalDerSingle_into(unsigned i, const 
                                    gsMatrix<T> & u, int n, 
                                    gsMatrix<T>& result) const
{ GISMO_NO_IMPLEMENTATION }


template<class T>
typename gsBasis<T>::uPtr gsBasis<T>::create() const
{ GISMO_NO_IMPLEMENTATION }


template<class T>
typename gsBasis<T>::uPtr gsBasis<T>::tensorize(const gsBasis & other) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
typename gsBasis<T>::domainIter
gsBasis<T>::makeDomainIterator() const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
typename gsBasis<T>::domainIter
gsBasis<T>::makeDomainIterator(const boxSide & s) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
int gsBasis<T>::numElements() const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
int gsBasis<T>::numElements(boxSide const & s) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
int gsBasis<T>::elementIndex(const gsVector<T> & u) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
const gsBasis<T>& gsBasis<T>::component(unsigned i) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
gsBasis<T>& gsBasis<T>::component(unsigned i)
{ return const_cast<gsBasis<T>&>(component(i));}

template<class T>
void gsBasis<T>::refine(gsMatrix<T> const & boxes, int refExt)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::refineElements(std::vector<unsigned> const & boxes)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::refineElements_withCoefs(gsMatrix<T> & coefs,std::vector<unsigned> const & boxes)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::uniformRefine(int numKnots, int mul)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::uniformRefine_withCoefs(gsMatrix<T>& coefs, int numKnots, int mul)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::uniformRefine_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, 
                                            int numKnots, int mul)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::degreeElevate(int const & i, int dir)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::degreeReduce(int const & i)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::degreeIncrease(int const & i, int dir)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::degreeDecrease(int const & i, int dir)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::setDegree(int const& i)
{ 
    //TODO: If the degree is not the same in all directions, then this does not what is specified
    const int p = maxDegree();
    if ( i > p )
    {
        degreeElevate(i-p); 
    }
    else if  ( i < p )
    {
        degreeReduce(p-i); 
    }
}

template<class T>
void gsBasis<T>::setDegreePreservingMultiplicity(int const& i)
{ 
    for ( index_t d = 0; d < dim(); ++ d )
    {
        if ( i > degree(d) )
            degreeIncrease(i-degree(d),d);
        else if ( i < degree(d) )
            degreeDecrease(-i+degree(d),d);
    }
}


template<class T>
void gsBasis<T>::reduceContinuity(int const & i)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
gsDomain<T> * gsBasis<T>::domain() const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
int gsBasis<T>::maxDegree() const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
int gsBasis<T>::minDegree() const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
int gsBasis<T>::totalDegree() const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
int gsBasis<T>::degree(int i) const 
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::reverse()
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::matchWith(const boundaryInterface & bi, const gsBasis<T> & other,
               gsMatrix<unsigned> & bndThis, gsMatrix<unsigned> & bndOther) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
T gsBasis<T>::getMinCellLength() const
{
    const domainIter it = this->makeDomainIterator();
    T h = 0;
    for (; it->good(); it->next() )
    {
        const T sz = it->getMinCellLength();
        if ( sz < h || h == 0 ) h = sz;
    }
    return h;
}

template<class T>
T gsBasis<T>::getMaxCellLength() const
{
    const domainIter it = this->makeDomainIterator();
    T h = 0;
    for (; it->good(); it->next() )
    {
        const T sz = it->getMaxCellLength();
        if ( sz > h ) h = sz;
    }
    return h;
}

// gsBasis<T>::linearComb(active, evals, m_tmpCoefs, result);
// gsBasis<T>::jacobianFromGradients(active, grads, m_tmpCoefs, result);

/*
template<class T>
void gsBasis<T>::linearComb(const gsMatrix<unsigned>  & actives, 
                            const gsMatrix<T>         & basisVals,
                            const gsMatrix<T>         & coefs, 
                            gsMatrix<T>&                result )
{
    // basisVals.rows()==1 (or else basisVals.rows() == coefs.cols() and .cwiseProd)
	result.resize(coefs.cols(), basisVals.cols()) ;

	for ( index_t j=0; j!=basisVals.cols() ; j++ ) // for all basis function values
    {
        //todo grab result.col(j)
		result.col(j) =  basisVals(0,j) * coefs.row( actives(0,j) ) ;//transpose ?
		for ( index_t i=1; i< actives.rows() ; i++ )
		    result.col(j) += basisVals(i,j) * coefs.row( actives(i,j) ) ; 
    }
}
*/

}; // namespace gismo
