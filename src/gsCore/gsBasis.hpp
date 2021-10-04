/** @file gsBasis.hpp

    @brief Provides implementation of Basis default operations.

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
gsBasis<T>::gsBasis(const gsBasis& other) : Base(other)
{ }

template<class T>
gsBasis<T>::~gsBasis()
{ }

template<class T>
gsBasisFun<T> gsBasis<T>::function(index_t i) const
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
    gsMatrix<index_t> actives;

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
    gsMatrix<index_t> ind;

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
    gsMatrix<index_t> actives;

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
    gsMatrix<index_t> actives;

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
    gsMatrix<index_t> actives;

    this->evalAllDers_into(u,n,B);
    this->active_into(u,actives);

    // for derivatives 0 to n, evaluate the function by linear combination
    // of coefficients with the respective function values/derivatives
    for( unsigned i = 0; i <= n; i++)
        linearCombination_into( coefs, actives, B[i], result[i] );
}


template<class T>
void gsBasis<T>::linearCombination_into(const gsMatrix<T> & coefs,
                                        const gsMatrix<index_t> & actives,
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
        for ( index_t i = 0; i < actives.rows(); ++i )  // for all nonzero basis functions
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
    gsMatrix<index_t> act;

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

    // gsDebug <<"gs Interpolate error : " << solver.error() << std::"\n";
    // gsDebug <<"gs Interpolate iters : " << solver.iterations() << std::"\n";
    // gsDebug <<"intpl sol : " << x.transpose() << std::"\n";

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
void gsBasis<T>::anchors_into(gsMatrix<T>&) const
{ GISMO_NO_IMPLEMENTATION }

template<class T> inline
void gsBasis<T>::anchor_into(index_t, gsMatrix<T>&) const
{ GISMO_NO_IMPLEMENTATION }


template<class T> inline
void gsBasis<T>::connectivityAtAnchors(gsMesh<T> & mesh) const
{
    gsMatrix<T> nodes = anchors();
    nodes.transposeInPlace();// coefficient vectors have ctrl points at rows
    connectivity(nodes, mesh);
}

template<class T>
void gsBasis<T>::connectivity(const gsMatrix<T> &, gsMesh<T> &) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::active_into(const gsMatrix<T> &, gsMatrix<index_t>&) const
{ GISMO_NO_IMPLEMENTATION }

template <class T>
bool gsBasis<T>::isActive(const index_t, const gsVector<T> &) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::numActive_into(const gsMatrix<T> &, gsVector<index_t>&) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::activeCoefs_into(const gsVector<T> &, const gsMatrix<T> &,
                                  gsMatrix<T>&) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
gsMatrix<index_t>
gsBasis<T>::allBoundary() const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
gsMatrix<index_t>
gsBasis<T>::boundaryOffset(boxSide const &,index_t) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
index_t
gsBasis<T>::functionAtCorner(boxCorner const &) const
{ GISMO_NO_IMPLEMENTATION }

/// @cond
template<class T>
gsBasis<T>* gsBasis<T>::boundaryBasis_impl(boxSide const &) const
{ GISMO_NO_IMPLEMENTATION }
/// @endcond

template<class T>
typename gsBasis<T>::uPtr gsBasis<T>::componentBasis(boxComponent b) const
{
    GISMO_ASSERT( b.totalDim() == this->dim(), "The dimensions do not agree." );

    const short_t dim = this->dim();

    uPtr result;
    short_t d=0;
    for (short_t i=0; i<dim; ++i)
    {
        boxComponent::location loc = b.locationForDirection(i);
        if (loc)
        {
            if (result)
                result = result->boundaryBasis( boxSide(loc+2*d) );
            else
                result =   this->boundaryBasis( boxSide(loc+2*d) );
        }
        else
            ++d;
    }

    if (!result)
        result = clone();

    return result;
}

template<class T>
typename gsBasis<T>::uPtr gsBasis<T>::componentBasis_withIndices(boxComponent b, gsMatrix<index_t>& indices, bool noBoundary) const
{
    GISMO_ASSERT( b.totalDim() == this->dim(), "The dimensions do not agree." );
    const short_t dim = this->dim();

    uPtr result;
    short_t d=0;
    for (short_t i=0; i<dim; ++i)
    {
        boxComponent::location loc = b.locationForDirection(i);
        if (loc)
        {
            if (result)
            {
                gsMatrix<index_t> tmp = result->boundary( boxSide(loc+2*d) );
                for (index_t j=0; j<tmp.size(); ++j)
                    tmp(j,0) = indices(tmp(j,0),0);
                tmp.swap(indices);
                result = result->boundaryBasis( boxSide(loc+2*d) );
            }
            else
            {
                indices = this->boundary( boxSide(loc+2*d) );
                result = this->boundaryBasis( boxSide(loc+2*d) );
            }
        }
        else
            ++d;
    }

    if (!result)
    {
        result = clone();
        const index_t sz = this->size();
        indices.resize(sz,1);
        for (index_t i=0;i<sz;++i)
            indices(i,0) = i;
    }

    if (noBoundary && d > 0)
    {

        gsMatrix<index_t> bdy_indices = result->allBoundary();

        const index_t indices_sz = indices.rows();
        const index_t bdy_indices_sz = bdy_indices.rows();

        // Copy all entries from indices to indices_cleaned except
        // those with indices in bdy_indices

        gsMatrix<index_t> indices_cleaned(indices_sz - bdy_indices_sz, 1);
        index_t j = 0, t = 0;
        for (index_t i = 0; i < indices_sz; ++i)
        {
            if (util::greater(i, bdy_indices(j, 0)) && j < bdy_indices_sz)
                ++j;
            if (util::less(i, bdy_indices(j, 0)) || j == bdy_indices_sz)
            {
                indices_cleaned(t, 0) = indices(i, 0);
                ++t;
            }
        }
        GISMO_ASSERT(t == indices_cleaned.rows(), "Internal error.");
        indices.swap(indices_cleaned);
    }

    return result;

}

template<class T>
gsMatrix<T> gsBasis<T>::support() const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
gsMatrix<T> gsBasis<T>::support(const index_t &) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
gsMatrix<T> gsBasis<T>::supportInterval(index_t dir) const
{ return support().row(dir); }

template<class T>
void gsBasis<T>::eval_into(const gsMatrix<T> &, gsMatrix<T>&) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::evalSingle_into(index_t, const gsMatrix<T> &, gsMatrix<T>&) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::deriv_into(const gsMatrix<T> &, gsMatrix<T>&) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::derivSingle_into(index_t,
                                  const gsMatrix<T> &,
                                  gsMatrix<T>&) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::deriv2_into(const gsMatrix<T> &, gsMatrix<T>&) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::deriv2Single_into(index_t,
                                   const gsMatrix<T> &,
                                   gsMatrix<T>&) const
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
        GISMO_ERROR("evalAllDers implemented for order up to 2<"<<n<< " for "<<*this);
        break;
    }
}

template<class T>
void gsBasis<T>::evalAllDersSingle_into(index_t, const gsMatrix<T> &,
                                        int, gsMatrix<T>&) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::evalDerSingle_into(index_t, const
                                    gsMatrix<T> &, int,
                                    gsMatrix<T>&) const
{ GISMO_NO_IMPLEMENTATION }


template<class T>
typename gsBasis<T>::uPtr gsBasis<T>::create() const
{ GISMO_NO_IMPLEMENTATION }


template<class T>
typename gsBasis<T>::uPtr gsBasis<T>::tensorize(const gsBasis &) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
typename gsBasis<T>::domainIter
gsBasis<T>::makeDomainIterator() const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
typename gsBasis<T>::domainIter
gsBasis<T>::makeDomainIterator(const boxSide &) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
size_t gsBasis<T>::numElements() const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
size_t gsBasis<T>::numElements(boxSide const &) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
size_t gsBasis<T>::elementIndex(const gsVector<T> &) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
gsMatrix<T> gsBasis<T>::elementInSupportOf(index_t j) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
const gsBasis<T>& gsBasis<T>::component(short_t) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
gsBasis<T>& gsBasis<T>::component(short_t i)
{ return const_cast<gsBasis<T>&>(const_cast<const gsBasis<T>*>(this)->component(i));}

template<class T>
std::vector<index_t> gsBasis<T>::asElements(gsMatrix<T> const &, int) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::refine(gsMatrix<T> const &, int)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::refineElements(std::vector<index_t> const &)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::refineElements_withCoefs(gsMatrix<T> &,std::vector<index_t> const &)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::uniformRefine(int, int)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::uniformRefine_withCoefs(gsMatrix<T>& , int , int )
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::uniformRefine_withTransfer(gsSparseMatrix<T,RowMajor> &,
                                            int, int)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::uniformCoarsen(int)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::uniformCoarsen_withTransfer(gsSparseMatrix<T,RowMajor> &,
                                            int)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::degreeElevate(short_t const &, short_t const)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::degreeReduce(short_t const &, short_t const)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::degreeIncrease(short_t const &, short_t const)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::degreeDecrease(short_t const &, short_t const)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::setDegree(short_t const& i)
{
    const short_t dm = this->dim();
    for (short_t k = 0; k!=dm; ++k)
    {
        const short_t p = this->degree(k);
        if ( i > p )
        {
            this->degreeElevate(i-p, k);
        }
        else if  ( i < p )
        {
            this->degreeReduce(p-i, k);
        }
    }
}

template<class T>
void gsBasis<T>::setDegreePreservingMultiplicity(short_t const& i)
{
    for ( short_t d = 0; d < dim(); ++ d )
    {
        if ( i > degree(d) )
            degreeIncrease(i-degree(d),d);
        else if ( i < degree(d) )
            degreeDecrease(-i+degree(d),d);
    }
}

template<class T>
void gsBasis<T>::elevateContinuity(int const &)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::reduceContinuity(int const &)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
gsDomain<T> * gsBasis<T>::domain() const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
short_t gsBasis<T>::maxDegree() const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
short_t gsBasis<T>::minDegree() const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
short_t gsBasis<T>::totalDegree() const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
short_t gsBasis<T>::degree(short_t) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::reverse()
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::matchWith(const boundaryInterface &, const gsBasis<T> &,
               gsMatrix<index_t> &, gsMatrix<index_t> &) const
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
void gsBasis<T>::linearComb(const gsMatrix<index_t>  & actives,
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
