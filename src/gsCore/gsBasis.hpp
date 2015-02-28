/** @file gsBasis.hpp

    @brief Provides implementation of Basis default operatiions.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsUtils/gsCollocationMatrix.h>

namespace gismo
{

// Evaluates a linear combination of basis functions (default implementation)
template<class T>
void gsBasis<T>::eval_into(const gsMatrix<T> &u, const gsMatrix<T> & coefs, gsMatrix<T>& result) const {

    result.resize(coefs.cols(), u.cols()) ;
    gsMatrix<T> B ;
    gsMatrix<unsigned> ind;

    this->eval_into(u,B)    ; // col j = nonzero basis functions at column point u(..,j)
    this->active_into(u,ind);  // col j = indices of active functions at column point u(..,j)

    for ( index_t j=0; j< u.cols() ; j++ ) // for all points (columns of u)
    {
        result.col(j) =  coefs.row( ind(0,j) ) * B(0,j)  ;
        for ( index_t i=1; i< ind.rows() ; i++ ) // for all nonzero basis functions
            result.col(j)  +=   coefs.row( ind(i,j) ) * B(i,j)  ;
    }

}


// Evaluates the Jacobian of the function given by coefs (default implementation)
template<class T>
void gsBasis<T>::deriv_into(const gsMatrix<T> &u, const gsMatrix<T> & coefs, gsMatrix<T>& result) const {  
    unsigned n = coefs.cols();
    unsigned numPts = u.cols();       // at how many points to evaluate the gradients
    int pardim = this->dim();

    result.setZero( n, numPts * pardim );
    gsMatrix<T> B;
    gsMatrix<unsigned> ind;

    this->deriv_into(u,B);     // col j = nonzero derivatives at column point u(..,j)
    this->active_into(u,ind);  // col j = indices of active functions at column point u(..,j)
  
    for (unsigned j = 0; j < numPts; ++j)
        for (unsigned k=0; k<n; ++k ) // for all rows of the jacobian
            for ( index_t i=0; i< ind.rows() ; i++ ) // for all nonzero basis functions)
            {
                result.block(k,j*pardim,  1, pardim).noalias() +=  
                    coefs(ind(i,j), k) * B.block(i*pardim, j, pardim, 1).transpose(); 
            }
}

// Evaluates the second derivatives of the function given by coefs (default implementation)
template<class T>
void gsBasis<T>::deriv2_into(const gsMatrix<T> &u, const gsMatrix<T> & coefs, gsMatrix<T>& result) const 
{  
    unsigned n = coefs.cols();
    unsigned numPts = u.cols();  // at how many points to evaluate the gradients

    gsMatrix<T> B;
    gsMatrix<unsigned> ind;
    this->deriv2_into(u,B);    // col j = nonzero 2nd derivatives at column point u(..,j)
    this->active_into(u,ind);  // col j = indices of active functions at column point u(..,j)

    const unsigned stride = B.rows() / ind.rows() ;

    result.setZero( n*stride, numPts);  
    for (unsigned j = 0; j < numPts; ++j) // For column point u(..,j)
    {
        for ( index_t i=0; i< ind.rows() ; i++ )  // for all nonzero basis functions
            for (unsigned k=0; k<n; ++k )         // for all components of the geometry
                result.block (stride * k, j, stride, 1).noalias() +=  
                    coefs(ind(i,j), k) * B.block( stride * i, j, stride, 1); 
    }
}


template<class T>
inline gsMatrix<T> * gsBasis<T>::laplacian(const gsMatrix<T> & u ) const 
{
    gsMatrix<T> tmp;
    this->deriv2_into(u,tmp);
    gsMatrix<T> * res = new gsMatrix<T>(tmp.colwise().sum());
    return res;
}

template<class T> inline
void gsBasis<T>::collocationMatrix(const gsMatrix<T> & u, gsSparseMatrix<T> & result) const 
{
    gsCollocationMatrix_into (*this, u, result);
}

template<class T> inline
gsGeometry<T> * gsBasis<T>::interpolate( gsMatrix<T> const& vals,
                                         gsMatrix<T> const& pts) const
{
    GISMO_ASSERT (dim()  == pts.rows() , "Wrong dimension of the points("<<
                  pts.rows()<<", expected "<<dim() <<").");
    GISMO_ASSERT (size() == pts.cols() , "Expecting as many points as the basis functions." );
    GISMO_ASSERT (size() == vals.cols(), "Expecting as many values as the number of points." );

    gsSparseMatrix<T>  Cmat;
    gsCollocationMatrix_into( *this , pts,Cmat );
    gsMatrix<T> x ( size(), vals.rows());

    //Eigen::ConjugateGradient< gsSparseMatrix<T> > solver(Cmat);// for symmetric - does not apply here
    //Eigen::BiCGSTAB<gsSparseMatrix<T>, Eigen::DiagonalPreconditioner<T> > solver;
    //Eigen::BiCGSTAB<gsSparseMatrix<T>, Eigen::IdentityPreconditioner > solver;
    Eigen::BiCGSTAB< gsSparseMatrix<T>,  Eigen::IncompleteLUT<T> > solver( Cmat );

    // Solves for many right hand side  columns
    x =  solver.solve( vals.transpose() );
  
    // gsInfo <<"gs Interpolate error : " << solver.error() << std::"\n";
    // gsInfo <<"gs Interpolate iters : " << solver.iterations() << std::"\n";
    // gsInfo <<"intpl sol : " << x.transpose() << std::"\n";
  
    return makeGeometry( give(x) );
}

template<class T> inline
gsGeometry<T> * gsBasis<T>::interpolate(gsMatrix<T> const & vals) const
{
    GISMO_ASSERT (this->size() == vals.cols(), 
                  "Expecting as many values as the number of basis functions." );
    gsMatrix<T> pts;
    anchors_into(pts);
    return interpolate(vals, pts);
}


template<class T> inline
void gsBasis<T>::anchors_into(gsMatrix<T>& result) const
{ GISMO_NO_IMPLEMENTATION }

template<class T> inline
void gsBasis<T>::anchor_into(unsigned i, gsMatrix<T>& result) const
{ GISMO_NO_IMPLEMENTATION }


template<class T> inline
void gsBasis<T>::connectivity(gsMesh<T> & mesh) const
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

template<class T>
void gsBasis<T>::numActive_into(const gsMatrix<T> & u, gsVector<unsigned>& result) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::activeCoefs_into(const gsVector<T> & u, const gsMatrix<T> & coefs, 
                                  gsMatrix<T>& result) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
gsMatrix<unsigned> *
gsBasis<T>::boundary( ) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
gsMatrix<unsigned> *
gsBasis<T>::boundary(boxSide const & s ) const
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
void gsBasis<T>::evalAllDers_into(const gsMatrix<T> & u, int n, gsMatrix<T>& result) const
{ GISMO_NO_IMPLEMENTATION }

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
gsBasis<T> * gsBasis<T>::create() const
{ GISMO_NO_IMPLEMENTATION }


template<class T>
gsBasis<T> * gsBasis<T>::tensorize(const gsBasis & other) const 
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
{ 
    return this->component( s.direction() ).numElements();
}

template<class T>
int gsBasis<T>::elementIndex(const gsVector<T> & u) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
gsBasis<T>& gsBasis<T>::component(unsigned i) const
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::refine(gsMatrix<T> const & boxes)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::refine(gsMatrix<T> const & boxes, int refExt)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::refineElements(std::vector<unsigned> const & boxes)
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
void gsBasis<T>::degreeElevate(int const & i)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::degreeElevateComponent(unsigned dir, int const & i)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::degreeReduce(int const & i)
{ GISMO_NO_IMPLEMENTATION }

template<class T>
void gsBasis<T>::setDegree(int const& i)
{ GISMO_NO_IMPLEMENTATION }

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


//template<class T> inline
// void gsBasis<T>::getLinearCombination(
// const gsMatrix<T>         & scalars,
// const gsMatrix<T> * const & coefs, 
// const gsMatrix<unsigned>  & indices, 
// gsMatrix<T>&                result )
// {
// 	result.resize(coefs->cols(), scalars.cols()) ;

// 	for ( index_t j=0; j< scalars.cols() ; j++ ) // for all values
// 	    {
// 		result.col(j).noalias() =  coefs.row( indices(0,j) ) * scalars(0,j) ;
// 		for ( index_t i=1; i< indices.rows() ; i++ ) // TO DO: numActive
// 		    result.col(j)  +=   coefs.row( indices(i,j) ) * scalars(i,j)  ; 
// 	    }
// }


}; // namespace gismo
