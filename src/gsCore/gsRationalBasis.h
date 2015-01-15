/** @file gsRationalBasis.h

    @brief Provides declaration of RationalBasis class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/
 
#pragma once

#include <gsCore/gsBasis.h>

namespace gismo
{

  /** \brief
      Class that creates a rational counterpart for a given basis.

      A rational basis holds an inner basis of the given source type,
      and a matrix of coefficients defining a weight function in terms
      of the source basis.

      If \em w_i is the i-th weight coefficient and \em b_i the i-th basis
      function of the source basis, then the i-th basis function of the
      resulting rational basis is given by
        r_i(u) = b_i(u) w_i / (sum(b_j(u) w_j, j = 1, ..., size))

      If the weights are constant, the rational basis is identical
      to the source basis.

      Example: the rational version of a B-spline basis is a NURBS basis.

      \tparam SrcT the basis of which to create a rational version

      \ingroup basis
      \ingroup Core
  */

template<class SrcT>
class gsRationalBasis : public gsBasis< typename SrcT::Scalar_t >
{
public:
    
    typedef typename SrcT::Scalar_t Scalar_t;
    typedef Scalar_t T;
    
    // Obtain dimension from the source basis
    static const int Dim = SrcT::Dim;
    
    static const bool IsRational = true;
    
    typedef gsBasis<T> Base;
    
    /// Associated source basis type
    typedef SrcT SourceBasis;
    
    /// Associated geometry type
    typedef typename
    gsTraits<SrcT,Dim>::RationalGeometryType  GeometryType;
    
    /// Associated boundary basis type
    typedef typename
    gsTraits<SrcT,Dim>::RationalBoundaryType  BoundaryBasisType;
    
public:

    /// Default empty constructor
    gsRationalBasis() :  Base() { }
    
    /// Construct a rational counterpart of basis
    gsRationalBasis(SrcT* basis)
        : m_src(basis)
    { 
        m_weights.setOnes(basis->size(), 1);
    }
    
    /// Construct a rational counterpart of basis
    gsRationalBasis(const SrcT & basis)
        : m_src(basis.clone() )
    { 
        m_weights.setOnes(basis.size(), 1);
    }
    
    /// Construct a rational counterpart of basis
    gsRationalBasis(SrcT * basis, const gsMatrix<T> & w)
        : m_src (basis), m_weights(w)
    { 
        GISMO_ASSERT(m_weights.rows() == m_src->size(), "Invalid basis/weights");
    }
    
    /// Construct a rational counterpart of basis
    gsRationalBasis(SrcT * basis, gsMovable< gsMatrix<T> > w)
        : m_src(basis), m_weights(w)
    { 
        GISMO_ASSERT(m_weights.rows() == m_src->size(), "Invalid basis/weights");
    }
    
    /// Copy Constructor
    gsRationalBasis(const gsRationalBasis & o) : Base(o)
    { 
        m_src= o.source().clone();
        m_weights = o.weights()  ; 
    }
    
    /// Assignment operator
    gsRationalBasis& operator=( const gsRationalBasis & o)
    {
        if ( this == &o )
            return *this;
        return gsRationalBasis(o);
    }
    
    // Destructor
    ~gsRationalBasis()
    { 
        delete m_src; 
    }
    
    /// Check the rational basis for consistency
    bool check() const
    { 
        return ( 
                m_weights.size() == m_src->size()
               ); 
    }
    
    gsBasis<T> * makeNonRational() const
    { return m_src->clone(); }
    
public:
    
///////////////////////////////////////////////////
// Virtual member functions overriding source basis
///////////////////////////////////////////////////

    int dim() const { return Dim; };
    
    int size() const { return m_src->size(); }

    int size(int const& k) const{ return m_src->size(k); }

    int numElements() const { return m_src->numElements(); }
    
    void active_into(const gsMatrix<T> & u, gsMatrix<unsigned>& result) const 
    { m_src->active_into(u, result); }
    
    gsMatrix<unsigned> * boundary( ) const {return m_src->boundary(); }
    
    gsMatrix<unsigned> * boundary(boxSide const & s ) const 
    {return m_src->boundary(s); };
    
    // Look at gsBasis class for a description
    int degree(int i = 0) const {return m_src->degree(i); }

    // Look at gsBasis class for a description
    int maxDegree()   const   {return m_src->maxDegree(); }

    // Look at gsBasis class for a description
    int minDegree()   const    {return m_src->minDegree(); }

    // Look at gsBasis class for a description
    int totalDegree() const     {return m_src->totalDegree(); }

    void uniformRefine(int numKnots = 1)
    {
        m_src->uniformRefine_withCoefs(m_weights, numKnots);
    }
    
    void uniformRefine_withCoefs(gsMatrix<T>& coefs, int numKnots = 1)
    {
        assert( coefs.rows() == this->size() && m_weights.rows() == this->size() );
        gsSparseMatrix<T, RowMajor> transfer;
        m_src->uniformRefine_withTransfer(transfer, numKnots);
      
        for (int i = 0; i < coefs.rows(); ++i) // transform to projective 
            coefs.row(i) *= m_weights(i);
        
        coefs      = transfer * coefs;
        m_weights = transfer * m_weights;
        // Alternative way
        // gsBasis<T> * tmp = m_src->clone();
        // tmp->uniformRefine_withCoefs(coefs, numKnots);
        // delete tmp;
        // m_src->uniformRefine_withCoefs(m_weights, numKnots);
        
        for (int i = 0; i < coefs.rows(); ++i) // back to affine coefs
            coefs.row(i) /= (m_weights)(i);
    }
    
    
    void uniformRefine_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, int numKnots = 1)
    {
        assert( m_weights.rows() == this->size() );
        
        // 1. Get source transfer matrix (while refining m_src)
        m_src->uniformRefine_withTransfer(transfer, numKnots);
        
        // 2. Compute rational basis transfer matrix
        // To be applied on affine coefficients, as usual.
        // Transfer matrix for rational bases. Formula:
        //
        //   ( (T*m_weights)'*I )^{-1} * T*W 
        //
        // Where ' denotes transpose, W is a diagonal matrix with
        // diagonal=m_weights, I is the identity and T the source
        // transfer matrix.
        // i.e. apply weight transform ( to compute weighted coefs),
        // apply T, return to affine by the inverse weight transform
        // (T*W)'*I )^{-1}.
        // In Eigen this could be something like
        // (transfer * (*m_weights) ).asDiagonal().inverse() * transfer * m_weights.asDiagonal() ; 
        // but that has troubles with the sparse/diagonal expressions etc.
        // So we do it by using a temporary.
        
        gsVector<T> tmp = m_weights ;
        m_weights.noalias() = transfer * tmp;  // Refine the weights as well
        
        transfer =  m_weights.cwiseInverse().asDiagonal() * transfer *  tmp.asDiagonal(); 
    }


    void degreeElevate(int const& i = 1) 
    {
        GISMO_NO_IMPLEMENTATION
    }
    
    virtual gsBasis<T>& component(unsigned i) const        { return m_src->component(i); }
    
    gsBasis<T> * boundaryBasis(boxSide const & s ) const   
    { 
        typename SrcT::BoundaryBasisType * bb = m_src->boundaryBasis(s);
        gsMatrix<unsigned> * ind = m_src->boundary(s);
      
        gsMatrix<T> ww( ind->size(),1);
        for ( index_t i=0; i<ind->size(); ++i)
            ww(i,0) = m_weights( (*ind)(i,0), 0);
        
        delete ind;
        return new BoundaryBasisType(bb, give(ww));
    }

    gsDomain<T> * domain() const { return m_src->domain(); }

    void anchors_into(gsMatrix<T> & result) const { return m_src->anchors_into(result) ;}

    // Look at gsBasis class for documentation 
    void connectivity(const gsMatrix<T> & nodes, gsMesh<T> & mesh) const
    { return m_src->connectivity(nodes, mesh); }
    
    gsMatrix<T> support() const {return m_src->support(); }
    
    gsMatrix<T> support(const unsigned & i) const {return m_src->support(i); }
    
    void eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const;
    
    void evalSingle_into(unsigned i, const gsMatrix<T> & u, gsMatrix<T>& result) const ;

    void eval_into(const gsMatrix<T> & u, const gsMatrix<T> & coefs, gsMatrix<T>& result) const ;

    void evalAllDers_into(const gsMatrix<T> & u, int n, gsMatrix<T>& result) const;
    
    void deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const ;
    
    void deriv2_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const;
    
//////////////////////////////////////////////////
// Additional members for rational bases
//////////////////////////////////////////////////

    /// Returns the source basis of the rational basis
    const SrcT & source () const
    { return *m_src; }
    
    SrcT & source ()
    { return *m_src; }
    
    /// Returns the weights of the rational basis
    const gsMatrix<T> & weights() const { return m_weights; }
    
    /// Returns the weights of the rational basis
    gsMatrix<T> & weights()  { return m_weights; }
    
    /// Access to i-th weight
    T & weight(int i)             { return m_weights(i); }
    const T & weight(int i) const { return m_weights(i); }
    
    /// Set weights
    void setWeights(gsMatrix<T> const & w) 
    { 
        GISMO_ASSERT( w.cols() == 1, "Weights should be scalars" ) ;
        m_weights = w;
    }
    
    typename gsBasis<T>::domainIter makeDomainIterator() const
    {
        gsWarn<< "rational domain iterator with evaluate the source.\n";
        return m_src->makeDomainIterator();
    }

    typename gsBasis<T>::domainIter makeDomainIterator(const boxSide & s) const
    {
        gsWarn<< "rational domain boundary iterator with evaluate the source.\n";
        return m_src->makeDomainIterator(s);
    }
    
// Data members
protected:
    
    // Source basis
    SrcT * m_src;
    
    // Weight vector of size: size() x 1
    gsMatrix<T> m_weights;
    
}; // class gsRationalBasis
    
    
template<class SrcT>
void gsRationalBasis<SrcT>::evalSingle_into(unsigned i, const gsMatrix<T> & u, gsMatrix<T>& result) const
{ 
    m_src->evalSingle_into(i, u, result);  
    result.array() *= m_weights.at(i,0);
    gsMatrix<T> denom;
    m_src->eval_into(u, m_weights, denom); 
    
    result.array() /= denom.array();
}
    
    
template<class SrcT>
void gsRationalBasis<SrcT>::eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const
{ 
    m_src->eval_into(u, result);
    const gsMatrix<unsigned> act = m_src->active(u);
    
    gsMatrix<T> denom;
    m_src->eval_into(u, m_weights, denom); 
    
    for ( index_t j=0; j< act.cols(); ++j)
    {
        result.col(j) /= denom(j);
        for ( index_t i=0; i< act.rows(); ++i)    
            result(i,j) *=  m_weights( act(i,j) ) ;
    }
}
  
  
/// For non-specialized version (e.g. tensor product)
template<class SrcT>
void gsRationalBasis<SrcT>::eval_into(const gsMatrix<T> & u, const gsMatrix<T> & coefs, gsMatrix<T>& result) const
{ 
    assert( coefs.rows() == m_weights.rows() ) ;
  
    // Compute projective coefficients
    gsMatrix<T> tmp = m_weights.asDiagonal() * coefs;
  
    gsMatrix<T> denom;
  
    // Evaluate the projective numerator and the denominator
    m_src->eval_into(u, tmp, result);
    m_src->eval_into(u, m_weights, denom);
  
    // Divide numerator by denominator
    for ( index_t j=0; j < u.cols() ; j++ ) // for all points (columns of u)
        result.col(j) /= denom(0,j) ;
}
    
    
template<class SrcT>
void gsRationalBasis<SrcT>::evalAllDers_into(const gsMatrix<T> & u, int n, gsMatrix<T>& result) const
{
    if (n == 0)
    {
        eval_into(u, result);
        return;
    }
    else if (n > 1)
    {
        GISMO_ERROR("gsRationalBasis::evalAllDers_into not implemented for n > 1");
        return;
    }
    assert( n == 1 );
    static const int d = Dim;

    // find active basis functions
    gsMatrix<unsigned> act;
    m_src->active_into(u,act);
    const int numAct = act.rows();

    // initialize result matrix with values and derivatives of source
    // basis
    m_src->evalAllDers_into(u, n, result);
    typename gsMatrix<T>::Block ev  = result.topRows(numAct);
    typename gsMatrix<T>::Block der = result.bottomRows(result.rows() - numAct);

    // evaluate weights and their derivatives
    gsMatrix<T> Wval, Wder;
    m_src->eval_into(u,  m_weights, Wval);
    m_src->deriv_into(u, m_weights, Wder);

    for (index_t i = 0; i < u.cols(); ++i)// for all parametric points
    {
        const T & Wval_i = Wval(0,i);
        // compute derivatives first since they depend on the function
        // values of the source basis
        der.col(i) *= Wval_i;
      
        for (index_t k = 0; k < numAct; ++k)
        {
            der.template block<d,1>(k*d,i).noalias() -=  
                ev(k,i) * Wder.block(0,i*d,1,d).transpose();
            
            der.template block<d,1>(k*d,i) *=
                m_weights( act(k,i), 0 ) / (Wval_i*Wval_i);
        }
        
        // compute function values
        ev.col(i) /= Wval_i;
      
        for (index_t k = 0; k < numAct; ++k)    
            ev(k,i) *= m_weights( act(k,i), 0 ) ;
    }
}
    

template<class SrcT>
void gsRationalBasis<SrcT>::deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result) const
{ 
  // Formula:
  // R_i' = (w_i N_i / W)' = w_i ( N_i'W - N_i W' ) / W^2 

  static const int d = Dim;
  //const index_t numPts = u.cols();// at how many points to evaluate the gradients
  gsMatrix<T> ev, Wval, Wder;
  
  m_src->deriv_into(u, result         );
  m_src->eval_into (u, ev             ); 
  m_src->eval_into (u, m_weights, Wval); 
  m_src->deriv_into(u, m_weights, Wder); 

  //std::cout<<"Src Deriv:\n"<< result <<"\n";
  //std::cout<<"Weights:\n"<< *m_weights <<"\n";
  //std::cout<<"Wval:\n"<< Wval <<"\n";
  //std::cout<<"Wder:\n"<< Wder <<"\n";

  gsMatrix<unsigned> act;
  m_src->active_into(u,act);
  
  for ( index_t i=0; i!= result.cols(); ++i )// for all parametric points
  {  
    result.col(i).array()  *=  Wval(0,i) ;  // N_i'W
    for ( index_t k=0; k<ev.rows() ; ++k ) // for all basis funct. values
    {
      result.template block<d,1>(k*d,i).noalias() -=  ev(k,i)  * // - N_i W'
	Wder.block(0,i*d,1,d).transpose();
      result.template block<d,1>(k*d,i) *= m_weights.at( act(k,i), 0) ;
    }
    result.col(i) /= Wval(0,i) * Wval(0,i);
  } 
};


// TODO: is this version faster? try it for Dim=1
//
//~ /// For specialized version (e.g. gsBSplineBasis )
//~ template<class T, template<class T> class source_basis>
//~ gsMatrix<T> * gsRationalBasis<SrcT>::deriv(const gsMatrix<T> & u ) const
//~ { 
//~     gsMatrix<T> * res   = m_src->deriv(u); 
//~     gsMatrix<T> * ev    = m_src->eval(u); 
//~     gsMatrix<T> * Wval  = m_src->eval(u,m_weights); 
//~ 
//~     // TO DO Should not use this, if we want to avoid infinite loops !!!
//~     gsMatrix<T> * Wder  = m_src->deriv(u,m_weights); 
//~     gsMatrix<unsigned> * act   = m_src->active(u);
//~ 
//~     for ( index_t i=0; i!= res->cols(); ++i )
//~     {
//~         res->col(i)  *= Wval->at(0,i) ;
//~         res->col(i)  -= ev->col(i)  * Wder->at(0,i) ;
//~ 
//~         for ( index_t j=0; j!= res->rows(); ++j )
//~         {
//~             res->at(j,i)  *= m_weights.at( act->at(j,i), 0) ;
//~             res->at(j,i)  /= Wval->at(0,i) * Wval->at(0,i);
//~         }
//~         
//~     }
//~ 
//~     delete ev;
//~     delete act;
//~     delete Wval;
//~     delete Wder;
//~     
//~     return res;
//~ };


template<class SrcT>
void gsRationalBasis<SrcT>::deriv2_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const
{   
    static const int d = Dim;

    GISMO_ASSERT(d==1,"gsRationalBasis::deriv2_into is currently not working correctly for dim>1.");

    gsMatrix<T> ev, Wval;

    //std::cout<<"TP rational deriv \n";
    m_src->deriv2_into(u,result); 
    m_src->eval_into(u, ev); 
    m_src->eval_into(u, m_weights, Wval); 

    gsMatrix<T> Wder;
    m_src->deriv_into(u, m_weights, Wder); 
    gsMatrix<unsigned> act;
    m_src->active_into(u,act);

    for ( index_t i=0; i!= result.cols(); ++i )
    {
        result.col(i).array()  *=  Wval(0,i) ;

        for ( index_t k=0; k < ev.rows() ; ++k )
        {
	  result.template block<d,1>(k*d,i).array() -=  ev(k,i)  * Wder(0,i) ;// TO DO: correct as deriv -- this is currently not right
            result.template block<d,1>(k*d,i) *= m_weights( act(k,i), 0) ;
            result.template block<d,1>(k*d,i) /= Wval(0,i) * Wval(0,i);
        }
    }
}


}; // namespace gismo
