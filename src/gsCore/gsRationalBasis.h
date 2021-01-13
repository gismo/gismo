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
#include <gsCore/gsBoundary.h>


namespace gismo
{

/** \brief
    Class that creates a rational counterpart for a given basis.

    A rational basis holds an inner (referred to as "source") basis of
    the given source type, and a matrix of coefficients defining a
    weight function in terms of the source basis.

    If \em w_i is the i-th weight coefficient and \em b_i the i-th basis
    function of the source basis, then the i-th basis function of the
    resulting rational basis is given by
    r_i(u) = b_i(u) w_i / (sum(b_j(u) w_j, j = 1, ..., size))

    If the weights are all equal to one (or all equal to a constant), the
    rational basis is identical (or a scalar multiple) to the source basis.

    Example: the rational version of a B-spline basis is a NURBS basis.

    \tparam SrcT the basis of which to create a rational version

    \ingroup basis
    \ingroup Core
*/

template<class SrcT>
class gsRationalBasis : public gsBasis< typename SrcT::Scalar_t >
{
public:
    typedef memory::shared_ptr< gsRationalBasis > Ptr;
    typedef memory::unique_ptr< gsRationalBasis > uPtr;

    typedef typename SrcT::Scalar_t Scalar_t;
    typedef Scalar_t T;

    // Obtain dimension from the source basis
    static const int Dim = SrcT::Dim;

    static const bool IsRational = true;

    typedef gsBasis<T> Base;

    /// Associated source basis type
    typedef SrcT SourceBasis;

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
    : m_src(basis.clone().release())
    {
        m_weights.setOnes(basis.size(), 1);
    }

    /// Construct a rational counterpart of basis
    gsRationalBasis(SrcT * basis, gsMatrix<T> w)
    : m_src (basis), m_weights(give(w))
    {
        GISMO_ASSERT(m_weights.rows() == m_src->size(),
                     "Invalid basis/weights ("<<m_weights.rows()<<"/"<<m_src->size());
    }

    /// Copy Constructor
    gsRationalBasis(const gsRationalBasis & o) : Base(o)
    {
        m_src= o.source().clone().release();
        m_weights = o.weights()  ;
    }

    /// Assignment operator
    gsRationalBasis& operator=( const gsRationalBasis & o)
    {
        if ( this != &o )
        {
            delete m_src;
            m_src = o.source().clone();
            m_weights = o.weights()  ;
        }
        return *this;
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

    memory::unique_ptr<gsBasis<T> > makeNonRational() const
    { return m_src->clone(); }

public:

// ***********************************************
// Virtual member functions overriding source basis
// ***********************************************

    short_t domainDim() const { return Dim; }

    index_t size() const { return m_src->size(); }

    int size(int const& k) const{ return m_src->size(k); }

    size_t numElements() const { return m_src->numElements(); }
    using Base::numElements; //unhide

    void active_into(const gsMatrix<T> & u, gsMatrix<index_t>& result) const
    { m_src->active_into(u, result); }
    
    virtual const gsBasis<T> & component(short_t i) const { return m_src->component(i); }
    using Base::component;

    gsMatrix<index_t> allBoundary( ) const {return m_src->allBoundary(); }
    
    gsMatrix<index_t> boundaryOffset(boxSide const & s, index_t offset ) const
    { return m_src->boundaryOffset(s,offset); }

    virtual index_t functionAtCorner(boxCorner const & c) const
    { return m_src->functionAtCorner(c); }

    // Look at gsBasis class for a description
    short_t degree(short_t i = 0) const {return m_src->degree(i); }

    // Look at gsBasis class for a description
    short_t maxDegree()   const   {return m_src->maxDegree(); }

    // Look at gsBasis class for a description
    short_t minDegree()   const    {return m_src->minDegree(); }

    // Look at gsBasis class for a description
    short_t totalDegree() const     {return m_src->totalDegree(); }

    void uniformRefine(int numKnots = 1, int mul=1)
    {
        m_src->uniformRefine_withCoefs(m_weights, numKnots, mul);
    }

    void uniformRefine_withCoefs(gsMatrix<T>& coefs, int numKnots = 1,  int mul=1);

    void uniformRefine_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, int numKnots = 1, int mul=1);

    /**
     * @brief Refines specified areas or boxes, depending on underlying basis.
     *
     * @param boxes See the function gsBasis::refineElements() of the underlying
     * basis for syntax.
     */
    void refineElements( std::vector<index_t> const & boxes)
    {
        // call the refineElements_withCoefs-function of the underlying
        // basis, where the weights are used as coefficients
        m_src->refineElements_withCoefs( m_weights, boxes );
    }

    /**
     * @brief Refines specified areas or boxes, depending on underlying basis.
     *
     * @param coefs Coefficients, given as gsMatrix of size \f$ n \times d\f$,
     * where \f$n\f$ is the number of basis functions and \f$d\f$ is the target
     * dimension.
     * @param boxes See the function gsBasis::refineElements() of the underlying
     * basis for syntax.
     */
    void refineElements_withCoefs(gsMatrix<T> & coefs,std::vector<index_t> const & boxes);

    void degreeElevate(short_t const& i = 1, short_t const dir = -1)
    {
        typename SourceBasis::GeometryType tmp(*m_src,give(m_weights));
        tmp.degreeElevate(i,dir);
        tmp.coefs().swap(m_weights);
        std::swap(*m_src, tmp.basis() );
    }

    void degreeReduce(short_t const& i = 1, short_t const dir = -1)
    {
        typename SourceBasis::GeometryType tmp(*m_src, give(m_weights));
        tmp.degreeReduce(i,dir);
        tmp.coefs().swap(m_weights);
        std::swap(*m_src, tmp.basis() );
    }

    /* if ever be reused, change to actual and current GISMO_UPTR_FUNCTION stuff und uPtr
      GISMO_UPTR_FUNCTION_DEF(gsBasis<T>, boundaryBasis, boxSide const &)
      {
      typename SrcT::BoundaryBasisType * bb = m_src->boundaryBasis(s);
      gsMatrix<index_t> ind = m_src->boundary(s);
      gsMatrix<T> ww( ind.size(),1);
      for ( index_t i=0; i<ind.size(); ++i)
      ww(i,0) = m_weights( (*ind)(i,0), 0);

      return new BoundaryBasisType(bb, give(ww));
      return 0;
      }
    */

    gsDomain<T> * domain() const { return m_src->domain(); }

    void anchors_into(gsMatrix<T> & result) const
    { return m_src->anchors_into(result); }

    void anchor_into(index_t i, gsMatrix<T> & result) const
    { return m_src->anchor_into(i,result); }

    // Look at gsBasis class for documentation
    void connectivity(const gsMatrix<T> & nodes, gsMesh<T> & mesh) const
    { return m_src->connectivity(nodes, mesh); }

    gsMatrix<T> support() const {return m_src->support(); }
    
    gsMatrix<T> support(const index_t & i) const {return m_src->support(i); }
    
    void eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const;
    
    void evalSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result) const ;

    void evalFunc_into(const gsMatrix<T> & u, const gsMatrix<T> & coefs, gsMatrix<T>& result) const;

    //void evalAllDers_into(const gsMatrix<T> & u, int n,
    //                      std::vector<gsMatrix<T> >& result) const;

    void deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const ;

    void deriv2_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const;

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

    /// Const access to i-th weight
    const T & weight(int i) const { return m_weights(i); }

    /// Set weights
    void setWeights(gsMatrix<T> const & w)
    {
        GISMO_ASSERT( w.cols() == 1, "Weights should be scalars" ) ;
        m_weights = w;
    }

    virtual void matchWith(const boundaryInterface & bi, const gsBasis<T> & other,
                           gsMatrix<index_t> & bndThis, gsMatrix<index_t> & bndOther) const
    { 
        if ( const gsRationalBasis * _other = dynamic_cast<const gsRationalBasis*>(&other) )
            m_src->matchWith(bi,*_other->m_src,bndThis,bndOther);
        else
            m_src->matchWith(bi,other,bndThis,bndOther);
    }

    /// Returns a matrix of projective coefficients. The input \a
    /// coefs are affine coefficients for this basis
    gsMatrix<T> projectiveCoefs(const gsMatrix<T> & coefs) const
    { return projectiveCoefs(coefs, m_weights); }

    /// Returns a matrix of projective coefficients. The input \a
    /// coefs are affine coefficients and weights
    static gsMatrix<T> projectiveCoefs(const gsMatrix<T> & coefs, const gsMatrix<T> & weights)
    {
        GISMO_ASSERT(coefs.rows() == weights.rows(),
                     "Invalid basis/coefficients ("<<coefs.rows()<<"/"<<weights.rows());
        const index_t n = coefs.cols();
        gsMatrix<T> rvo(coefs.rows(), n + 1);
        // switch from control points (n-dimensional) to
        // "projective control points" ((n+1)-dimensional),
        // where the last coordinate is the weight.
        rvo.leftCols(n).noalias() = weights.asDiagonal() * coefs;
        rvo.col(n)                = weights;
        return rvo;
    }

    /// Sets the weights and the \a coefs to be the affine
    /// coefficients corresponding to the projective coefficients
    /// \a pr_coefs
    static void setFromProjectiveCoefs(const gsMatrix<T> & pr_coefs,
                                       gsMatrix<T> & coefs, gsMatrix<T> & weights)
    {
        const index_t n = pr_coefs.cols() - 1;
        weights = pr_coefs.col( n );
        coefs   = pr_coefs.leftCols(n).array().colwise() / weights.col(0).array();
        // equiv: coefs = pr_coefs.leftCols(n).array() / weights.replicate(1,n).array();
    }


    typename gsBasis<T>::domainIter makeDomainIterator() const
    {
//        gsWarn<< "rational domain iterator with evaluate the source.\n";
        return m_src->makeDomainIterator();
    }

    typename gsBasis<T>::domainIter makeDomainIterator(const boxSide & s) const
    {
//        gsWarn<< "rational domain boundary iterator with evaluate the source.\n";
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
void gsRationalBasis<SrcT>::evalSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result) const
{ 
    m_src->evalSingle_into(i, u, result);  
    result.array() *= m_weights.at(i);
    gsMatrix<T> denom;
    m_src->evalFunc_into(u, m_weights, denom);

    result.array() /= denom.array();
}


template<class SrcT>
void gsRationalBasis<SrcT>::eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const
{
    m_src->eval_into(u, result);
    const gsMatrix<index_t> act = m_src->active(u);

    gsMatrix<T> denom;
    m_src->evalFunc_into(u, m_weights, denom);

    for ( index_t j=0; j< act.cols(); ++j)
    {
        result.col(j) /= denom(j);
        for ( index_t i=0; i< act.rows(); ++i)
            result(i,j) *=  m_weights( act(i,j) ) ;
    }
}


// For non-specialized version (e.g. tensor product)
template<class SrcT>
void gsRationalBasis<SrcT>::evalFunc_into(const gsMatrix<T> & u, const gsMatrix<T> & coefs, gsMatrix<T>& result) const
{
    assert( coefs.rows() == m_weights.rows() ) ;

    // Compute projective coefficients
    const gsMatrix<T> tmp = m_weights.asDiagonal() * coefs;

    gsMatrix<T> denom;

    // Evaluate the projective numerator and the denominator
    m_src->evalFunc_into( u, tmp, result);
    m_src->evalFunc_into( u, m_weights, denom);

    // Divide numerator by denominator
    result.array().rowwise() /= denom.row(0).array();
    // equivalent:
    //for ( index_t j=0; j < u.cols(); j++ ) // for all points (columns of u)
    //    result.col(j) /= denom(0,j);
}


/* TODO
template<class SrcT>
void gsRationalBasis<SrcT>::evalAllDers_into(const gsMatrix<T> & u, int n,
                                             std::vector<gsMatrix<T> >& result) const
{
    result.resize(n+1);

    std::vector<gsMatrix<T> > ev(n+1);

    m_src->evalAllDers_into(u, n, ev);

    // find active basis functions
    gsMatrix<index_t> act;
    m_src->active_into(u,act);

    const int numAct = act.rows();

    // evaluate weights and their derivatives
    gsMatrix<T> Wval, Wder;
    m_src->evalFunc_into (u,  m_weights, Wval);
    m_src->derivFunc_into(u, m_weights, Wder);

    for (index_t i = 0; i < u.cols(); ++i)// for all parametric points
    {
        const T & Wval_i = Wval(0,i);
        // compute derivatives first since they depend on the function
        // values of the source basis
        der.col(i) *= Wval_i;

        for (index_t k = 0; k < numAct; ++k)
        {
            der.template block<Dim,1>(k*Dim,i).noalias() -=
                ev(k,i) * Wder.col(i);

            der.template block<Dim,1>(k*Dim,i) *=
                m_weights( act(k,i), 0 ) / (Wval_i*Wval_i);
        }

        // compute function values
        ev.col(i) /= Wval_i;

        for (index_t k = 0; k < numAct; ++k)
            ev(k,i) *= m_weights( act(k,i), 0 ) ;
    }
}
//*/

template<class SrcT>
void gsRationalBasis<SrcT>::deriv_into(const gsMatrix<T> & u,
                                       gsMatrix<T>& result) const
{
    // Formula:
    // R_i' = (w_i N_i / W)' = w_i ( N_i'W - N_i W' ) / W^2

    gsMatrix<index_t> act;
    m_src->active_into(u,act);

    result.resize( act.rows()*Dim, u.cols() );

    std::vector<gsMatrix<T> > ev(2);
    m_src->evalAllDers_into(u, 1, ev);

    T W;
    gsMatrix<T> dW(Dim,1);

    for ( index_t i = 0; i!= result.cols(); ++i )// for all parametric points
    {
        // Start weight function computation
        W = 0.0;
        dW.setZero();
        for ( index_t k = 0; k != act.rows(); ++k ) // for all basis functions
        {
            const T curw = m_weights.at(act(k,i));
            W  += curw * ev[0](k,i);
            dW += curw * ev[1].template block<Dim,1>(k*Dim,i);
        }
        // End weight function computation

        result.col(i) = W * ev[1].col(i);  // N_i'W

        for ( index_t k = 0; k != act.rows() ; ++k) // for all basis functions
        {
            const index_t kd = k * Dim;

            result.template block<Dim,1>(kd,i).noalias() -=
                ev[0](k,i) * dW; // - N_i W'

            result.template block<Dim,1>(kd,i) *= m_weights.at( act(k,i) ) / (W * W);
        }
    }
}

template<class SrcT>
void gsRationalBasis<SrcT>::deriv2_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const
{
    // Formula:
    // ( W^2 / w_k) * R_k'' = ( N_k'' W - N_k W'' ) - 2 N_k' W' + 2 N_k (W')^2 / W
    // ( W^2 / w_k) * d_ud_vR_k = ( d_ud_vN_k W - N_k d_ud_vW )
    //                          - d_uN_k d_vW - d_vN_k d_uW + 2 N_k d_uW d_vW / W

    static const int str = Dim * (Dim+1) / 2;

    gsMatrix<index_t> act;
    m_src->active_into(u,act);

    result.resize( act.rows()*str, u.cols() );

    T W;
    gsMatrix<T> dW(Dim,1), ddW(str,1);
    std::vector<gsMatrix<T> > ev(3);
    m_src->evalAllDers_into(u, 2, ev);

    for ( index_t i = 0; i!= u.cols(); ++i ) // for all points
    {
        // Start weight function computation
        W = 0.0;
        dW .setZero();
        ddW.setZero();
        for ( index_t k = 0; k != act.rows(); ++k ) // for all basis functions
        {
            //to do with lweights (local weights)
            const T curw = m_weights.at(act(k,i));
            W   += curw * ev[0](k,i);
            dW  += curw * ev[1].template block<Dim,1>(k*Dim,i);
            ddW += curw * ev[2].template block<str,1>(k*str,i);
        }
        // End weight function computation

        result.col(i) = W * ev[2].col(i); // N_k'' W

        for ( index_t k=0; k != act.rows() ; ++k ) // for all basis functions
        {
            const index_t kstr = k * str;
            const index_t kd   = k * Dim;

            result.template block<str,1>(kstr,i) -=
                ev[0](k,i) * ddW; // - N_k * W''

            result.template block<Dim,1>(kstr,i) +=
                // - 2 N_k' W' + 2 N_k (W')^2 / W
                ( 2 * ev[0](k,i) / W ) * dW.cwiseProduct(dW)
                - 2 * ev[1].template block<Dim,1>(kd,i).cwiseProduct(dW);

            int m = Dim;
            for ( int _u=0; _u != Dim; ++_u ) // for all mixed derivatives
                for ( int _v=_u+1; _v != Dim; ++_v )
                {
                    result(kstr + m++, i) +=
                        - ev[1](kd+_u,i) * dW.at(_v) // - du N_k * dv W
                        - ev[1](kd+_v,i) * dW.at(_u) // - dv N_k * du W
                        // + 2 * N_k * du W * dv W / W
                        + 2 * ev[0](k,i) * dW.at(_u) * dW.at(_v) / W;
                }

            result.template block<str,1>(kstr,i) *=
                m_weights.at( act(k,i) ) / (W*W); // * (w_k / W^2)
        }
    }
}


template<class SrcT>
void gsRationalBasis<SrcT>::uniformRefine_withCoefs(gsMatrix<T>& coefs, int numKnots,  int mul)
{
    GISMO_ASSERT( coefs.rows() == this->size() && m_weights.rows() == this->size(),
                  "Invalid dimensions" );
    gsSparseMatrix<T, RowMajor> transfer;
    m_src->uniformRefine_withTransfer(transfer, numKnots, mul);

    coefs     = transfer * ( m_weights.asDiagonal() * coefs);
    m_weights = transfer * m_weights;
    // Alternative way
    // gsBasis<T> * tmp = m_src->clone();
    // tmp->uniformRefine_withCoefs(coefs, numKnots);
    // delete tmp;
    // m_src->uniformRefine_withCoefs(m_weights, numKnots);

    // back to affine coefs
    coefs.array().colwise() /= m_weights.col(0).array();
    // equiv:
    // for (int i = 0; i < coefs.rows(); ++i)
    //    coefs.row(i) /= m_weights.at(i);
}

template<class SrcT>
void gsRationalBasis<SrcT>::uniformRefine_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, int numKnots, int mul)
{
    // 1. Get source transfer matrix (while refining m_src)
    m_src->uniformRefine_withTransfer(transfer, numKnots, mul);

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
    // (transfer * m_weights).asDiagonal().inverse() * transfer * m_weights.asDiagonal() ;
    // but that has troubles with the sparse/diagonal expressions etc.
    // So we do it by using a temporary.

    const gsVector<T> tmp = m_weights ;
    m_weights.noalias() = transfer * tmp;  // Refine the weights as well

    transfer =  m_weights.cwiseInverse().asDiagonal() * transfer *  tmp.asDiagonal();
}


template<class SrcT>
void gsRationalBasis<SrcT>::refineElements_withCoefs(gsMatrix<T> & coefs,
                                                     std::vector<index_t> const & boxes)
{
    // switch from control points (n-dimensional) to
    // "projective control points" ((n+1)-dimensional),
    // where the last coordinate is the weight.
    gsMatrix<T> rw = projectiveCoefs(coefs, m_weights);

    // refine with these projective control points as coefficients.
    m_src->refineElements_withCoefs( rw, boxes );

    // Regain the new n-dimensional control points and the new
    // weights.
    setFromProjectiveCoefs(rw, coefs, m_weights);
}


} // namespace gismo
