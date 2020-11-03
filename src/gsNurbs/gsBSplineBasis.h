/** @file gsBSplineBasis.h

    @brief Provides declaration of BSplineBasis class

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris, D. Mokris
*/


#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsConstantBasis.h>

#include <gsTensor/gsTensorBasis.h>
#include <gsTensor/gsTensorDomainIterator.h>
#include <gsTensor/gsTensorDomainBoundaryIterator.h>

#include <gsNurbs/gsKnotVector.h>

namespace gismo
{

/// @brief Traits for BSplineBasis in more dimensions
template<short_t d, class T>
struct gsBSplineTraits
{
    typedef gsKnotVector<T> KnotVectorType;

    typedef gsTensorBSplineBasis<d,T> Basis;
    typedef gsTensorNurbsBasis<d,T>   RatBasis;
    typedef gsTensorBSpline<d,T>      Geometry;
    typedef gsTensorNurbs<d,T>        RatGeometry;
};
template<class T>
struct gsBSplineTraits<1,T>
{
    typedef gsKnotVector<T> KnotVectorType;
    typedef gsBSplineBasis<T>         Basis;
    typedef gsNurbsBasis<T>           RatBasis;
    typedef gsBSpline<T>              Geometry;
    typedef gsNurbs<T>                RatGeometry;
};
template<class T>
struct gsBSplineTraits<0,T>
{
    typedef gsKnotVector<T> KnotVectorType;
    typedef gsConstantBasis<T>                       Basis;
    typedef gsConstantBasis<T>                       RatBasis;
    typedef gsConstantFunction<T>                    Geometry;
    typedef gsConstantFunction<T>                    RatGeometry;
};

/** \brief
    A univariate B-spline basis.

    \tparam T coefficient type
    \tparam KnotVectorType the type of knot vector to use

    \ingroup basis
    \ingroup Nurbs
*/
template<class T>
class gsTensorBSplineBasis<1,T> : public gsTensorBasis<1,T>
{
public:
    typedef gsKnotVector<T> KnotVectorType;

    typedef gsTensorBasis<1,T> Base;

    typedef gsBSplineBasis<T> Self_t;

    typedef gsTensorBSplineBasis TensorSelf_t;
    typedef gsTensorBSplineBasis<1,T> TensorSelf_tt;

    /// @brief Coefficient type
    typedef T Scalar_t;

    /// @brief Associated geometry type
    typedef typename gsBSplineTraits<1,T>::Geometry GeometryType;

    /// @brief Associated Boundary basis type
    typedef typename gsBSplineTraits<0,T>::Basis BoundaryBasisType;

    /// @brief Dimension of the parameter domain
    static const short_t Dim = 1;

    /// @brief Smart pointer for gsTensorBSplineBasis
    typedef memory::shared_ptr< Self_t > Ptr;

    /// @brief Smart pointer for gsTensorBSplineBasis
    typedef memory::unique_ptr< Self_t > uPtr;
  
public:

    // Look at gsBasis class for a description
    // Note: Specializing pointer type at return
    //GISMO_UPTR_FUNCTION_PURE(TensorSelf_tt, clone)
    private: virtual gsTensorBSplineBasis * clone_impl() const = 0;
    public: uPtr clone() const { return uPtr(dynamic_cast<Self_t*>(clone_impl())); }
    
    // gsTensorBSplineBasis( const Base & o)
    // { 
    //     const gsTensorBSplineBasis * a;
    //     if ( ( a = dynamic_cast<const gsTensorBSplineBasis *>( &o )) )
    //     {
    //         m_p        = a->degree() ;
    //         m_knots    = KnotVectorType( a->knots() );
    //         m_periodic = a->m_periodic;
    //     }
    //     else
    //         GISMO_ERROR("Cannot convert "<<o<<" to gsTensorBSplineBasis\n");
    // }

    static Self_t * New(std::vector<gsBasis<T>*> & bb );

    static Self_t * New(std::vector<Self_t*> & bb )
    { 
        return new Self_t(*bb.front());
    }

    static uPtr make(std::vector<gsBasis<T>*> & bb )
    {
        return uPtr(New(bb));
    }

    static uPtr make(std::vector<Self_t*> & bb )
    { 
        return uPtr(new Self_t(*bb.front()));
    }

    static uPtr make( const KnotVectorType & KV )
    {
        return uPtr(new Self_t(KV));
    }
    
    // Note: these casts can be dangerous
    // operator Self_t &() { return dynamic_cast<Self_t&>(*this);}
    // operator const Self_t &() const { return dynamic_cast<const Self_t&>(*this);}

public:

    void swap(gsTensorBSplineBasis& other)
    {
        std::swap(m_p, other.m_p);
        std::swap(m_periodic, other.m_periodic);
        m_knots.swap(other.m_knots);
    }

/* Virtual member functions required by the base class */

    // Look at gsBasis class for a description
    short_t domainDim() const { return Dim; }

    // Unhide/forward gsTensorBasis<1,T>::size(k), since the overload
    // with size() automatically hides it in this class
    using Base::size;

    // Look at gsBasis class for a description
    index_t size() const { return m_knots.size() - m_p - 1 - m_periodic; }

    // Look at gsBasis class for a description
    size_t numElements() const { return m_knots.numElements(); }
    using Base::numElements; //unhide

    // Look at gsBasis class for a description
    size_t elementIndex(const gsVector<T> & u ) const;

    // Same as gsBasis::elementIndex but argument is a value instead of a vector
    size_t elementIndex(T u ) const;

    /// @brief Returns span (element) indices of the beginning and end
    /// of the support of the i-th basis function.
    void elementSupport_into(const index_t i, gsMatrix<index_t,1,2>& result) const
    {
        gsMatrix<index_t> tmp_vec;
        m_knots.supportIndex_into(i, tmp_vec);
        result = tmp_vec;
    }

    // Look at gsBasis class for a description
    const TensorSelf_t & component(short_t i) const = 0;

    // Look at gsBasis class for a description
    TensorSelf_t & component(short_t i) = 0;

    /// @brief Returns the anchors (greville points) of the basis
    void anchors_into(gsMatrix<T> & result) const 
    { 
        m_knots.greville_into(result); 
    }

    /// @brief Returns the anchors (greville points) of the basis
    void anchor_into(index_t i, gsMatrix<T> & result) const
    { 
        result.resize(1,1);
        result(0,0) = m_knots.greville(i);
    }

    // Look at gsBasis class for a description
    void connectivity(const gsMatrix<T> & nodes, 
                      gsMesh<T> & mesh) const;

    // Look at gsBasis class for a description
    void active_into(const gsMatrix<T> & u, gsMatrix<index_t>& result) const;

    // Look at gsBasis class for a description
    bool isActive(const index_t i, const gsVector<T> & u) const;

    // Look at gsBasis class for a description
    gsMatrix<index_t> allBoundary( ) const ;

    // Look at gsBasis class for a description
    gsMatrix<index_t> boundaryOffset(boxSide const & s,index_t offset) const;

#ifdef __DOXYGEN__
    /// @brief Gives back the boundary basis at boxSide s
    typename BoundaryBasisType::uPtr boundaryBasis(boxSide const & s);
#endif
    GISMO_UPTR_FUNCTION_DEC(gsConstantBasis<T>, boundaryBasis, boxSide const &)

    // Look at gsBasis class for a description
    gsMatrix<T> support() const ;

    // Look at gsBasis class for a description
    gsMatrix<T> support(const index_t & i ) const ;

    /// @brief Only meaningfull for periodic basis: For basis members that have
    /// a twin, this function returns the other twin index, otherwise it
    /// returns the same index as the argument
    index_t twin(index_t i) const ;

    // Look at gsBasis class for a description
    virtual void eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const;

    // Look at gsBasis class for a description
    virtual void evalSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result) const;

    // Look at gsBasis class for a description
    virtual void evalFunc_into(const gsMatrix<T> & u, const gsMatrix<T> & coefs, gsMatrix<T>& result) const;

    // Look at gsBasis class for a description
    void deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const ;

    // Look at gsBasis class for a description
    void derivSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result ) const ;

    // Look at gsBasis class for a description
    void deriv_into(const gsMatrix<T> & u, const gsMatrix<T> & coefs, gsMatrix<T>& result ) const ;

    // Look at gsBasis class for a description
    void deriv2_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const ;

    // Look at gsBasis class for a description
    void deriv2Single_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result ) const ;

    // Look at gsBasis class for a description
    void deriv2_into(const gsMatrix<T> & u, const gsMatrix<T> & coefs, gsMatrix<T>& result ) const ;

    // Look at gsBasis class for a description
    gsMatrix<T> laplacian(const gsMatrix<T> & u ) const ;

    // Look at gsBasis class for a description
    typename gsBasis<T>::uPtr tensorize(const gsBasis<T> & other) const;
    
    /// @brief Check the BSplineBasis for consistency
    bool check() const
    { 
        if ( m_periodic > 0 )
        {
            // Periodicity check wrt knot values
            return ( 
                m_knots.degree()      == m_p &&
                (int)m_knots.size()   >  2*m_p+1
                );           
        }
        else
        {
            return ( 
                m_knots.degree()      == m_p &&
                (int)m_knots.size()   >  2*m_p+1
                ); 
        }
    }
  
    /// @brief Prints the object as a string.
    std::ostream &print(std::ostream &os) const;

    /// @brief Return a string with detailed information on the basis.
    std::string detail() const;

    // Look at gsBasis class for a description
    virtual void evalDerSingle_into(index_t i, const gsMatrix<T> & u,
                                    int n, gsMatrix<T>& result) const;

    // Look at gsBasis class for a description
    virtual void evalAllDers_into(const gsMatrix<T> & u, int n,
                                  std::vector<gsMatrix<T> >& result) const;

    // Look at gsBasis class for a description
    virtual void evalAllDersSingle_into(index_t i, const gsMatrix<T> & u,
                                        int n, gsMatrix<T>& result) const;

    // Look at gsBasis class for a description
    short_t degree(short_t i) const
    {
        GISMO_UNUSED(i);
        GISMO_ASSERT(i==0,"Asked for degree(i) in 1D basis.");
        return m_p; 
    }

    short_t degree() const {return m_p;}

    // Look at gsBasis class for a description
    short_t maxDegree()   const { return m_p; }

    // Look at gsBasis class for a description
    short_t minDegree()   const { return m_p; }

    // Look at gsBasis class for a description
    short_t totalDegree() const { return m_p; }
 
    /// @brief Returns the order of the B-spline  basis
    inline unsigned order() const { return m_p+1; }

    /// @brief True iff the point \a pp is in the domain of the basis
    inline bool inDomain(T const & pp) const 
    { return ( (pp >= *(m_knots.begin()+m_p)) &&  (pp <= *(m_knots.end()-m_p-1) ) ); }

    /// @brief Returns the starting value of the domain of the basis
    T domainStart() const { return *(m_knots.begin()+m_p); }

    /// @brief Returns the ending value of the domain of the basis
    T domainEnd() const { return *(m_knots.end()-m_p-1); }

    /// @brief Returns length of the ``active" part of the knot vector.
    T _activeLength() const { return domainEnd() - domainStart(); }

    /// @brief Returns the index of the first active (ie. non-zero) basis function at point u
    /// Takes into account non-clamped knots.
    inline index_t firstActive(T u) const {
        return ( inDomain(u) ? (m_knots.iFind(u)-m_knots.begin()) - m_p : 0 );
    }

    // Number of active functions at any point of the domain
    inline index_t numActive() const { return m_p + 1; }

    // Look at gsBasis class for a description
    gsDomain<T> * domain() const { return const_cast<KnotVectorType *>(&m_knots); }

    /// @brief Returns the knot vector of the basis
    const KnotVectorType & knots (int i  = 0) const 
    { 
        GISMO_ENSURE(i==0, "Invalid knots requested");
        return m_knots;
    }
    KnotVectorType & knots (int i  = 0)
    {
        GISMO_ENSURE(i==0, "Invalid knots requested");
        return m_knots;
    }

    T knot(index_t const & i) const { return m_knots[i];}

    /// @brief Inserts the knot \em knot in the underlying knot vector.
    void insertKnot(T knot, int mult=1)
    { m_knots.insert( knot, mult); }

    // compatibility with tensor-bsplines
    void insertKnots(const std::vector< std::vector<T> >& refineKnots)
    {
        GISMO_ASSERT( refineKnots.size() == 1, "refineKnots vector has wrong size" );
        this->knots().insert(refineKnots.front());
    }

    // Look at gsBasis class for a description
    void refineElements(std::vector<index_t> const & elements)
    { m_knots.refineSpans(elements); }

    // Look at gsBasis class for a description
    void uniformRefine(int numKnots = 1, int mul=1)
    { m_knots.uniformRefine(numKnots,mul); }

    // Look at gsBasis class for a description
    void uniformRefine_withCoefs(gsMatrix<T>& coefs, int numKnots = 1, int mul=1);

    // Look at gsBasis class for a description
    void uniformRefine_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, int numKnots = 1, int mul=1);
    
    // Look at gsBasis class for a description
    void uniformCoarsen(int numKnots = 1)
    { m_knots.coarsen(numKnots); }

    // Look at gsBasis class for a description
    void uniformCoarsen_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, int numKnots = 1);

    /// @brief Refine the basis by inserting the given knots and perform knot
    /// refinement for the given coefficient matrix.
    void refine_withCoefs(gsMatrix<T>& coefs, const std::vector<T>& knots);

    // compatibility with tensor-bsplines
    void refine_withCoefs(gsMatrix<T> & coefs, 
                          const std::vector< std::vector<T> >& refineKnots)
    {
        refine_withCoefs(coefs,refineKnots.front() );
    }

    /// @brief Refine the basis by inserting the given knots and produce a sparse matrix which maps coarse coefficient vectors to refined ones.
    void refine_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, const std::vector<T>& knots);

    void refine_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, 
                             const std::vector<std::vector<T> >& knots)
    {
        GISMO_ASSERT(knots.size()==1, "the knots you want to insert do not have the right dimension");
        refine_withTransfer(transfer,knots.front());
        //GISMO_NO_IMPLEMENTATION
    }

    /// @brief Increases the degree without adjusting the smoothness at inner
    /// knots, except from the knot values in \a knots (constrained
    /// knots of initial geometry)
    /// 
    /// This type of refinement is known as k-refinement.  Note that
    /// this type of refinement is ment to be performed after one (or
    /// more) h-refinement steps the parent mesh (\a other),
    /// otherwise this function is equivalent to p-refinement of \a other.
    ///
    /// @param other parent/reference mesh determining the
    /// smoothness at the inner knots.
    /// @param i number of k-refinement steps to perform
    ///
    /// @remarks Not tested yet!
    void refine_k(const TensorSelf_t & other, int const & i = 1)
    { 
        GISMO_ASSERT( m_p >= other.m_p, "Degree of other knot-vector should be lower.");
        //for (typename std::vector<T>::iterator it = 
        //         m_knots.begin(); it != m_knots.end(); ++it)
        //    GISMO_ASSERT( has(*it), "Knot "<< *it<<" is not in the knot vector.");

        // grab unique knots
        const std::vector<T> knots = other.m_knots.unique();

        // Increase the degree without adjusting any knot
        m_p += i; 
        m_knots.set_degree(m_p);
        // Adjust (reduce) smoothness to satisfy initial constraint knots
        m_knots.insert(knots,i);
    }

    /// @brief p-refinement (essentially degree elevation)
    void refine_p(short_t const & i = 1)
    { degreeElevate(i); }

    /// @brief Uniform h-refinement (placing \a i new knots inside each knot-span
    void refine_h(short_t const & i = 1)
    { uniformRefine(i); }
  
    /// @brief Elevate the degree of the basis and preserve the smoothness
    void degreeElevate(short_t const & i = 1, short_t const dir = -1)
    {
        GISMO_UNUSED(dir);
        GISMO_ASSERT( dir == -1 || dir == 0, "Invalid direction");
        m_p+=i; m_knots.degreeElevate(i); 
    }

    // Look at gsBasis for documentation
    void degreeReduce (short_t const & i = 1, short_t const dir = -1)
    {
        GISMO_UNUSED(dir);
        GISMO_ASSERT( dir == -1 || dir == 0, "Invalid direction");
        GISMO_ASSERT( i<=m_p, "Cannot reduce degree to negative");
        m_p-=i; m_knots.degreeReduce(i);
        //m_periodic =
    }

    // Look at gsBasis for documentation
    void degreeIncrease(short_t const & i = 1, short_t const dir = -1)
    {
        GISMO_UNUSED(dir);
        GISMO_ASSERT( dir == -1 || dir == 0, "Invalid direction");
        m_p+=i; m_knots.degreeIncrease(i);
    }

    // Look at gsBasis for documentation
    void degreeDecrease(short_t const & i = 1, short_t const dir = -1)
    {
        GISMO_UNUSED(dir);
        GISMO_ASSERT( dir == -1 || dir == 0, "Invalid direction");
        m_p-=i; m_knots.degreeDecrease(i);
    }

    /// @brief Elevate spline continuity at interior knots by \a i
    void elevateContinuity(int const & i = 1) 
    { 
        GISMO_ASSERT( i>=0 && ( m_knots.size()>static_cast<size_t>(2*(m_p+1)) || i<=m_p ), 
                      "Cannot achieve continuity less than C^{-1} at interior knots.");
        m_knots.reduceMultiplicity(i);
    }

    /// @brief Reduces spline continuity at interior knots by \a i
    void reduceContinuity(int const & i = 1) 
    { 
        GISMO_ASSERT( i>=0 && ( m_knots.size()>static_cast<size_t>(2*(m_p+1)) || i<=m_p ), 
                      "Cannot achieve continuity less than C^{-1} at interior knots.");
        // TODO check: max interior mult + i <= m_p+1
        m_knots.increaseMultiplicity(i);
    }

    /// @brief Tells, whether the basis is periodic.
    bool isPeriodic() const
    {
        return (m_periodic > 0);
    }

    // Look at gsBasis class for a description
    index_t functionAtCorner(boxCorner const & c) const;

    /// @brief Returns number of functions crossing the boundary of the knot vector.
    int numCrossingFunctions () const
    {
        return m_periodic;
    }

    /// @brief Checks, if both endknots have multiplicity m_p + 1.
    bool isClamped() const
    {
        if( m_knots[0] != m_knots[m_p])
            return false;

        else if( m_knots[m_knots.size() - m_p -1] != m_knots[m_knots.size()-1])
            return false;

        else
            return true;
    }

    /// @brief If flag is true, tries to convert the basis to periodic
    /// (succeeds only if the knot vector is suitable).
    void setPeriodic(bool flag = true)
    {
        if ( flag )
            _convertToPeriodic();
        else
            m_periodic = 0;
    }

    // Compatible with tensor B-spline basis
    void setPeriodic(int dir)
    {
        GISMO_UNUSED(dir);
        GISMO_ASSERT(dir==0, "Invalid direction");
            _convertToPeriodic();
    }

    /// @brief Returns the multiplicity of the first ``significant" knot
    /// (i.e., the m_p+1st). If it is different from the multiplicity
    /// of the corresponding knot at the end, returns zero.
    int borderKnotMult() const;

    typename gsBasis<T>::domainIter makeDomainIterator() const
    {
        return typename gsBasis<T>::domainIter(new gsTensorDomainIterator<T,1>(*this));
    }

    typename gsBasis<T>::domainIter makeDomainIterator(const boxSide & s) const
    {
        return ( s == boundary::none ? 
                 typename gsBasis<T>::domainIter(new gsTensorDomainIterator<T,1>(*this)) :
                 typename gsBasis<T>::domainIter(new gsTensorDomainBoundaryIterator<T,1>(*this, s))
                );
    }

    /// @brief Moves the knot vectors to enforce periodicity.
    void enforceOuterKnotsPeriodic();

    // Look at gsBasis class for a description
    void reverse() { m_knots.reverse(); }

    void matchWith(const boundaryInterface & bi,
                   const gsBasis<T> & other,
                   gsMatrix<index_t> & bndThis,
                   gsMatrix<index_t> & bndOther) const;

protected:

    /// @brief Tries to convert the basis into periodic
    void _convertToPeriodic();

    /// @brief Adjusts endknots so that the knot vector can be made periodic.
    void _stretchEndKnots();

public:

    /// @brief Helper function for evaluation with periodic basis.
    ///
    /// @param coefs coefficients (control points, one per row) before
    /// converting the basis into periodic.  
    ///
    /// @return copy of coefs with the first m_periodic rows copied to
    /// the last m_periodic rows.
    gsMatrix<T> perCoefs( const gsMatrix<T>& coefs ) const
    {
        gsMatrix<T> per_coefs = coefs;
        per_coefs.bottomRows( m_periodic ) = coefs.topRows( m_periodic );
        return per_coefs;
    }

    gsMatrix<T> perCoefs(const gsMatrix<T>& coefs, int dir) const
    {
        GISMO_UNUSED(dir);
        GISMO_ASSERT(dir==0, "Error");
        return perCoefs(coefs);
    }

    /// @brief Helper function for transforming periodic coefficients
    /// to full coefficients
    void expandCoefs(gsMatrix<T> & coefs) const
    {
        const index_t sz = coefs.rows();
        coefs.conservativeResize(sz+m_periodic, Eigen::NoChange);
        coefs.bottomRows( m_periodic ) = coefs.topRows( m_periodic );
    }

    /// @brief Helper function for transforming full coefficients to
    /// periodic coefficients
    void trimCoefs(gsMatrix<T> & coefs) const
    {
        const index_t sz = coefs.rows();
        coefs.conservativeResize(sz-m_periodic, Eigen::NoChange);
    }

    /// @brief Returns the size of the basis ignoring the bureaucratic way of
    /// turning the basis into periodic.
    int trueSize() const
    { return this->size() + m_periodic; }
 
// Data members
protected:

    /// @brief Degree
    short_t m_p;

    /// @brief Knot vector
    KnotVectorType m_knots;
    
    /// @brief Denotes whether the basis is periodic, ( 0 -- non-periodic, >0 -- number of ``crossing" functions)
    int m_periodic;

    /*/// @brief Multiplicity of the p+1st knot from the beginning and from the end.
      int m_bordKnotMulti;*/

}; // class gsTensorBSplineBasis<1>


//Using C++11 alias:
// template<class T, class KnotVectorType>
// using gsBSplineBasis = gsTensorBSplineBasis<1,T>

/** \brief
    A univariate B-spline basis.

    \tparam T coefficient type
    \tparam KnotVectorType the type of knot vector to use

    \ingroup basis
    \ingroup Nurbs
*/
template<class T>
class gsBSplineBasis : public gsTensorBSplineBasis<1,T>
{
public:
    typedef memory::shared_ptr< gsBSplineBasis > Ptr;
    typedef memory::unique_ptr< gsBSplineBasis > uPtr;

    typedef gsKnotVector<T> KnotVectorType;
    typedef gsTensorBSplineBasis<1,T> Base;
    typedef gsBSplineBasis<T> Self_t;

    /// @brief Associated geometry type
    typedef typename gsBSplineTraits<1,T>::Geometry GeometryType;

    /// @brief Associated Boundary basis type
    typedef typename gsBSplineTraits<0,T>::Basis BoundaryBasisType;

public:

    // Default empty constructor
    explicit gsBSplineBasis(const bool periodic = false )
    { 
        m_p = 0;
        m_knots.initClamped(0);
        m_periodic = 0;

        if( periodic )
            this->_convertToPeriodic();

        if( ! this->check() )
            gsWarn << "Warning: Inconsistent "<< *this<< "\n";
    }

    /// @brief Construct BSpline basis of a knot vector
    explicit gsBSplineBasis(KnotVectorType KV, const bool periodic = false)
    { 
        m_p        = KV.degree();
        m_knots.swap(KV);
        m_periodic = 0;

        if( periodic )
            this->_convertToPeriodic();

        if( ! this->check() )
            gsWarn << "Warning: Insconsistent "<< *this<< "\n";
    }

    /// @brief Compatibility constructor with input an std::vector containing
    /// a single knotvector
    explicit gsBSplineBasis(std::vector<KnotVectorType> KV)
    { 
        GISMO_ASSERT(1 == KV.size(), "Expecting a single knotvector." );

        m_p        = KV.front().degree();
        m_knots.swap(KV.front());
        m_periodic = 0;
    }

    /// @brief Construct a BSpline basis
    /// @param u0 starting parameter
    /// @param u1 end parameter parameter
    /// @param interior number of interior knots
    /// @param degree degree of the spline space
    /// @param mult_interior multiplicity at the interior knots
    /// @param periodic specifies if basis is periodic or not
    gsBSplineBasis(const T u0, const T u1, const unsigned interior, 
                   const int degree, const unsigned mult_interior=1,
                   const bool periodic = false )
    { 
        m_p = degree;
        m_knots.initUniform(u0, u1, interior, m_p+1, mult_interior, m_p);
        m_periodic = 0;

        if( periodic )
            this->_convertToPeriodic();
        
        if( ! this->check() )
            gsWarn << "Warning: Insconsistent "<< *this<< "\n";
    }

/*
    /// @brief Copy Constructor
    gsBSplineBasis( const gsBSplineBasis & o)
    : Base(o)
    { }
*/

    // Look at gsBasis class for a description
    GISMO_CLONE_FUNCTION(gsBSplineBasis)

    // Look at gsBasis class for a description
    Self_t & component(short_t i);
    
    // Look at gsBasis class for a description
    const Self_t & component(short_t i) const;

    memory::unique_ptr<gsGeometry<T> > makeGeometry( gsMatrix<T> coefs ) const;
            
private:
    
    using Base::m_p;
    using Base::m_knots;
    using Base::m_periodic;
};


} // namespace gismo


// *****************************************************************
#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsBSplineBasis.hpp)
#else
#ifdef gsBSplineBasis_EXPORT
#include GISMO_HPP_HEADER(gsBSplineBasis.hpp) // needed on thebrain?
#undef  EXTERN_CLASS_TEMPLATE
#define EXTERN_CLASS_TEMPLATE CLASS_TEMPLATE_INST
#endif
namespace gismo
{
template<class T> const short_t gsTensorBSplineBasis<1,T>::Dim; //-O3 (SLE11) fix
EXTERN_CLASS_TEMPLATE gsTensorBSplineBasis<1,real_t>;
EXTERN_CLASS_TEMPLATE gsBSplineBasis<real_t>;
}
#endif
// *****************************************************************
