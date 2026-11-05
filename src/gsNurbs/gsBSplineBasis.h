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
#include <gsDomain/gsTensorDomainIterator.h>
#include <gsDomain/gsTensorDomainBoundaryIterator.h>

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

        // Default empty constructor
    explicit gsTensorBSplineBasis(const bool periodic = false )
    {
        m_p = 0;
        m_knots.initClamped(0);

        if( periodic )
            m_knots.setPeriodic(true);

        if( ! this->check() )
            gsWarn << "Warning: Inconsistent "<< *this<< "\n";
    }

    /// @brief Construct BSpline basis of a knot vector
    explicit gsTensorBSplineBasis(KnotVectorType KV, const bool periodic = false)
    {
        m_p        = KV.degree();
        m_knots.swap(KV);

        if( periodic )
            m_knots.setPeriodic(true);

        if( ! this->check() )
            gsWarn << "Warning: Insconsistent "<< *this<< "\n";
    }

    /// @brief Compatibility constructor with input an std::vector containing
    /// a single knotvector
    explicit gsTensorBSplineBasis(std::vector<KnotVectorType> KV)
    {
        GISMO_ASSERT(1 == KV.size(), "Expecting a single knotvector." );

        m_p        = KV.front().degree();
        m_knots.swap(KV.front());
    }

    /// @brief Construct a BSpline basis
    /// @param u0 starting parameter
    /// @param u1 end parameter parameter
    /// @param interior number of interior knots
    /// @param degree degree of the spline space
    /// @param mult_interior multiplicity at the interior knots
    /// @param periodic specifies if basis is periodic or not
    ///
    /// If \a periodic, the space is the periodic spline space with
    /// interior+1 breakpoints per period, all of multiplicity
    /// \a mult_interior (also at the seam), i.e. C^{degree-mult_interior}
    /// continuity, with size() == (interior+1)*mult_interior functions and
    /// numCrossingFunctions() == degree+1-mult_interior.
    /// Refused with GISMO_ERROR unless 1 <= mult_interior <= degree and
    /// (interior+1)*mult_interior >= degree+1.
    gsTensorBSplineBasis(const T u0, const T u1, const unsigned interior,
                   const int degree, const unsigned mult_interior=1,
                   const bool periodic = false )
    {
        if ( periodic )
        {
            GISMO_ENSURE( mult_interior >= 1,
                "Periodic basis: the multiplicity must be at least 1." );
            GISMO_ENSURE( degree >= 1 && mult_interior <= static_cast<unsigned>(degree),
                "Periodic basis: the multiplicity must be at most the degree." );
            GISMO_ENSURE( (interior+1)*static_cast<index_t>(mult_interior) >= degree+1,
                "Periodic basis: the number of functions must be at least degree+1." );

            m_p = degree;
            KnotVectorType kv( _periodicKnots(u0, u1, interior, degree, mult_interior), degree );
            kv.setPeriodic(true);
            m_knots.swap(kv);
        }
        else
        {
            m_p = degree;
            m_knots.initUniform(u0, u1, interior, m_p+1, mult_interior, m_p);
        }

        if( ! this->check() )
            gsWarn << "Warning: Insconsistent "<< *this<< "\n";
    }

    using gsBasis<T>::create; //unhide from gsBasis

    static typename gsBasis<T>::uPtr create(KnotVectorType KV, short_t dim)
    {
        typedef typename gsBasis<T>::uPtr basisPtr;

        if ( dim > 1 && KV.isPeriodic() )
            GISMO_ERROR("gsBSplineBasis::create: a periodic knot vector would be periodic in every direction; "
                        "periodicity in one direction of a tensor basis is not supported.");

        switch (dim)
        {
        case 1:
            return basisPtr(new gsBSplineBasis<T>(give(KV)));
            break;
        case 2:
            return basisPtr(new gsTensorBSplineBasis<2,T>(KV,KV));
            break;
        case 3:
            return basisPtr(new gsTensorBSplineBasis<3,T>(KV,KV,KV));
            break;
        case 4:
            return basisPtr(new gsTensorBSplineBasis<4,T>(KV,KV,KV,KV));
            break;
        }
        GISMO_ERROR("Dimension should be between 1 and 4.");
    }

    GISMO_CLONE_FUNCTION(gsTensorBSplineBasis)

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

    static typename gsBasis<T>::uPtr create(std::vector<KnotVectorType> cKV);

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
        m_knots.swap(other.m_knots);
    }

/* Virtual member functions required by the base class */

    // Look at gsBasis class for a description
    short_t domainDim() const { return Dim; }

    // Unhide/forward gsTensorBasis<1,T>::size(k), since the overload
    // with size() automatically hides it in this class
    using Base::size;

    // Look at gsBasis class for a description
    index_t size() const { return m_knots.size() - m_p - 1 - m_knots.numCrossingFunctions(); }

    // Look at gsBasis class for a description
    size_t numElements(boxSide const & s = 0) const { return (0==s?m_knots.numElements():1); }

    // Look at gsBasis class for a description
    size_t elementIndex(const gsVector<T> & u ) const;

    // Same as gsBasis::elementIndex but argument is a value instead of a vector
    size_t elementIndex(T u ) const;

    /// @brief Returns an element (as a 1x2 matrix of parameter values) inside
    /// the support of the j-th basis function.
    ///
    /// For a crossing function of a periodic basis this is the first element of
    /// the domain.
    gsMatrix<T> elementInSupportOf(index_t j) const;

    /// @brief Returns span (element) indices of the beginning and end
    /// of the support of the i-th basis function.
    ///
    /// On a periodic basis a crossing function (i < numCrossingFunctions())
    /// has two disjoint pieces of support; the result is then the hull,
    /// the range of all elements [0, numElements()], which is a superset of
    /// the true support.
    void elementSupport_into(const index_t i, gsMatrix<index_t,1,2>& result) const
    {
        if ( isPeriodic() && i < static_cast<index_t>(m_knots.numCrossingFunctions()) )
        {
            result.resize(1,2);
            result << 0, static_cast<index_t>(m_knots.numElements());
            return;
        }
        gsMatrix<index_t> tmp_vec;
        m_knots.supportIndex_into(i, tmp_vec);
        result = tmp_vec.cwiseMax(0).cwiseMin(m_knots.numElements());
    }

    template <int _Rows>
    gsMatrix<T> elementDom(const gsMatrix<index_t,_Rows,2> & box) const
    {
        gsMatrix<T> rvo(1,2);
        rvo.at(0) = m_knots.uValue(box.at(0));
        rvo.at(1) = m_knots.uValue(box.at(1));
        return rvo;
    }

    // Look at gsBasis class for a description
    const TensorSelf_t & component(short_t i) const override;

    // Look at gsBasis class for a description
    TensorSelf_t & component(short_t i) override;

    memory::unique_ptr<gsGeometry<T> > makeGeometry( gsMatrix<T> coefs ) const;
    
    /// @brief Returns the anchors (greville points) of the basis
    ///
    /// On a periodic basis the result holds n = size() anchors: the Greville
    /// points of the stretched functions 0..n-1, reduced into [a,b) modulo the
    /// period b-a. They are not sorted. An interior knot of multiplicity p+1
    /// gives near-duplicate anchors (greville_into perturbs duplicates by
    /// 1e-9), so a collocation matrix built on them is ill-conditioned rather
    /// than singular.
    void anchors_into(gsMatrix<T> & result) const override
    {
        if ( isPeriodic() )
        {
            gsMatrix<T> g;
            m_knots.greville_into(g, false); // clamping would give Greville points of a clamped basis
            result = g.leftCols(size());
            for ( index_t k = 0; k != result.cols(); ++k )
                result(0,k) = _reduceToPeriod(result(0,k));
        }
        else
            m_knots.greville_into(result);
    }

    /// @brief Returns the anchor (greville point) of the i-th basis function
    ///
    /// On a periodic basis this is column i of anchors_into(), i.e. reduced
    /// into [a,b) modulo the period.
    void anchor_into(index_t i, gsMatrix<T> & result) const override
    {
        result.resize(1,1);
        if ( isPeriodic() )
        {
            GISMO_ASSERT( i >= 0 && i < size(), "Invalid index of basis function." );
            gsMatrix<T> g;
            m_knots.greville_into(g, false);
            result(0,0) = _reduceToPeriod(g(0,i));
        }
        else
            result(0,0) = m_knots.greville(i);
    }

    /// @brief Control net of the spline with coefficient rows \a nodes: one
    /// vertex per row and the edges (i-1,i). \a nodes has size() rows (n
    /// rows on a periodic basis, not n+m). On a periodic basis the control
    /// polygon of the closed curve is closed by the edge (n-1,0).
    void connectivity(const gsMatrix<T> & nodes,
                      gsMesh<T> & mesh) const override;

    // Look at gsBasis class for a description
    void active_into(const gsMatrix<T> & u, gsMatrix<index_t>& result) const override;

    /// @brief True iff u lies in the (closed) support of the i-th function.
    ///
    /// On a periodic basis all stretched copies of function i are tested, and
    /// the result is false outside [a,b].
    bool isActive(const index_t i, const gsVector<T> & u) const;

    /// @brief Indices of the boundary functions. Empty on a periodic basis.
    gsMatrix<index_t> allBoundary( ) const ;

    /// @brief Indices of the boundary functions at side s and offset.
    /// Throws (GISMO_ERROR) on a periodic basis, which has no boundary.
    gsMatrix<index_t> boundaryOffset(boxSide const & s,index_t offset) const override;

#ifdef __DOXYGEN__
    /// @brief Gives back the boundary basis at boxSide s
    ///
    /// Throws (GISMO_ERROR) on a periodic basis, which has no boundary.
    typename BoundaryBasisType::uPtr boundaryBasis(boxSide const & s);
#endif
    GISMO_UPTR_FUNCTION_DEC(gsConstantBasis<T>, boundaryBasis, boxSide const &)

    // Look at gsBasis class for a description
    gsMatrix<T> support() const ;

    /// @brief Support of the i-th basis function.
    ///
    /// On a periodic basis a crossing function (i < numCrossingFunctions()) has
    /// two pieces of support in the domain, [a, t_{i+p+1}] and [t_{i+n}, b];
    /// the hull [a,b], a superset of the true support, is returned.
    gsMatrix<T> support(const index_t & i ) const ;

    /// @brief Only meaningfull for periodic basis: For basis members that have
    /// a twin, this function returns the other twin index, otherwise it
    /// returns the same index as the argument
    index_t twin(index_t i) const ;

    // Look at gsBasis class for a description
    virtual void eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const;

    /// @brief Evaluates the i-th basis function at the points u.
    ///
    /// On a periodic basis this is the folded function i (the sum of its
    /// stretched copies): zero outside [a,b], and the left limit at b. Cost
    /// O(p^2) per point.
    virtual void evalSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result) const;

    /// @brief Evaluates the spline with coefficient rows \a coefs at the points u.
    ///
    /// \a coefs has size() rows; on a periodic basis these are the n
    /// periodic coefficients (not n+m): they are expanded internally to the
    /// n+m stretched coefficients for De Boor's algorithm. Cost O(n d)
    /// for the expansion plus O(p^2) per point.
    virtual void evalFunc_into(const gsMatrix<T> & u, const gsMatrix<T> & coefs, gsMatrix<T>& result) const;

    // Look at gsBasis class for a description
    void deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const ;

    /// @brief First derivative of the i-th basis function at the points u.
    ///
    /// On a periodic basis this is the folded function i (the sum of its
    /// stretched copies): zero outside [a,b], and the left limit at b. Cost
    /// O(p^2) per point.
    void derivSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result ) const ;

    /// @brief First derivative of the spline with coefficient rows \a coefs.
    /// \a coefs has size() rows (n rows on a periodic basis, not n+m).
    void deriv_into(const gsMatrix<T> & u, const gsMatrix<T> & coefs, gsMatrix<T>& result ) const ;

    // Look at gsBasis class for a description
    void deriv2_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const ;

    /// @brief Second derivative of the i-th basis function at the points u.
    ///
    /// On a periodic basis this is the folded function i (the sum of its
    /// stretched copies): zero outside [a,b], and the left limit at b. Cost
    /// O(p^2) per point.
    void deriv2Single_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result ) const ;

    /// @brief Second derivative of the spline with coefficient rows \a coefs.
    /// \a coefs has size() rows (n rows on a periodic basis, not n+m).
    void deriv2_into(const gsMatrix<T> & u, const gsMatrix<T> & coefs, gsMatrix<T>& result ) const ;

    // Look at gsBasis class for a description
    gsMatrix<T> laplacian(const gsMatrix<T> & u ) const ;

    // Look at gsBasis class for a description
    typename gsBasis<T>::uPtr tensorize(const gsBasis<T> & other) const;

    /// @brief Check the BSplineBasis for consistency
    ///
    /// For a periodic basis it checks that the degrees agree, that
    /// 1 <= m <= p and n >= p+1, that the seam multiplicity is exactly
    /// p+1-m, and that the ghost knots satisfy t_{k+n} = t_k + (b-a)
    /// for k = 0..m+p.
    bool check() const
    {
        if ( m_knots.isPeriodic() )
        {
            if ( m_knots.degree() != m_p )
                return false;
            const index_t K = m_knots.size();
            const index_t m = m_knots.numCrossingFunctions();
            const index_t n = K - m_p - 1 - m;
            if ( m < 1 || m > m_p || n < m_p + 1 )
                return false;
            const T a = m_knots[m_p];
            const T b = m_knots[K-m_p-1];
            const T L = b - a;
            if ( !( m_knots[m-1] < a ) )
                return false;
            const T tol = 1e-12 * (std::max)( (std::max)(math::abs(a), math::abs(b)), math::abs(L) );
            for ( index_t k = 0; k <= m + m_p; ++k )
                if ( math::abs( m_knots[k+n] - m_knots[k] - L ) > tol )
                    return false;
            return true;
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

    /// @brief The n-th derivative of the i-th basis function at the points u.
    ///
    /// On a periodic basis this is the folded function i (the sum of its
    /// stretched copies): zero outside [a,b], and the left limit at b. Cost
    /// O(p^2) per point (O(n p) more for derivative order n).
    virtual void evalDerSingle_into(index_t i, const gsMatrix<T> & u,
                                    int n, gsMatrix<T>& result) const;

    // Look at gsBasis class for a description
    virtual void evalAllDers_into(const gsMatrix<T> & u, int n,
                                  std::vector<gsMatrix<T> >& result,
                                  bool sameElement = false) const;

    using Base::evalAllDersSingle_into;
    /// @brief Values and derivatives up to order n of the i-th basis function.
    ///
    /// On a periodic basis this is the folded function i (the sum of its
    /// stretched copies): zero outside [a,b], and the left limit at b. Cost
    /// O(p^2) per point (O(n p) more for derivative order n).
    GISMO_DEPRECATED
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

    /// @brief Reduces \a x into [a,b) modulo the period b-a.
    T _reduceToPeriod(T x) const
    {
        const T a = domainStart(), b = domainEnd(), L = b - a;
        x -= L * math::floor((x - a) / L);
        // x is now within rounding of [a,b]; a value still outside [a,b) is
        // the seam point, whose representative is a.
        if ( x >= b || x < a ) x = a;
        return x;
    }

    /// @brief Returns the index of the first active (ie. non-zero) basis function at point u
    /// Takes into account non-clamped knots.
    inline index_t firstActive(T u) const {
        return ( inDomain(u) ? (m_knots.iFind(u)-m_knots.begin()) - m_p : 0 );
    }

    // Number of active functions at any point of the domain
    inline index_t numActive() const { return m_p + 1; }

    // Look at gsBasis class for a description
    memory::shared_ptr<gsDomain<T> > domain() const { return memory::make_shared_not_owned(&m_knots); }

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

    /// @brief Removes the knot \em knot in the underlying knot vector.
    void removeKnot(T knot, int mult=1)
    { m_knots.remove( knot, mult); }

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
    void uniformRefine(int numKnots = 1, int mul=1, int dir=-1)
    { GISMO_UNUSED(dir); m_knots.uniformRefine(numKnots,mul); }

    /// Look at gsBasis class for a description.
    ///
    /// On a periodic basis the knots are inserted in each element of the
    /// domain [a,b] only (ghost spans contribute none), the same knots as
    /// uniformRefine(), through refine_withCoefs(); \a coefs must hold
    /// the n = size() reduced coefficients and are replaced by the n'
    /// refined ones, T*coefs with the exact n'xn transfer T of
    /// uniformRefine_withTransfer(). Throws std::runtime_error, leaving the basis and
    /// \a coefs unchanged, if the knots are refused (multiplicity above
    /// degree+1) or if \a coefs does not have n rows. Cost as refine_withCoefs().
    void uniformRefine_withCoefs(gsMatrix<T>& coefs, int numKnots = 1, int mul = 1, short_t const dir = -1);

    /// Look at gsBasis class for a description.
    ///
    /// On a periodic basis the knots are inserted in each element of the
    /// domain [a,b] only, as in uniformRefine(), through
    /// refine_withTransfer(); the result is the exact n'xn transfer of
    /// reduced coefficients. Throws std::runtime_error, leaving the basis
    /// and \a transfer unchanged, if the knots are refused. Cost as refine_withTransfer().
    void uniformRefine_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, int numKnots = 1, int mul = 1);

    // Look at gsBasis class for a description
    void uniformCoarsen(int numKnots = 1)
    { m_knots.coarsen(numKnots); }

    /// Look at gsBasis class for a description.
    ///
    /// On a periodic basis the removed period knots are re-inserted into
    /// the coarsened basis by the periodic refine_withTransfer(), so the
    /// transfer is the exact n x n_c map from the n_c reduced coarse
    /// coefficients to the n reduced fine ones. Throws std::runtime_error,
    /// leaving the basis and \a transfer unchanged, if the knot vector
    /// refuses the coarsening. Cost as refine_withTransfer().
    void uniformCoarsen_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, int numKnots = 1);

    /// @brief Refine the basis by inserting the given knots and perform knot
    /// refinement for the given coefficient matrix.
    ///
    /// On a periodic basis the coefficients are the n = size() reduced ones
    /// (one per row) and are replaced by T*coefs, with T the exact n'xn
    /// transfer of refine_withTransfer(); it is a GISMO_ENSURE that
    /// coefs.rows() == size(). A knot in (a,b) is inserted, a knot in
    /// (a-L,a) or (b,b+L) is inserted as knot+-L, with L = b-a; a knot on
    /// the seam or outside (a-L,b+L), or one raising an interior
    /// multiplicity above degree+1, is refused. Every refusal throws
    /// std::runtime_error and leaves the basis and \a coefs unchanged.
    /// Cost O((n+p)p + rp) for r inserted knots.
    void refine_withCoefs(gsMatrix<T>& coefs, const std::vector<T>& knots);
    void refine_withCoefs(gsMatrix<T>& /* coefs */, const gsMatrix<T> & /* boxes */)
    {GISMO_NO_IMPLEMENTATION}

    // compatibility with tensor-bsplines
    void refine_withCoefs(gsMatrix<T> & coefs,
                          const std::vector< std::vector<T> >& refineKnots)
    {
        refine_withCoefs(coefs,refineKnots.front() );
    }

    /// @brief Refine the basis by inserting the given knots and produce a sparse matrix which maps coarse coefficient vectors to refined ones.
    ///
    /// On a periodic basis with n = size() functions the result is the
    /// exact n'xn matrix T, n' = n + (number of inserted knots), mapping
    /// reduced coefficients to reduced coefficients, c' = T*c, and the
    /// basis then holds the regenerated periodic knot vector. Knots are
    /// accepted and refused as in refine_withCoefs(); a refusal throws
    /// std::runtime_error and leaves the basis and \a transfer unchanged.
    /// The matrix has O(n'p) nonzeros and costs O((n+p)p + rp) for r
    /// inserted knots. See _periodicRefineTransfer() for the construction.
    void refine_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, const std::vector<T>& knots);
    void refine_withTransfer(gsMatrix<T>& /* coefs */, const gsMatrix<T> & /* boxes */)
    {GISMO_NO_IMPLEMENTATION}
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
        GISMO_ENSURE( !isPeriodic() && !other.isPeriodic(),
                      "refine_k: not defined for periodic bases (the reference knots contain ghosts and the seam).");
        GISMO_ASSERT( m_p >= other.m_p, "Degree of other knot-vector should be lower.");
        //for (typename std::vector<T>::iterator it =
        //         m_knots.begin(); it != m_knots.end(); ++it)
        //    GISMO_ASSERT( has(*it), "Knot "<< *it<<" is not in the knot vector.");

        // grab unique knots
        const std::vector<T> knots = other.m_knots.unique();

        // Work on a copy: a refusal by the knot vector leaves *this unchanged
        KnotVectorType kv(m_knots);
        // Increase the degree without adjusting any knot
        kv.set_degree(kv.degree() + i);
        // Adjust (reduce) smoothness to satisfy initial constraint knots
        kv.insert(knots,i);
        m_knots.swap(kv);
        m_p = m_knots.degree();
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
        m_knots.degreeElevate(i);
        m_p = m_knots.degree();
    }

    // Look at gsBasis for documentation
    void degreeReduce (short_t const & i = 1, short_t const dir = -1)
    {
        GISMO_UNUSED(dir);
        GISMO_ASSERT( dir == -1 || dir == 0, "Invalid direction");
        GISMO_ASSERT( i<=m_p, "Cannot reduce degree to negative");
        m_knots.degreeReduce(i);
        m_p = m_knots.degree();
    }

    // Look at gsBasis for documentation
    void degreeIncrease(short_t const & i = 1, short_t const dir = -1)
    {
        GISMO_UNUSED(dir);
        GISMO_ASSERT( dir == -1 || dir == 0, "Invalid direction");
        m_knots.degreeIncrease(i);
        m_p = m_knots.degree();
    }

    // Look at gsBasis for documentation
    void degreeDecrease(short_t const & i = 1, short_t const dir = -1)
    {
        GISMO_UNUSED(dir);
        GISMO_ASSERT( dir == -1 || dir == 0, "Invalid direction");
        m_knots.degreeDecrease(i);
        m_p = m_knots.degree();
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
        return m_knots.isPeriodic();
    }

    /// @brief Index of the function at corner c. Throws (GISMO_ERROR) on a
    /// periodic basis, which has no boundary.
    index_t functionAtCorner(boxCorner const & c) const;

    /// @brief Returns number of functions crossing the boundary of the knot vector.
    int numCrossingFunctions () const
    {
        return static_cast<int>(m_knots.numCrossingFunctions());
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
    /// (GISMO_ERROR, leaving the basis unchanged, if the knot vector is
    /// not suitable); if false, drops the periodic property.
    void setPeriodic(bool flag = true)
    {
        m_knots.setPeriodic(flag);
    }

    // Compatible with tensor B-spline basis
    void setPeriodic(int dir)
    {
        GISMO_UNUSED(dir);
        GISMO_ASSERT(dir==0, "Invalid direction");
        m_knots.setPeriodic(true);
    }

    GISMO_DEPRECATED
    typename gsBasis<T>::domainIter makeDomainIterator() const
    {
        return m_knots.beginAll();
    }

    GISMO_DEPRECATED
    typename gsBasis<T>::domainIter makeDomainIterator(const boxSide & s) const
    {
        return m_knots.beginBdr(s);
    }

    // Look at gsBasis class for a description
    void reverse() { m_knots.reverse(); }

    /// @brief Not available on a periodic basis (throws, a periodic basis has
    /// no boundary); a no-op otherwise.
    void matchWith(const boundaryInterface &, const gsBasis<T> &,
                   gsMatrix<index_t> &, gsMatrix<index_t> &, index_t) const
    {
        if ( isPeriodic() )
            GISMO_ERROR("gsBSplineBasis: a periodic basis has no boundary.");
    }


    /// @brief Matches the boundary functions across an interface. Throws
    /// (GISMO_ERROR) if this basis or \a other is periodic.
    void matchWith(const boundaryInterface & bi,
                   const gsBasis<T> & other,
                   gsMatrix<index_t> & bndThis,
                   gsMatrix<index_t> & bndOther) const;

protected:

    /// @brief Periodic knot insertion by pad-and-fold. Does not modify
    /// *this. Returns in \a transfer the exact n'xn map of reduced
    /// coefficients and in \a kvNew the refined periodic knot vector.
    ///
    /// Let t_0..t_{K-1} be the stretched knots, [a,b] the domain and L = b-a.
    /// The knot vector is continued periodically by P = 2p+1 knots on
    /// each side, giving a non-periodic basis whose domain strictly
    /// contains [a,b] (no knot has multiplicity above p+1, so
    /// t_{-p-1} < t_0 and t_{K-1} < t_{K+p}). The inserted knots and all
    /// their translates by multiples of L inside the open padded domain
    /// are inserted there by Boehm's algorithm.
    ///
    /// The periodic curve sum_j c_j P_j equals sum_{i<n+m} c_{i mod n} B_i
    /// on the stretched basis, i.e. the same function on the padded
    /// basis with zero padding coefficients (pad functions vanish on
    /// (a,b)). Knot insertion reproduces it exactly on the padded domain,
    /// where the window of n'+m' refined functions starting m' functions
    /// before the first one whose first knot equals a (padded index of
    /// that knot minus m') is linearly independent and the only one active
    /// on [a,b]; hence its coefficients are those of the curve in the
    /// new stretched basis, and by periodicity they repeat with period
    /// n'. The transfer is T = S_new * T_pad * E_old: E_old embeds the n
    /// reduced coefficients (stretched, then padded), T_pad is the padded
    /// Boehm transfer and S_new selects the first n' window rows. In
    /// the fold, columns i of the stretched block are summed into column
    /// i mod n (setFrom sums duplicate triplets).
    ///
    /// Cost O((n+p)p) for the transfer, O(rp) for r translates, O(K)
    /// for the copies.
    void _periodicRefineTransfer(const std::vector<T> & knots,
                                 gsSparseMatrix<T,RowMajor> & transfer,
                                 KnotVectorType & kvNew) const;

    /// @brief Knots of uniform refinement (\a numKnots new knots of
    /// multiplicity \a mul per element) of the domain elements of a
    /// periodic basis; none for ghost spans.
    void _domainUniformRefinementKnots(int numKnots, int mul,
                                       std::vector<T> & newKnots) const;

    /// @brief Stretched knot vector (ghost knots included) of the periodic
    /// space with \a interior+1 breakpoints per period on [\a u0, \a u1],
    /// each of multiplicity \a mult, i.e. t_k = s_{(k-m) mod n} + L*floor((k-m)/n)
    /// for k = 0..n+m+p, with n = (interior+1)*mult, m = degree+1-mult, L = u1-u0.
    static std::vector<T> _periodicKnots(const T u0, const T u1,
                                         const unsigned interior,
                                         const int degree, const unsigned mult);

public:

    /// @brief Helper for the n+m stretched coefficients of a periodic basis.
    ///
    /// @param coefs the n+m = trueSize() stretched coefficient rows (control
    /// points, one per row). It must not be called with the n rows of a
    /// periodic curve: it would overwrite rows n-m..n-1.
    ///
    /// @return copy of coefs with the first numCrossingFunctions() rows copied to
    /// the last numCrossingFunctions() rows.
    gsMatrix<T> perCoefs( const gsMatrix<T>& coefs ) const
    {
        const index_t m = numCrossingFunctions();
        gsMatrix<T> per_coefs = coefs;
        per_coefs.bottomRows( m ) = coefs.topRows( m );
        return per_coefs;
    }

    gsMatrix<T> perCoefs(const gsMatrix<T>& coefs, int dir) const
    {
        GISMO_UNUSED(dir);
        GISMO_ASSERT(dir==0, "Error");
        return perCoefs(coefs);
    }

    /// @brief Expands the n periodic coefficient rows to the n+m stretched
    /// rows (m = numCrossingFunctions()): the first m rows are appended,
    /// so that stretched function i takes row i mod n (m <= degree+1 <= n).
    void expandCoefs(gsMatrix<T> & coefs) const
    {
        const index_t m = numCrossingFunctions();
        const index_t sz = coefs.rows();
        coefs.conservativeResize(sz+m, gsEigen::NoChange);
        coefs.bottomRows( m ) = coefs.topRows( m );
    }

    /// @brief Drops the last numCrossingFunctions() rows, turning the n+m
    /// stretched coefficient rows into n periodic rows. The ghost rows are
    /// not checked against their twins.
    void trimCoefs(gsMatrix<T> & coefs) const
    {
        const index_t sz = coefs.rows();
        coefs.conservativeResize(sz-numCrossingFunctions(), gsEigen::NoChange);
    }

    /// @brief Returns the size of the basis ignoring the bureaucratic way of
    /// turning the basis into periodic.
    int trueSize() const
    { return this->size() + numCrossingFunctions(); }

// Data members
protected:

    /// @brief Degree
    short_t m_p;

    /// @brief Knot vector
    KnotVectorType m_knots;

    /*/// @brief Multiplicity of the p+1st knot from the beginning and from the end.
      int m_bordKnotMulti;*/

}; // class gsTensorBSplineBasis<1>

#ifdef GISMO_WITH_PYBIND11

  /**
   * @brief Initializes the Python wrapper for the class: gsBSplineBasis
   */
  void pybind11_init_gsBSplineBasis(pybind11::module &m);

#endif // GISMO_WITH_PYBIND11

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
}
#endif
// *****************************************************************
