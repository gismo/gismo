/** @file gsBSpline.h

    @brief Represents a B-spline curve/function with one parameter

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <limits>

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsGeometry.h>

#include <gsNurbs/gsBSplineBasis.h>
#include <gsNurbs/gsBoehm.h>
#include <gsNurbs/gsBSplineSolver.h>

#include <gsUtils/gsPointGrid.h>

namespace gismo {
namespace internal {

template<class T>
struct gsCurveIntersectionResult;

} // namespace internal
} // namespace gismo

namespace gismo
{

/** \brief
    A B-spline function of one argument, with arbitrary target dimension.

    This is the geometry type associated with gsBSplineBasis.

    \tparam T coefficient type

    # Periodic curves #

    A curve on a periodic basis is a closed curve over [a,b] (see
    gsKnotVector for p, n, m, mu and the period L). It stores reduced
    coefficients, exactly \c basis().size() = n rows, one per periodic
    function (I3). Ghost rows of the stretched representation exist only
    inside an algorithm.

    Construction routes:
     + the from-scratch constructor with \c periodic = true;
     + setPeriodic() on a clamped curve whose first and last control points
       coincide. The curve is kept exactly, with a C^0 seam (the n+m rows of
       the clamped curve pass the twin check);
     + any constructor, or setPeriodic(), given n+m rows: accepted when ghost
       row n+j equals row j within 100 eps max|coefs|, then compressed to n
       rows. Any other row count is refused.

    setPeriodic(false) unwraps exactly onto the stretched knots.

    Exact operations: knot insertion, degreeElevate() (exact in exact
    arithmetic) and reverse() (keeps [a,b], c'_j = c_{(m-1-j) mod n}).
    Approximations: degreeIncrease(), degreeDecrease() and degreeReduce()
    interpolate the curve at the periodic anchors of the new basis with a
    direct sparse LU. They are refused when the backward error exceeds
    1000 eps (||C|| ||x|| + ||v||), where ||.|| is the largest absolute
    entry (also for C).

    Refused (std::runtime_error, curve unchanged): merge(), segmentFromTo(),
    splitAt(), toBezier(), intersect(), the corner functions
    (isPatchCorner(), findCorner(), setOriginCorner(), setFurthestCorner()),
    and knots on the seam or farther than one period away.

    Evaluation outside [a,b] returns 0 (no wrapping), as for the basis.

    \ingroup geometry
    \ingroup Nurbs
*/


template<class T>
class gsBSpline : public gsGeoTraits<1,T>::GeometryBase
{
public:
    typedef gsKnotVector<T> KnotVectorType;

    typedef typename gsGeoTraits<1,T>::GeometryBase Base;

    typedef gsBSplineBasis<T> Basis;

    /// Shared pointer for gsBSpline
    typedef memory::shared_ptr< gsBSpline > Ptr;

    /// Unique pointer for gsBSpline
    typedef memory::unique_ptr< gsBSpline > uPtr;

public:

    /// Default empty constructor.
    gsBSpline() { }

    // enable swap
    using gsGeometry<T>::swap;

    /// @brief Construct B-Spline by basis and coefficient matrix.
    ///
    /// A curve on a periodic basis stores exactly basis.size() = n
    /// coefficient rows (one per periodic function). \a coefs may have n
    /// rows, or n+m rows (m = basis.numCrossingFunctions()) whose ghost
    /// rows n+j equal their twins j within 100*eps*max|coefs|; the ghost
    /// rows are then dropped. Any other row count, or a ghost row that
    /// differs from its twin, raises GISMO_ERROR.
    gsBSpline( const Basis & basis, gsMatrix<T> coefs )
    : Base( basis, _reducePeriodicCoefs(basis, give(coefs)) )
    { }

    /// @brief Construct B-Spline by a knot vector and coefficient matrix.
    ///
    /// If \a periodic is true the knot vector is converted by
    /// gsKnotVector::setPeriodic() (GISMO_ERROR if it is not suitable). A
    /// knot vector that is already periodic stays periodic whatever the value of
    /// \a periodic. On a periodic basis the curve stores n = basis size rows;
    /// \a coefs may have n or n+m rows, as in gsBSpline(const Basis&,
    /// gsMatrix<T>). Nothing is stored if a check fails.
    gsBSpline(KnotVectorType KV, gsMatrix<T> coefs, bool periodic = false )
    {
        Basis b(give(KV));
        if( periodic )
            b.setPeriodic(true);

        if( b.isPeriodic() )
            coefs = _reducePeriodicCoefs(b, give(coefs));
        else if( coefs.rows() + b.degree() + 1 != static_cast<index_t>( b.knots().size() ) )
            gsWarn << "gsBSpline Warning: #Knots="<< b.knots().size()<< ", #coefs="<< coefs.rows() <<"\n";

        this->m_basis = new Basis(give(b));
        m_coefs.swap(coefs);
    }


    /// @brief Construct a B-spline from an interval and knot vector specification.
    /// \param u0 starting parameter
    /// \param u1 end parameter
    /// \param interior number of interior knots
    /// \param degree degree of the spline space
    /// \param coefs coefficients of the spline space
    /// \param mult_interior multiplicity at the interior knots
    /// \param periodic specifies whether the B-spline is periodic
    ///
    /// A periodic curve on the basis with (interior+1)*mult_interior
    /// functions stores that many rows; \a coefs may also carry the ghost rows
    /// (n+m rows, see gsBSpline(const Basis&, gsMatrix<T>)). Other row counts
    /// are refused with GISMO_ERROR.
    ///
    /// \ingroup Nurbs
    gsBSpline(T u0, T u1, unsigned interior, int degree,
              gsMatrix<T> coefs, unsigned mult_interior=1, bool periodic = false)
    {
        this->m_basis = new Basis(u0, u1, interior, degree, mult_interior, periodic );
        if( periodic )
        {
            this->m_coefs = _reducePeriodicCoefs(this->basis(), give(coefs));
        }
        else // non-periodic
        {
            this->m_coefs.swap(coefs);
            GISMO_ASSERT( this->m_coefs.rows() == this->basis().size(),
                          "Number of coefficients does not match the size of the basis.");
        }
    }

    GISMO_CLONE_FUNCTION(gsBSpline)

    GISMO_BASIS_ACCESSORS

public:

    /// Prints the object as a string. On a periodic curve the control points
    /// shown are the first and the last of the n stored rows.
    std::ostream &print(std::ostream &os) const override
    {
        os << "BSpline curve "<< "of degree "<< this->basis().degree()<< ", "<<  this->basis().knots() <<".\n";
        os << "with control points "<< this->m_coefs.row(0)<<" ... "<<this->m_coefs.bottomRows(1) << ".\n";
        if( this->basis().isPeriodic() )
            os << "Periodic with overlay " << this->basis().numCrossingFunctions() << ".\n";
        return os;
    }

    /// Returns the starting value of the domain of the basis
    T domainStart() const { return *this->basis().knots().domainBegin(); }

    /// Returns the end value of the domain of the basis
    T domainEnd() const { return *this->basis().knots().domainEnd(); }

    /// Returns a reference to the knot vector
    KnotVectorType & knots(const int i = 0)
    {
        GISMO_UNUSED(i);
        GISMO_ASSERT( i==0, "Requested knots of invalid direction "<< i );
        return this->basis().knots();
    }

    /// Returns a (const )reference to the knot vector
    const KnotVectorType & knots(const int i = 0) const
    {
        GISMO_UNUSED(i);
        GISMO_ASSERT( i==0, "Requested knots of invalid direction "<< i );
        return this->basis().knots();
    }

    /// Returns the degree of the B-spline
    short_t degree(short_t i = 0) const
    {
        GISMO_UNUSED(i);
        GISMO_ASSERT( i==0, "Requested knots of invalid direction "<< i );
        return this->basis().degree();
    }

    // compatible curves: same degree, same first/last p+1 knots
    void isCompatible( gsGeometry<T> * )
    { GISMO_NO_IMPLEMENTATION }

    // compatible curves: same degree, same first/last p+1 knots
    void makeCompatible( gsGeometry<T> * )
    { GISMO_NO_IMPLEMENTATION }


    /// Merge other B-spline into this one. Refused (GISMO_ERROR, curve
    /// unchanged) if this curve or \a other is periodic: a periodic curve has no
    /// end point to join.
    void merge( gsGeometry<T> * otherG ) override;

    /// Segment this BSpline curve between u0 and u1. Either of these values
    /// can be outside the parameter range of the curve, but u1 must be
    /// greater than u0. The degree of the curve is not modified.
    /// Refused (GISMO_ERROR) on a periodic curve; use setPeriodic(false) first.
    /// \param u0 starting parameter
    /// \param u1 end parameter
    /// \param tolerance proximity of the segment boundaries and B-spline knots to treat them as equal
    gsBSpline<T> segmentFromTo(T u0, T u1, T tolerance=1e-15) const;

    /// Split this BSpline curve at u0. This value must be inside the
    /// parameter range of the curve. Refused (GISMO_ERROR) on a periodic curve.
    /// \param u0 split parameter
    /// \param left curve segment from the left of parameter domain to u0
    /// \param right curve segment from u0 to the right of parameter domain
    /// \param tolerance proximity of the segment boundaries and B-spline knots to treat them as equal
    void splitAt(T u0, gsBSpline<T>& left,  gsBSpline<T>& right, T tolerance=1e-15) const;

    /// Convert BSpline curve into Bezier segments. Refused (GISMO_ERROR) on a
    /// periodic curve.
    /// \param tolerance proximity of the segment boundaries and B-spline knots to treat them as equal
    /// @note: parameter range do NOT change, extra knots().affineTransform() needed
    gsMultiPatch<T> toBezier(T tolerance=1e-15) const;

    /// Compute all intersections with another B-Spline curve. Refused
    /// (GISMO_ERROR) if either curve is periodic.
    /// \param other the other B-Spline curve
    /// \param tolerance the tolerance for computing the intersection
    /// \param curvatureTolerance the curvature tolerance for the splitting stage
    /// @note curvatureTolerance should be slightly greater than 1.0, e.g., 1+1e-6
    /// \return a vector of intersection results
    std::vector<internal::gsCurveIntersectionResult<T>> intersect(const gsBSpline<T>& other, T tolerance = 1e-5, T curvatureTolerance=1+1e-6) const;

    /// calculates curvature approximation for a B-spline curve by dividing
    /// the total length of the polyline formed by its control points by
    /// the direct distance between the first and last control point.
    /// On a periodic curve this measures the open polygon c_0..c_{n-1}
    /// of the n stored rows.
    T pseudoCurvature() const;

    /// Insert the given new knot (multiplicity \a i) without changing
    /// the curve.
    ///
    /// On a periodic curve the knot is handled by
    /// gsBSplineBasis::refine_withCoefs(): a value in (a,b) is inserted, a value
    /// in (a-L,a) or (b,b+L) is inserted as its translate by +-L (L=b-a), and
    /// the seam values a, b, anything farther, or a multiplicity above
    /// degree+1 raise GISMO_ERROR. The curve keeps its n rows and is
    /// unchanged if the insertion is refused. Cost O(K + n p d) for K
    /// stretched knots and d coordinates.
    void insertKnot( T knot, index_t i = 1);
private:
    // Resolve hidden overload w.r.t. gsGeometry
    virtual void insertKnot( T knot, index_t dir, index_t i = 1) override;
public:

    /// Insert the given new knots in the range \a [\em inBegin .. \em inEend )
    /// without changing the curve.
    ///
    /// On a periodic curve the knots are accepted and refused as in
    /// insertKnot(T,index_t), jointly: nothing changes if any of them is
    /// refused.
    /// \param inBegin iterator pointing to the first knot to be inserted
    /// \param inEnd iterator pointing to one position after the last
    /// knot to be inserted
    template <class It>
    void insertKnots( It inBegin, It inEnd)
    {
        if( this->basis().isPeriodic() )
        {
            const std::vector<T> v(inBegin, inEnd);
            Basis nb(this->basis());
            gsMatrix<T> c = m_coefs;
            nb.refine_withCoefs(c, v);
            this->basis().swap(nb);
            m_coefs.swap(c);
        }
        else // non-periodic
            gsBoehmRefine( this->basis().knots(), this->coefs(), this->degree(), inBegin, inEnd);
    }

    /// @brief Elevates the degree by \a i, keeping the curve.
    ///
    /// On a periodic curve the construction is exact in exact arithmetic
    /// (no solver): the stretched curve is unwrapped onto a clamped
    /// non-periodic pad, elevated by the non-periodic algorithm and the
    /// window of functions active on [a,b] is selected (see the definition
    /// for the argument). The floating-point accuracy is that of
    /// bspline::degreeElevateBSpline: round-off of order
    /// (p+i)^2 eps max|c| on uniform and quasi-uniform knot vectors,
    /// independent of the number of knots. The curve
    /// keeps n' = basis size rows. Cost
    /// O((n+m+p) p i d).
    void degreeElevate(short_t const i = 1, short_t const dir = -1) override;

    /// @brief Increases the degree by \a i, in general only approximately.
    ///
    /// Non-periodic: as gsGeometry::degreeIncrease(). Periodic: the new
    /// curve interpolates the old one at the n' Greville anchors of the
    /// new basis (a direct sparse LU of the n'xn' collocation matrix,
    /// checked by info() and by a backward-error residual) and is in
    /// general not equal to it. Refused (GISMO_ERROR, curve unchanged) if
    /// the basis cannot be increased, e.g. n < p+i+1, or if the
    /// system is numerically singular. Cost O(n' p^2) for the LU.
    void degreeIncrease(short_t const i = 1, short_t const dir = -1) override;

    /// @brief Decreases the degree by \a i, in general only approximately.
    ///
    /// Non-periodic: as gsGeometry::degreeDecrease(). Periodic: anchor
    /// interpolation as in degreeIncrease(), with the same complexity
    /// (the LU of the n'xn' collocation matrix, O(n' p^2)); refused
    /// (GISMO_ERROR, curve unchanged) if the number of crossing
    /// functions m <= i.
    void degreeDecrease(short_t const i = 1, short_t const dir = -1) override;

    /// @brief Reduces the degree by \a i, in general only approximately.
    ///
    /// Non-periodic: as gsGeometry::degreeReduce(). Periodic: anchor
    /// interpolation as in degreeIncrease(), with the same complexity
    /// (the LU of the n'xn' collocation matrix, O(n' p^2)); refused
    /// (GISMO_ERROR, curve unchanged) if the reduced basis would have
    /// fewer than degree+1 functions.
    void degreeReduce(short_t const i = 1, short_t const dir = -1) override;

    /// @brief Reverses the parameterization, keeping the domain [a,b].
    ///
    /// Non-periodic: the coefficient rows are reversed. Periodic: the curve
    /// keeps its n rows, with c'_j = c_{(m-1-j) mod n}.
    void reverse() override;

    /// @brief Returns true iff the point p is contained (approximately) on
    /// the curve, with the given tolerance.
    // Under Construction..( TO DO )
    bool contains( gsMatrix<T> const & p, T const & tol = 1e-6 )
        {
            assert( p.cols()==1 );
            gsBSplineSolver<T> slv;
            std::vector<T> roots;
            short_t dim = this->geoDim();
            gsMatrix<T> ev, tmp(1,1);
            int i(1);

            slv.allRoots( *this, roots, 0, p(0,0), 1e-6, 100);
            for ( typename std::vector<T>::const_iterator it=roots.begin();
                  it != roots.end(); ++it)
            {
                tmp(0,0)= *it;
                this->eval_into(tmp,ev);
                for (i=1; i!=dim; ++i )
                    if ( math::abs( ev(i,0) - p(i,0) ) > tol )
                        break;
                if (i==dim)
                    return true;
            }
            return false;
        }

    /// Sample \a npoints uniformly distributed (in parameter domain) points on the curve.
    gsMatrix<T> sample(int npoints = 50) const
    {
        gsMatrix<T> images;
        gsMatrix<T> interval = this->parameterRange();
        const gsMatrix<T> pts = gsPointGrid(interval, npoints );
        this->eval_into( pts, images );
        return images;
    }

    /// @brief Converts the curve to periodic (\a flag true) or back.
    ///
    /// flag == true: the knot vector is converted by
    /// gsKnotVector::setPeriodic() and the n+m coefficient rows are
    /// compressed to n; the ghost rows n+j must equal their twins j
    /// within 100*eps*max|coef|, otherwise GISMO_ERROR. A refused
    /// conversion leaves the curve unchanged. No-op if already periodic.
    ///
    /// flag == false: exact unwrapping onto the stretched knots with the
    /// n+m expanded rows (c_{i mod n} for stretched function i), the same
    /// curve in a non-periodic basis. No-op if not periodic.
    void setPeriodic(bool flag = true)
    {
        if( flag == this->basis().isPeriodic() )
            return;

        if( flag )
        {
            Basis nb(this->basis());
            nb.setPeriodic(true);
            gsMatrix<T> c = _reducePeriodicCoefs(nb, m_coefs);
            this->basis().swap(nb);
            m_coefs.swap(c);
        }
        else
        {
            gsMatrix<T> c = m_coefs;
            this->basis().expandCoefs(c);
            this->basis().setPeriodic(false);
            m_coefs.swap(c);
        }
    }

    /// Number of coefficient rows (the basis size)
    inline int numCoefs() const { return this->m_coefs.rows(); }

    inline bool isPeriodic() const { return this->basis().isPeriodic(); }


    /// Returns true if this curve is closed
    bool isClosed(T tol = 1e-10) const
    {
        return ( this->basis().isPeriodic() ||
                 (m_coefs.row(0) - m_coefs.row(m_coefs.rows()-1)).squaredNorm()<tol );
    }

    /// \brief Return true if point \a u is on the curve with
    /// tolerance \a tol
    bool isOn(gsMatrix<T> const &u, T tol = 1e-3) const override;

    /// \brief Return true if point \a u is an endpoint (corner) of
    /// the curve with tolerance \a tol. Refused (GISMO_ERROR) on a periodic
    /// curve, which has no end points.
    bool isPatchCorner(gsMatrix<T> const &v, T tol = 1e-3) const;

    /// \brief returns the tensor-index \a curr of the corner control
    /// point \a v, or an invalid index if the corner is not found
    /// within the tolerance \a tol. Refused (GISMO_ERROR) on a periodic
    /// curve, which has no corners.
    void findCorner(const gsMatrix<T>   & v,
                    gsVector<index_t,1> & curr,
                    T tol = 1e-3);

    /// \brief Modifies the parameterization such that the point \a v
    /// is the starting point of the curve. Assumes that \a v is
    /// either the starting or the ending point of the curve. Refused
    /// (GISMO_ERROR) on a periodic curve.
    void setOriginCorner(gsMatrix<T> const &v);

    /// \brief Modifies the parameterization such that the point \a v
    /// is the ending point of the curve. Assumes that \a v is
    /// either the starting or the ending point of the curve. Refused
    /// (GISMO_ERROR) on a periodic curve.
    void setFurthestCorner(gsMatrix<T> const &v);

    void swapDirections(const unsigned i, const unsigned j);

private:

    // Avoid hidden overloads w.r.t. gsGeometry
    using Base::insertKnot;

    /// @brief Coefficients of a curve on basis \a b, reduced to b.size() rows
    /// if \a b is periodic.
    ///
    /// Non-periodic \a b, or b.size() rows: \a coefs is returned unchanged.
    /// b.size()+m rows (m = b.numCrossingFunctions()): the ghost rows n+j
    /// must equal their twins j within 100*eps*max|coefs| and are then
    /// dropped; GISMO_ERROR otherwise. Any other row count: GISMO_ERROR.
    /// Cost O((n+m) d).
    static gsMatrix<T> _reducePeriodicCoefs(const Basis & b, gsMatrix<T> coefs)
    {
        if( !b.isPeriodic() )
            return coefs;

        const index_t n = b.size();
        const index_t m = b.numCrossingFunctions();
        if( coefs.rows() == n )
            return coefs;

        GISMO_ENSURE( coefs.rows() == n + m,
                      "gsBSpline: a periodic basis with n = " << n << " functions and m = "
                      << m << " crossing functions needs n or n+m coefficient rows, got "
                      << coefs.rows() << "." );

        const T twinTol = 100 * std::numeric_limits<T>::epsilon();
        const T scale = coefs.size() == 0 ? T(0) : coefs.cwiseAbs().maxCoeff();
        for( index_t j = 0; j < m; ++j )
        {
            const T diff = (coefs.row(n+j) - coefs.row(j)).cwiseAbs().maxCoeff();
            GISMO_ENSURE( diff <= twinTol * scale,
                          "gsBSpline: ghost coefficient row " << n+j << " differs from its twin row "
                          << j << " by " << diff << "." );
        }
        coefs.conservativeResize(n, gsEigen::NoChange);
        return coefs;
    }

    /// Replaces the basis by the periodic basis \a nb and the coefficients by
    /// the interpolant of the current curve at the n' anchors of \a nb.
    /// Everything is computed on locals; a throw leaves the curve unchanged.
    void _periodicDegreeChange(Basis & nb);

protected:

    using Base::m_coefs;

    // TODO Check function
    // check function: check the coefficient number, degree, knot vector ...
}; // class gsBSpline


#ifdef GISMO_WITH_PYBIND11

  /**
   * @brief Initializes the Python wrapper for the class: gsBSpline
   */
  void pybind11_init_gsBSpline(pybind11::module &m);

#endif // GISMO_WITH_PYBIND11


/*
// Product of spline functions
template<class T>
gsBSpline<T> operator*(const gsBSpline<T> & lhs, const gsBSpline<T> & rhs)
{
    // Note: quick-fix implementation. Better algorithms  exist
    gsKnotVector<T> kv = lhs.knots().knotUnion( lhs.knots() );
    kv.degreeElevate( lhs.degree() + rhs.degree() - kv.degree() );

    gsMatrix<T> pts;
    kv.greville_into(pts);
    const gsMatrix<T> ev = (lhs.eval(pts)).array() * (rhs.eval(pts)).array();

    // fixme: avoid temporaries here
    return *(static_cast<gsBSpline<T>*>(gsBSplineBasis<T>(kv).interpolateData(ev,pts)));
}
*/

} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsBSpline.hpp)
#endif
