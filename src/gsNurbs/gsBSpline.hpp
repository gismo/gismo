/** @file gsBSpline.hpp

    @brief Implementation of a B-spline curve/function with one parameter

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/


#pragma once 

#include <gsNurbs/gsBSplineAlgorithms.h>
#include <gsIO/gsXml.h>
#include <gsIO/gsXmlGenericUtils.hpp>

#include <gsCore/gsMultiPatch.h>
#include <gsNurbs/gsCurveCurveIntersection.h>

namespace gismo
{


template<class T>
void gsBSpline<T>::merge( gsGeometry<T> * otherG )
{
    GISMO_ENSURE( !this->basis().isPeriodic(),
                  "gsBSpline: merge is not defined on a periodic curve (it has no end point to join)." );
    // See also gsNurbs::merge().
    // check geometric dimension
    GISMO_ASSERT(this->geoDim()==otherG->geoDim(),
                 "gsNurbs: cannot merge curves in different spaces ( R^"
                 << this->geoDim() << ", R^" << otherG->geoDim() << " ).");

    // check if the type of other is BSpline
    gsBSpline *  other = dynamic_cast<gsBSpline *>( otherG );
    GISMO_ASSERT( other!=NULL, "Can only merge with B-spline curves.");
    GISMO_ENSURE( !other->basis().isPeriodic(),
                  "gsBSpline: cannot merge a periodic curve (it has no end point to join)." );
    other= other->clone().release();

    // check degree
    const short_t mDeg = this ->basis().degree();
    const short_t oDeg = other->basis().degree();
    const short_t deg  = math::max(mDeg,oDeg);

    other->gsBSpline::degreeElevate( deg - oDeg ); // degreeElevate(0) does nothing (and very quickly)
    this ->gsBSpline::degreeElevate( deg - mDeg );

    // check whether the resulting curve will be continuous
    // TODO: ideally, the tolerance should be a parameter of the function
    T tol = 1e-8;
    gsMatrix<T> mValue = this ->eval(this ->support().col(1));
    gsMatrix<T> oValue = other->eval(other->support().col(0));
    bool continuous = gsAllCloseAbsolute(mValue,oValue,tol);

    // merge knot vectors.
    KnotVectorType& mKnots = this ->basis().knots();
    KnotVectorType& oKnots = other->basis().knots();
    T lastKnot = mKnots.last();
    if (continuous) // reduce knot multiplicity
    {
        // TODO check for clamped knot vectors otherwise
        // we should do knot insertion beforehands
        mKnots.remove(lastKnot);
    }// else there is a knot of multiplicity deg + 1 at the discontinuity

    oKnots.addConstant(lastKnot-oKnots.first());
    mKnots.append( oKnots.begin()+deg+1, oKnots.end());

    // merge coefficients
    unsigned n= this->coefsSize();
    index_t skip = continuous ? 1 : 0;
    this->m_coefs.conservativeResize( n + other->coefsSize() -skip, gsEigen::NoChange ) ;

    this->m_coefs.block( n,0,other->coefsSize()-skip,other->geoDim() ) =
        other->m_coefs.block( skip,0,other->coefsSize()-skip,other->geoDim() ) ;

    delete other;
}

template<class T>
gsBSpline<T> gsBSpline<T>::segmentFromTo(T u0, T u1, T tolerance) const
{
  GISMO_ENSURE(!this->basis().isPeriodic(),
               "gsBSpline: segmentFromTo is not defined on a periodic curve (it has no end points); use setPeriodic(false) first.");

  GISMO_ASSERT(domainStart()-tolerance < u0 && u0 < domainEnd()+tolerance,
               "starting point "<< u0 <<" not in the knot vector");
  GISMO_ASSERT(domainStart()-tolerance < u1 && u1 < domainEnd()+tolerance,
               "end point "<< u1 <<" not in the knot vector");

  //First make a copy of the actual object, to allow const
  gsBSpline<T> copy(*this);

  // Extract a reference to the knots, the basis and coefs of the copy
  KnotVectorType & knots = copy.basis().knots();

  // some constants
  const short_t p = basis().degree();                       // degree
  const index_t multStart = p + 1 - knots.multiplicity(u0); // multiplicity
  const index_t multEnd   = p + 1 - knots.multiplicity(u1); // multiplicity

  // insert the knot, such that its multiplicity is p+1
  if (multStart>0)
    copy.insertKnot(u0, 0, multStart);
  if (multEnd>0)
    copy.insertKnot(u1, 0, multEnd  );

  gsMatrix<T>& coefs = copy.coefs();
  const index_t tDim  = coefs.cols();

  // Split the coefficients
  // find the number of coefs left from u0
  const index_t nL  = knots.uFind(u0).firstAppearance();
  // find the number of coefs left from u1
  index_t nL2 = knots.uFind(u1).firstAppearance();
  bool isEnd = math::abs(u1 - this->domainEnd()) < tolerance;
//  if ( isEnd ) { nL2 += 1; }       // Adjust for end parameter
  if ( isEnd )
    nL2 = copy.numCoefs();       // Adjust for end parameter

  // Prepare control points for new geometry
  gsMatrix<T> coefRes = coefs.block(nL, 0, nL2-nL, tDim);

  // Prepare new knot vector
  typename KnotVectorType::iterator itStart = knots.iFind(u0);
  typename KnotVectorType::iterator itEnd = knots.iFind(u1) + (isEnd ? p + 1 : 0);
  typename KnotVectorType::knotContainer matRes(itStart-p, itEnd+1);
  KnotVectorType knotsRes(give(matRes), p);

  return gsBSpline<T>(Basis(give(knotsRes)), give(coefRes));
}

template<class T>
gsMultiPatch<T> gsBSpline<T>::toBezier(T tolerance) const
{
    GISMO_ENSURE(!this->basis().isPeriodic(),
                 "gsBSpline: toBezier is not defined on a periodic curve; use setPeriodic(false) first.");
    gsMultiPatch<T> bezierSegments;

    gsBSpline<T> currentSegment(*this);
    gsBSpline<T> leftPart;

    for (auto iter = this->knots().ubegin() + 1; iter != this->knots().uend() - 1; ++iter)
    {
        currentSegment.splitAt(*iter, leftPart, currentSegment, tolerance);
        bezierSegments.addPatch(leftPart);
    }

    bezierSegments.addPatch(currentSegment); // Add the last segment
    return bezierSegments;
}

template<class T>
std::vector<internal::gsCurveIntersectionResult<T>> gsBSpline<T>::intersect(const gsBSpline<T>& other,
                                                                            T tolerance,
                                                                            T curvatureTolerance) const
{
    GISMO_ENSURE(!this->basis().isPeriodic() && !other.basis().isPeriodic(),
                 "gsBSpline: intersect is not defined on a periodic curve; use setPeriodic(false) first.");

    std::vector<internal::gsBoundingBoxPair<T>> hulls = internal::getPotentialIntersectionRanges<T>(*this, other, curvatureTolerance);

    std::vector<internal::gsCurveIntersectionResult<T>> results;
    for (typename std::vector<internal::gsBoundingBoxPair<T>>::const_iterator hull = hulls.begin(); hull!=hulls.end(); hull++)
    {
        gsBSpline<T> crv1 = this->segmentFromTo(hull->b1.getRange().getMin(), hull->b1.getRange().getMax());
        gsBSpline<T> crv2 = other.segmentFromTo(hull->b2.getRange().getMin(), hull->b2.getRange().getMax());

        internal::gsCurveCurveDistanceSystem<T> obj(crv1, crv2);
		gsVector<T,2> uv;
        uv(0) = 0.5 * (crv1.domainStart() + crv1.domainEnd());
        uv(1) = 0.5 * (crv2.domainStart() + crv2.domainEnd());
        T distance = obj.compute(uv, tolerance);

        if (distance < math::max( (T)1e-10, tolerance))
        {
            gsMatrix<T> uCrv1(1,1), uCrv2(1,1);
			uCrv1(0,0) = uv(0);
			uCrv2(0,0) = uv(1);
		  	gsMatrix<T> pt = 0.5 * (crv1.eval(uCrv1) + crv2.eval(uCrv2));
            internal::gsCurveIntersectionResult<T> result(uv(0), uv(1), pt);
            results.push_back(result);
        }
    }

  return results;
}

template<class T>
T gsBSpline<T>::pseudoCurvature() const
{
    index_t coefsSize = m_coefs.rows();

    T len = (m_coefs.row(0)-m_coefs.row(coefsSize-1)).norm();
    T total = 0.0;
    for (index_t ipt = 0; ipt != coefsSize - 1; ++ipt)
    {
        T dist = (m_coefs.row(ipt) - m_coefs.row(ipt + 1)).norm();
        total += dist;
    }
    return total / len;
}

template<class T>
void gsBSpline<T>::insertKnot( T knot, index_t dir, index_t i)
{
    GISMO_UNUSED(dir);
    if (i==0) return;
    //if ( i==1)
    //single knot insertion: Boehm's algorithm
    //else
    //knot with multiplicity:   Oslo algorithm
    if( this->basis().isPeriodic() )
    {
        // Pad-and-fold refinement on copies: a refused insertion leaves the curve unchanged.
        Basis nb(this->basis());
        gsMatrix<T> c = m_coefs;
        nb.refine_withCoefs(c, std::vector<T>(i, knot));
        this->basis().swap(nb);
        m_coefs.swap(c);
    }
    else // non-periodic
        gsBoehm( this->basis().knots(), this->coefs() , knot, i);
}

template<class T>
void gsBSpline<T>::insertKnot( T knot, index_t i)
{
    insertKnot(knot,0,i);
}

template<class T>
void gsBSpline<T>::splitAt(T u0, gsBSpline<T>& left,  gsBSpline<T>& right, T tolerance) const
{
  GISMO_ENSURE(!this->basis().isPeriodic(),
               "gsBSpline: splitAt is not defined on a periodic curve (it has no end points); use setPeriodic(false) first.");
  GISMO_ASSERT(domainStart()-tolerance < u0 && u0 < domainEnd()+tolerance,
               "splitting point "<< u0 <<" not in the knot vector");

  left  = segmentFromTo(this->domainStart(), u0, tolerance);
  right = segmentFromTo(u0, this->domainEnd(), tolerance);
}

template<class T>    
bool gsBSpline<T>::isOn(gsMatrix<T> const &u, T tol) const
{
    GISMO_ASSERT( u.cols() == 1, "Expecting single point.");
    gsBSplineSolver<T> slv;
    gsMatrix<T> e;
    
    for ( index_t k = 0; k != u.rows(); ++k)
    {        
        std::vector<T> roots;
        slv.allRoots(*this, roots, k, u(k,0) ) ;

        if( roots.size()!=0 )
        {
            gsAsMatrix<T> xx(roots); 
            this->eval_into( xx, e ); 
            
            //if( math::abs( e(!k,0)-u(!k,0) ) < tol )
            for(index_t j=0; j!=e.cols(); j++)
                if( ( e.col(j)-u ).norm() < tol )
                    return true;
        }
    }
    return false;
}

template<class T>    
bool gsBSpline<T>::isPatchCorner(gsMatrix<T> const &v, T tol) const
{
    GISMO_ENSURE(!this->basis().isPeriodic(),
                 "gsBSpline: a periodic curve has no corners (no end points).");
    return (( v - m_coefs.row(0)       ).squaredNorm() < tol ||
            ( v - m_coefs.bottomRows(1)).squaredNorm() < tol );
}

template<class T>    
void gsBSpline<T>::findCorner(const gsMatrix<T> & v, 
                                             gsVector<index_t,1> & curr,
                                             T tol)
{
    GISMO_ENSURE(!this->basis().isPeriodic(),
                 "gsBSpline: a periodic curve has no corners (no end points).");
    if ((v - m_coefs.row(0)).squaredNorm() < tol)
        curr[0] = 0;
    else if ((v - m_coefs.bottomRows(1)).squaredNorm() < tol)
        curr[0] = m_coefs.rows()-1;
    else
    {
        curr[0] = m_coefs.rows(); // invalidate result
        gsWarn<<"Point "<< v <<" is not an corner of the patch. (Call isPatchCorner() first!).\n";
    }
}


template<class T>    
void gsBSpline<T>::setOriginCorner(gsMatrix<T> const &v)
{
    GISMO_ENSURE(!this->basis().isPeriodic(),
                 "gsBSpline: a periodic curve has no corners (no end points).");
    if ((v - m_coefs.row(0)).squaredNorm() < (T)(1e-3))
        return;
    else if ((v - m_coefs.bottomRows(1)).squaredNorm() < (T)(1e-3))
        this->reverse();
    else
        gsWarn<<"Point "<< v <<" is not an endpoint of the curve.\n";
}

template<class T>    
void gsBSpline<T>::setFurthestCorner(gsMatrix<T> const &v)
{
    GISMO_ENSURE(!this->basis().isPeriodic(),
                 "gsBSpline: a periodic curve has no corners (no end points).");
    if ((v - m_coefs.bottomRows(1)).squaredNorm() < (T)(1e-3))
        return;
    else if ((v - m_coefs.row(0)).squaredNorm() < (T)(1e-3))
        this->reverse();
    else
        gsWarn<<"Point "<< v <<" is not an endpoint of the curve.\n";
}

template <class T>
void gsBSpline<T>::swapDirections(const unsigned i, const unsigned j)
{
    GISMO_UNUSED(i);
    GISMO_UNUSED(j);
    GISMO_ASSERT( static_cast<int>(i) == 0 && static_cast<int>(j) == 0,
                  "Invalid basis components "<<i<<" and "<<j<<" requested" );
}


template<class T>
void gsBSpline<T>::degreeElevate(short_t const i, short_t const dir)
{
    GISMO_UNUSED(dir);
    GISMO_ASSERT( (dir == -1) || (dir == 0),
                  "Invalid basis component "<< dir <<" requested for degree elevation" );

    if (this->basis().isPeriodic())
    {
        if (i == 0)
            return;

        // Periodic degree elevation by unwrap, elevate, select. The construction
        // is exact in exact arithmetic; in floating point its accuracy is that
        // of bspline::degreeElevateBSpline, which loses accuracy as the number
        // of knots grows.
        //
        // Let t_0..t_{K-1} be the stretched knots, p the degree, [a,b] the domain,
        // L = b-a, and c_{i mod n} the coefficient of stretched function i
        // (0 <= i < n+m). The curve is the spline sum_i c_{i mod n} B_i on the
        // stretched knots.
        //  1. Pad: the clamped knot vector tau = [c^(p+1), t_0..t_{K-1}, e^(p+1)]
        //     (c < t_0, e > t_{K-1}) contains every stretched B-spline B_i as
        //     the function p+1+i of tau, on the same knots t_i..t_{i+p+1}. The
        //     p+1 functions added on each side get coefficient zero; they are
        //     supported in [c,a] (left pad j <= p ends at t_j <= t_p = a) or in
        //     [b,e] (right pads start at t_{n+m+r} >= b). Hence the padded
        //     curve equals the periodic curve on [a,b].
        //  2. Elevate the padded clamped curve by the non-periodic algorithm,
        //     which represents the same curve exactly in degree p+i.
        //  3. Select: the new periodic basis has stretched knots t'_0..t'_{K'-1}
        //     that raise the multiplicity of every distinct knot t by i, which
        //     is invariant under translation by L. The functions of the
        //     elevated padded basis that are nonzero on [a,b] are exactly the n'+m'
        //     functions o..o+n'+m'-1, o = f-m', f the index of the first knot
        //     equal to a (the first span after the last copy of a has
        //     active functions f+mu'-1-p' = f-m'). Their knots are the new
        //     stretched knots. They are linearly independent on [a,b], so the
        //     coefficients of the curve in them are unique. The periodic curve
        //     written in the elevated periodic space, which contains it, is one
        //     such representation; the selected rows are therefore the new
        //     periodic coefficients and the ghost rows equal their twins (in
        //     exact arithmetic).
        const short_t p = this->basis().degree();
        const index_t n = this->basis().size();
        const index_t m = this->basis().numCrossingFunctions();
        const index_t d = m_coefs.cols();
        const T a = this->basis().domainStart();
        const T L = this->basis().domainEnd() - a;

        Basis nb(this->basis());
        nb.degreeElevate(i);
        const index_t n1 = nb.size();
        const index_t m1 = nb.numCrossingFunctions();
        GISMO_ASSERT( m1 == m, "Degree elevation changed the number of crossing functions." );

        const KnotVectorType & kv = this->basis().knots();
        std::vector<T> tau;
        tau.reserve(kv.size() + 2*(p+1));
        tau.insert(tau.end(), p+1, *kv.begin() - L);
        tau.insert(tau.end(), kv.begin(), kv.end());
        tau.insert(tau.end(), p+1, *(kv.end()-1) + L);
        Basis pb( KnotVectorType(p, tau.begin(), tau.end()) );

        gsMatrix<T> ex = m_coefs;
        this->basis().expandCoefs(ex);
        gsMatrix<T> P = gsMatrix<T>::Zero(n+m+2*p+2, d);
        P.middleRows(p+1, n+m) = ex;

        bspline::degreeElevateBSpline(pb, P, i);

        const KnotVectorType & kp = pb.knots();
        const index_t f = std::lower_bound(kp.begin(), kp.end(), a) - kp.begin();
        const index_t o = f - m1;
        GISMO_ASSERT( o >= 0 && o + n1 + m1 <= P.rows(),
                      "Degree elevation window out of range: offset " << o << ", rows " << P.rows() << "." );

#ifndef NDEBUG
        {
            const T b = a + L;
            const T knotTol = 50 * std::numeric_limits<T>::epsilon() * math::max(math::max(math::abs(a), math::abs(b)), L);
            const index_t K1 = n1 + m1 + p + i + 1;
            GISMO_ASSERT( static_cast<index_t>(nb.knots().size()) == K1 && o + K1 <= static_cast<index_t>(kp.size()),
                          "Degree elevation: unexpected knot count." );
            for (index_t k = 0; k < K1; ++k)
                GISMO_ASSERT( math::abs( *(kp.begin()+o+k) - *(nb.knots().begin()+k) ) <= knotTol,
                              "Degree elevation: window knot " << k << " differs from the new stretched knot." );
        }
#endif

        gsMatrix<T> Cnew = P.middleRows(o, n1);
        this->basis().swap(nb);
        m_coefs.swap(Cnew);
        return;
    }

    if (knots().numLeftGhosts() != 0 || knots().numRightGhosts() != 0)
        gsGeometry<T>::degreeElevate(i);
    else
        bspline::degreeElevateBSpline(this->basis(), this->m_coefs, i);
}

template<class T>
void gsBSpline<T>::_periodicDegreeChange(Basis & nb)
{
    // The new curve is the interpolant, in the new periodic space, of the
    // current curve at the n' Greville anchors of the new basis. It equals the
    // old curve only if the old curve lies in the new space (it does not when
    // the degree goes down or the continuity at the seam goes up); in general
    // it is an approximation that agrees with the old curve at the anchors.
    //
    // A direct sparse LU is used, not gsBasis::interpolateData (an iterative
    // solver without convergence check). The LU does not fail on a nearly
    // singular system, so both the factorization status and the residual are
    // checked. Periodic Greville interpolation is well posed for simple
    // uniform knots (Schoenberg's cardinal interpolation) but is not proven
    // nonsingular in general or for repeated knots; anchors that coincide are
    // only perturbed by about 1e-9 by gsKnotVector::greville_into(), which makes
    // the matrix ill-conditioned rather than singular.
    //
    // The collocation matrix is cyclic banded, not banded: an anchor near b
    // sees the crossing functions, whose stretched indices n..n+m-1 fold to
    // the periodic indices 0..m-1, so its row has entries in the first
    // columns (and symmetrically for an anchor near a).
    //
    // Cost O(n' p) to assemble, O(n' p^2) for the LU (fill in the wrap-around corners).
    gsMatrix<T> A;
    nb.anchors_into(A);
    gsMatrix<T> v;
    this->eval_into(A, v);
    gsSparseMatrix<T> C = nb.collocationMatrix(A);
    C.makeCompressed();

    typename gsSparseSolver<T>::LU solver;
    solver.compute(C);
    GISMO_ENSURE( solver.succeed(),
                  "gsBSpline: LU factorization of the anchor interpolation system failed (info "
                  << solver.info() << ")." );
    gsMatrix<T> X = solver.solve(v.transpose());

    // Backward-error test: a direct LU guarantees ||C x - v|| = O(eps (||C|| ||x|| + ||v||)),
    // not a small residual relative to ||v|| alone.
    const T resTol = 1000 * std::numeric_limits<T>::epsilon();
    const T res = (C * X - v.transpose()).cwiseAbs().maxCoeff();
    const T normC = C.coeffs().cwiseAbs().maxCoeff();
    const T bound = resTol * (normC * X.cwiseAbs().maxCoeff() + v.cwiseAbs().maxCoeff());
    GISMO_ENSURE( res <= bound,
                  "gsBSpline: anchor interpolation system is numerically singular, residual "
                  << res << " exceeds " << bound << "." );

    this->basis().swap(nb);
    m_coefs.swap(X);
}

template<class T>
void gsBSpline<T>::degreeIncrease(short_t const i, short_t const dir)
{
    GISMO_ASSERT( (dir == -1) || (dir == 0),
                  "Invalid basis component "<< dir <<" requested for degree increase" );
    if (!this->basis().isPeriodic())
    {
        gsGeometry<T>::degreeIncrease(i, dir);
        return;
    }
    Basis nb(this->basis());
    nb.degreeIncrease(i);
    _periodicDegreeChange(nb);
}

template<class T>
void gsBSpline<T>::degreeDecrease(short_t const i, short_t const dir)
{
    GISMO_ASSERT( (dir == -1) || (dir == 0),
                  "Invalid basis component "<< dir <<" requested for degree decrease" );
    if (!this->basis().isPeriodic())
    {
        gsGeometry<T>::degreeDecrease(i, dir);
        return;
    }
    Basis nb(this->basis());
    nb.degreeDecrease(i);
    _periodicDegreeChange(nb);
}

template<class T>
void gsBSpline<T>::degreeReduce(short_t const i, short_t const dir)
{
    GISMO_ASSERT( (dir == -1) || (dir == 0),
                  "Invalid basis component "<< dir <<" requested for degree reduction" );
    if (!this->basis().isPeriodic())
    {
        gsGeometry<T>::degreeReduce(i, dir);
        return;
    }
    Basis nb(this->basis());
    nb.degreeReduce(i);
    _periodicDegreeChange(nb);
}

template<class T>
void gsBSpline<T>::reverse()
{
    if (!this->basis().isPeriodic())
    {
        gsCurve<T>::reverse();
        return;
    }

    // The reversed stretched function j is the old stretched function
    // n+m-1-j, whose coefficient is c_{(n+m-1-j) mod n} = c_{(m-1-j) mod n}.
    // A plain row reversal would give the right curve cyclically shifted by m.
    const index_t n = this->basis().size();
    const index_t m = this->basis().numCrossingFunctions();
    gsMatrix<T> c(n, m_coefs.cols());
    for (index_t j = 0; j < n; ++j)
        c.row(j) = m_coefs.row( ((m-1-j) % n + n) % n );
    this->basis().reverse();
    m_coefs.swap(c);
}

namespace internal
{

/// Get a BSpline from XML data
template<class T>
class gsXml< gsBSpline<T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsBSpline<T>);
    static std::string tag () { return "Geometry"; }
    static std::string type () { return "BSpline"; }

    static gsBSpline<T> * get (gsXmlNode * node)
    {
        return getGeometryFromXml< gsBSpline<T> >(node);
    }

    static gsXmlNode * put (const gsBSpline<T> & obj,
                            gsXmlTree & data )
    {
        return putGeometryToXml< gsBSpline<T> >(obj,data);
    }
};


}// namespace internal

} // namespace gismo
