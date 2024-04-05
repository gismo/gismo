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
    // See also gsNurbs::merge().
    // check geometric dimension
    GISMO_ASSERT(this->geoDim()==otherG->geoDim(),
                 "gsNurbs: cannot merge curves in different spaces ( R^"
                 << this->geoDim() << ", R^" << otherG->geoDim() << " ).");

    // check if the type of other is BSpline
    gsBSpline *  other = dynamic_cast<gsBSpline *>( otherG );
    GISMO_ASSERT( other!=NULL, "Can only merge with B-spline curves.");
    other= other->clone().release();

    GISMO_ASSERT( this ->basis().isPeriodic() == false &&
                  other->basis().isPeriodic() == false,
                  "Cannot merge a closed curve with anything." );

    // check degree
    const int mDeg = this ->basis().degree();
    const int oDeg = other->basis().degree();
    const int deg  = math::max(mDeg,oDeg);

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
    int n= this->coefsSize();
    int skip = continuous ? 1 : 0;
    this->m_coefs.conservativeResize( n + other->coefsSize() -skip, gsEigen::NoChange ) ;

    this->m_coefs.block( n,0,other->coefsSize()-skip,other->geoDim() ) =
        other->m_coefs.block( 1,0,other->coefsSize()-skip,other->geoDim() ) ;

    delete other;
}

template<class T>
gsBSpline<T> gsBSpline<T>::segmentFromTo(T u0, T u1, T tolerance) const
{

  GISMO_ASSERT(domainStart()-tolerance < u0 && u0 < domainEnd()+tolerance,
               "starting point "<< u0 <<" not in the knot vector");
  GISMO_ASSERT(domainStart()-tolerance < u1 && u1 < domainEnd()+tolerance,
               "end point "<< u1 <<" not in the knot vector");

  //First make a copy of the actual object, to allow const
  gsBSpline<T> copy(*this);

  // Extract a reference to the knots, the basis and coefs of the copy
  KnotVectorType & knots = copy.basis().knots();

  // some constants
  const int p = basis().degree();                         // degree
  const index_t multStart = p + 1 - knots.multiplicity(u0); // multiplicity
  const index_t multEnd   = p + 1 - knots.multiplicity(u1);   // multiplicity

  // insert the knot, such that its multiplicity is p+1
  if (multStart>0) { copy.insertKnot(u0, 0, multStart); }
  if (multEnd>0)   { copy.insertKnot(u1, 0, multEnd  ); }

  gsMatrix<T>& coefs = copy.coefs();
  const index_t tDim  = coefs.cols();

  // Split the coefficients
  // find the number of coefs left from u0
  const index_t nL  = knots.uFind(u0).firstAppearance();
  // find the number of coefs left from u1
  index_t nL2 = knots.uFind(u1).firstAppearance();
  bool isEnd = math::abs(u1 - this->domainEnd()) < tolerance;
//  if ( isEnd ) { nL2 += 1; }       // Adjust for end parameter
  if ( isEnd ) { nL2 = copy.numCoefs(); }       // Adjust for end parameter

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
gsMultiPatch<T> gsBSpline<T>::toBezier(T tolerance) const {
  gsMultiPatch<T> bezierSegments;

  gsBSpline<T> currentSegment(*this);
  gsBSpline<T> leftPart;

  for (auto iter = this->knots().ubegin() + 1; iter != this->knots().uend() - 1; ++iter) {
    currentSegment.splitAt(*iter, leftPart, currentSegment, tolerance);
    bezierSegments.addPatch(leftPart);
  }

  bezierSegments.addPatch(currentSegment); // Add the last segment
  return bezierSegments;
}

template<class T>
std::vector<internal::gsCurveIntersectionResult<T>> gsBSpline<T>::intersect(const gsBSpline<T>& other,
                                                    T tolerance) const {
  std::vector<internal::gsBoundingBoxPair<T>> hulls = internal::getPotentialIntersectionRanges<T>(*this, other);

  std::vector<internal::gsCurveIntersectionResult<T>> results;
  for (const auto &hull : hulls) {
    gsBSpline<T> crv1 = this->segmentFromTo(hull.b1.getRange().getMin(), hull.b1.getRange().getMax());
    gsBSpline<T> crv2 = other.segmentFromTo(hull.b2.getRange().getMin(), hull.b2.getRange().getMax());

    internal::gsCurveCurveDistanceSystem<T> obj(crv1, crv2);
    gsMatrix<T, 2, 1> uv;
    uv(0, 0) = 0.5 * (crv1.domainStart() + crv1.domainEnd());
    uv(1, 0) = 0.5 * (crv2.domainStart() + crv2.domainEnd());
    T distance = obj.compute(uv, tolerance);

    if (distance < math::max((T)1e-10, tolerance)) {
      internal::gsCurveIntersectionResult<T> result(uv(0), uv(1), 0.5 * (crv1.eval(uv.row(0)) + crv2.eval(uv.row(1))));
      results.push_back(result);
    }
  }

  return results;
}

template<class T>
T gsBSpline<T>::pseudoCurvature() const {
  int coefsSize = m_coefs.rows();

  T len = (m_coefs.row(0)-m_coefs.row(coefsSize-1)).norm();
  T total = 0.0;
  for (int ipt = 0; ipt != coefsSize - 1; ++ipt) {
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
        int borderKnotMult = this->basis().borderKnotMult();
        KnotVectorType & knots = this->knots();
        unsigned deg = this->basis().degree();

        GISMO_ASSERT( knot != knots[deg] && knot != knots[knots.size() - deg - 1],
                      "You are trying to increase the multiplicity of the p+1st knot but the code is not ready for that.\n");

        // If we would be inserting to "passive" regions, we
        // rather insert the knot into the mirrored part.
        // Adjustment of the mirrored knots is then desirable.
        if( knot < knots[deg - borderKnotMult + 1] )
        {
            knot += this->basis()._activeLength();
        }
        else if( knot > knots[knots.size() - deg + borderKnotMult - 2] )
        {
            knot -= this->basis()._activeLength();
        }
        // If necessary, we update the mirrored part of the knot vector.
        if((knot < knots[2*deg + 1 - borderKnotMult]) || (knot >= knots[knots.size() - 2*deg - 2 + borderKnotMult]))
            this->basis().enforceOuterKnotsPeriodic();

        // We copy some of the control points to pretend for a while
        // that the basis is not periodic.

        //gsMatrix<T> trueCoefs = this->basis().perCoefs( this->coefs() );
        gsBoehm( this->basis().knots(), this->coefs(), knot, i );
        //this->coefs() = trueCoefs;
        //this->coefs().conservativeResize( this->basis().size(), this->coefs().cols() );
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
    return (( v - m_coefs.row(0)       ).squaredNorm() < tol ||
            ( v - m_coefs.bottomRows(1)).squaredNorm() < tol );
}

template<class T>    
void gsBSpline<T>::findCorner(const gsMatrix<T> & v, 
                                             gsVector<index_t,1> & curr,
                                             T tol)
{
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
    
    bspline::degreeElevateBSpline(this->basis(), this->m_coefs, i);
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
