/** @file gsCurveCurveIntersection.h

    @brief Implementation of B Spline Curve/Curve intersection

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Ye Ji
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>

namespace gismo{

namespace internal {

#define EPSILON_CCI ( 100*std::numeric_limits<T>::epsilon() )

template<class T=real_t>
struct gsPoint2d{
  gsPoint2d(T xCoord, T yCoord) {
    pt.x() = xCoord;
    pt.y() = yCoord;
  }

  inline T x() const { return pt.x();}
  inline T y() const { return pt.y();}

  gsPoint<2,T> pt;
};

// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are collinear
// 1 --> Clockwise
// 2 --> Counterclockwise
template<class T>
int orientationxx(const gsPoint2d<T> &p, const gsPoint2d<T> &q, const gsPoint2d<T> &r) {
  // See https://www.geeksforgeeks.org/orientation-3-ordered-points/
  // for details of below formula.
  double val = (q.y() - p.y()) * (r.x() - q.x()) - (q.x() - p.x()) * (r.y() - q.y());

  if (val == 0) return 0;  // collinear

  return (val > 0) ? 1 : 2; // clock or counterclock wise
}

// Given three collinear points p, q, r, the function checks if
// point q lies on the line segment 'pr'
template<class T>
bool onSegment(const gsPoint2d<T> &p, const gsPoint2d<T> &q, const gsPoint2d<T> &r) {
  if (q.x() <= math::max(p.x(), r.x()) && q.x() >= math::min(p.x(), r.x()) &&
      q.y() <= math::max(p.y(), r.y()) && q.y() >= math::min(p.y(), r.y()))
    return true;

  return false;
}

// The main function that returns true if line segment 'p1q1'
// and 'p2q2' intersect.
template<class T>
bool doIntersect(const gsPoint2d<T> &p1, const gsPoint2d<T> &q1,
                 const gsPoint2d<T> &p2, const gsPoint2d<T> &q2) {
  // Find the four orientations needed for general and special cases
  int o1 = orientationxx(p1, q1, p2);
  int o2 = orientationxx(p1, q1, q2);
  int o3 = orientationxx(p2, q2, p1);
  int o4 = orientationxx(p2, q2, q1);

  // General case
  if (o1 != o2 && o3 != o4)
    return true;

  // Special Cases
  // p1, q1 and p2 are collinear and p2 lies on segment p1q1
  if (o1 == 0 && onSegment(p1, p2, q1)) return true;

  // p1, q1 and q2 are collinear and q2 lies on segment p1q1
  if (o2 == 0 && onSegment(p1, q2, q1)) return true;

  // p2, q2 and p1 are collinear and p1 lies on segment p2q2
  if (o3 == 0 && onSegment(p2, p1, q2)) return true;

  // p2, q2 and q1 are collinear and q1 lies on segment p2q2
  if (o4 == 0 && onSegment(p2, q1, q2)) return true;

  return false; // Doesn't fall in any of the above cases
}

/** \brief
    Curve/Curve intersection result class.

    This is the geometry type associated with gsBSplineBasis.

    \tparam T coefficient type

    \param m_paramOnCurve1 intersection parameter on Curve1
    \param m_paramOnCurve2 intersection parameter on Curve2
    \param m_point physical coordinate of intersection point
*/
template<class T=real_t>
struct gsCurveIntersectionResult {
  // Explicit constructor for direct initialization
  explicit gsCurveIntersectionResult(T paramOnCurve1, T paramOnCurve2, gsVector<T> pt)
      : m_paramOnCurve1(paramOnCurve1), m_paramOnCurve2(paramOnCurve2), m_point(std::move(pt)) {}

  // Equality operators
  bool operator==(const gsCurveIntersectionResult& other) const {
    return m_paramOnCurve1 == other.paramOnCurve1 &&
        m_paramOnCurve2 == other.paramOnCurve2 &&
        m_point == other.point;
  }

  bool operator!=(const gsCurveIntersectionResult& other) const {
    return !(*this == other);
  }

  T getParamOnCurve1() const { return m_paramOnCurve1; }
  T getParamOnCurve2() const { return m_paramOnCurve2; }
  gsVector<T> getPoint() const { return m_point; }

 private:
  T m_paramOnCurve1{};
  T m_paramOnCurve2{};
  gsVector<T> m_point{};
};

/// Managing an interval defined by a minimum and maximum value
template<class T=real_t>
class gsInterval {
 public:
  gsInterval(T min, T max) : m_min(min), m_max(max) {}

  bool almostEqual(T a, T b, T tolerance = EPSILON_CCI) const {
    return std::abs(a - b) <= tolerance;
  }

  bool operator==(const gsInterval<T> &other) const {
    return almostEqual(m_min, other.getMin()) && almostEqual(m_max, other.getMax());
  }

  T getMin() const { return m_min; }
  T getMax() const { return m_max; }

  void setMin(T mmin) { m_min = mmin; }
  void setMax(T mmax) { m_min = mmax; }

 private:
  T m_min{0}, m_max{0};
};

/// representing and manipulating bounding boxes for gsBSpline object
template<class T=real_t>
class gsCurveBoundingBox {
 public:
  explicit gsCurveBoundingBox(const gsBSpline<T> &curve)
      : range(curve.domainStart(), curve.domainEnd()) {
    low.resize( curve.geoDim() );
    high.resize(curve.geoDim() );
    for (int i = 0; i != curve.geoDim(); ++i) {
      low[i] = std::numeric_limits<T>::max();
      high[i] = std::numeric_limits<T>::lowest();
    }
    // compute min / max from control points
    for (int ipt = 0; ipt != curve.coefsSize(); ++ipt) {
      const auto &p = curve.coef(ipt); // Access once per iteration
      for (int i = 0; i != curve.geoDim(); ++i) {
        high[i] = math::max(high[i], p(i));
        low[i] = math::min(low[i], p(i));
      }
    }
  }

  [[nodiscard]] bool intersect(const gsCurveBoundingBox<T> &other,
                                T eps = EPSILON_CCI) const {
    assert( this->getLow().size() == other.getLow().size() );
    for (int i = 0; i != other.getLow().size(); ++i) {
      if (math::max(low[i], other.low[i]) > math::min(high[i], other.high[i]) + eps) {
        return false;
      }
    }
    return true;
  }

  gsVector<T> getLow() const { return low; }
  gsVector<T> getHigh() const { return high; }
  gsInterval<T> getRange() const { return range; }

 private:
  gsVector<T> low{}, high{};
  gsInterval<T> range{};
};

template<class T=real_t>
struct gsBoundingBoxPair {
  gsBoundingBoxPair(const gsCurveBoundingBox<T> &i1, const gsCurveBoundingBox<T> &i2)
      : b1(i1), b2(i2) {}

  gsCurveBoundingBox<T> b1{};
  gsCurveBoundingBox<T> b2{};
};

/// compute potential intersection range by recursively subdividing curves
/// and bounding box intersection detection
template<class T=real_t>
std::vector<gsBoundingBoxPair<T>> getPotentialIntersectionRanges(const gsBSpline<T> &curve1,
                                                                 const gsBSpline<T> &curve2) {
  gsCurveBoundingBox<T> h1(curve1);
  gsCurveBoundingBox<T> h2(curve2);

  // Bounding boxes do not intersect. No intersection possible
  if (!h1.intersect(h2)) {
    return {};
  }

  T crv1Curvature = curve1.pseudoCurvature();
  T crv2Curvature = curve2.pseudoCurvature();
  constexpr T MAX_CURVATURE = 1.0 + 5e-6;

  // Check for intersection between endpoint line segments if curves are linear enough
  if (crv1Curvature <= MAX_CURVATURE && crv2Curvature <= MAX_CURVATURE) {
    if ( curve1.geoDim() == 2 ) {
      // line segment intersection check for the planar case
      gsPoint2d<T> pt1(curve1.coef(0).x(), curve1.coef(0).y());
      gsPoint2d<T> pt2(curve1.coef(curve1.coefsSize() - 1).x(), curve1.coef(curve1.coefsSize() - 1).y());
      gsPoint2d<T> pt3(curve2.coef(0).x(), curve2.coef(0).y());
      gsPoint2d<T> pt4(curve2.coef(curve2.coefsSize() - 1).x(), curve2.coef(curve2.coefsSize() - 1).y());
      if (doIntersect(pt1, pt2, pt3, pt4)) {
        return {gsBoundingBoxPair<T>(h1, h2)};
      }
    } else {
      return {gsBoundingBoxPair<T>(h1, h2)};
    }
    return {};
  }

  // Recursive subdividing for curves with high curvature
  if (crv1Curvature > MAX_CURVATURE && crv2Curvature > MAX_CURVATURE) {
    // Subdivide both curves by splitting them in the parametric center
    T curve1MidParm = 0.5 * (curve1.domainStart() + curve1.domainEnd());
    T curve2MidParm = 0.5 * (curve2.domainStart() + curve2.domainEnd());

    gsBSpline<T> c11, c12, c21, c22;
    curve1.splitAt(curve1MidParm, c11, c12);
    curve2.splitAt(curve2MidParm, c21, c22);

    auto result1 = getPotentialIntersectionRanges<T>(c11, c21);
    auto result2 = getPotentialIntersectionRanges<T>(c11, c22);
    auto result3 = getPotentialIntersectionRanges<T>(c12, c21);
    auto result4 = getPotentialIntersectionRanges<T>(c12, c22);

    // append all results
    result1.insert(result1.end(), result2.begin(), result2.end());
    result1.insert(result1.end(), result3.begin(), result3.end());
    result1.insert(result1.end(), result4.begin(), result4.end());

    return result1;
  } else if (crv1Curvature <= MAX_CURVATURE && MAX_CURVATURE < crv2Curvature) {
    // Subdivide only curve 2
    T curve2MidParm = 0.5 * (curve2.domainStart() + curve2.domainEnd());
    gsBSpline<T> c21, c22;
    curve2.splitAt(curve2MidParm, c21, c22);

    auto result1 = getPotentialIntersectionRanges<T>(curve1, c21);
    auto result2 = getPotentialIntersectionRanges<T>(curve1, c22);

    result1.insert(result1.end(), result2.begin(), result2.end());
    return result1;
  } else if (crv2Curvature <= MAX_CURVATURE && MAX_CURVATURE < crv1Curvature) {
    // Subdivide only curve 1
    T curve1MidParm = 0.5 * (curve1.domainStart() + curve1.domainEnd());
    gsBSpline<T> c11, c12;
    curve1.splitAt(curve1MidParm, c11, c12);

    auto result1 = getPotentialIntersectionRanges<T>(c11, curve2);
    auto result2 = getPotentialIntersectionRanges<T>(c12, curve2);

    result1.insert(result1.end(), result2.begin(), result2.end());
    return result1;
  }

  return {};
}

/// Refine the intersection by using Projected Gauss-Newton method
template<class T=real_t>
class gsCurveCurveDistanceSystem {
 public:
  gsCurveCurveDistanceSystem(const gsBSpline<T> &crv1, const gsBSpline<T> &crv2)
      : m_crv1(crv1), m_crv2(crv2) {
    assert(crv1.geoDim() == crv2.geoDim());
  }

  void Values(const gsMatrix<T, 2, 1> &uv,
              gsMatrix<T> &funcVal,
              gsMatrix<T> &invJac) const {
    std::vector<gsMatrix<>> crv1Der = m_crv1.evalAllDers(uv.row(0), 1);
    std::vector<gsMatrix<>> crv2Der = m_crv2.evalAllDers(uv.row(1), 1);

    funcVal = crv1Der[0] - crv2Der[0];

    gsMatrix<T> jac(m_crv1.geoDim(), 2);
    jac.col(0) =  crv1Der[1];
    jac.col(1) = -crv2Der[1];
    invJac = jac.completeOrthogonalDecomposition().pseudoInverse();
  }

  /// Projected Gauss-Newton method
  T compute(gsMatrix<T,2,1> &uv, T tolerance = 1e-5, size_t maxIter = 20) const
  {
    gsMatrix<T> funVal;
    gsMatrix<T> invJac;
    Values(uv, funVal, invJac); // Assumes Values correctly inverts the Jacobian
    T currentResidual = funVal.norm();

    gsMatrix<T, 2, 1> deltaUv = invJac * funVal;
    gsMatrix<T, 2, 1> newUv;
    for (size_t iter = 1; iter != maxIter; ++iter) {
      uv -= deltaUv; // update uv

      // project uv into its parameter range
      newUv(0, 0) = math::min(math::max(uv(0, 0), m_crv1.domainStart()), m_crv1.domainEnd());
      newUv(1, 0) = math::min(math::max(uv(1, 0), m_crv2.domainStart()), m_crv2.domainEnd());
      // using C++17 standard
//      newUv(0, 0) = std::clamp(newUv(0, 0), m_crv1.domainStart(), m_crv1.domainEnd());
//      newUv(1, 0) = std::clamp(newUv(1, 0), m_crv2.domainStart(), m_crv2.domainEnd());

      T delta = (uv - newUv).norm();
      uv = newUv;
      if (delta > EPSILON_CCI) {
        Values(uv, funVal, invJac); // Re-evaluate funVal and invJac
        currentResidual = funVal.norm();
        break;
      }

      Values(uv, funVal, invJac); // Re-evaluate funVal and invJac
      currentResidual = funVal.norm();
      if (currentResidual < tolerance) break;

      deltaUv = invJac * funVal;
      if (deltaUv.norm() < 1e-3*tolerance) break;

//      gsInfo << "iter. = " << iter << ", residual = " << currentResidual << "\n";
    }
    return currentResidual;
  }

 private:
  const gsBSpline<T> m_crv1{}, m_crv2{};
};

} // namespace internal

} // namespace gismo

