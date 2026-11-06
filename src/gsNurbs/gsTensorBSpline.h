/** @file gsTensorBSpline.h

    @brief Represents a tensor-product B-spline patch

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsGeometry.h>
#include <gsNurbs/gsTensorBSplineBasis.h>

namespace gismo
{

/// Helper function for the slice function
/// selects the row of coefficients from coefficients of geo that are suitable
/// for the isoparametric slice in \a dir_fixed with \a par.
/// Note that geo has to have already C^0 continuity at \a par in direction \a dir.
template <short_t d, class T>
void constructCoefsForSlice(index_t dir_fixed, const index_t index,
                            const gsMatrix<T> & fullCoefs,
                            const gsVector<index_t, d> & sizes,
                            gsMatrix<T>& result);

/** \brief
    A tensor product of \em d B-spline functions, with arbitrary target dimension.

    This is the geometry type associated with gsTensorBSplineBasis.
    
    \tparam d the parametric dimension of the tensor product
    \tparam T coefficient type

    # Periodic patches #

    A patch whose basis is periodic in direction pd stores n reduced
    coefficient slices in pd (see gsKnotVector for the notation). That is
    exactly \c basis().size() rows (I3), before and after every operation.
    setPeriodic(dir) compresses n+m ghost slices to n.

    Operations along pd reshape the coefficients along pd, apply the 1D
    gsBSpline operation, and commit only on success, so a refusal leaves the
    patch unchanged:
     + insertKnot() (exact);
     + degreeElevate() (exact in exact arithmetic);
     + degreeReduce(), degreeIncrease(), degreeDecrease() (anchor
       interpolation, approximate);
     + reverse(k) (keeps [a,b]).

    slice() works at a fixed parameter in either direction. At a fixed
    periodic parameter it is exact, with no knot insertion.

    Refused (std::runtime_error): splitAt(), uniformSplit() and
    splitAtMult() in pd; the corner functions; the four-corner constructor
    with periodic knots.

    \ingroup geometry
    \ingroup Nurbs
*/
template<short_t d, class T>
class gsTensorBSpline GISMO_FINAL : public gsGeoTraits<d,T>::GeometryBase
{
public: 
    typedef gsKnotVector<T> KnotVectorType;

    typedef T Scalar_t;

    typedef typename gsBSplineTraits<d,T>::Basis Basis;

    typedef typename gsGeoTraits<d,T>::GeometryBase Base;

    /// Family type
    typedef gsBSplineBasis<T>  Family_t;

    /// Shared pointer for gsTensorBSpline
    typedef memory::shared_ptr< gsTensorBSpline > Ptr;

    /// Unique pointer for gsTensorBSpline
    typedef memory::unique_ptr< gsTensorBSpline > uPtr;

    /// Associated Boundary basis type
    typedef typename gsBSplineTraits<static_cast<short_t>(d-1),T>::Geometry BoundaryGeometryType;

    /// Associated Boundary basis type
    typedef typename gsBSplineTraits<static_cast<short_t>(d-1),T>::Basis BoundaryBasisType;

public:

    /// Default empty constructor
    gsTensorBSpline() { }

    using gsGeometry<T>::swap;

#if !EIGEN_HAS_RVALUE_REFERENCES
    gsTensorBSpline & operator=(gsTensorBSpline other)
    { this->swap(other); return *this;}
#endif

    // Construct B-Spline by basis functions and coefficient matrix
    //gsTensorBSpline( const gsConstantBasis<T> & basis, const gsMatrix<T> & coefs )
    //{ }

    /// Construct B-Spline by basis functions and coefficient matrix
    gsTensorBSpline( const Basis & basis, gsMatrix<T> coefs )
    : Base( basis, give(coefs) )
    { }
        
    /// Construct 2D tensor B-Spline by knot vectors, degrees and
    /// coefficient matrix (copying coefficient matrix)
    template<typename U>
    gsTensorBSpline( KnotVectorType KV1, gsKnotVector<U> KV2, 
                     gsMatrix<T> tcoefs,
                     typename util::enable_if<d==2,U>::type * = NULL)
    {
        GISMO_ASSERT(d==2, "Wrong dimension: tried to make a "
                     << d<<"D tensor B-spline using 2 knot-vectors.");

        std::vector<Family_t*> cbases;
        cbases.push_back(new gsBSplineBasis<T>(give(KV1)) );
        cbases.push_back(new gsBSplineBasis<T>(give(KV2)) );
        Basis * tbasis = Basis::New(cbases); //d==2
        
        GISMO_ASSERT(tbasis->size()==tcoefs.rows(),
                     "Coefficient matrix for the tensor B-spline does not have "
                     "the expected number of control points (rows)." );
        
        this->m_basis = tbasis;
        this->m_coefs = give(tcoefs);
    }
    
    /// Construct 2D tensor B-Spline by knot vectors, degrees and 4 corner vertices,
    /// the four boundary curves are linear interpolations of the corners,
    /// the size of the matrix *corner* is 4 by 3.
    /// Refused (GISMO_ENSURE) if a knot vector is periodic: a patch given by
    /// its four corners has no periodic meaning.
    gsTensorBSpline(gsMatrix<T> const & corner,
                    KnotVectorType KV1,
                    KnotVectorType KV2);

    /// Construct 3D tensor B-Spline by knot vectors, degrees and
    /// coefficient matrix (copying coefficient matrix)
    gsTensorBSpline( KnotVectorType KV1,
                     KnotVectorType KV2,
                     KnotVectorType KV3,
                     gsMatrix<T> tcoefs )
    {
        GISMO_ASSERT(d==3, "Wrong dimension: tried to make a "
                     << d<<"D tensor B-spline using 3 knot-vectors.");

        std::vector<Family_t*> cbases;
        cbases.push_back(new gsBSplineBasis<T>(give(KV1)) );
        cbases.push_back(new gsBSplineBasis<T>(give(KV2)) );
        cbases.push_back(new gsBSplineBasis<T>(give(KV3)) );
        Basis * tbasis = Basis::New(cbases); //d==3
        
        GISMO_ASSERT(tbasis->size()==tcoefs.rows(),
                     "Coefficient matrix for the tensor B-spline does not have "
                     "the expected number of control points (rows)." );
        
        this->m_basis = tbasis;
        this->m_coefs = give(tcoefs);
    }

    /// Construct 4D tensor B-Spline by knot vectors, degrees and
    /// coefficient matrix (copying coefficient matrix)
    gsTensorBSpline( KnotVectorType KV1,
                     KnotVectorType KV2,
                     KnotVectorType KV3,
                     KnotVectorType KV4,
                     gsMatrix<T> tcoefs )
    {
        GISMO_ASSERT(d==4, "Wrong dimension: tried to make a "
                     << d<<"D tensor B-spline using 4 knot-vectors.");

        std::vector<Family_t*> cbases;
        cbases.reserve(4);
        cbases.push_back(new gsBSplineBasis<T>(give(KV1)) );
        cbases.push_back(new gsBSplineBasis<T>(give(KV2)) );
        cbases.push_back(new gsBSplineBasis<T>(give(KV3)) );
        cbases.push_back(new gsBSplineBasis<T>(give(KV4)) );
        Basis * tbasis = Basis::New(cbases); //d==4
        
        GISMO_ASSERT(tbasis->size()==tcoefs.rows(),
                     "Coefficient matrix for the tensor B-spline does not have "
                     "the expected number of control points (rows)." );
        
        this->m_basis = tbasis;
        this->m_coefs = give(tcoefs);
    }
    
    GISMO_CLONE_FUNCTION(gsTensorBSpline)

    GISMO_BASIS_ACCESSORS

public:

    /// @brief Elevates the degree by \a i in direction \a dir (all directions
    /// if \a dir = -1), keeping the geometry.
    ///
    /// On a geometry that is periodic in \a dir the elevation is exact in
    /// exact arithmetic (no solver): the 1D periodic elevation of
    /// gsBSpline::degreeElevate is applied to every iso-curve along \a dir.
    /// The periodic direction keeps its n slices.
    void degreeElevate(short_t const i = 1, short_t const dir = -1) override;

    /// @brief Reduces the degree by \a i in direction \a dir (all
    /// directions if \a dir = -1).
    ///
    /// Non-periodic geometry: as gsGeometry::degreeReduce(). Periodic
    /// geometry: an approximation, not the same geometry in general. The
    /// 1D periodic anchor interpolation of gsBSpline::degreeReduce is applied
    /// along \a dir, so the result interpolates the old geometry on the
    /// iso-lines through the new anchors. A refusal of the 1D operation
    /// (GISMO_ERROR) leaves the geometry unchanged.
    void degreeReduce(short_t const i = 1, short_t const dir = -1) override;

    /// @brief Increases the degree by \a i in direction \a dir (all
    /// directions if \a dir = -1), in general only approximately on a periodic
    /// geometry; see degreeReduce().
    void degreeIncrease(short_t const i = 1, short_t const dir = -1) override;

    /// @brief Decreases the degree by \a i in direction \a dir (all
    /// directions if \a dir = -1), in general only approximately on a periodic
    /// geometry; see degreeReduce().
    void degreeDecrease(short_t const i = 1, short_t const dir = -1) override;

    /// @brief Inserts knot \a knot at direction \a dir, \a i times.
    ///
    /// Exact on a periodic direction, by the 1D periodic transfer of
    /// gsBSpline::insertKnot (the knot rules and refusals are those of the
    /// curve; the geometry is unchanged if refused).
    void insertKnot( T knot, int dir, int i = 1) override;

    /// Returns a reference to the knot vector in direction \a i
    KnotVectorType & knots(const int i) { return this->basis().knots(i); }

    /// Returns a reference to the knot vector \a i
    const KnotVectorType & knots(const int i) const { return this->basis().knots(i); }

    /*** Virtual member functions required by the base class ***/

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const override;

    /*** Additional members for tensor B-Splines ***/

    /// Returns the degree of the basis wrt direction i
    short_t degree(const unsigned & i) const
    { return this->basis().component(i).degree(); }

    /// @brief Toggle orientation wrt coordinate k.
    ///
    /// In a periodic direction the domain [a,b] is kept and the n slices are
    /// permuted: new slice j is old slice (m-1-j) mod n, m being the number of
    /// crossing functions.
    /// \todo use flipTensor to generalize to any dimension
    void reverse(unsigned k);

    /// Swap directions. A periodic direction moves with its component.
    void swapDirections(const unsigned i, const unsigned j);

    void toggleOrientation() override;

    /// \brief Return true if point \a u is a corner of
    /// the patch with tolerance \a tol.
    /// Refused (GISMO_ENSURE) on a periodic geometry: every corner touches a
    /// side of the periodic direction. The same holds for findCorner(),
    /// setOriginCorner() and setFurthestCorner().
    bool isPatchCorner(gsMatrix<T> const &v, T tol = 1e-3) const;

    /// \brief returns the tensor-index \a curr of the corner control
    /// point \a v, or an invalid index if the corner is not found
    /// within the tolerance \a tol
    void findCorner(const gsMatrix<T>   & v,
                    gsVector<index_t,d> & curr,
                    T tol = 1e-3);

    /// \brief Modifies the parameterization such that the point \a v
    /// is the origin of the parametrization of the patch. Assumes
    /// that \a v is either input is indeed a corner of this patch
    void setOriginCorner(gsMatrix<T> const &v);

    /// \brief Modifies the parameterization such that the point \a v
    /// is the ending corner of the parametrization of the
    /// patch. Assumes that \a v is either input is indeed a corner of
    /// this patch
    void setFurthestCorner(gsMatrix<T> const &v);
    
    /// @brief Makes the geometry periodic in direction \a dir.
    ///
    /// The coefficients are given with n+m slices in \a dir (m = number of
    /// crossing functions) and are compressed to the n slices of the reduced
    /// storage, which every operation of a periodic geometry keeps. Ghost
    /// slice n+j must equal its twin slice j (the tolerance is that of
    /// gsBSpline::setPeriodic). Nothing happens if \a dir is already periodic.
    /// Refused (GISMO_ENSURE, geometry unchanged) for \a dir out of range, if
    /// another direction is periodic, if the knots do not admit a periodic
    /// basis, or if a ghost slice differs from its twin.
    void setPeriodic( int dir );
    
    /// Returns a local representation of the geometry in the cell
    /// containing the point \a u
    ///
    /// In a periodic direction the local coefficients are read from the n
    /// reduced slices through the folded indices of the 1D basis
    /// (active_into()); the local knot vectors are not periodic.
    typename gsGeometry<T>::uPtr localRep(const gsMatrix<T> & u) const;

public:

    /// Constructs an isoparametric slice of this tensorBSpline by fixing
    /// \a par in direction \a dir_fixed. The resulting tensorBSpline has
    /// one less dimension and is given back in \a result.
    ///
    /// On a periodic geometry, a fixed non-periodic parameter gives a result
    /// that is periodic with n rows in the periodic direction. A fixed periodic
    /// parameter is supported too: the coefficients are the contraction of the
    /// coefficient tensor with the p+1 active periodic basis functions at
    /// \a par; exact, with no knot insertion.
    void slice(index_t dir_fixed,T par,BoundaryGeometryType & result) const;

    /// Splits the geometry either two parts in direction \a dir, or if \a dir = -1
    /// in 2^d parts, by calling splitAt() for each direction.
    /// The function automatically searches for the midpoint the corresponding knot vector.
    /// Refused (GISMO_ENSURE) in a periodic direction, see splitAt().
    std::vector<gsGeometry<T>* > uniformSplit(index_t dir = -1) const override;

    /// Split the patch into smaller patches at the position of all
    /// knots with multiplicity at least \a minMult.
    /// Refused (GISMO_ENSURE) in a periodic direction, see splitAt().
    std::vector<gsGeometry<T>* > splitAtMult(index_t minMult = 1, index_t dir = -1) const;

    /// Splits the geometry into two pieces (\a left, \a right) along direction \a dir at \a xi. The splitting
    /// is performed by increasing the multiplicity of knot \a xi to p+1, or if \a xi does not exist as knot,
    /// it is inserted p+1 times.
    /// Refused (GISMO_ENSURE) if \a dir is a periodic direction (a closed
    /// direction has no end points to split at); \a left and \a right are
    /// unchanged then. Splitting in a non-periodic direction keeps the
    /// periodicity of the other direction.
    void splitAt( index_t dir,T xi, gsTensorBSpline<d,T>& left,  gsTensorBSpline<d,T>& right) const;

    typename gsGeometry<T>::uPtr iface(const boundaryInterface & bi,
                                       const gsGeometry<T> & other) const override;

protected:
    // TODO Check function
    // check function: check the coefficient number, degree, knot vector ...

    using Base::m_basis;
    using Base::m_coefs;

private:

    /// Index of the periodic direction, or -1.
    short_t _periodicDir() const
    {
        for (short_t j = 0; j < d; ++j)
            if (knots(j).isPeriodic())
                return j;
        return -1;
    }

    /// @brief Applies a 1D gsBSpline operation along direction \a dir.
    ///
    /// The coefficient tensor is reshaped so that row r is index r along
    /// \a dir (size(dir) rows) and each column is one multi-index of the
    /// other directions and one geometric coordinate, wrapped with a copy of
    /// component(dir) in a temporary curve, and op(curve) is called. As the
    /// 1D operations are linear in the rows, this is the tensor-product
    /// operation (identity in the other directions). If op throws,
    /// the geometry is unchanged. Cost: a copy of the coefficients,
    /// O(N geoDim), plus the 1D operation.
    template<class Op> void _alongDirection(short_t dir, Op op);

}; // class gsTensorBSpline

#ifdef GISMO_WITH_PYBIND11

  /**
   * @brief Initializes the Python wrapper for the class: gsTensorBSpline
   */
  template <short_t d>
  void pybind11_init_gsTensorBSpline(pybind11::module &m);

    /**
     * @brief Initializes the Python wrapper for the factory constructor of gsTensorBSpline
     */
    void pybind11_init_gsTensorBSpline_factory(pybind11::module &m);

#endif // GISMO_WITH_PYBIND11

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTensorBSpline.hpp)
#endif
