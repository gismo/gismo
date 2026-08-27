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
    /// the size of the matrix *corner* is 4 by 3
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

    // Look at gsGeometry class for a description
    void degreeElevate(short_t const i = 1, short_t const dir = -1) override;

    /// Inserts knot \a knot at direction \a dir, \a i times
    void insertKnot( T knot, int dir, int i = 1) override;

    /// Removes knot \a knot from direction \a dir up to \a i times.
    /// A removal is accepted when the feasibility residual of Tiller's
    /// algorithm is below \a tol RELATIVE to the largest control-point entry
    /// of the coefficient fiber; see gsKnotRemoveSingle().
    /// Returns the number of times the knot was successfully removed.
    index_t removeKnot( T knot, short_t dir, short_t i = 1, T tol = (T)(1e-10));

    // -----------------------------------------------------------------------
    //  Calculus methods
    // -----------------------------------------------------------------------

    /// Returns a new gsTensorBSpline in full Bézier (C^{-1}) form by inserting
    /// each interior knot to multiplicity degree+1 in every direction.
    gsTensorBSpline<d,T> toBezier() const;

    /// Returns the pointwise product A·B as a new gsTensorBSpline.
    /// Both arguments must be scalar (targetDim==1) and live on the same
    /// parametric domain; their degrees and knot vectors may differ.
    ///
    /// Both factors are brought to Bézier form on the union of their
    /// breakpoints and multiplied element by element with the Bernstein
    /// product formula
    /// \f[ (AB)_k = \sum_{i+j=k}
    ///        \frac{\binom{p_A}{i}\binom{p_B}{j}}{\binom{p_A+p_B}{k}}
    ///        a_i b_j \f]
    /// applied per direction, so the result has degree \f$p_A+p_B\f$.
    ///
    /// If \a keepBezier is false (default), interior knots are then removed to
    /// reach the minimal space: where A is \f$C^{a}\f$ and B is \f$C^{b}\f$
    /// the product is \f$C^{\min(a,b)}\f$, i.e. interior multiplicity
    /// \f$p_A+p_B-\min(a,b)\f$.  Whether that removal actually succeeds is
    /// decided by the tolerance of removeKnot().
    /// If \a keepBezier is true, the raw Bézier (C^{-1}) result is returned.
    ///
    /// Complexity: O( nel · (p_A+p_B+1)^d · (p_A+1)^d ) for the product.
    static gsTensorBSpline<d,T> multiply(const gsTensorBSpline<d,T> & A,
                                         const gsTensorBSpline<d,T> & B,
                                         bool keepBezier = false);

    /// Brings \a splines to a common tensor-product space, in place, so that
    /// their coefficient matrices may be combined entry by entry.
    /// Each spline is first degree-elevated to the per-direction maximum
    /// degree, then knots are inserted to the per-direction union of the
    /// interior knots.  Elevation must come first because it raises interior
    /// multiplicities in order to preserve continuity.
    /// All splines must live on the same parametric domain.
    static void makeCompatible(std::vector< gsTensorBSpline<d,T> > & splines);

    /// Returns \a a·A + \a b·B, evaluated in the common space produced by
    /// makeCompatible().  A and B need not share a basis; they must have the
    /// same target dimension and parametric domain.
    static gsTensorBSpline<d,T> linearCombination(
                                    T a, const gsTensorBSpline<d,T> & A,
                                    T b, const gsTensorBSpline<d,T> & B);

    /// Returns c² as a new gsTensorBSpline (only for scalar targetDim==1).
    /// Equivalent to multiply(*this, *this, \a keepBezier).
    gsTensorBSpline<d,T> squared(bool keepBezier = false) const;

    /// Returns c³ as a new gsTensorBSpline (only for scalar targetDim==1).
    /// Equivalent to multiply(squared(), *this, \a keepBezier).
    gsTensorBSpline<d,T> cubed(bool keepBezier = false) const;

    /// Returns the partial derivative in direction \a dir as a gsTensorBSpline.
    /// The result has degree p-1 in direction \a dir and p in all other directions.
    gsTensorBSpline<d,T> grad(short_t dir) const;

    /// Returns all d partial derivatives as a vector of gsTensorBSpline.
    std::vector< gsTensorBSpline<d,T> > grad() const;

    /// Returns the divergence: for targetDim==d, sum_k ∂f_k/∂x_k.
    /// Result is a scalar gsTensorBSpline.
    gsTensorBSpline<d,T> div() const;

    /// Returns the Laplacian: for scalar targetDim==1, sum_k ∂²c/∂x_k².
    /// Every term is degree-elevated back to p in the direction it was
    /// differentiated in, then all terms are summed in their common space.
    /// The result therefore has degree p and interior multiplicity m+2, where
    /// m is the input multiplicity — two orders less continuous than the input.
    ///
    /// \warning \a keepBezier does NOT reach Bézier form, despite its name.
    /// It only inserts each interior knot twice more in the directions that
    /// were not differentiated, which the knot union already supplies, so for
    /// simple interior knots it returns exactly the same space as the default.
    /// FIXME: either make it return toBezier() of the result, or drop the
    /// parameter; see gsSplineFactory_test:lapl_keepBezier_versus_minimal.
    gsTensorBSpline<d,T> lapl(bool keepBezier = false) const;

    /// Returns the Hessian of a scalar (targetDim==1) spline as a single
    /// gsTensorBSpline with targetDim() == d*d, where column i*d+j holds
    /// ∂²c/∂x_i∂x_j.  Requires degree >= 2 in every direction.
    /// All d² entries are brought to the common space S(p, C^{p-3}) by
    /// makeCompatible(), so the matrix is stored on one basis; the symmetric
    /// off-diagonal columns i*d+j and j*d+i are therefore identical.
    gsTensorBSpline<d,T> hess() const;

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

    /// Toggle orientation wrt coordinate k
    /// \todo use flipTensor to generalize to any dimension
    void reverse(unsigned k);

    /// Swap directions
    void swapDirections(const unsigned i, const unsigned j);

    void toggleOrientation() override;

    /// \brief Return true if point \a u is a corner of
    /// the patch with tolerance \a tol
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
    
    /// Sets the resulting BSpline to be periodic in direction \em dir.
    /// \param dir
    inline void setPeriodic( int dir )
    {
        this->m_coefs = this->basis().perCoefs( this->m_coefs, dir );
        this->basis().setPeriodic( dir );
    }
    
    /// Returns a local representation of the geometry in the cell
    /// containing the point \a u
    typename gsGeometry<T>::uPtr localRep(const gsMatrix<T> & u) const;

public:

    /// Constructs an isoparametric slice of this tensorBSpline by fixing
    /// \a par in direction \a dir_fixed. The resulting tensorBSpline has
    /// one less dimension and is given back in \a result.
    void slice(index_t dir_fixed,T par,BoundaryGeometryType & result) const;

    /// Splits the geometry either two parts in direction \a dir, or if \a dir = -1
    /// in 2^d parts, by calling splitAt() for each direction.
    /// The function automatically searches for the midpoint the corresponding knot vector.
    std::vector<gsGeometry<T>* > uniformSplit(index_t dir = -1) const override;

    /// Split the patch into smaller patches at the position of all
    /// knots with multiplicity at least \a minMult
    std::vector<gsGeometry<T>* > splitAtMult(index_t minMult = 1, index_t dir = -1) const;

    /// Splits the geometry into two pieces (\a left, \a right) along direction \a dir at \a xi. The splitting
    /// is performed by increasing the multiplicity of knot \a xi to p+1, or if \a xi does not exist as knot,
    /// it is inserted p+1 times.
    void splitAt( index_t dir,T xi, gsTensorBSpline<d,T>& left,  gsTensorBSpline<d,T>& right) const;

    typename gsGeometry<T>::uPtr iface(const boundaryInterface & bi,
                                       const gsGeometry<T> & other) const override;

protected:
    // TODO Check function
    // check function: check the coefficient number, degree, knot vector ...

    using Base::m_basis;
    using Base::m_coefs;

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
