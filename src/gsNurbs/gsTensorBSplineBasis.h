/** @file gsTensorBSplineBasis.h

    @brief Provides declaration of TensorBSplineBasis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/


#pragma once

#include <gsTensor/gsTensorBasis.h>
#include <gsNurbs/gsBSplineBasis.h>
#include <gsNurbs/gsTensorBSpline.h>

namespace gismo
{

/**
    @brief A tensor product B-spline basis.

    \tparam T coefficient type
    \tparam d dimension of the parameter domain
    \tparam KnotVectorType type of the knot-vector

    \ingroup basis
    \ingroup Nurbs
*/

template<short_t d, class T>
class gsTensorBSplineBasis : public gsTensorBasis<d,T>
{
public:
    typedef gsKnotVector<T> KnotVectorType;

    /// Base type
    typedef gsTensorBasis<d,T> Base;

    /// Family type
    typedef gsBSplineBasis<T>  Family_t;

    typedef gsTensorBSplineBasis Self_t;

    /// Coordinate basis type
    typedef gsBSplineBasis<T> CoordinateBasis;
    typedef CoordinateBasis                  Basis_t;

    /// Coefficient type
    typedef T Scalar_t;

    /// Associated Boundary basis type
    typedef typename gsBSplineTraits<d,T>::Geometry GeometryType;

    /// Associated Boundary basis type
    typedef typename gsBSplineTraits<d-1,T>::Basis BoundaryBasisType;

    typedef typename Base::iterator        iterator;
    typedef typename Base::const_iterator  const_iterator;

    /// Smart pointer for gsTensorBSplineBasis
    typedef memory::shared_ptr< Self_t > Ptr;

    /// Smart pointer for gsTensorBSplineBasis
    typedef memory::unique_ptr< Self_t > uPtr;

public:

    /// \brief Default constructor
    gsTensorBSplineBasis() : Base()
    {
        for(short_t i = 0; i!=d; ++i)
            this->m_bases[i] = new Basis_t();
    }

    void swap(gsTensorBSplineBasis & other)
    {
        this->Base::swap(static_cast<Base&>(other));
    }

#if !EIGEN_HAS_RVALUE_REFERENCES
    gsTensorBSplineBasis & operator=(gsTensorBSplineBasis other)
    { gsTensorBSplineBasis::swap(other); return *this;}
// #else // defined implicitly
//     gsTensorBSplineBasis(gsTensorBSplineBasis && other) : Base(give(other)) { }
//     gsTensorBSplineBasis(const gsTensorBSplineBasis &  other) : Base(other) { }
//     gsTensorBSplineBasis & operator=(gsTensorBSplineBasis&& other)
//     { return (gsTensorBSplineBasis &)Base::operator=(give(other)); }
//     gsTensorBSplineBasis & operator=(const gsTensorBSplineBasis& other)
//     { return (gsTensorBSplineBasis &)Base::operator=(other); }
#endif

    /**
        \brief Constructs a d-dimensional tensor product B-spline basis
               from exactly d knot-vectors.
    */
    template<typename... U,
             typename util::enable_if<(sizeof...(U)==d), int>::type = 0>
    gsTensorBSplineBasis( gsKnotVector<U>... KV )
    {
        std::vector<KnotVectorType> knotVectors{ KnotVectorType(give(KV))... };
        initFromKnotVectors(give(knotVectors));
    }

    explicit gsTensorBSplineBasis(std::vector<KnotVectorType> KV)
    { initFromKnotVectors(give(KV)); }

    template<typename... B,
             typename util::enable_if<(sizeof...(B)==d), int>::type = 0>
    gsTensorBSplineBasis( B*... bb ) : Base()
    {
        Basis_t * bases[d] = { bb... };
        for(short_t i = 0; i != d; ++i)
            this->m_bases[i] = bases[i];
    }

    explicit gsTensorBSplineBasis(std::vector< gsBasis<T>*> & bb );

    explicit gsTensorBSplineBasis(std::vector< Basis_t*> & bb );

#ifdef __DOXYGEN__
    /// \brief Returns the boundary basis for side s.
    typename gsBSplineTraits<static_cast<short_t>(d-1),T>::Basis::uPtr boundaryBasis(boxSide const & s);
#endif
    /// On a periodic basis the sides of the periodic direction are not
    /// boundaries: asking for one raises std::runtime_error.
    GISMO_UPTR_FUNCTION_DEF(BoundaryBasisType, boundaryBasis, boxSide const &) override
    {
        if (static_cast<int>(n1.direction()) == periodicDirection())
            GISMO_ERROR("gsTensorBSplineBasis::boundaryBasis: the periodic direction has no boundary.");
        std::vector<gsBasis<T>*> rr;
        this->getComponentsForSide(n1,rr);
        return BoundaryBasisType::New(rr);
    }

    GISMO_CLONE_FUNCTION(gsTensorBSplineBasis)

    static Self_t * New(std::vector<gsBasis<T>*> & bb )
    { return new Self_t(bb); }

    static Self_t * New(std::vector<Basis_t*> & bb )
    { return new Self_t(bb); }

    static uPtr make(std::vector<gsBasis<T>*> & bb )
    { return uPtr( new Self_t(bb) ); }

    static uPtr make(std::vector<Basis_t*> & bb )
    { return uPtr( new Self_t(bb) ); }

public:

    KnotVectorType & knots (int i)
    { return Self_t::component(i).knots(); }

    const KnotVectorType & knots (int i) const
    { return Self_t::component(i).knots(); }

    // knot \a k of direction \a i
    T knot(int i, int k) const
    { return Self_t::component(i).knots()[k]; }


    const Basis_t & component(short_t dir) const override
    {
        return static_cast<const Basis_t &>(Base::component(dir));
    }

    Basis_t & component(short_t dir) override
    {
        return static_cast<Basis_t &>(Base::component(dir));
    }

    /// \brief Active functions at the points \a u.
    ///
    /// On a periodic basis the index of the periodic direction is folded
    /// modulo size(pd): function j of that direction is the sum of its
    /// stretched copies j, j+n, ..., so the indices are unique and lie in
    /// [0, size()). Cost O(npts (p+1)^d) as for a non-periodic basis, which
    /// does not evaluate a modulo.
    void active_into(const gsMatrix<T> & u, gsMatrix<index_t>& result) const override;

    /// Returns a box with the coordinate-wise active functions
    ///
    /// On a periodic basis the box is in raw stretched indices of the
    /// periodic direction pd (the upper corner may exceed size(pd)-1),
    /// because folded indices do not form a box; callers that index
    /// coefficients must reduce modulo size(pd) themselves.
    /// \param u evaluation points
    /// \param low lower left corner of the box
    /// \param upp upper right corner of the box
    void active_cwise(const gsMatrix<T> & u, gsVector<index_t,d>& low,
                      gsVector<index_t,d>& upp ) const;

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const override
    {
        os << "TensorBSplineBasis: dim=" << this->dim()<< ", size="<< this->size() <<".";
        const int pd = periodicDirection();
        if( pd != -1 )
            os << "Periodic in " << pd << "-th direction.\n";
        for ( short_t i = 0; i!=d; ++i )
            os << "\n  Direction "<< i <<": "<< Self_t::component(i).knots() <<" ";
        os << "\n";
        return os;
    }

    //
    // param other parent/reference mesh determining the smoothness at the inner knots.
    // param i number of k-refinement steps to perform

    /// \brief Perform k-refinement coordinate-wise, in all directions.
    ///
    /// \copydetails gsBSplineBasis::refine_k
    void k_refine(Self_t & other, int const & i = 1)
    {
        for (short_t j = 0; j < d; ++j)
            Self_t::component(j).refine_k(other.component(j), i);
    }

    /// \brief p-refinement (essentially degree elevation in all
    /// directions)
    void refine_p(int const & i = 1)
    {
        for (short_t j = 0; j < d; ++j)
            Self_t::component(j).refine_p(i);
    }

    /// \brief Uniform h-refinement (placing \a i new knots inside
    /// each knot-span, for all directions
    void refine_h(int const & i = 1)
    {
        for (short_t j = 0; j < d; ++j)
            Self_t::component(j).refine_h(i);
    }

    /**
     * \brief Takes a vector of coordinate wise knot values and
     * inserts these values to the basis.
     *
     * Also constructs and returns
     * the transfer matrix that transfers coefficients to the new
     * basis.
     *
     * \param[out]  transfer Transfer matrix
     * \param[in] refineKnots Coordinate-wise knot values to be inserted
     */
    void refine_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, const std::vector< std::vector<T> >& refineKnots);
    void refine_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, gsMatrix<T> const & boxes)
    {
        this->refine_withTransfer(transfer,this->_boxToKnots(boxes));
    }

    /**
     * \brief Takes a vector of coordinate wise knot values and
     * inserts these values to the basis.
     *
     * Also takes the old coefficients and changes them to reflect the new coefficients.
     *
     * \param[out]  coefs new coefficients
     * \param refineKnots Coordinate-wise knot values to be inserted
     *
     * \todo rename to insertKnots_withCoefs
     */
    void refine_withCoefs(gsMatrix<T> & coefs,const std::vector< std::vector<T> >& refineKnots);
    void refine_withCoefs(gsMatrix<T> & coefs, gsMatrix<T> const & boxes)
    {
        this->refine_withCoefs(coefs,this->_boxToKnots(boxes));
    }

    /// Inserts the knot \em knot with multiplicity \em mult in the knot
    /// vector of direction \a dir.
    void insertKnot(T knot, index_t dir, int mult=1)
    { Self_t::component(dir).insertKnot( knot, mult); }

    /// Removes the knot \em knot with multiplicity \em mult in the knot
    /// vector of direction \a dir.
    void removeKnot(T knot, index_t dir, int mult=1)
    { Self_t::component(dir).removeKnot( knot, mult); }

    /**
     * \brief Takes a vector of coordinate wise knot values and
     * inserts these values to the basis.
     *
     * \param refineKnots Coordinate-wise knot values to be inserted
     */
    void insertKnots(const std::vector< std::vector<T> >& refineKnots)
    {
        GISMO_ASSERT( refineKnots.size() == d, "refineKnots vector has wrong size" );
        for (short_t j = 0; j < d; ++j) // refine basis in each direction
            this->knots(j).insert(refineKnots[j]);
    }

    /** \brief Refinement of the tensor basis on the area defined by \em boxes.
     *
     * Applies "local" refinement within the tensor-product structure of
     * gsTensorBSplineBasis. The areas for refinement are specified in \em boxes.\n
     * \em boxes is a gsMatrix of size <em>d</em> x <em>(2*N)</em>, where\n
     * \em d is the dimension of the parameter domain and\n
     * \em N is the number of refinement boxes.\n
     * \n
     * Every two successive columns in \em boxes correspond to the coordinates
     * of the lower and upper corners of one refinement box, respectively (see example below).
     * Note that a new knot will be inserted in every knot span contained in this area.
     * If some of the given boxes overlap, the refinement will only be done once.
     *
     * <b>Example</b>, let
     * \verbatim
     d = 2
     knotvector1 = knotvector2 = [ 0 0 0  0.25  0.5  0.75  1 1 1 ]

     boxes = [ 0.25  0.75  0     0.5 ]
     [ 0     0.25  0.75  1   ]
     \endverbatim
     * The areas
     * <em>[ 0.25, 0.75 ] x [ 0, 0.25 ]</em> and
     * <em>[ 0, 0.5 ] x [ 0.75, 1 ]</em> will be refined. The knots \em 0.125, \em 0.375, and \em 0.625 will
     * be inserted in \em knotvector1, and the knots \em 0.125 and \em 0.875 in \em knotvector2.
     *
     * \param[in] boxes gsMatrix of size <em>d</em> x <em>(2*N)</em>;
     * specifies areas for refinement.\n
     * See above for details and format.
     * \param[in] refExt is ignored
     *
     * \remarks NOTE This function directly modifies the basis (by inserting
     * knots in the underlying univariate B-spline bases).
     *
     * \ingroup Nurbs
     *
     *    \note: the \a refExt parameter is ignored in this implementation
     */
    void refine( gsMatrix<T> const & boxes, int refExt = 0) override;

    GISMO_MAKE_GEOMETRY_NEW

    /// \brief Reduces spline continuity (in all directions) at
    /// interior knots by \a i
    void reduceContinuity(int const & i = 1) override
    {
        for (short_t j = 0; j < d; ++j)
            Self_t::component(j).reduceContinuity(i);
    }

    /// \brief Returns span (element) indices of the beginning and end
    /// of the support of the i-th basis function.
    ///
    /// Each direction uses the univariate rule, so along a periodic
    /// direction the result is clipped to the domain [a,b] (the hull of the
    /// support of a crossing function) and never contains ghost spans.
    template <int _Rows>
    void elementSupport_into(const index_t i, gsMatrix<index_t, _Rows, 2> & result) const
    {
        gsMatrix<index_t,1,2> tmp;
        const gsVector<index_t, d> ti = this->tensorIndex(i);

        for (short_t dim = 0; dim < d; ++dim)
        {
            Self_t::component(dim).elementSupport_into(ti[dim], tmp);
            result.row(dim) = tmp;
        }
    }

    /// \brief Returns span (element) indices of the beginning and end
    /// of the support of the i-th basis function.
    gsMatrix<index_t, d, 2> elementSupport(const index_t & i) const
    {
        gsMatrix<index_t, d, 2> result(d, 2);
        elementSupport_into(i, result);
        return result;
    }

    /// \brief Returns the indices of active basis functions in the
    /// given input element box
    ///
    /// On a periodic basis the element indices are domain-relative and the
    /// indices of the periodic direction are folded modulo size(pd).
    template <int _Rows>
    void elementActive_into(const gsMatrix<index_t,_Rows,2> & box,
                             gsMatrix<index_t> & result) const
    {
        GISMO_ASSERT( box.rows() == static_cast<index_t>(d), "Invalid input box");
        gsMatrix<index_t,d,2> tmp;

        gsVector<index_t,d> str;
        this->stride_cwise(str);

        const int pd = periodicDirection();
        if (-1 == pd)
        {
            for (short_t dm = 0; dm != d; ++dm)
            {
                tmp(dm,0) = Self_t::component(dm).knots().lastKnotIndex (box(dm,0)) - this->degree(dm);
                tmp(dm,1) = Self_t::component(dm).knots().firstKnotIndex(box(dm,1)) - 1;
            }

            const gsVector<index_t,d> sz = tmp.col(1)- tmp.col(0) + gsVector<index_t,d>::Ones();

            result = gsVector<index_t>::LinSpaced(sz[0], tmp(0,0), tmp(0,1));
            for (short_t dm = 1; dm != d; ++dm)
                result = result.replicate(1,sz[dm]) +
                    gsVector<index_t>::Constant(result.rows(), str[dm] )
                    * gsVector<index_t>::LinSpaced(sz[dm], tmp(dm,0), tmp(dm,1))
                    .transpose();
            return;
        }

        // One index vector per direction. Element indices count from the
        // domain start a, the knot indices from the first ghost knot.
        std::vector< gsVector<index_t> > idx(d);
        for (short_t dm = 0; dm != d; ++dm)
        {
            const gsKnotVector<T> & kv = Self_t::component(dm).knots();
            if (dm != pd)
            {
                const index_t lo = kv.lastKnotIndex (box(dm,0)) - this->degree(dm);
                const index_t hi = kv.firstKnotIndex(box(dm,1)) - 1;
                idx[dm] = gsVector<index_t>::LinSpaced(hi-lo+1, lo, hi);
            }
            else
            {
                const index_t sh = kv.domainUBegin() - kv.ubegin();
                const index_t lo = kv.lastKnotIndex (box(dm,0)+sh) - this->degree(dm);
                const index_t hi = kv.firstKnotIndex(box(dm,1)+sh) - 1;
                const index_t n  = this->size(dm);
                if (hi-lo+1 >= n)
                    idx[dm] = gsVector<index_t>::LinSpaced(n, 0, n-1);
                else
                {
                    idx[dm].resize(hi-lo+1);
                    for (index_t k = 0; k != idx[dm].size(); ++k)
                        idx[dm][k] = (lo+k) % n;
                }
            }
        }

        result = idx[0];
        for (short_t dm = 1; dm != d; ++dm)
            result = result.replicate(1,idx[dm].size()) +
                gsVector<index_t>::Constant(result.rows(), str[dm] )
                * idx[dm].transpose();
    }

    template <int _Rows>
    gsMatrix<T> elementDom(const gsMatrix<index_t,_Rows,2> & box) const
    {
        GISMO_ASSERT( box.rows() == static_cast<index_t>(d), "Invalid input box");
        gsMatrix<T> rvo;
        rvo.resize(d,2);
        for (short_t dm = 0; dm != d; ++dm)
        {
            rvo(dm,0) = Self_t::component(dm).knots().uValue(box(dm,0));
            rvo(dm,1) = Self_t::component(dm).knots().uValue(box(dm,1));
        }
        return rvo;
    }

    /// Index of the periodic direction, or -1 if no component is periodic.
    ///
    /// Derived from the components' knot vectors on every call (O(d)), so it
    /// follows component edits and swapDirections. At most one periodic
    /// direction is supported; two raise GISMO_ENSURE.
    int periodicDirection() const
    {
        int pd = -1;
        for (short_t i = 0; i != d; ++i)
            if (Self_t::component(i).isPeriodic())
            {
                GISMO_ENSURE(-1 == pd, "gsTensorBSplineBasis: periodicity in more than one direction is not supported.");
                pd = i;
            }
        return pd;
    }

    /// Tells, whether there is a coordinate direction in which the basis is periodic.
    bool isPeriodic() const { return -1 != periodicDirection(); }

    /// Makes the \a dir -th component periodic.
    ///
    /// Returns without change if \a dir is already periodic. Throws
    /// std::runtime_error, leaving the basis unchanged, if another
    /// direction is periodic (checked here) or if the knots of \a dir do
    /// not admit a periodic basis (refused by the 1D component's
    /// setPeriodic).
    void setPeriodic( const int dir )
    {
        GISMO_ENSURE(dir >= 0 && dir < d, "Invalid direction");
        const int pd = periodicDirection();
        if (pd == dir) return;
        GISMO_ENSURE(-1 == pd, "gsTensorBSplineBasis::setPeriodic: the basis is already periodic in another direction.");
        Self_t::component(dir).setPeriodic(true);
    }

    /// Sets the coefficients so that the resulting TensorBSpline is periodic in direction dir.
    ///
    /// Ghost-layout helper: input and output have trueSize() = n+m slices in
    /// \a dir (m = numCrossingFunctions()), and ghost slice n+j takes slice j.
    /// Reduced storage has n slices and needs no such copy. Cost O(rows*cols).
    gsMatrix<T> perCoefs( const gsMatrix<T>& originalCoefs, short_t dir ) const
    {
        gsMatrix<T> result = originalCoefs;
        const index_t m = Self_t::component(dir).numCrossingFunctions();
        if (0 == m) return result;
        gsVector<index_t,d> sz;
        this->size_cwise(sz);
        const index_t n = sz[dir];
        sz[dir] = n + m;
        GISMO_ENSURE(originalCoefs.rows() == sz.prod(), "perCoefs expects trueSize() slices in direction dir");
        index_t str = 1;
        for (short_t k = 0; k < dir; ++k) str *= sz[k];   // stride of dir in the expanded layout
        for (index_t r = 0; r != result.rows(); ++r)
            if ((r / str) % (n + m) >= n)
                result.row(r) = originalCoefs.row(r - n * str);
        return result;
    }

    /// Indices of the boundary functions. The sides of a periodic direction are not boundaries
    /// and are omitted.
    gsMatrix<index_t> allBoundary() const override;

    /// Raises std::runtime_error for a side of the periodic direction.
    gsMatrix<index_t> boundaryOffset(boxSide const & s, index_t offset) const override;

    /// Raises std::runtime_error on a periodic basis, where every corner touches a periodic side.
    index_t functionAtCorner(boxCorner const & c) const override;

private:

    void initFromKnotVectors(std::vector<KnotVectorType> KV)
    {
        GISMO_ENSURE(d == KV.size(), "Invalid number of knot-vectors given." );
        for(short_t i = 0; i!=d; ++i)
            this->m_bases[i] = new Basis_t( give(KV[i]) );
    }

    /// For each knot span, check if its midpoint is
    /// contained in any of the refinement boxes and returns
    /// The mid-point of that span
    std::vector<std::vector<T>> _boxToKnots(gsMatrix<T> const & boxes);

};

#ifdef GISMO_WITH_PYBIND11

    /**
     * @brief Initializes the Python wrapper for the class: gsTensorBSplineBasis
     */
    template <short_t d>
    void pybind11_init_gsTensorBSplineBasis(pybind11::module &m);

    /**
     * @brief Initializes the Python wrapper for the factory constructor of gsTensorBSplineBasis
     */
    void pybind11_init_gsTensorBSplineBasis_factory(pybind11::module &m);


#endif // GISMO_WITH_PYBIND11

} // namespace gismo


// *****************************************************************
#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTensorBSplineBasis.hpp)
/*
#else
#ifdef gsTensorBSplineBasis_EXPORT
#include GISMO_HPP_HEADER(gsTensorBSplineBasis.hpp)
#undef  EXTERN_CLASS_TEMPLATE
#define EXTERN_CLASS_TEMPLATE CLASS_TEMPLATE_INST
#endif
namespace gismo
{
EXTERN_CLASS_TEMPLATE gsTensorBSplineBasis<2,real_t>;
EXTERN_CLASS_TEMPLATE gsTensorBSplineBasis<3,real_t>;
EXTERN_CLASS_TEMPLATE gsTensorBSplineBasis<4,real_t>;
}
*/
#endif
// *****************************************************************
