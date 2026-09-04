/** @file dnk_expr.h

    @brief Defines the dnk expression

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
               H.M. Verhelst
*/

#pragma once

#include <gsUtils/gsCombinatorics.h>
#include <gsTensor/gsTensorBasis.h>

namespace gismo
{
namespace expr
{

/// Number of scalar entries per active basis function in the k-th block of
/// gsFuncData::values, for a parametric domain of dimension \a d: the number
/// of distinct partial derivatives of total order \a k, i.e. the number of
/// weak compositions of \a k into \a d non-negative parts.
///
/// values[1] and values[2] are laid out by gsTensorBasis::evalAllDers_into
/// (pure derivatives, then mixed ones) rather than by the general composition
/// walk used for k>=3, hence the two special cases.
static inline index_t dnk_blockSize(const index_t k, const short_t d)
{
    return 1==k ? d : (2==k ? d*(d+1)/2 : (index_t)numCompositions(k,d));
}

/// Row offset, within the k-th block of gsFuncData::values, of the pure
/// partial derivative @f$ \partial^k / \partial \xi_{dir}^k @f$.
///
/// For k=1,2 the pure derivatives occupy the first \a d rows of their block
/// (see dnk_blockSize()), so the offset is simply \a dir. For k>=3 the
/// composition-walk order used by gsTensorBasis::evalAllDers_into is not a
/// closed form in \a dir, so the offset is found by replicating the same walk
/// (firstComposition/nextComposition, gsUtils/gsCombinatorics.h) and counting
/// steps until the pure composition (0,..,k,..,0) with the \a k in slot
/// \a dir is reached.
static inline index_t dnk_pureOffset(const index_t k, const short_t d, const short_t dir)
{
    if (k <= 2) return dir;
    gsVector<unsigned> cc;
    firstComposition(static_cast<unsigned>(k), d, cc);
    index_t off = 0;
    do { if ( cc[dir] == static_cast<unsigned>(k) ) return off; ++off; }
    while ( nextComposition(cc) );
    GISMO_ERROR("dnk: pure derivative of order "<<k<<" in direction "<<dir<<" not found");
}

/// Largest absolute off-diagonal entry of a square matrix; used to assert
/// that a geometry map's Jacobian is (numerically) diagonal, i.e. that the
/// map is an axis-aligned affine scaling in the face-normal direction.
template<class Mat>
static inline typename Mat::Scalar dnk_offDiagMax(const Mat & J)
{
    return (J - J.diagonal().asDiagonal().toDenseMatrix()).cwiseAbs().maxCoeff();
}

/**
 * @brief Expression for the pure k-th derivative of a finite element variable
 *        in the +\f$\xi_{dir}\f$ parametric direction, scaled to a physical
 *        derivative under an axis-aligned affine geometry map.
 *
 * Building block of the Hoang et al. (2019, CMAME 344) ghost/skeleton
 * penalty, which needs the jump @f$[\![\partial^k_n u]\!]@f$ on background-
 * mesh faces: for a @f$C^{k-1}@f$ spline of degree @f$k@f$, only the k-th
 * normal derivative can jump across an element face, all lower derivatives
 * being continuous by construction. Concretely,
 * @f[
 *   \mathrm{dnk}(u,G,k)_r = \frac{\partial^k u_r/\partial\xi_{dir}^k}{J(dir,dir)^k},
 *   \qquad dir = \text{G.data().side.direction()},\ J = \text{G.data().jacobian()},
 * @f]
 * i.e. the derivative always taken in the *increasing* \f$\xi_{dir}\f$
 * direction, regardless of which of the two opposite sides of the face is
 * being evaluated. This is deliberate, not an oversight: it is exactly what
 * makes \f$\mathrm{dnk}_{\text{left}} - \mathrm{dnk}_{\text{right}}\f$ the
 * correct face jump for two elements sharing a face (their \c boxSide%s have
 * opposite \c parameter()). Consequently \c dnk only equals the physical
 * outward-normal derivative \f$\partial^k u/\partial n^k\f$ up to a sign: on
 * a side with \c parameter()\c ==false (e.g. west/south) the outward normal
 * points in the *decreasing* \f$\xi_{dir}\f$ direction, so
 * \f$\partial^k u/\partial n^k = (-1)^k\,\mathrm{dnk}(u,G,k)\f$ there, while
 * on a \c parameter()\c ==true side (east/north) the two agree exactly; for
 * even \a k the sign cancels and the two always agree.
 * The axis-aligned-affine assumption on \a G is checked by a debug-only
 * \c GISMO_ASSERT (compiled out under \c -DNDEBUG, the configuration this
 * expression is normally built in), so it is not enforced in release builds.
 *
 * For a scalar source (\c dim()==1, e.g. a pressure space or a plain
 * \c gsFeVariable) the result is the \f$na\times 1\f$ column of per-active-
 * function derivatives. For a \a d-component space (\c dim()==d, e.g. a
 * velocity space built via \c getSpace(basis,d)) the result is instead the
 * block-diagonal \f$(na\cdot d)\times d\f$ matrix with the same \f$na\times 1\f$
 * column repeated once per component in column \a c, rows
 * \f$[c\cdot na,(c+1)\cdot na)\f$, zero elsewhere -- the same layout
 * \c symbol_expr::eval() uses for a space's value expression, so that
 * \c dnk(v,...)\ *\ dnk(v,...).tr() assembles the block-diagonal-per-component
 * matrix a d-component ghost penalty needs.
 *
 * @ingroup Expressions
 * @tparam E The expression type
 */
template<class E>
class dnk_expr : public _expr<dnk_expr<E> >
{
    typename E::Nested_t _u;
    typename gsGeometryMap<typename E::Scalar>::Nested_t _G;
    index_t _k;

public:
    typedef typename E::Scalar Scalar;
    enum {Space = E::Space, ScalarValued = 0, ColBlocks = 0};

    dnk_expr(const E & u, const gsGeometryMap<Scalar> & G, const index_t k)
    : _u(u), _G(G), _k(k)
    { GISMO_ENSURE(k >= 1, "dnk: the derivative order must be >= 1 (got "<<k<<")."); }

    mutable gsMatrix<Scalar> res;
    const gsMatrix<Scalar> & eval(const index_t pt) const
    {
        const gsFuncData<Scalar> & fd = _u.data();
        GISMO_ENSURE( 0 != _G.data().side.index(),
                      "dnk: the geometry map has no side set; dnk is only valid "
                      "inside a boundary or face loop." );
        GISMO_ENSURE( _k < (index_t)fd.values.size() && 0 != fd.values[_k].rows(),
                      "dnk: derivatives of order "<<_k<<" were not computed "
                      "(values.size()="<<fd.values.size()<<"); the source's "
                      "compute()/evalAllDers_into does not deliver this order "
                      "(gsGeometry-backed variables fill values[0..2] only)." );

        const short_t d   = _u.parDim();
        const short_t dir = _G.data().side.direction();
        const index_t bsz = dnk_blockSize(_k, d);
        const index_t off = dnk_pureOffset(_k, d, dir);
        const index_t na  = fd.actives.rows();

        // Materialise: gsFuncData::jacobian() returns an Eigen Transpose that nests
        // the temporary gsAsConstMatrix BY REFERENCE (gsFuncData.h:368-373).
        // Storing the view by value would dangle; in -O3 -DNDEBUG that reads
        // garbage instead of crashing.
        const gsMatrix<Scalar> J = _G.data().jacobian(pt);
        GISMO_ASSERT( dnk_offDiagMax(J) <= 1e-10 * (1.0 + J.diagonal().cwiseAbs().maxCoeff()),
                      "dnk assumes an axis-aligned affine geometry map; "
                      "the Jacobian has non-negligible off-diagonal entries:\n"<<J );

        const Scalar s = math::pow(J(dir,dir), (Scalar)_k);
        const index_t dm = _u.dim();
        res.setZero(na*dm, dm);
        for (index_t c = 0; c != dm; ++c)
            for (index_t i = 0; i != na; ++i)
                res(c*na + i, c) = fd.values[_k](i*bsz + off, pt) / s;
        return res;
    }

    index_t rows() const { return _u.data().actives.rows() * _u.dim(); }
    index_t cols() const { return _u.dim(); }

    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        GISMO_ENSURE( 0 != (int)Space || 1 == _u.dim(),
                      "dnk: a vector-valued source is supported for spaces only; "
                      "a function-backed variable carries its components in the rows "
                      "of gsFuncData::values (actives.rows()==1), which the "
                      "block-diagonal layout cannot express." );

        // Fail-fast, patch-0-only sanity check: hierarchical/mapped bases do
        // not implement evalAllDers_into beyond order 2 (the generic
        // gsFunctionSet::evalAllDers_into GISMO_ERRORs for order>2), so a
        // request for k>=3 on such a source would otherwise fail deep inside
        // evaluation instead of here. This inspects only piece(0), so it is a
        // nicety, not a guarantee across all patches; the values.size()/rows()
        // check in eval() is what actually protects every call.
        const gsFunctionSet<Scalar> & src = _u.source();
        const gsFunctionSet<Scalar> & p0  = src.piece(0);
        // Template-argument commas would be mis-parsed as macro-argument
        // separators inside GISMO_ENSURE(...); name the types first.
        typedef gsTensorBasis<1,Scalar> dnk_TB1;
        typedef gsTensorBasis<2,Scalar> dnk_TB2;
        typedef gsTensorBasis<3,Scalar> dnk_TB3;
        typedef gsTensorBasis<4,Scalar> dnk_TB4;
        GISMO_ENSURE( _k <= 2 ||
                      nullptr != dynamic_cast<const dnk_TB1*>(&p0) ||
                      nullptr != dynamic_cast<const dnk_TB2*>(&p0) ||
                      nullptr != dynamic_cast<const dnk_TB3*>(&p0) ||
                      nullptr != dynamic_cast<const dnk_TB4*>(&p0),
                      "dnk(u,G,k) with k>=3 requires a tensor B-spline basis; "
                      "hierarchical/mapped bases do not implement evalAllDers_into "
                      "beyond order 2." );

        evList.add(_u);
        evList.add(_G);
        _u.data().flags |= NEED_ACTIVE | NEED_VALUE | NEED_DERIV_N;
        _u.data().derivOrder = math::max(_u.data().derivOrder, _k);
        _G.data().flags |= NEED_DERIV;
    }

    const gsFeSpace<Scalar> & rowVar() const {return _u.rowVar();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const
    { os << "dnk("; _u.print(os); os << "," << _k << ")"; }
};

/**
 * @brief Expression for the pure k-th normal derivative of a finite element
 *        solution; see dnk_expr for the definition.
 * @ingroup Expressions
 * @tparam T The expression type
 */
template<class T>
class dnk_expr<gsFeSolution<T> > : public _expr<dnk_expr<gsFeSolution<T> > >
{
protected:
    const gsFeSolution<T> _u;
    typename gsGeometryMap<T>::Nested_t _G;
    index_t _k;

public:
    typedef T Scalar;
    enum {Space = 0, ScalarValued = 0, ColBlocks = 0};

    dnk_expr(const gsFeSolution<T> & u, const gsGeometryMap<T> & G, const index_t k)
    : _u(u), _G(G), _k(k)
    { GISMO_ENSURE(k >= 1, "dnk: the derivative order must be >= 1 (got "<<k<<")."); }

    mutable gsMatrix<T> res;
    const gsMatrix<T> & eval(const index_t pt) const
    {
        GISMO_ASSERT(_u.check(), "Invalid state in gsFeSolution");
        const gsFuncData<T> & fd = _u.data();
        GISMO_ENSURE( 0 != _G.data().side.index(),
                      "dnk: the geometry map has no side set; dnk is only valid "
                      "inside a boundary or face loop." );
        GISMO_ENSURE( _k < (index_t)fd.values.size() && 0 != fd.values[_k].rows(),
                      "dnk: derivatives of order "<<_k<<" were not computed "
                      "(values.size()="<<fd.values.size()<<"); the source's "
                      "compute()/evalAllDers_into does not deliver this order "
                      "(gsGeometry-backed variables fill values[0..2] only)." );

        const gsDofMapper & map = _u.mapper();
        res.setZero(_u.dim(), 1); // scalar, but per component

        const short_t d   = _u.parDim();
        const short_t dir = _G.data().side.direction();
        const index_t bsz = dnk_blockSize(_k, d);
        const index_t off = dnk_pureOffset(_k, d, dir);

        const gsMatrix<T> J = _G.data().jacobian(pt);
        GISMO_ASSERT( dnk_offDiagMax(J) <= 1e-10 * (1.0 + J.diagonal().cwiseAbs().maxCoeff()),
                      "dnk assumes an axis-aligned affine geometry map; "
                      "the Jacobian has non-negligible off-diagonal entries:\n"<<J );
        const T s = math::pow(J(dir,dir), (T)_k);

        auto & act = fd.actives.col(1 == fd.actives.cols() ? 0:pt );
        for (index_t c = 0; c!= _u.dim(); c++)
            for (index_t i = 0; i!=fd.actives.rows(); ++i)
            {
                const index_t ii = map.index(act[i], fd.patchId, c);
                const T val = fd.values[_k](i*bsz + off, pt) / s;
                if ( map.is_free_index(ii) ) // DoF value is in the solVector
                    res.at(c) += _u.coefs().at(ii) * val;
                else
                    res.at(c) += _u.fixedPart().at( map.global_to_bindex(ii) ) * val;
            }
        return res;
    }

    index_t rows() const { return _u.dim(); }
    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_u.space());
        evList.add(_G);
        _u.data().flags |= NEED_ACTIVE | NEED_VALUE | NEED_DERIV_N;
        _u.data().derivOrder = math::max(_u.data().derivOrder, _k);
        _G.data().flags |= NEED_DERIV;
    }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<T>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<T>::get();}

    void print(std::ostream &os) const { os << "dnk(s," << _k << ")"; }
};

/**
 * @brief Returns the pure k-th derivative of an expression in the
 *        +\f$\xi_{dir}\f$ parametric direction (\a dir taken from \a G's
 *        side), side-independent -- see dnk_expr for the sign relation to
 *        the physical outward-normal derivative.
 * @ingroup Expressions
 * @param u The expression
 * @param G The geometry map (its side determines the direction)
 * @param k The derivative order
 */
template<class E> EIGEN_STRONG_INLINE
dnk_expr<E> dnk(const symbol_expr<E> & u,
                const gsGeometryMap<typename E::Scalar> & G, const index_t k)
{ return dnk_expr<E>(u, G, k); }

/**
 * @brief Returns the pure k-th derivative of a finite element solution in the
 *        +\f$\xi_{dir}\f$ parametric direction (\a dir taken from \a G's
 *        side), side-independent -- see dnk_expr for the sign relation to
 *        the physical outward-normal derivative.
 * @ingroup Expressions
 * @param u The solution
 * @param G The geometry map (its side determines the direction)
 * @param k The derivative order
 */
template<class T> EIGEN_STRONG_INLINE
dnk_expr<gsFeSolution<T> > dnk(const gsFeSolution<T> & u,
                               const gsGeometryMap<T> & G, const index_t k)
{ return dnk_expr<gsFeSolution<T> >(u, G, k); }

}// namespace expr
}// namespace gismo
