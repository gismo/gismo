/** @file cross_expr.h

    @brief Defines the cross expression

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
               H.M. Verhelst
*/

#pragma once

namespace gismo
{
namespace expr
{

/**
 * @brief Expression for the cross product of two expressions
 * @ingroup Expressions
 * @tparam E1 The first expression type
 * @tparam E2 The second expression type
 *
 * This class is a partial specialization for the case where the first
 * expression is not a column block.
 * B x [A1 A2 A3] = [BxA1  BxA2  BxA3]
 */
template <typename E1, typename E2>
class cross_expr<E1,E2,false> : public _expr<cross_expr<E1, E2, false> >
{
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;

public:
    enum {ScalarValued = E1::ScalarValued && E2::ScalarValued,
        ColBlocks = E2::ColBlocks};
    enum {Space = (int)E1::Space + (int)E2::Space };

    typedef typename E1::Scalar Scalar;

    cross_expr(_expr<E1> const& u,
              _expr<E2> const& v)
    : _u(u), _v(v) { }

    mutable Temporary_t tmp;
    const Temporary_t & eval(const index_t k) const
    {
        GISMO_ASSERT(_u.rows() == 3 || _u.cols() == 3, "cross(.) requires 3D variable.");
        GISMO_ASSERT(_v.rows() == 3 || _v.cols() == 3, "cross(.) requires 3D variable.");

        tmp = _u.eval(k).cross(_v.eval(k));
        return tmp; // assumes result is not scalarvalued
    }

    typename E1::Nested_t const & left() const { return _u; }
    typename E2::Nested_t const & right() const { return _v; }

    index_t rows() const { return _u.rows(); }
    index_t cols() const { return _u.cols(); }
    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); _v.parse(evList); }


    index_t cardinality_impl() const
    { return 0==E1::Space ? _v.cardinality(): _u.cardinality(); }

    const gsFeSpace<Scalar> & rowVar() const
    { return 0==E1::Space ? _v.rowVar() : _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const
    { return 0==E2::Space ? _u.colVar() : _v.colVar(); }

    void print(std::ostream &os) const { _u.print(os); os<<"x"; _v.print(os); }
};

/*
  Expression for cross-product operation (second version)

  First argument E1 has ColBlocks = true

  Partial specialization for (right) blockwise cross-product
  [A1 A2 A3] x B = [A1xB  A2xB  A3xB]

  as well as

  both are ColBlocks: [A1 A2 A3] x [B1 B2 B3] = [A1xB1  A2xB2  A3xB3]
                                                [A2xB1 ..           ]
                                                [                   ]
*/

/**
 * @brief Expression for the cross product of two expressions
 * @ingroup Expressions
 * @tparam E1 The first expression type
 * @tparam E2 The second expression type
 *
 * This class is a partial specialization for the case where the first
 * expression is a column block.
 *
 * Partial specialization for (right) blockwise cross-product
 * [A1 A2 A3] x B = [A1xB  A2xB  A3xB]
 *
 * as well as
 *
 * both are ColBlocks: [A1 A2 A3] x [B1 B2 B3] = [A1xB1  A2xB2  A3xB3]
 *                                               [A2xB1 ..           ]
 *                                               [                   ]
 */
template <typename E1, typename E2>
class cross_expr<E1, E2, true> : public _expr<cross_expr<E1, E2, true> >
{
public:
    typedef typename E2::Scalar Scalar;
private:
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;

    mutable gsMatrix<Scalar> res;
public:
    enum {ScalarValued = 0, ColBlocks = E1::ColBlocks}; //(!)
    enum {Space = (int)E1::Space + (int)E2::Space };

    cross_expr(_expr<E1> const& u,
              _expr<E2> const& v)
    : _u(u), _v(v)
    {
        GISMO_ASSERT(3==u.rows() || 3==u.cols(), "cross(.) requires 3D variable.");
        GISMO_ASSERT(3==v.rows() || 3==v.cols(), "cross(.) requires 3D variable.");
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        const index_t uc = _u.cols();
        const index_t ur = _u.rows();
        const index_t nb = _u.cardinality();
        const auto tmpA = _u.eval(k);
        const auto tmpB = _v.eval(k);

        const index_t vc = _v.cols();

        // either _v.cardinality()==1 or _v.cardinality()==_u.cardinality()
        if (  1 == _v.cardinality() ) //second is not ColBlocks
        {
            res.resize(ur, vc*nb);
            GISMO_ASSERT(tmpA.cols()==uc*nb, "Dimension error.. "<< tmpA.cols()<<"!="<<uc*nb );
            //gsInfo<<"cols = "<<res.cols()<<"; rows = "<<res.rows()<<"\n";
            for (index_t i = 0; i!=nb; ++i)
                res.middleCols(i*vc,vc).noalias()
                    = tmpA.middleCols(i*uc,uc).cross(tmpB);
        }
        // both are ColBlocks: [A1 A2 A3] x [B1 B2 B3] = [A1xB1  A2xB2  A3xB3]
        //                                               [A2xB1 ..           ]
        //                                               [                   ]
        else
        {
            const index_t nbv = _v.cardinality();
            res.resize(ur*nb, vc*nbv);
            for (index_t i = 0; i!=nb; ++i)
                for (index_t j = 0; j!=nbv; ++j)
                {
                    res.block(i*ur,j*vc,ur,vc).noalias() =
                        tmpA.middleCols(i*uc,uc).cross(tmpB.middleCols(j*vc,vc));
                    // res.middleCols(i*vc,vc).noalias()
                    //     = tmpA.middleCols(i*uc,uc) * tmpB.middleCols(i*vc,vc);
                }
        }
        return res;
    }

    typename E1::Nested_t const & left() const { return _u; }
    typename E2::Nested_t const & right() const { return _v; }

    index_t rows() const {
        return _u.rows();
    }
    index_t cols() const {
        return _v.cols();
    }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); _v.parse(evList); }


    index_t cardinality_impl() const { return  _u.cardinality(); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const
    {
        if ( 1 == _v.cardinality() )
            return _u.colVar();
        else
            return _v.colVar();
    }

    void print(std::ostream &os) const
    { os << "("; _u.print(os);os <<"x";_v.print(os);os << ")"; }
};

}// namespace expr
}// namespace gismo