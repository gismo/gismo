/** @file mult_expr.h

    @brief Defines the mult expression

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
 * @brief Expression for the multiplication operation (first version)
 *        First argument E1 has ColBlocks = false
 *        Partial specialization for (right) blockwise multiplication:
 *        B * [A1 A2 A3] = [B*A1  B*A2  B*A3]
 * @ingroup Expressions
 * @tparam E1 the type of the first expression
 * @tparam E2 the type of the second expression
 */
template <typename E1, typename E2>
class mult_expr<E1,E2,false> : public _expr<mult_expr<E1, E2, false> >
{
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;

public:
    enum {ScalarValued = E1::ScalarValued && E2::ScalarValued,
        ColBlocks = E2::ColBlocks};
    enum {Space = (int)E1::Space + (int)E2::Space };

    typedef typename E1::Scalar Scalar;

    mult_expr(_expr<E1> const& u,
              _expr<E2> const& v)
    : _u(u), _v(v) { }

    mutable Temporary_t tmp;
    const Temporary_t & eval(const index_t k) const
    {
        GISMO_ASSERT(0==_u.cols()*_v.rows() || _u.cols() == _v.rows(),
                     "Wrong dimensions "<<_u.cols()<<"!="<<_v.rows()<<" in * operation:\n"
                     << _u <<" times \n" << _v );

        // Note: a * b * c --> (a*b).eval()*c
        tmp = _u.eval(k) * _v.eval(k);
        return tmp; // assumes result is not scalarvalued
    }

    typename E1::Nested_t const & first() const { return _u; }
    typename E2::Nested_t const & second() const { return _v; }

    index_t rows() const { return E1::ScalarValued ? _v.rows()  : _u.rows(); }
    index_t cols() const { return E2::ScalarValued ? _u.cols()  : _v.cols(); }
    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); _v.parse(evList); }


    index_t cardinality_impl() const
    { return 0==E1::Space ? _v.cardinality(): _u.cardinality(); }

    const gsFeSpace<Scalar> & rowVar() const
    { return 0==E1::Space ? _v.rowVar() : _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const
    { return 0==E2::Space ? _u.colVar() : _v.colVar(); }

    void print(std::ostream &os) const { _u.print(os); os<<"*"; _v.print(os); }
};

/**
 * @brief Expression for the multiplication operation (second version)
 *        First argument E1 has ColBlocks = true
 *        Partial specialization for (right) blockwise multiplication:
 *        [A1 A2 A3] * B = [A1*B  A2*B  A3*B]
 *        As well as for when both are ColBlocks:
 *        [A1 A2 A3] * [B1 B2 B3] = [A1*B1  A2*B2  A3*B3]
 *                                   [A2*B1  ...        ]
 *                                   [                  ]
 * @ingroup Expressions
 * @tparam E1 the type of the first expression
 * @tparam E2 the type of the second expression
 */
template <typename E1, typename E2>
class mult_expr<E1, E2, true> : public _expr<mult_expr<E1, E2, true> >
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

    mult_expr(_expr<E1> const& u,
              _expr<E2> const& v)
    : _u(u), _v(v)
    {

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
                    = tmpA.middleCols(i*uc,uc) * tmpB;
        }
        // both are ColBlocks: [A1 A2 A3] * [B1 B2 B3] = [A1*B1  A2*B2  A3*B3]
        //                                               [A2*B1 ..           ]
        //                                               [                   ]
        else
        {
            const index_t nbv = _v.cardinality();
            res.resize(ur*nb, vc*nbv);
            for (index_t i = 0; i!=nb; ++i)
                for (index_t j = 0; j!=nbv; ++j)
                {
                    res.block(i*ur,j*vc,ur,vc).noalias() =
                        tmpA.middleCols(i*uc,uc) * tmpB.middleCols(j*vc,vc);
                    // res.middleCols(i*vc,vc).noalias()
                    //     = tmpA.middleCols(i*uc,uc) * tmpB.middleCols(i*vc,vc);
                }
        }
        return res;
    }

    typename E1::Nested_t const & first() const { return _u; }
    typename E2::Nested_t const & second() const { return _v; }

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
    { os << "("; _u.print(os);os <<"*";_v.print(os);os << ")"; }
};

/**
 * @brief Expression for multiplication operation (third version)
 *        Scalar version
 * @ingroup Expressions
 * @tparam E2 the type of the expression
 */
template <typename E2>
class mult_expr<_expr<typename E2::Scalar,true>, E2, false>
    : public _expr<mult_expr<_expr<typename E2::Scalar,true>, E2, false> >
// template <typename E> class scmult_expr : public _expr<scmult_expr<E> >
{
public:
    typedef typename E2::Scalar Scalar;
private:
    _expr<Scalar,true>    _c;
    typename E2::Nested_t _v;

    //mult_expr(const mult_expr&);
public:
    enum {ScalarValued = E2::ScalarValued, ColBlocks = E2::ColBlocks};
    enum {Space = E2::Space};

    mult_expr(Scalar const & c, _expr<E2> const& v)
    : _c(c), _v(v) { }

    mult_expr(_expr<Scalar,true> const & c, _expr<E2> const& v)
    : _c(c), _v(v) { }

    EIGEN_STRONG_INLINE AutoReturn_t eval(const index_t k) const
    {
        return ( _c.eval(0) * _v.eval(k) );
    }

    _expr<Scalar,true>    const & first()  const { return _c; }
    typename E2::Nested_t const & second() const { return _v; }

    index_t rows() const { return _v.rows(); }
    index_t cols() const { return _v.cols(); }

    void parse(gsExprHelper<Scalar> & evList) const
    { _v.parse(evList); }

    index_t cardinality_impl() const
    { return _v.cardinality(); }

    const gsFeSpace<Scalar> & rowVar() const { return _v.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _v.colVar(); }

    void print(std::ostream &os) const { os << _c <<"*";_v.print(os); }
};

/// Multiplication operator for expressions

/**
 * @brief Multiplication operator for expressions
 * @ingroup Expressions
 * @param u The first expression
 * @param v The second expression
 */
template <typename E1, typename E2> EIGEN_STRONG_INLINE
mult_expr<E1,E2> const operator*(_expr<E1> const& u, _expr<E2> const& v)
{ return mult_expr<E1, E2>(u, v); }

/**
 * @brief Multiplication operator for expressions
 * @ingroup Expressions
 * @param u The first expression
 * @param v The second expression
 */
template <typename E2> EIGEN_STRONG_INLINE
mult_expr<_expr<typename E2::Scalar,true>,E2,false> const
operator*(typename E2::Scalar const& u, _expr<E2> const& v)
{ return mult_expr<_expr<typename E2::Scalar,true>, E2, false>(u, v); }

/**
 * @brief Multiplication operator for expressions
 * @ingroup Expressions
 * @param u The first expression
 * @param v The second expression
 */
template <typename E2> EIGEN_STRONG_INLINE
mult_expr<_expr<typename E2::Scalar,true>,E2,false> const
operator*(_expr<typename E2::Scalar,true> const& u, _expr<E2> const& v)
{ return mult_expr<_expr<typename E2::Scalar,true>, E2, false>(u, v); }

/**
 * @brief Multiplication operator for expressions
 * @ingroup Expressions
 * @param u The first expression
 * @param v The second expression
 */
template <typename E1> EIGEN_STRONG_INLINE
mult_expr<_expr<typename E1::Scalar,true>,E1,false> const
operator*(_expr<E1> const& v, typename E1::Scalar const& u)
{ return mult_expr<_expr<typename E1::Scalar,true>,E1, false>(u, v); }

/**
 * @brief Multiplication operator for expressions
 * @ingroup Expressions
 * @param u The first expression
 * @param v The second expression
 */
template <typename E1> EIGEN_STRONG_INLINE
mult_expr<_expr<typename E1::Scalar,true>,E1,false> const
operator*(_expr<E1> const& v, _expr<typename E1::Scalar,true> const& u)
{ return mult_expr<_expr<typename E1::Scalar,true>,E1, false>(u, v); }

/**
 * @brief Multiplication operator for expressions
 * @ingroup Expressions
 * @param u The first expression
 * @param v The second expression
 */
template <typename E1> EIGEN_STRONG_INLINE
mult_expr<_expr<typename E1::Scalar,true>,E1,false> const
operator-(_expr<E1> const& u)
{ return mult_expr<_expr<typename E1::Scalar,true>,E1, false>(-1, u); }

/**
 * @brief Multiplication operator for expressions
 * @ingroup Expressions
 * @param u The first expression
 * @param v The second expression
 */
template <typename E> mult_expr<constMat_expr, E> const
operator*( gsMatrix<typename E::Scalar> const& u, _expr<E> const& v)
{ return mult_expr<constMat_expr, E>(mat(u), v); }

/**
 * @brief Multiplication operator for expressions
 * @ingroup Expressions
 * @param u The first expression
 * @param v The second expression
 */
template <typename E> mult_expr<E, constMat_expr> const
operator*(_expr<E> const& u, gsMatrix<typename E::Scalar> const& v)
{ return mult_expr<E, constMat_expr>(u, mat(v) ); }

}// namespace expr
}// namespace gismo