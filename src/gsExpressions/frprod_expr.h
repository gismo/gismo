/** @file frprod_expr.h

    @brief Defines the frprod expression

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

/*
  Expression for the Frobenius matrix (or double dot) product (first
  version) Also block-wise

  [A1 A2 A3] . [B1 B2 B3]
  =
  [ A1.B1  A1.B2  A1.B3 ]
  [ A2.B1  A2.B2  A2.B3 ]
  [ A3.B1  A3.B2  A3.B3 ]
*/

/**
 * @brief Expression for the Frobenius product of two expressions (first version)
 *
 *        This expression works block-wise
 *          [A1 A2 A3] . [B1 B2 B3]
 *          =
 *          [ A1.B1  A1.B2  A1.B3 ]
 *          [ A2.B1  A2.B2  A2.B3 ]
 *          [ A3.B1  A3.B2  A3.B3 ]
 */
template <typename E1, typename E2, bool = E2::ColBlocks>
class frprod_expr : public _expr<frprod_expr<E1, E2> >
{
public:
    typedef typename E1::Scalar Scalar;
    enum {ScalarValued = 0, ColBlocks=E2::ColBlocks};
    enum { Space = (int)E1::Space + (int)E2::Space };
    // E1 E2 this (16 cases..)
    // 0  0  0
    // 1  1
    // 2  2
    // 3  3

private:
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;

    mutable gsMatrix<Scalar> res;

public:

    frprod_expr(_expr<E1> const& u, _expr<E2> const& v)
    : _u(u), _v(v)
    {
        //todo: add check() functions, which will evaluate expressions on an empty matrix (no points) to setup initial dimensions ???
        //GISMO_ASSERT(_u.rows() == _v.rows(),
        //             "Wrong dimensions "<<_u.rows()<<"!="<<_v.rows()<<" in % operation");
        //GISMO_ASSERT(_u.cols() == _v.cols(),
        //             "Wrong dimensions "<<_u.cols()<<"!="<<_v.cols()<<" in % operation");
    }

    const gsMatrix<Scalar> & eval(const index_t k) const //todo: specialize for nb==1
    {
        // The blocks are cols() wide, which differs from rows() unless
        // they are square (e.g. jac(u) of a space whose dim() is not the
        // parametric dimension), and the two sides may have different
        // numbers of blocks (test and trial spaces)
        const index_t cb = _u.cols();
        const index_t nb = _u.cardinality();
        const index_t nbv = _v.cardinality();
        auto A = _u.eval(k);
        auto B = _v.eval(k);
        GISMO_ASSERT(A.cols()==cb*nb && B.cols()==_v.cols()*nbv,
                     "Dimension error in % operation");
        GISMO_ASSERT(A.rows()==B.rows() && cb==_v.cols(), "Wrong block dimensions "
                     <<A.rows()<<"x"<<cb<<" and "<<B.rows()<<"x"<<_v.cols()<<" in % operation");
        res.resize(nb, nbv);
        for (index_t i = 0; i!=nb; ++i) // all with all
            for (index_t j = 0; j!=nbv; ++j)
                res(i,j) =
                    (A.middleCols(i*cb,cb).array() * B.middleCols(j*cb,cb).array()).sum();
        return res;
    }

    // One scalar per pair of basis functions, and a single scalar if
    // neither side is a space (eval() is 1 x 1 then)
    index_t rows() const { return 1; }
    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); _v.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _v.colVar(); }

    void print(std::ostream &os) const
    { os << "("; _u.print(os); os<<" % "; _v.print(os); os<<")";}
};

/**
 * @brief Expression for the Frobenius product of two expressions (second version)
 *        When left hand only side is block-wise
 *        [A1 A2 A3] : B = [A1:B  A2:B  A3:B]
 * @ingroup Expressions
 * @tparam E1 The first expression type
 * @tparam E2 The second expression type
 */
template <typename E1, typename E2>
class frprod_expr<E1,E2,false> : public _expr<frprod_expr<E1, E2,false> >
{
public:
    typedef typename E1::Scalar Scalar;
    enum {ScalarValued = 0, Space = E1::Space, ColBlocks= E1::ColBlocks};

private:
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;

    mutable gsMatrix<Scalar> res;

public:

    frprod_expr(_expr<E1> const& u, _expr<E2> const& v)
    : _u(u), _v(v)
    {
        // gsInfo << "expression is space ? "<<E1::Space <<"\n"; _u.print(gsInfo);
        // GISMO_ASSERT(_u.rows() == _v.rows(),
        //              "Wrong dimensions "<<_u.rows()<<"!="<<_v.rows()<<" in % operation");
        // GISMO_ASSERT(_u.cols() == _v.cols(),
        //              "Wrong dimensions "<<_u.cols()<<"!="<<_v.cols()<<" in % operation");
    }

    const gsMatrix<Scalar> & eval(const index_t k) const //todo: specialize for nb==1
    {
        auto A = _u.eval(k);
        auto B = _v.eval(k);
        const index_t cb = _u.cols(); // block width, see above
        const index_t nb = _u.cardinality();
        GISMO_ASSERT(A.cols()==cb*nb, "Dimension error in % operation");
        GISMO_ASSERT(A.rows()==B.rows() && cb==B.cols(), "Wrong block dimensions "
                     <<A.rows()<<"x"<<cb<<" and "<<B.rows()<<"x"<<B.cols()<<" in % operation");
        res.resize(nb, 1);
        for (index_t i = 0; i!=nb; ++i) // all with all
            res(i,0) =
                (A.middleCols(i*cb,cb).array() * B.array()).sum();
        return res;
    }

    // One scalar per basis function, or a single scalar (see above)
    index_t rows() const { return 1; }
    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); _v.parse(evList); }


    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _v.rowVar(); }

    void print(std::ostream &os) const
    { os << "("; _u.print(os); os<<" % "; _v.print(os); os<<")";}
};

/**
 * @brief Returns the Frobenius product of two expressions
 * @param u The first expression
 * @param v The second expression
 * @ingroup Expressions
 */
template <typename E1, typename E2> EIGEN_STRONG_INLINE
frprod_expr<E1,E2> const  operator%(_expr<E1> const& u, _expr<E2> const& v)
{ return frprod_expr<E1, E2>(u, v); }

}// namespace expr
}// namespace gismo