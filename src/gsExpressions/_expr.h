/** @file _expr.h

    @brief Defines the expression for a constant value

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
    \brief Base class for all expressions
    \tparam E The expression type
    \ingroup Expressions
*/
template <typename E>
class _expr<E, false>
{
protected://private:
    _expr(){}
    _expr(const _expr&) { }
public:
    // Defined in derived classes: enum { Space, ScalarValued, ColBlocks }
    // - ScalarValued: 0 is a scalar (must have Space=0),1 one denotes gsMatrix
    // - ColBlocks: the expression stacks matrices per basis function
    // - Space: 0: not a trial nor a test object (eg. normal vector, force function)
    //          1: a test object  (essentially a right-hand side vector expression)
    //          2: a trial object
    //          3: a trial+trial object (essentially a matrix expression)

    typedef typename expr_traits<E>::Nested_t Nested_t;
    typedef typename expr_traits<E>::Scalar   Scalar;

    /// Prints the expression as a string to \a os
    void print(std::ostream &os) const
    {
        //gsInfo<<"\n Space="<<E::Space<<", ScV="<<E::ScalarValued<<", ColBlocks="<<E::ColBlocks<<"\n";
        static_cast<E const&>(*this).print(os);
        os<<"\n";
        /*
          std::string tmp(__PRETTY_FUNCTION__);
          tmp.erase(0,74);
          tmp.erase(tmp.size()-42,42);
          size_t pos = 0;
          while((pos=tmp.find(", false",0))!=std::string::npos) tmp.erase(pos,7);
          while((pos=tmp.find(", true",0))!=std::string::npos) tmp.erase(pos,6);
          while((pos=tmp.find("gismo::expr::",0))!=std::string::npos) tmp.erase(pos,13);
          while((pos=tmp.find("_expr",0))!=std::string::npos) tmp.erase(pos,5);
          while((pos=tmp.find("<double>",0))!=std::string::npos) tmp.erase(pos,8);
          // while((pos=tmp.find("<long double>",0))!=std::string::npos) tmp.erase(pos,13);
          // while((pos=tmp.find("<float>",0))!=std::string::npos) tmp.erase(pos,7);
          tmp.erase(std::remove_if(tmp.begin(),tmp.end(),::isspace),tmp.end());
          os<<tmp<<"\n";
        */
    }

    std::ostream & printDetail(std::ostream &os) const
    {
        os << (isVectorTr() ? "VectorTr " :
               (isVector() ? "Vector " :
                (isMatrix() ? "Matrix " :
                 "Scalar ") ) )
           <<"expression of size "<< rows() // bug: this might be invalid if unparsed
           << " x "<<cols()<<"\n";
        print(os);
        return os;
    }

    /// Evaluates the expression at evaluation point indexed by \a k
    MatExprType eval(const index_t k) const
    { return static_cast<E const&>(*this).eval(k); }

    /// Returns the transpose of the expression
    transpose_expr<E> tr() const
    { return transpose_expr<E,false>(static_cast<E const&>(*this)); }

    /// Returns the coordinate-wise transpose of the expression
    transpose_expr<E,true> cwisetr() const
    { return transpose_expr<E,true>(static_cast<E const&>(*this)); }

    /// Returns the puts the expression to colBlocks
    colBlocks_expr<E> cb() const
    { return colBlocks_expr<E>(static_cast<E const&>(*this)); }

    /// Returns the sign of the expression
    sign_expr<E> sgn(Scalar tolerance=0) const
    { return sign_expr<E>(static_cast<E const&>(*this), tolerance); }

    /// Returns exp(expression)
    exp_expr<E> exp() const
    { return exp_expr<E>(static_cast<E const&>(*this)); }

    /// Returns the expression's positive part
    ppart_expr<E> ppart() const
    { return ppart_expr<E>(static_cast<E const&>(*this)); }
    ppartval_expr<E> ppartval() const
    { return ppartval_expr<E>(static_cast<E const&>(*this)); }

    /// Returns the expression's negative part
    mult_expr<real_t, ppart_expr<mult_expr<double,E,false>> , false>
    npart() const { return -1* ( -(*this) ).ppart() ; }

    /// Returns an evaluation of the (sub-)expression in temporary memory
    temp_expr<E> temp() const
    { return temp_expr<E>(static_cast<E const&>(*this)); }

    /// Returns the inverse of the expression (for matrix-valued expressions)
    inv_expr<E> const inv() const
    { return inv_expr<E>(static_cast<E const&>(*this)); }

    /// Returns the trace of the expression (for matrix-valued expressions)
    trace_expr<E> trace() const
    { return trace_expr<E>(static_cast<E const&>(*this)); }

    /// Returns the adjugate of the expression (for matrix-valued expressions)
    adjugate_expr<E> adj() const
    { return adjugate_expr<E>(static_cast<E const&>(*this)); }

    /// Returns the Euclidean norm of the expression
    norm_expr<E> norm() const
    { return norm_expr<E>(static_cast<E const&>(*this)); }

    /// Returns the vector normalized to unit length
    normalized_expr<E> normalized() const
    { return normalized_expr<E>(static_cast<E const&>(*this)); }

    /// Returns the determinant of the expression
    det_expr<E> det() const
    { return det_expr<E>(static_cast<E const&>(*this)); }

    /// Returns the squared Euclidean norm of the expression
    sqNorm_expr<E> sqNorm() const
    { return sqNorm_expr<E>(static_cast<E const&>(*this)); }

    /// Returns the square root of the expression (component-wise)
    mult_expr<E,E,0> (sqr)() const { return (*this)*(*this); }

    symm_expr<E> symm() const
    { return symm_expr<E>(static_cast<E const&>(*this)); }

    symmetrize_expr<E> symmetrize() const
    { return symmetrize_expr<E>(static_cast<E const&>(*this)); }

    /// For matrix-valued expressions which are actually 1x1 matrix,
    /// returns a scalar valued expression
    value_expr<E> val() const
    { return value_expr<E>(static_cast<E const&>(*this)); }

    /// Returns a diagonal matrix expression of the vector expression
    asdiag_expr<E> asDiag() const
    { return asdiag_expr<E>(static_cast<E const&>(*this)); }

    /// Returns the rowSum of a matrix
    max_expr<E> max() const
    { return max_expr<E>(static_cast<E const&>(*this)); }

    /// Returns the rowSum of a matrix
    rowsum_expr<E> rowSum() const
    { return rowsum_expr<E>(static_cast<E const&>(*this)); }

    /// Returns the colSum of a matrix
    colsum_expr<E> colSum() const
    { return colsum_expr<E>(static_cast<E const&>(*this)); }

    col_expr<E> operator[](const index_t i) const
    { return col_expr<E>(static_cast<E const&>(*this),i); }

    /// Returns the row-size of the expression
    index_t rows() const
    { return static_cast<E const&>(*this).rows(); }

    /// Returns the column-size of the expression
    index_t cols() const
    { return static_cast<E const&>(*this).cols(); }

    index_t cardinality() const
    { return static_cast<E const&>(*this).cardinality_impl(); }

    static index_t cardinality_impl() { return 1; }

    ///\brief Returns true iff the expression is scalar-valued.
    /// \note This is a runtime check, for compile-time check use E::ScalarValued
    bool isScalar() const { return rows()*cols()<=1; } //!rowSpan && !colSpan

    static constexpr bool isVector  () { return 1==E::Space; }
    static constexpr bool isVectorTr() { return 2==E::Space; }
    static constexpr bool isMatrix  () { return 3==E::Space; }

    ///\brief Parse the expression and discover the list of evaluation
    ///sources, also sets the required evaluation flags
    void parse(gsExprHelper<Scalar> & evList) const
    { static_cast<E const&>(*this).parse(evList); }

    template<class op> void apply(op & _op) const
    { static_cast<E const&>(*this).apply(_op); }

    /// Returns the space that is found on the left-most of the
    /// expression
    const gsFeSpace<Scalar> & rowVar() const
    {
        // assert ValueType!=0
        return static_cast<E const&>(*this).rowVar();
    }

    /// Returns the space that is found on the right-most of
    /// the expression
    const gsFeSpace<Scalar> & colVar() const
    {
        // assert ValueType==2
        return static_cast<E const&>(*this).colVar();
    }

    // Overload conversions, eg. converts _expr<mult_expr> to
    // mult_expr.
    operator E&()             { return static_cast<      E&>(*this); }
    operator E const&() const { return static_cast<const E&>(*this); }

    E const & derived() const { return static_cast<const E&>(*this); }
};

/// Stream operator for expressions
template <typename E>
std::ostream &operator<<(std::ostream &os, const _expr<E> & b)
{b.print(os); return os; }

/*
    \brief Expression for a constant value
    \ingroup Expressions
    \tparam T The type of the constant value
*/
template<class T>
class _expr<T, true> : public _expr<_expr<T> >
{
    const T _c;
public:
    typedef T Scalar;
    typedef const _expr<T> Nested_t;

    explicit _expr(Scalar c) : _c(give(c)) { }

public:
    enum {Space = 0, ScalarValued = 1, ColBlocks= 0};

    inline Scalar eval(const index_t = 0) const { return _c; }

    inline const _expr<T, true> & val() const { return *this; }
    index_t rows() const { return 0; }
    index_t cols() const { return 0; }
    void parse(gsExprHelper<Scalar> &) const { }

    const gsFeSpace<T> & rowVar() const { return gsNullExpr<T>::get(); }
    const gsFeSpace<T> & colVar() const { return gsNullExpr<T>::get(); }

    void print(std::ostream &os) const { os<<_c; }
};

}// namespace expr
}// namespace gismo
