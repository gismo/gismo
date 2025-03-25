/** @file gsExpressions.h

    @brief Defines different expressions

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsFuncData.h>
#include <gsMSplines/gsMappedBasis.h>

namespace gismo
{

// Adaptor to compute Hessian
template <typename Derived>
void secDerToHessian(const gsEigen::DenseBase<Derived> &  secDers,
                     const index_t dim,
                     gsMatrix<typename Derived::Scalar> & hessian)
{
    const index_t sz = dim*(dim+1)/2;
    const gsAsConstMatrix<typename Derived::Scalar>
        ders(secDers.derived().data(), sz, secDers.size() / sz );
    hessian.resize(dim*dim, ders.cols() );

    switch ( dim )
    {
    case 1:
        hessian = secDers.transpose(); //==ders
        break;
    case 2:
        hessian.row(0)=ders.row(0);//0,0
        hessian.row(1)=//1,0
            hessian.row(2)=ders.row(2);//0,1
        hessian.row(3)=ders.row(1);//1,1
        break;
    case 3:
        hessian.row(0)=ders.row(0);//0,0
        hessian.row(3)=//0,1
            hessian.row(1)=ders.row(3);//1,0
        hessian.row(6)=//0,2
            hessian.row(2)=ders.row(4);//2,0
        hessian.row(4)=ders.row(1);//1,1
        hessian.row(7)=//1,2
            hessian.row(5)=ders.row(5);//2,1
        hessian.row(8)=ders.row(2);//2,2
        break;
    default:
        break;
    }
}



// Forward declaration in gismo namespace
template<class T> class gsExprHelper;

/** @namespace gismo::expr

    @brief
    This namespace contains expressions used for FE computations

    \ingroup Assembler
*/
namespace expr
{

template <class E> struct is_arithmetic{enum{value=0};};
template <> struct is_arithmetic<real_t>{enum{value=1};};
template <typename E, bool = is_arithmetic<E>::value >
class _expr {using E::GISMO_ERROR_expr;};

template<class E> class symbol_expr;
template<class T> class gsNullExpr;
template<class T> class gsFeSpace;
template<class T> class gsFeElement;
template<class T> class gsFeVariable;
template<class T> class gsComposition;
template<class T> class gsFeSolution;
template<class E> class symm_expr;
template<class E> class symmetrize_expr;
template<class E> class normalized_expr;
template<class E> class trace_expr;
template<class E> class integral_expr;
template<class E> class adjugate_expr;
template<class E> class norm_expr;
template<class E> class sqNorm_expr;
template<class E> class det_expr;
template<class E> class value_expr;
template<class E> class asdiag_expr;
template<class E> class max_expr;
template<class E> class rowsum_expr;
template<class E> class colsum_expr;
template<class E> class col_expr;
template<class T> class meas_expr;
template<class E> class inv_expr;
template<class E, bool cw = false> class transpose_expr;
template<class E> class colBlocks_expr;
template<class E> class abs_expr;
template<class E> class pow_expr;
template<class E> class sign_expr;
template<class E> class ppart_expr;
template<class E> class exp_expr;
template<class E> class ppartval_expr;
template<class T> class cdiam_expr;
template<class E> class temp_expr;
template<class E1, class E2, bool = E1::ColBlocks && !E1::ScalarValued && !E2::ScalarValued> class mult_expr
{using E1::GISMO_ERROR_mult_expr_has_invalid_template_arguments;};

// Call as pow(a,b)
template<class E> pow_expr<E>
pow(_expr<E> const& u, real_t q) { return pow_expr<E>(u,q); }

/*
  Traits class for expressions
*/
template <typename E> struct expr_traits
{
public:
//    typedef typename E::Scalar Scalar;
    typedef real_t Scalar;//todo
    typedef const E Nested_t;
};

#  define Temporary_t typename util::conditional<ScalarValued,Scalar,   \
        typename gsMatrix<Scalar>::Base >::type
#if __cplusplus >= 201402L || _MSVC_LANG >= 201402L // c++14
#  define MatExprType  auto
#  define AutoReturn_t auto
//note: in c++11 auto-return requires -> decltype(.)
#else // 199711L, 201103L
#  define MatExprType typename gsMatrix<real_t>::constRef
#  define AutoReturn_t typename util::conditional<ScalarValued,real_t,MatExprType>::type
#endif

/**
   \brief Base class for all expressions
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

}
}

#include <gsExpressions/precomputed_expr.h>

#include <gsExpressions/gsFeSpace.h>
#include <gsExpressions/gsFeVariable.h>
#include <gsExpressions/gsFeSolution.h>
#include <gsExpressions/gsFeElement.h>
#include <gsExpressions/gsComposition.h>

#include <gsExpressions/gsNullExpr.h>
#include <gsExpressions/gsGeometryMap.h>

#include <gsExpressions/symbol_expr.h>


// A
#include <gsExpressions/abs_expr.h>
#include <gsExpressions/add_expr.h>
#include <gsExpressions/adjugate_expr.h>
#include <gsExpressions/asdiag_expr.h>
// B
// C
#include <gsExpressions/col_expr.h>
#include <gsExpressions/colsum_expr.h>
#include <gsExpressions/constMat_expr.h>
#include <gsExpressions/colBlocks_expr.h>
#include <gsExpressions/curl_expr.h>
// D
#include <gsExpressions/diag_expr.h>
// #include <gsExpressions/div_expr.h> //todo
#include <gsExpressions/divide_expr.h>
// E
#include <gsExpressions/exp_expr.h>
// F
#include <gsExpressions/flat_expr.h>
#include <gsExpressions/frprod_expr.h>
// G
#include <gsExpressions/grad_expr.h>
// H
// I
#include <gsExpressions/idMat_expr.h>
#include <gsExpressions/integral_expr.h>
// J
// K
// L
#include <gsExpressions/lapl_expr.h>
// M
#include <gsExpressions/max_expr.h>
#include <gsExpressions/meas_expr.h>
#include <gsExpressions/mult_expr.h>
// N
#include <gsExpressions/nabla_expr.h>
#include <gsExpressions/nabla2_expr.h>
#include <gsExpressions/normal_expr.h>
// O
#include <gsExpressions/onormal_expr.h>
// P
#include <gsExpressions/pow_expr.h>
// Q
// R
#include <gsExpressions/reshape_expr.h>
#include <gsExpressions/replicate_expr.h>
#include <gsExpressions/rowsum_expr.h>
// S
#include <gsExpressions/sign_expr.h>
#include <gsExpressions/sub_expr.h>
// T
#include <gsExpressions/tangent_expr.h>
#include <gsExpressions/temp_expr.h>
#include <gsExpressions/trace_expr.h>
#include <gsExpressions/transpose_expr.h>
// U
// V
#include <gsExpressions/value_expr.h>
#include <gsExpressions/voigt_expr.h> // @hverhelst todo: add and replace flat_expr
// W
// X
// Y
// Z







namespace gismo
{
namespace expr
{

/*
  Expression for a constant value
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

    inline Scalar eval(const index_t ) const { return _c; }

    inline _expr val() const { return *this; }
    index_t rows() const { return 0; }
    index_t cols() const { return 0; }
    void parse(gsExprHelper<Scalar> &) const { }

    const gsFeSpace<T> & rowVar() const { return gsNullExpr<T>::get(); }
    const gsFeSpace<T> & colVar() const { return gsNullExpr<T>::get(); }

    void print(std::ostream &os) const { os<<_c; }
};


#define GISMO_EXPR_VECTOR_EXPRESSION(name, mname, isSv)                 \
    template<class E> class name##_##expr  : public _expr<name##_##expr<E> > { \
        typename E::Nested_t _u;                                        \
    public:                                                             \
    typedef typename E::Scalar Scalar;                                  \
    enum {Space= E::Space, ScalarValued= isSv, ColBlocks= E::ColBlocks}; \
    name##_##expr(_expr<E> const& u) : _u(u) { }                        \
    mutable Temporary_t tmp;                                            \
    const Temporary_t & eval(const index_t k) const {                   \
        tmp = _u.eval(k).mname(); return tmp; }                         \
    index_t rows() const { return isSv ? 0 : _u.rows(); }               \
    index_t cols() const { return isSv ? 0 : _u.cols(); }               \
    void parse(gsExprHelper<Scalar> & evList) const { _u.parse(evList); } \
    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();} \
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();} \
    void print(std::ostream &os) const                                  \
        { os << #name <<"("; _u.print(os); os <<")"; }                  \
    };

/// Eucledian Norm
GISMO_EXPR_VECTOR_EXPRESSION(norm,norm,1);
/// Squared Eucledian Norm
GISMO_EXPR_VECTOR_EXPRESSION(sqNorm,squaredNorm,1);
/// Normalization of a vector to unit measure
GISMO_EXPR_VECTOR_EXPRESSION(normalized,normalized,0);
/// Inverse of a matrix expression
GISMO_EXPR_VECTOR_EXPRESSION(inv,cramerInverse,0);
// GISMO_EXPR_VECTOR_EXPRESSION(cwSqr,array().square,0)
// GISMO_EXPR_VECTOR_EXPRESSION(sum,array().sum,1)
// GISMO_EXPR_VECTOR_EXPRESSION(sqrt,array().sqrt,0)
//GISMO_EXPR_VECTOR_EXPRESSION(abs,array().abs,0)

//Determinant
GISMO_EXPR_VECTOR_EXPRESSION(det,determinant,1);

//GISMO_EXPR_VECTOR_EXPRESSION(replicate,replicate,0);

#undef GISMO_EXPR_VECTOR_EXPRESSION

/**
   Expression for the component-wise positive part
*/
template<class E>
class ppart_expr : public _expr<ppart_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;
    enum {ScalarValued = E::ScalarValued, Space = E::Space, ColBlocks= E::ColBlocks};
private:
    typename E::Nested_t _u;
    mutable gsMatrix<Scalar> res;
public:

    ppart_expr(_expr<E> const& u) : _u(u) { }

    const gsMatrix<Scalar> & eval(index_t k) const
    {
        res = _u.eval(k).cwiseMax(0.0); // component-wise maximum with zero
        return res;
    }


    index_t rows() const { return _u.rows(); }
    index_t cols() const { return _u.cols(); }

    void parse(gsExprHelper<Scalar> & el) const
    { _u.parse(el); }

    const gsFeSpace<Scalar> & rowVar() const {return _u.rowVar();}
    const gsFeSpace<Scalar> & colVar() const {return _u.colVar();}

    void print(std::ostream &os) const { os<<"posPart("; _u.print(os); os <<")"; }
};

/**
   Expression for the positive part of a given expression
*/
template<class E>
class ppartval_expr : public _expr<ppartval_expr<E> >
{
  typename E::Nested_t _u;
 public:
  typedef typename E::Scalar Scalar;
  enum {ScalarValued = 1, Space = 0, ColBlocks= 0};
  mutable Scalar res;
 public:

  ppartval_expr(_expr<E> const& u) : _u(u) { }

  Scalar & eval(index_t k) const
  {
    res = std::max(0.0,_u.eval(k));
    return res; // component-wise maximum with zero
  }

  index_t rows() const { return 0; }
  index_t cols() const { return 0; }

  void parse(gsExprHelper<Scalar> & evList) const
  { _u.parse(evList); }

  const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
  const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

  void print(std::ostream &os) const { os<<"posPart("; _u.print(os); os <<")"; }
};



/**
   computes outer products of a matrix by a space of dimension > 1
   [Jg Jg Jg] * Jb ..
   (d x d^2)  * (d^2 x N*d)  --> (d x N*d)
*/
template <typename E1, typename E2>
class matrix_by_space_expr  : public _expr<matrix_by_space_expr<E1,E2> >
{
public:
    typedef typename E1::Scalar Scalar;
    enum {ScalarValued = 0, ColBlocks = 1};
    enum {Space = E2::Space};

private:
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;
    mutable gsMatrix<Scalar> res;

public:
    matrix_by_space_expr(E1 const& u, E2 const& v) : _u(u), _v(v) { }


    // choose if ColBlocks
    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        const index_t r   = _u.rows();
        const index_t N  = _v.cols() / (r*r);

        const auto uEv  = _u.eval(k);
        const auto vEv  = _v.eval(k);

        res.resize(r, N*r*r);
        // gsDebugVar(res.cols());
        for (index_t s = 0; s!=r; ++s)
            for (index_t i = 0; i!=N; ++i)
            {
                res.middleCols((s*N + i)*r,r).noalias() =
                    uEv.col(s) * vEv.middleCols((s*N + i)*r,r).row(s);
                //uEv*vEv.middleCols((s*N + i)*r,r);
            }
        //meaning: [Jg Jg Jg] * Jb ..
        return res;
    }

    index_t rows() const { return _u.cols(); }
    index_t cols() const { return _v.cols(); }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); _v.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _v.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _v.colVar(); }

    void print(std::ostream &os) const { os << "matrix_by_space("; _u.print(os); os<<")"; }
};

/**
   computes outer products of a matrix by a space of dimension > 1
   [Jg Jg Jg] * Jb ..
   (d x d^2)  * (d^2 x N*d)  --> (d x N*d)
*/
template <typename E1, typename E2>
class matrix_by_space_expr_tr  : public _expr<matrix_by_space_expr_tr<E1,E2> >
{
public:
    typedef typename E1::Scalar Scalar;
    enum {ScalarValued = 0, ColBlocks = 1};
    enum {Space = E2::Space};

private:
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;
    mutable gsMatrix<Scalar> res;

public:
    matrix_by_space_expr_tr(E1 const& u, E2 const& v) : _u(u), _v(v) { }


    // choose if ColBlocks
    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        const index_t r  = _u.rows();
        const index_t N  = _v.cols() / (r*r);

        const auto uEv  = _u.eval(k);
        const auto vEv  = _v.eval(k);

        res.resize(r, N*r*r);
        for (index_t s = 0; s!=r; ++s)
            for (index_t i = 0; i!=N; ++i)
            {
                res.middleCols((s*N + i)*r,r).noalias() =
                    uEv.transpose()*vEv.middleCols((s*N + i)*r,r).transpose();
            }
        //meaning: [Jg Jg Jg] * Jb ..
        return res;
    }

    index_t rows() const { return _u.cols(); }
    index_t cols() const { return _v.cols(); }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); _v.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _v.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _v.colVar(); }

    void print(std::ostream &os) const { os << "matrix_by_space_tr("; _u.print(os); os<<")"; }
};


/*
  Expression for the derivative of the jacobian of a spline geometry map,
  with respect to the coordinate c.

  It returns a matrix with the gradient of u in row d.
*/
template<class E>
class dJacdc_expr : public _expr<dJacdc_expr<E> >
{
    typename E::Nested_t _u;
public:
    enum{ Space = E::Space, ScalarValued = 0, ColBlocks = (1==E::Space?1:0)};

    typedef typename E::Scalar Scalar;

    mutable gsMatrix<Scalar> res;
    index_t _c;

    dJacdc_expr(const E & u, index_t c) : _u(u), _c(c)
    { GISMO_ASSERT(1==u.dim(),"grad(.) requires 1D variable, use jac(.) instead.");}

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        index_t dd = _u.source().domainDim();
        index_t n = _u.rows();
        res.setZero(dd, dd*n);

        gsMatrix<Scalar> grad = _u.data().values[1].reshapeCol(k, dd, n);
        for(index_t i = 0; i < n; i++){
            res.row(_c).segment(i*dd,dd) = grad.col(i);
        }
        return res;
    }

    index_t rows() const { return _u.source().domainDim(); }

    index_t cols() const { return _u.source().domainDim()*_u.rows(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_u);
        _u.data().flags |= NEED_GRAD;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const
    {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "dJacdc("; _u.print(os); os <<")"; }
};






/// The nabla2 (\f$\nabla^2\f$) of a finite element variable
template<class T>
nabla2_expr<T> nabla2(const gsFeVariable<T> & u) { return nabla2_expr<T>(u); }
// #define lapl(x) nabla2(x).sum() // assume tarDim==1





/*
  Expression for the (precomputed) second fundamental form of a surface
*/
template<class T>
class fform2nd_expr  : public _expr<fform2nd_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;
public:
    typedef T Scalar;
    enum {Space = 0, ScalarValued = 0, ColBlocks = 0};

    fform2nd_expr(const gsGeometryMap<T> & G) : _G(G) { }

    const gsAsConstMatrix<Scalar> eval(const index_t k) const
    {
        return gsAsConstMatrix<Scalar>(_G.data().fundForms.col(k).data(),rows(),cols());
    }

    index_t rows() const { return _G.source().domainDim() ; }
    index_t cols() const { return _G.source().domainDim() ; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_G);
        _G.data().flags |= NEED_2ND_FFORM;
    }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<T>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<T>::get();}

    void print(std::ostream &os) const { os << "fform2nd("; _G.print(os); os <<")"; }
};

/*
  Expression for the (precomputed) inverse of the Jacobian matrix of
  a geometry map
*/
template<class T>
class jacInv_expr  : public _expr<jacInv_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;
public:
    typedef T Scalar;
    enum {Space = 0, ScalarValued = 0, ColBlocks = 0};

    jacInv_expr(const gsGeometryMap<T> & G) : _G(G)
    {
        // Note: for non-square Jacobian matrices, generalized inverse, i.e.: (J^t J)^{-t} J^t
        //GISMO_ASSERT(rows() == cols(), "The Jacobian matrix is not square");
    }

    MatExprType eval(const index_t k) const { return _G.data().jacInvTr.reshapeCol(k,cols(),rows()).transpose(); }

    index_t rows() const { return _G.source().domainDim(); }
    index_t cols() const { return _G.source().targetDim(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_G);
        _G.data().flags |= NEED_GRAD_TRANSFORM;
    }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<T>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<T>::get();}

    // todo mat_expr ?
    // tr() const --> _G.data().fundForm(k)

    void print(std::ostream &os) const { os << "jacInv("; _G.print(os); os <<")"; }
};

/*
  Expression for the Jacobian matrix of a FE variable
*/
template<class E>
class jac_expr : public _expr<jac_expr<E> >
{
    typename E::Nested_t _u;
public:
    enum {ColBlocks = (1==E::Space?1:0) };
    enum {Space = E::Space, ScalarValued= 0 };

    typedef typename E::Scalar Scalar;

    mutable gsMatrix<Scalar> res;

    jac_expr(const E & _u_)
    : _u(_u_) { }

    MatExprType eval(const index_t k) const
    {
        if (0!=Space)
        {
            // Dim x (numActive*Dim)
            res = _u.data().values[1].col(k).transpose().blockDiag(_u.dim());
        }
        else
        {
            res = _u.data().values[1]
                .reshapeCol(k, _u.parDim(), _u.targetDim()).transpose()
                .blockDiag(_u.dim());
        }
        return res;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    //const gsFeSpace<Scalar> & rowVar() const { return rowVar_impl<E>(); }
    const gsFeSpace<Scalar> & colVar() const { return gsNullExpr<Scalar>::get(); }

    index_t rows() const { return rows_impl(_u); }
    index_t cols() const { return cols_impl(_u); }

    // index_t rows() const { return _u.dim(); }
    // index_t cols() const { return _u.source().domainDim(); }

    index_t cardinality_impl() const
    {
        return _u.dim() * _u.data().actives.rows();
    }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_u);
        _u.data().flags |= NEED_DERIV|NEED_ACTIVE;
        //note: cardinality() depends on actives
    }

    void print(std::ostream &os) const { os << "\u2207("; _u.print(os);os <<")"; }

private:

    // The jacobian is different for gsFeVariable, gsFeSolution and gsFeSpace
    // gsFeSolution: Does not work
    // gsFeVariable: dim()=1 and source().targetDim()=d
    // gsFeSpace: dim()=d and source().targetDim()=1
    template<class U> inline
    typename util::enable_if<!(util::is_same<U,gsFeSpace<Scalar> >::value), index_t >::type  // What about solution??
    rows_impl(const U & u)  const
    {
        return u.source().targetDim();
    }

    template<class U> inline
    typename util::enable_if< (util::is_same<U,gsFeSpace<Scalar> >::value), index_t >::type
    rows_impl(const U & u) const
    {
        return u.dim();
    }

    template<class U> inline
    typename util::enable_if<!(util::is_same<U,gsFeSpace<Scalar> >::value), index_t >::type
    cols_impl(const U & u)  const
    {
        return u.source().domainDim();
    }

    template<class U> inline
    typename util::enable_if< (util::is_same<U,gsFeSpace<Scalar> >::value), index_t >::type
    cols_impl(const U & u) const
    {
        return u.source().domainDim();
    }

    // The jacobian is different for gsFeVariable, gsFeSolution and gsFeSpace
    // gsFeSolution: Does not work
    // gsFeVariable: rowVar = NULL
    // gsFeSpace: rowVar = u
    template<class U> inline
    typename util::enable_if<!(util::is_same<U,gsFeSpace<Scalar> >::value), const gsFeSpace<Scalar> &  >::type
    rowVar_impl() const
    {
        return gsNullExpr<Scalar>::get();
    }

    template<class U> inline
    typename util::enable_if<(util::is_same<U,gsFeSpace<Scalar> >::value), const gsFeSpace<Scalar> &  >::type
    rowVar_impl() const
    {
        return _u;
    }
};

/*
  Expression for the Jacobian matrix of a geometry map
*/
template<class T>
class jac_expr<gsGeometryMap<T> > : public _expr<jac_expr<gsGeometryMap<T> > >
{
    typename gsGeometryMap<T>::Nested_t _G;

public:
    typedef T Scalar;
    enum {Space = 0, ScalarValued= 0, ColBlocks= 0};

    jac_expr(const gsGeometryMap<T> & G) : _G(G) { }
    MatExprType eval(const index_t k) const
    {
        // TarDim x ParDim
        return _G.data().values[1]
            .reshapeCol(k, _G.data().dim.first, _G.data().dim.second).transpose();
    }

    index_t rows() const { return _G.source().targetDim(); }

    index_t cols() const { return _G.source().domainDim(); }

    static const gsFeSpace<Scalar> & rowVar() { return gsNullExpr<Scalar>::get(); }
    static const gsFeSpace<Scalar> & colVar() { return gsNullExpr<Scalar>::get(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_G);
        _G.data().flags |= NEED_DERIV;
    }

    meas_expr<T> absDet() const
    {
        GISMO_ASSERT(rows() == cols(), "The Jacobian matrix is not square");
        return meas_expr<T>(_G);
    }

    jacInv_expr<T> inv() const
    {
        GISMO_ASSERT(rows() == cols(), "The Jacobian matrix is not square");
        return jacInv_expr<T>(_G);
    }

    /// The generalized Jacobian matrix inverse, i.e.: (J^t J)^{-t} J^t
    jacInv_expr<T> ginv() const { return jacInv_expr<T>(_G); }

    void print(std::ostream &os) const { os << "\u2207("; _G.print(os); os <<")"; }
};

template<class E>
class hess_expr : public _expr<hess_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;
private:
    typename E::Nested_t _u;
    mutable gsMatrix<Scalar> res;
public:
    enum {ScalarValued = 0, ColBlocks = (1==E::Space?1:0) };
    enum {Space = E::Space };

public:
    hess_expr(const E & u) : _u(u)
    {
        //gsInfo << "\n-expression is space ? "<<E::Space <<"\n"; _u.print(gsInfo);
        //GISMO_ASSERT(1==_u.dim(),"hess(.) requires 1D variable");
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        const gsFuncData<Scalar> & dd = _u.data();
        const index_t sz = cardinality_impl();
        res.resize(dd.dim.first, sz*dd.dim.first);
        secDerToHessian(dd.values[2].col(k), dd.dim.first, res);
        res.resize(dd.dim.first, res.cols()*dd.dim.first);
        // Note: auto returns by value here,
        // in C++11 we may add in -> decltype(res) &
        return res;
    }

    index_t rows() const { return _u.data().dim.first; }
    index_t cols() const
    {   return rows();
        //return 2*_u.data().values[2].rows() / (1+_u.data().dim.first);
    }

    index_t cardinality_impl() const
    {
        return 2*_u.data().values[2].rows()/
            (_u.data().dim.first*(1+_u.data().dim.first));
        //gsDebugVar(_u.data().values.front().rows());//empty!
    }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_u);
        _u.data().flags |= NEED_2ND_DER;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return gsNullExpr<Scalar>::get(); }

    void print(std::ostream &os) const
    //    { os << "hess("; _u.print(os);os <<")"; }
    { os << "\u210D(U)"; }
};

template<class T>
class hess_expr<gsFeSolution<T> > : public _expr<hess_expr<gsFeSolution<T> > >
{
protected:
    const gsFeSolution<T> _u;

public:
    typedef T Scalar;
    enum{Space = 0, ScalarValued = 0, ColBlocks = 0 };

    hess_expr(const gsFeSolution<T> & u) : _u(u) { }

    mutable gsMatrix<T> deriv2, res;
    const gsMatrix<T> & eval(const index_t k) const
    {
        GISMO_ASSERT(_u.check(), "Invalid state in gsFeSolution");
        const gsDofMapper & map = _u.mapper();
        const index_t numActs = _u.data().values[0].rows();
        const index_t pdim = _u.parDim();
        index_t numDers = pdim*(pdim+1)/2;
        auto & act = _u.data().actives.col(1 == _u.data().actives.cols() ? 0:k );

        // In the scalar case, the hessian is returned as a pdim x pdim matrix
        if (1==_u.dim())
        {
            res.setZero(numDers,1);
            for (index_t i = 0; i!=numActs; ++i)
            {
                const index_t ii = map.index(act[i], _u.data().patchId, 0);
                deriv2 = _u.data().values[2].block(i*numDers,k,numDers,1);
                if ( map.is_free_index(ii) ) // DoF value is in the solVector
                    res += _u.coefs().at(ii) * deriv2;
                else
                    res +=_u.fixedPart().at( map.global_to_bindex(ii) ) * deriv2;
            }
            secDerToHessian(res, pdim, deriv2);
            res.swap(deriv2);
            res.resize(pdim,pdim);
        }
        // In the vector case, the hessian is returned as a matrix where each row corresponds to the component of the solution and contains the derivatives in the columns
        else
        {
            res.setZero(rows(), numDers);
            for (index_t c = 0; c != _u.dim(); c++)
                for (index_t i = 0; i != numActs; ++i)
                {
                    const index_t ii = map.index(act[i], _u.data().patchId, c);
                    deriv2 = _u.space().data().values[2].block(i * numDers, k, numDers,
                                                               1).transpose(); // start row, start col, rows, cols
                    if (map.is_free_index(ii)) // DoF value is in the solVector
                        res.row(c) += _u.coefs().at(ii) * deriv2;
                    else
                        res.row(c) += _u.fixedPart().at(map.global_to_bindex(ii)) * deriv2;
                }
        }
        return res;
    }

    index_t rows() const
    {
        if (1==_u.dim())
            return _u.parDim();
        else
            return _u.dim(); //  number of components
    }
    index_t cols() const
    {
        if (1==_u.dim())
            return _u.parDim();
        // second derivatives in the columns; i.e. [d11, d22, d33, d12, d13, d23]
        else
            return _u.parDim() * (_u.parDim() + 1) / 2;
    }

    const gsFeSpace<Scalar> & rowVar() const { return gsNullExpr<Scalar>::get(); }
    const gsFeSpace<Scalar> & colVar() const { return gsNullExpr<Scalar>::get(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        _u.parse(evList);                         // add symbol
        evList.add(_u.space());
        _u.data().flags |= NEED_ACTIVE | NEED_VALUE | NEED_DERIV2;
    }

    void print(std::ostream &os) const { os << "\u210D(s)"; }
};


/*
  Expression for the partial derivative (matrices) of the Jacobian
  matrix of the geometry map
*/
template<class T>
class dJacG_expr : public _expr<dJacG_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;

    mutable gsMatrix<T> res;
public:
    typedef T Scalar;

    dJacG_expr(const gsGeometryMap<T> & G) : _G(G) { }

    MatExprType eval(const index_t k) const
    {
        const index_t sz = _G.data().values[0].rows();
        const index_t s = _G.data().derivSize(); //dim.first*(_G.data().dim.first+1)/2;
        (void)s;
        res.resize(_G.data().dim.second, sz*_G.data().dim.first);
        res.setOnes();//_G.data().values[2].segment(i*k,k); // todo
        return res;
    }

    index_t rows() const { return _G.source().targetDim(); }
    index_t cols() const { return _G.source().domainDim(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_G);
        _G.data().flags |= NEED_2ND_DER;
    }
};

template <typename E1, typename E2>
class collapse_expr : public _expr<collapse_expr<E1, E2> >
{
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;

public:
    enum {ScalarValued = 0, ColBlocks = 0};
    enum { Space = (int)E1::Space + (int)E2::Space };

    typedef typename E1::Scalar Scalar;

    mutable gsMatrix<Scalar> res;

    collapse_expr(_expr<E1> const& u,
                  _expr<E2> const& v)
    : _u(u), _v(v) { }

    //EIGEN_STRONG_INLINE MatExprType
    const gsMatrix<Scalar> &
    eval(const index_t k) const
    {
        const index_t nb = rows();
        const auto tmpA = _u.eval(k);
        const auto tmpB = _v.eval(k);

        if (E1::ColBlocks)
        {
            const index_t ur = _v.rows();
            res.resize(nb, ur);
            for (index_t i = 0; i!=nb; ++i)
            {
                res.row(i).transpose().noalias() = tmpA.middleCols(i*ur,ur) * tmpB;
            }
        }
        else if (E2::ColBlocks)
        {
            const index_t ur = _u.cols();
            res.resize(nb, ur);
            for (index_t i = 0; i!=nb; ++i)
            {
                res.row(i).noalias() = tmpA * tmpB.middleCols(i*ur,ur);
            }
        }

        return res;
    }

    index_t rows() const { return E1::ColBlocks ? _u.cols() / _v.rows() : _v.cols() / _u.cols() ; }
    index_t cols() const { return E1::ColBlocks ? _v.rows()  : _u.cols(); }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); _v.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const
    { return E1::ColBlocks ? _u.rowVar() : _v.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const
    {
        GISMO_ERROR("none");
    }

    void print(std::ostream &os) const { _u.print(os); os<<"~"; _v.print(os); }
};

// Multi-matrix collapsed by a vector
template <typename E1, typename E2> //EIGEN_STRONG_INLINE
//collapse_expr<E1,E2> const  operator&(<E1> const& u, _expr<E2> const& v)
collapse_expr<E1,E2> collapse( _expr<E1> const& u, _expr<E2> const& v)
{ return collapse_expr<E1, E2>(u, v); }



/*
  lincom_expr (lc) ?
  Expression for (square) matrix summation operation

  M [r x r*k] is a list of matrices
  Summation is done over k,
  [M1 M2 .. Mk]

  u [s x k] is a list of vectors

  Computed quantity is of size [r x r*s] and contains
  [ ... sum_k(Mk * u(s,k) ) ... ]_s
*/
template <typename E1, typename E2>
class summ_expr : public _expr<summ_expr<E1,E2> >
{
public:
    typedef typename E1::Scalar Scalar;

    enum {Space = E1::Space, ScalarValued= 0, ColBlocks= E2::ColBlocks};

    summ_expr(E1 const& u, E2 const& M) : _u(u), _M(M) { }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        auto sl   = _u.eval(k);
        const index_t sr = sl.rows();
        auto ml   = _M.eval(k);
        const index_t mr = ml.rows();
        const index_t mb = _M.cardinality();

        GISMO_ASSERT(_M.cols()==_M.rows(),"Matrix must be square: "<< _M.rows()<<" x "<< _M.cols() << " expr: "<< _M );
        GISMO_ASSERT(mb==_u.cols(),"cardinality must match vector, but card(M)="<<_M.cardinality()<<" and cols(u)="<<_u.cols());

        res.setZero(mr, sr * mr);
        for (index_t i = 0; i!=sr; ++i)
            for (index_t t = 0; t!=mb; ++t) // lc
                res.middleCols(i*mr,mr) += sl(i,t) * ml.middleCols(t*mr,mr);
        return res;
    }

    index_t rows() const { return _M.rows(); }
    index_t cols() const { return _M.rows(); } //_u.rows()

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); _M.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return gsNullExpr<Scalar>::get(); }

    index_t cardinality_impl() const
    { GISMO_ERROR("Something went terribly wrong"); }

    void print(std::ostream &os) const
    { os << "sum("; _M.print(os); os<<","; _u.print(os); os<<")"; }

private:
    typename E1::Nested_t _u;
    typename E2::Nested_t _M;

    mutable gsMatrix<Scalar> res;
};

/*
  Expression for symmetrization operation
*/
template <typename E>
class symm_expr : public _expr<symm_expr<E> >
{
    typename E::Nested_t _u;

    mutable gsMatrix<typename E::Scalar> tmp;
public:
    typedef typename E::Scalar Scalar;

    enum { Space = (0==E::Space ? 0 : E::Space), ScalarValued= E::ScalarValued, ColBlocks= E::ColBlocks };

    symm_expr(_expr<E> const& u)
    : _u(u) { }

    MatExprType eval(const index_t k) const
    {
        //const MatExprType tmp = _u.eval(k);
        tmp = _u.eval(k);
        // todo: avoid temporary or not ?
        return tmp * tmp.transpose();
    }

    index_t rows() const { return _u.rows(); }
    index_t cols() const { return _u.rows(); }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.rowVar(); }

    void print(std::ostream &os) const { os << "symm("; _u.print(os); os <<")"; }
};

template <typename E>
class symmetrize_expr : public _expr<symmetrize_expr<E> >
{
    typename E::Nested_t _u;

    mutable gsMatrix<typename E::Scalar> tmp;
public:
    enum { Space = (0==E::Space ? 0 : 3), ScalarValued=E::ScalarValued, ColBlocks= E::ColBlocks };
    typedef typename E::Scalar Scalar;

    symmetrize_expr(_expr<E> const& u)
    : _u(u) { }

    MatExprType eval(const index_t k) const
    {
        //const MatExprType tmp = _u.eval(k);
        tmp = _u.eval(k);
        // todo: avoid temporary or not ?
        return tmp + tmp.transpose();
    }

    index_t rows() const { return _u.rows(); }
    index_t cols() const { return _u.rows(); }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.rowVar(); }
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "symmetrize("; _u.print(os); os <<")"; }
};

/* Symmetrization operation
   template <typename E> symm_expr<E> const
   symm(_expr<E> const& u) { return symm_expr<E>(u);}
*/

#undef MatExprType
#undef AutoReturn_t
//----------------------------------------------------------------------------------

/// The identity matrix of dimension \a dim
EIGEN_STRONG_INLINE idMat_expr id(const index_t dim) { return idMat_expr(dim); }

EIGEN_STRONG_INLINE constMat_expr ones(const index_t dim)
{
    gsMatrix<real_t> ones(dim, dim);
    ones.fill(1);
    return constMat_expr(ones);
}

EIGEN_STRONG_INLINE constMat_expr mat(const gsMatrix<real_t> mat) { return constMat_expr(mat); }

// Returns the unit as an expression
//EIGEN_STRONG_INLINE _expr<real_t> one() { return _expr<real_t,true>(1); }

/// Replicate an expression
template <typename E> EIGEN_STRONG_INLINE
replicate_expr<E> const replicate(E const & u, index_t n, index_t m = 1)
{ return replicate_expr<E>(u, n, m); }

/// Reshape an expression
template <class E> EIGEN_STRONG_INLINE
reshape_expr<E> const reshape(E const & u, index_t n, index_t m)
{ return reshape_expr<E>(u, n, m); }

/// Make a matrix 2x2 expression "flat"
template <typename E> EIGEN_STRONG_INLINE
flat_expr<E> const flat(E const & u)
{ return flat_expr<E>(u); }

/// Get diagonal elements of matrix as a vector
template <typename E> EIGEN_STRONG_INLINE
diag_expr<E> const diagonal(E const & u)
{ return diag_expr<E>(u); }

/// Absolute value
template<class E> EIGEN_STRONG_INLINE
abs_expr<E> abs(const E & u) { return abs_expr<E>(u); }

/// The gradient of a variable
template<class E> EIGEN_STRONG_INLINE
grad_expr<E> grad(const E & u) { return grad_expr<E>(u); }

/// The derivative of the jacobian of a geometry map with respect to a coordinate.
template<class E> EIGEN_STRONG_INLINE
dJacdc_expr<E> dJacdc(const E & u, index_t c) { return dJacdc_expr<E>(u,c); }

/// The curl of a finite element variable
template<class T> EIGEN_STRONG_INLINE
curl_expr<T> curl(const gsFeVariable<T> & u) { return curl_expr<T>(u); }

/// The nabla (\f$\nabla\f$) of a finite element variable
template<class T> EIGEN_STRONG_INLINE
nabla_expr<T> nabla(const gsFeVariable<T> & u) { return nabla_expr<T>(u); }

/// The (outer pointing) boundary normal of a geometry map
template<class T> EIGEN_STRONG_INLINE
onormal_expr<T> nv(const gsGeometryMap<T> & u) { return onormal_expr<T>(u); }

template<class T> EIGEN_STRONG_INLINE
normal_expr<T> sn(const gsGeometryMap<T> & u) { return normal_expr<T>(u); }

/// The tangent boundary vector of a geometry map in 2D
template<class T> EIGEN_STRONG_INLINE
tangent_expr<T> tv(const gsGeometryMap<T> & u) { return tangent_expr<T>(u); }

template<class E> EIGEN_STRONG_INLINE
lapl_expr<E> lapl(const symbol_expr<E> & u) { return lapl_expr<E>(u); }

template<class T> EIGEN_STRONG_INLINE
lapl_expr<gsFeSolution<T> > lapl(const gsFeSolution<T> & u)
{ return lapl_expr<gsFeSolution<T> >(u); }

/// The second fundamental form of \a G
template<class T> EIGEN_STRONG_INLINE fform2nd_expr<T> fform2nd(const gsGeometryMap<T> & G)
{ return fform2nd_expr<T>(G); }

/// The Jacobian matrix of a FE variable
template<class E> EIGEN_STRONG_INLINE
jac_expr<E> jac(const symbol_expr<E> & u) { return jac_expr<E>(u); }

/// The Jacobian matrix of a geometry map
template<class T> EIGEN_STRONG_INLINE
jac_expr<gsGeometryMap<T> > jac(const gsGeometryMap<T> & G) {return jac_expr<gsGeometryMap<T> >(G);}

/// Jacobian matrix for a solution expression
template<class T> EIGEN_STRONG_INLINE
grad_expr<gsFeSolution<T> > jac(const gsFeSolution<T> & s) {return grad_expr<gsFeSolution<T> >(s);}

template<class E> EIGEN_STRONG_INLINE
hess_expr<E> hess(const symbol_expr<E> & u) { return hess_expr<E>(u); }

/// The hessian of a geometry map
template<class T> EIGEN_STRONG_INLINE
hess_expr<gsGeometryMap<T> > hess(const gsGeometryMap<T> & u) { return hess_expr<gsGeometryMap<T> >(u); }

/// The hessian of a solution variable
template<class T> EIGEN_STRONG_INLINE
hess_expr<gsFeSolution<T> > hess(const gsFeSolution<T> & u) { return hess_expr<gsFeSolution<T> >(u); }

/// The partial derivatives of the Jacobian matrix of a geometry map
template<class T> EIGEN_STRONG_INLINE
dJacG_expr<T> dJac(const gsGeometryMap<T> & G) { return dJacG_expr<T>(G); }

/// The measure of a geometry map
template<class T> EIGEN_STRONG_INLINE
meas_expr<T> meas(const gsGeometryMap<T> & G) { return meas_expr<T>(G); }

/// Multiplication operator for expressions
template <typename E1, typename E2> EIGEN_STRONG_INLINE
mult_expr<E1,E2> const operator*(_expr<E1> const& u, _expr<E2> const& v)
{ return mult_expr<E1, E2>(u, v); }

template <typename E2> EIGEN_STRONG_INLINE
mult_expr<typename E2::Scalar,E2,false> const
operator*(typename E2::Scalar const& u, _expr<E2> const& v)
{ return mult_expr<typename E2::Scalar, E2, false>(u, v); }

template <typename E1> EIGEN_STRONG_INLINE
mult_expr<typename E1::Scalar,E1,false> const
operator*(_expr<E1> const& v, typename E1::Scalar const& u)
{ return mult_expr<typename E1::Scalar,E1, false>(u, v); }

template <typename E1> EIGEN_STRONG_INLINE
mult_expr<typename E1::Scalar,E1,false> const
operator-(_expr<E1> const& u)
{ return mult_expr<typename E1::Scalar,E1, false>(-1, u); }

template <typename E> mult_expr<constMat_expr, E> const
operator*( gsMatrix<typename E::Scalar> const& u, _expr<E> const& v)
{ return mult_expr<constMat_expr, E>(mat(u), v); }

template <typename E> mult_expr<E, constMat_expr> const
operator*(_expr<E> const& u, gsMatrix<typename E::Scalar> const& v)
{ return mult_expr<E, constMat_expr>(u, mat(v) ); }

/// Frobenious product (also known as double dot product) operator for expressions
template <typename E1, typename E2> EIGEN_STRONG_INLINE
frprod_expr<E1,E2> const  operator%(_expr<E1> const& u, _expr<E2> const& v)
{ return frprod_expr<E1, E2>(u, v); }

/// Scalar division operator for expressions
template <typename E1, typename E2> EIGEN_STRONG_INLINE
divide_expr<E1,E2> const operator/(_expr<E1> const& u, _expr<E2> const& v)
{ return divide_expr<E1,E2>(u, v); }

template <typename E> EIGEN_STRONG_INLINE
divide_expr<E,typename E::Scalar> const
operator/(_expr<E> const& u, const typename E::Scalar v)
{ return divide_expr<E,typename E::Scalar>(u, v); }

template <typename E> EIGEN_STRONG_INLINE
divide_expr<typename E::Scalar,E> const
operator/(const typename E::Scalar u, _expr<E> const& v)
{ return divide_expr<typename E::Scalar,E>(u, v); }

/// Addition operator for expressions
template <typename E1, typename E2> EIGEN_STRONG_INLINE
add_expr<E1,E2> const operator+(_expr<E1> const& u, _expr<E2> const& v)
{ return add_expr<E1, E2>(u, v); }

/// Addition operator for expressions and numbers
template <typename E> EIGEN_STRONG_INLINE
add_expr< E, _expr<typename E::Scalar, true> >
operator+(_expr<E> const& u, const typename E::Scalar v)
{ return add_expr<E,_expr<typename E::Scalar>>(u, _expr<typename E::Scalar,true>(v)); }

/// Addition operator for expressions and numbers
template <typename E> EIGEN_STRONG_INLINE
add_expr< E, _expr<typename E::Scalar, true> >
operator+(const typename E::Scalar v, _expr<E> const& u)
{ return add_expr<E,_expr<typename E::Scalar>>(u, _expr<typename E::Scalar,true>(v)); }

/// Matrix-summation operator for expressions
template <typename E1, typename E2> EIGEN_STRONG_INLINE
summ_expr<E1,E2> const summ(E1 const & u, E2 const& M)
{ return summ_expr<E1,E2>(u, M); }

/// Matrix by space TODO: find better name and/or description? And is this the best place?
/// [Jg Jg Jg] * Jb ..
template <typename E1, typename E2> EIGEN_STRONG_INLINE
matrix_by_space_expr<E1,E2> const matrix_by_space(E1 const & u, E2 const& v)
{ return matrix_by_space_expr<E1,E2>(u, v); }

/// Matrix by space TODO: find better name and/or description? And is this the best place?
/// [Jg Jg Jg] * Jb ..
template <typename E1, typename E2> EIGEN_STRONG_INLINE
matrix_by_space_expr_tr<E1,E2> const matrix_by_space_tr(E1 const & u, E2 const& v)
{ return matrix_by_space_expr_tr<E1,E2>(u, v); }

/// Subtraction operator for expressions
template <typename E1, typename E2> EIGEN_STRONG_INLINE
sub_expr<E1,E2> const operator-(_expr<E1> const& u, _expr<E2> const& v)
{ return sub_expr<E1, E2>(u, v); }

template <typename E2> EIGEN_STRONG_INLINE
sub_expr<_expr<typename E2::Scalar>,E2> const
operator-(typename E2::Scalar const& s, _expr<E2> const& v)
{
    // assert E2::ScalarValued
    return sub_expr<_expr<typename E2::Scalar>, E2>(_expr<typename E2::Scalar>(s), v);
}


// Shortcuts for common quantities, for instance function
// transformations by the geometry map \a G
#define GISMO_SHORTCUT_VAR_EXPRESSION(name,impl) template<class E> EIGEN_STRONG_INLINE \
    auto name(const E & u) -> decltype(impl) { return impl; }
#define GISMO_SHORTCUT_MAP_EXPRESSION(name,impl) template<class T> EIGEN_STRONG_INLINE \
    auto name(const gsGeometryMap<T> & G)  -> decltype(impl) { return impl; }
#define GISMO_SHORTCUT_PHY_EXPRESSION(name,impl) template<class E> EIGEN_STRONG_INLINE \
    auto name(const E & u, const gsGeometryMap<typename E::Scalar> & G)  -> decltype(impl) { return impl; }

// Divergence
GISMO_SHORTCUT_VAR_EXPRESSION(  div, jac(u).trace() )
GISMO_SHORTCUT_PHY_EXPRESSION( idiv, ijac(u,G).trace()    )

// The unit (normalized) boundary (outer pointing) normal
GISMO_SHORTCUT_MAP_EXPRESSION(unv, nv(G).normalized()   )
// The unit (normalized) boundary (surface) normal
GISMO_SHORTCUT_MAP_EXPRESSION(usn, sn(G).normalized()   )

GISMO_SHORTCUT_PHY_EXPRESSION(igrad, grad(u)*jac(G).ginv() ) // transpose() problem ??
GISMO_SHORTCUT_VAR_EXPRESSION(igrad, grad(u) ) // u is presumed to be defined over G

GISMO_SHORTCUT_PHY_EXPRESSION( ijac, jac(u) * jac(G).ginv())

// note and todo: does this work for non-scalar solutions?
GISMO_SHORTCUT_PHY_EXPRESSION(ihess,
                              jac(G).ginv().tr()*( hess(u) - summ(igrad(u,G),hess(G)) ) * jac(G).ginv() )
GISMO_SHORTCUT_VAR_EXPRESSION(ihess, hess(u) )

GISMO_SHORTCUT_PHY_EXPRESSION(ilapl, ihess(u,G).trace()   )
GISMO_SHORTCUT_VAR_EXPRESSION(ilapl, hess(u).trace() )

GISMO_SHORTCUT_VAR_EXPRESSION(fform, jac(u).tr()*jac(u) )
GISMO_SHORTCUT_VAR_EXPRESSION(shapeop, fform(u).inv() * fform2nd(u) )

#undef GISMO_SHORTCUT_PHY_EXPRESSION
#undef GISMO_SHORTCUT_VAR_EXPRESSION
#undef GISMO_SHORTCUT_MAP_EXPRESSION

} // namespace expr

} //namespace gismo
