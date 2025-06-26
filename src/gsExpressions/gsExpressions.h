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

// Forward declaration in gismo namespace
template<class T> class gsExprHelper;

/** @namespace gismo::expr

    @brief This namespace contains expressions used for FE computations

    \ingroup Expressions
*/
namespace expr
{
template <typename E, bool = util::is_arithmetic<E>::value >
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

#include <gsExpressions/gsVector_expr.h>

#include <gsExpressions/symbol_expr.h>

// Other
#include <gsExpressions/_expr.h>
#include <gsExpressions/_expr_macros.h>
// A
#include <gsExpressions/abs_expr.h>
#include <gsExpressions/add_expr.h>
#include <gsExpressions/adjugate_expr.h>
#include <gsExpressions/asdiag_expr.h>
// B
// C
#include <gsExpressions/col_expr.h>
#include <gsExpressions/collapse_expr.h>
#include <gsExpressions/colsum_expr.h>
#include <gsExpressions/constMat_expr.h>
#include <gsExpressions/colBlocks_expr.h>
#include <gsExpressions/curl_expr.h>
// D
#include <gsExpressions/diag_expr.h>
// #include <gsExpressions/div_expr.h> //todo
#include <gsExpressions/divide_expr.h>
#include <gsExpressions/dJacdc_expr.h>
#include <gsExpressions/dJacG_expr.h>
// E
#include <gsExpressions/exp_expr.h>
// F
#include <gsExpressions/fform2nd_expr.h>
#include <gsExpressions/flat_expr.h>
#include <gsExpressions/frprod_expr.h>
// G
#include <gsExpressions/grad_expr.h>
// H
#include <gsExpressions/hess_expr.h>
// I
#include <gsExpressions/idMat_expr.h>
#include <gsExpressions/integral_expr.h>
// J
#include <gsExpressions/jac_expr.h>
#include <gsExpressions/jacInv_expr.h>
// K
// L
#include <gsExpressions/lapl_expr.h>
// M
#include <gsExpressions/matrix_by_space_expr.h>
#include <gsExpressions/matrix_by_space_tr_expr.h>
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
#include <gsExpressions/ppart_expr.h>
#include <gsExpressions/ppartval_expr.h>
#include <gsExpressions/pow_expr.h>
// Q
// R
#include <gsExpressions/reshape_expr.h>
#include <gsExpressions/replicate_expr.h>
#include <gsExpressions/rowsum_expr.h>
// S
#include <gsExpressions/sign_expr.h>
#include <gsExpressions/sub_expr.h>
#include <gsExpressions/summ_expr.h>
#include <gsExpressions/symmetrize_expr.h>
#include <gsExpressions/symm_expr.h>
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

#undef MatExprType
#undef AutoReturn_t

