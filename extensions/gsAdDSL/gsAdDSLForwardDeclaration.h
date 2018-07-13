#pragma once

namespace adDSL {

template< typename Impl >
struct RealExpression;

template< typename ExprType_v > struct Expr_abs_Real;
template< typename ExprType_v > struct Expr_acos_Real;
template< typename ExprType_v > struct Expr_asin_Real;
template< typename ExprType_v1, typename ExprType_v2 > struct Expr_atan2_Real_Real;
template< typename ExprType_v > struct Expr_atan_Real;
template< typename ExprType_v > struct Expr_ceil_Real;
template< typename ExprType_v > struct Expr_cos_Real;
template< typename ExprType_v > struct Expr_cosh_Real;
template< typename ExprType_v > struct Expr_exp_Real;
template< typename ExprType_v > struct Expr_floor_Real;
template< typename ExprType_v > struct Expr_log10_Real;
template< typename ExprType_v > struct Expr_log_Real;
template< typename ExprType_v1, typename ExprType_v2 > struct Expr_max_Real_Real;
template< typename ExprType_v1, typename ExprType_v2 > struct Expr_min_Real_Real;
template< typename ExprType_v1, typename ExprType_v2 > struct Expr_pow_Real_Real;
template< typename ExprType_v > struct Expr_sin_Real;
template< typename ExprType_v > struct Expr_sinh_Real;
template< typename ExprType_v > struct Expr_sqrt_Real;
template< typename ExprType_v > struct Expr_tan_Real;
template< typename ExprType_v > struct Expr_tanh_Real;

template< typename ExprType_v > Expr_abs_Real< ExprType_v >   abs(const RealExpression< ExprType_v >& v);
template< typename ExprType_v > Expr_acos_Real< ExprType_v >  acos(const RealExpression< ExprType_v >& v);
template< typename ExprType_v > Expr_asin_Real< ExprType_v >  asin(const RealExpression< ExprType_v >& v);
template< typename ExprType_v1, typename ExprType_v2 > Expr_atan2_Real_Real< ExprType_v1, ExprType_v2 > atan2(const RealExpression< ExprType_v1 >& v1, const RealExpression< ExprType_v2 >& v2);
template< typename ExprType_v > Expr_atan_Real< ExprType_v >  atan(const RealExpression< ExprType_v >& v);
template< typename ExprType_v > Expr_ceil_Real< ExprType_v >  ceil(const RealExpression< ExprType_v >& v);
template< typename ExprType_v > Expr_cos_Real< ExprType_v >   cos(const RealExpression< ExprType_v >& v);
template< typename ExprType_v > Expr_cosh_Real< ExprType_v >  cosh(const RealExpression< ExprType_v >& v);
template< typename ExprType_v > Expr_exp_Real< ExprType_v >   exp(const RealExpression< ExprType_v >& v);
template< typename ExprType_v > Expr_floor_Real< ExprType_v > floor(const RealExpression< ExprType_v >& v);
template< typename ExprType_v > Expr_log10_Real< ExprType_v > log10(const RealExpression< ExprType_v >& v);
template< typename ExprType_v > Expr_log_Real< ExprType_v >   log(const RealExpression< ExprType_v >& v);
template< typename ExprType_v1, typename ExprType_v2 > Expr_max_Real_Real< ExprType_v1, ExprType_v2 >   max(const RealExpression< ExprType_v1 >& v1, const RealExpression< ExprType_v2 >& v2);
template< typename ExprType_v1, typename ExprType_v2 > Expr_min_Real_Real< ExprType_v1, ExprType_v2 >   min(const RealExpression< ExprType_v1 >& v1, const RealExpression< ExprType_v2 >& v2);
template< typename ExprType_v1, typename ExprType_v2 > Expr_pow_Real_Real< ExprType_v1, ExprType_v2 >   pow(const RealExpression< ExprType_v1 >& v1, const RealExpression< ExprType_v2 >& v2);
template< typename ExprType_v > Expr_cos_Real< ExprType_v >   cos(const RealExpression< ExprType_v >& v);
template< typename ExprType_v > Expr_sin_Real< ExprType_v >   sin(const RealExpression< ExprType_v >& v);
template< typename ExprType_v > Expr_sinh_Real< ExprType_v >  sinh(const RealExpression< ExprType_v >& v);
template< typename ExprType_v > Expr_sqrt_Real< ExprType_v >  sqrt(const RealExpression< ExprType_v >& v);
template< typename ExprType_v > Expr_tan_Real< ExprType_v >   tan(const RealExpression< ExprType_v >& v);
template< typename ExprType_v > Expr_tanh_Real< ExprType_v >  tanh(const RealExpression< ExprType_v >& v);

}
