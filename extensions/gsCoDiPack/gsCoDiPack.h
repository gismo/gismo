/** @file gsCoDiPack.h

    @brief Header for CoDiPack package

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, M. Moeller
*/

#pragma once

// Enable implicit conversion and disable compiler warnings
#define CODI_ImplicitConversion 1
#define CODI_ImplicitConversionWarning 0
#include <codi.hpp>

#include <gsCore/gsTemplateTools.h>

#ifdef GISMO_BUILD_LIB
#ifdef gsCoDiPack_EXPORTS
#undef  EXTERN_CLASS_TEMPLATE
#define EXTERN_CLASS_TEMPLATE CLASS_TEMPLATE_INST
#endif

namespace codi {
template<class T>
using codi_real = std::conditional<ExpressionTraits::IsExpression<T>::value, GISMO_COEFF_TYPE, T>;

// General forward AD type: RealForwardGen<real_t>
EXTERN_CLASS_TEMPLATE
ActiveType<ForwardEvaluation<real_t, real_t> >;

// General vector forward AD type: RealForwardVec<real_t, 4>
EXTERN_CLASS_TEMPLATE
ActiveType<ForwardEvaluation<real_t, Direction<real_t, 4> > >;  
  
// General reverse AD type: RealReverseGen<real_t>
// Jacobian taping approach with linear index handling.
EXTERN_CLASS_TEMPLATE
ActiveType<JacobianLinearTape<JacobianTapeTypes<real_t, real_t,
                                                LinearIndexManager<index_t>,
                                                DefaultChunkedData> > >;

// General vector reverse AD type: RealReverseVec<real_t, 4>
// Jacobian taping approach with linear index handling.
EXTERN_CLASS_TEMPLATE
ActiveType<JacobianLinearTape<JacobianTapeTypes<real_t, Direction<real_t, 4>,
                                                LinearIndexManager<index_t>,
                                                DefaultChunkedData> > >;

// General unchecked reverse AD type: RealReverseUncheckedGen<real_t>
// Requires preallocation of data. See DataManagementTapeInterface.
// Jacobian taping approach with linear index handling.
EXTERN_CLASS_TEMPLATE
ActiveType<JacobianLinearTape<JacobianTapeTypes<real_t, real_t,
                                                LinearIndexManager<index_t>,
                                                DefaultBlockData> > >;

// General reverse AD type: RealReverseIndexGen<real_t>
// Jacobian taping approach with reuse index handling.
EXTERN_CLASS_TEMPLATE
ActiveType<JacobianReuseTape<JacobianTapeTypes<real_t, real_t,
                                               MultiUseIndexManager<index_t>,
                                               DefaultChunkedData> > >;

// General reverse AD type: RealReversePrimalGen<real_t>
// Primal value taping approach with linear index handling.
EXTERN_CLASS_TEMPLATE
ActiveType<PrimalValueLinearTape<
             PrimalValueTapeTypes<real_t, real_t,
                                  LinearIndexManager<index_t>,
                                  InnerStatementEvaluator, DefaultChunkedData> > >;

// General reverse AD type: RealReversePrimalIndexGen<real_t>
// Primal value taping approach with reuse index handling.
EXTERN_CLASS_TEMPLATE
ActiveType<PrimalValueReuseTape<
             PrimalValueTapeTypes<real_t, real_t,
                                  MultiUseIndexManager<index_t>,
                                  InnerStatementEvaluator, DefaultChunkedData> > >;

} // namespace codi
#endif
