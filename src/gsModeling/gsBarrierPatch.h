/** @file gsBarrierPatch.h

    @brief This software facilitates the creation of analysis-suitable
    parameterizations from given boundary representations. Serving as a
    reference implementation, it embodies the methods and concepts detailed
    in Ye Ji's doctoral research.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Ye Ji, H.M. Verhelst
*/

#pragma once

#include <gismo.h>
#include <gsHLBFGS/gsHLBFGS.h>

using namespace gismo;

namespace gismo
{
/**
 * \brief Computes a patch parametrization given a set of boundary geometries.
 * Parametrization is not guaranteed to be non-singular. Works for planar surfaces and volumes.
 *
 * \tparam d domain dimension
 * \tparam T Coefficient type
 * \ingroup Modeling
 */
template<short_t d, typename T=real_t>
class gsBarrierPatch
{
 public:
  /// Constructs the object using a given multi-patch and a degree of freedom mapper.
  explicit gsBarrierPatch(const gsMultiPatch<T> &mp, const gsDofMapper &mapper);

  /// Constructs the object using a given multi-patch and an optional patch-wise flag (default true).
  explicit gsBarrierPatch(const gsMultiPatch<T> &mp, bool patchWise = true);

  /// Sets the mapper.
  void setMapper(const gsDofMapper &mapper) { m_mapper = mapper; };

  /// Computes analysis-suitable parameterizations using different methods.
  void compute();

  /// Returns the result in a multi-patch format.
  const gsMultiPatch<T> &result() const { return m_mp; };

  /// Returns the options list.
  gsOptionList &options() { return m_options; };

  /// Sets the default options.
  void defaultOptions();

 private:
  /// Creates a mapper.
  void _makeMapper();

  /// Creates a mapper for a single patch.
  gsDofMapper _makeMapperOnePatch(const gsGeometry<T>&currPatch) const;

  /// Creates a mapper for global patches.
  void _makeMapperGlobalPatches();

  /// Creates a mapper for local patches.
  void _makeMapperLocalPatches();

  /// Log information about the mapper.
  void logMapperInformation();

  mutable gsExprEvaluator<T> m_evaluator;
  mutable gsExprAssembler<T> m_assembler;
  gsDofMapper m_mapper;

  mutable gsMultiPatch<T> m_mp;
  gsMultiBasis<T> m_mb;
  gsMultiPatch<T> m_bRep;

  std::string m_filename;
  T m_boxsize = 1.0;

  gsVector<T, d> m_boundingBoxLeftBottomCorner;
  gsVector<T, d> m_scalingVec;

  size_t m_freeInterface = 1;
  gsOptionList m_options;
};
}// namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsBarrierPatch.hpp)
#endif
