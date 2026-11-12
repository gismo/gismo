/** @file gsDofMapperCreator.h

    @brief Provides free factory functions building a gsDofMapper from a
    function set, a topology and boundary conditions.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst, A. Mantzaflaris
*/

#pragma once

#include <vector>

#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsDofMapper.h>
#include <gsCore/gsBoxTopology.h>
#include <gsPde/gsBoundaryConditions.h>
#include <gsAssembler/gsAssemblerOptions.h> // dirichlet::strategy, iFace::strategy

namespace gismo
{

/** @brief Builds a gsDofMapper with a possibly different basis per component.

    \param bases      \a bases[c] is the basis of component \a c; \a bases.size()
                      is the number of components. The pointers are not owned
                      and must stay valid for the duration of the call.
    \param topology   interface topology between patches. It is valid to pass
                      an empty gsBoxTopology() for single-patch bases or when
                      no inter-patch conformity is required; in that case no
                      interface matching is performed.
    \param bc         boundary conditions; may be empty
    \param unk        unknown index the conditions are filtered by; -1 accepts all
    \param conforming if true, matching dofs across the interfaces of \a topology
    \param finalize   if true, gsDofMapper::finalize() is called before returning

    Requirements, all checked with GISMO_ENSURE: \a bases is non-empty, no entry
    is null, all entries have the same nPieces(), and gsMappedBasis entries are
    only accepted when all entries are the same pointer.

    Interfaces are matched per component, with that component's own basis.
    Each component's basis must itself be conforming across the interfaces when
    \a conforming is true. Vector fields whose components are exchanged across a
    rotated interface (Piola-type maps) are not supported.

    \note Interfaces of type interaction::contact are skipped by the conforming
    loop: matching dofs across a contact interface is wrong (the two sides are
    physically distinct). This mirrors gsMultiBasis<T>::repairInterfaces and the
    pre-existing gsFeSpace behaviour.

    \ingroup Assembler
*/
template<class T>
gsDofMapper createMapper(const std::vector<const gsFunctionSet<T>*> & bases,
                         const gsBoxTopology           & topology,
                         const gsBoundaryConditions<T> & bc,
                         index_t unk,
                         bool    conforming,
                         bool    finalize = false);

/** @brief Builds a gsDofMapper for the (possibly multi-patch) function set
    \a bases, with the same basis in each of the \a nComp components.

    Equivalent to the per-component overload called with \a nComp copies of
    \a &bases.

    \param bases      the function set; each piece is one patch
    \param topology   interface topology between patches. It is valid to pass
                      an empty gsBoxTopology() for single-patch bases or when
                      no inter-patch conformity is required; in that case no
                      interface matching is performed.
    \param bc         boundary conditions; may be empty
    \param nComp      number of components
    \param unk        unknown index the conditions are filtered by; -1 accepts all
    \param conforming if true, matching dofs across the interfaces of \a topology
    \param finalize   if true, gsDofMapper::finalize() is called before returning

    \note Interfaces of type interaction::contact are skipped by the conforming
    loop: matching dofs across a contact interface is wrong (the two sides are
    physically distinct). This mirrors gsMultiBasis<T>::repairInterfaces and the
    pre-existing gsFeSpace behaviour.

    \ingroup Assembler
*/
template<class T>
gsDofMapper createMapper(const gsFunctionSet<T>        & bases,
                         const gsBoxTopology           & topology,
                         const gsBoundaryConditions<T> & bc,
                         index_t nComp,
                         index_t unk,
                         bool    conforming,
                         bool    finalize = false);

/// \brief Convenience overload; the topology is taken from \a bases when it is
/// a gsMultiBasis, otherwise an empty topology is used.
/// \ingroup Assembler
template<class T>
gsDofMapper createMapper(const gsFunctionSet<T> & bases,
                         index_t nComp = 1, bool conforming = true,
                         bool finalize = false);

/// \brief Convenience overload with an explicit \a topology and no boundary
/// conditions.
/// \ingroup Assembler
template<class T>
gsDofMapper createMapper(const gsFunctionSet<T> & bases, const gsBoxTopology & topology,
                         index_t nComp = 1, bool conforming = true,
                         bool finalize = false);

/// \brief Convenience overload with boundary conditions; the topology is taken
/// from \a bases when it is a gsMultiBasis, otherwise an empty topology is used.
/// \ingroup Assembler
template<class T>
gsDofMapper createMapper(const gsFunctionSet<T> & bases, const gsBoundaryConditions<T> & bc,
                         index_t nComp = 1, index_t unk = 0, bool conforming = true,
                         bool finalize = false);

/** @brief Strategy-enum form; replaces the gsMultiBasis<T>::getMapper family.

    The boundary conditions are only taken into account when
    \a ds == dirichlet::elimination; interfaces are glued iff
    \a is == iFace::glue.

    \note On the \c multigrid branch the corresponding
    gsMultiBasis<T>::getMapper(dirichlet::strategy,iFace::strategy,...) path was
    broken: its no-BC branch forwarded to a constructor that always matched the
    interfaces, silently discarding the requested interface strategy (so a
    \c dg or \c smooth request was glued anyway). Honouring \a is here is a
    deliberate BUG FIX restoring the semantics of branch \c dev.

    \ingroup Assembler
*/
template<class T>
gsDofMapper createMapper(const gsFunctionSet<T> & bases, const gsBoundaryConditions<T> & bc,
                         dirichlet::strategy ds, iFace::strategy is,
                         index_t nComp = 1, index_t unk = 0, bool finalize = false);

#ifdef GISMO_WITH_PYBIND11
void pybind11_init_gsDofMapperCreator(pybind11::module &m);
#endif

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsDofMapperCreator.hpp)
#endif
