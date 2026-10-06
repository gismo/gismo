/** @file gsDofMapperCreator.hpp

    @brief implementation file for the gsDofMapper factory functions

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, C. Hofreither, A. Mantzaflaris
**/

#pragma once

#include <gsAssembler/gsDofMapperCreator.h>

#include <gsCore/gsFunctionSet.h>
#include <gsCore/gsBasis.h>
#include <gsCore/gsBoxTopology.h>
#include <gsCore/gsMultiBasis.h>

#include <gsMSplines/gsMappedBasis.h>

#include <gsPde/gsBoundaryConditions.h>

namespace gismo
{

namespace internal
{

// The creators take the function set of component c from bases[c].  The
// single-basis overloads repeat one set for every component, the
// per-component overloads supply one set each; everything below is written
// against that per-component view only, so that both paths extract
// interface traces, boundary dofs and corner functions the same way -- from
// the basis of the component being constrained.

// The domain dimension of \a bases if it is a gsMappedBasis, 0 otherwise.
// The dimensions share no common base class but gsFunctionSet, so each
// instantiated one is asked for in turn.
template<class T>
short_t mappedBasisDim(const gsFunctionSet<T> & bases)
{
    if (dynamic_cast<const gsMappedBasis<1,T>*>(&bases) != nullptr) return 1;
    if (dynamic_cast<const gsMappedBasis<2,T>*>(&bases) != nullptr) return 2;
    if (dynamic_cast<const gsMappedBasis<3,T>*>(&bases) != nullptr) return 3;
    return 0;
}

template<class T>
void checkPatchSide(const gsFunctionSet<T> & bases, const patchSide & ps,
                    const char * what)
{
    GISMO_ENSURE(0 <= ps.patch && ps.patch < bases.nPieces(),
                 "createMapper: "<<what<<" is set on patch "<<ps.patch
                 <<", but there are "<<bases.nPieces()<<" patches.");
    const index_t d = bases.basis(ps.patch).domainDim();
    GISMO_ENSURE(1 <= ps.index() && ps.index() <= 2*d,
                 "createMapper: "<<what<<" is set on side "<<ps.index()
                 <<" of patch "<<ps.patch<<", which has "<<2*d<<" sides.");
}

inline void checkComponentSelection(index_t cc, index_t nComp, const char * what)
{
    GISMO_ENSURE(-1 == cc || (0 <= cc && cc < nComp),
                 "createMapper: "<<what<<" selects component "<<cc
                 <<"; expected -1 (every component) or a component in [0,"<<nComp<<").");
}

// Patch and side references of the conditions that apply to unknown \a unk,
// checked before the mapper is touched.  A reference to a patch that does
// not exist would otherwise index past the function set.  Component
// selections are checked for every kind when \a allComponents is set.
// Otherwise only the Dirichlet selection is: the single-basis creator has
// always passed it on to the mapper, which requires a valid one, whereas it
// has always skipped a Clamped, Collapsed, coupled or corner condition whose
// component matches none, and it keeps doing so.
template<class T>
void checkConditions(const gsBoundaryConditions<T> & bc, index_t unk,
                     const gsFunctionSet<T> & bases, index_t nComp,
                     bool allComponents)
{
    static const char * const sideKinds[3] = {"Dirichlet", "Clamped", "Collapsed"};
    for (index_t t = 0; t != 3; ++t)
        for (typename gsBoundaryConditions<T>::const_iterator
             it = bc.begin(sideKinds[t]) ; it != bc.end(sideKinds[t]); ++it )
        {
            if (unk!=-1 && it->unknown() != unk) continue;
            const std::string what = std::string("a ") + sideKinds[t] + " condition";
            checkPatchSide(bases, it->ps, what.c_str());
            if (allComponents || 0 == t)
                checkComponentSelection(it->unkComponent(), nComp, what.c_str());
        }

    for (typename gsBoundaryConditions<T>::const_cpliterator
         it = bc.coupledBegin(); it != bc.coupledEnd(); ++it )
    {
        if (unk!=-1 && it->unknown!=-1 && it->unknown != unk) continue;
        checkPatchSide(bases, it->ifc.first() , "a coupled condition");
        checkPatchSide(bases, it->ifc.second(), "a coupled condition");
        if (allComponents)
            checkComponentSelection(it->component, nComp, "a coupled condition");
    }

    for (typename gsBoundaryConditions<T>::const_citerator
         it = bc.cornerBegin() ; it != bc.cornerEnd(); ++it )
    {
        if (unk!=-1 && it->unknown!=-1 && it->unknown != unk) continue;
        GISMO_ENSURE(0 <= it->patch && it->patch < bases.nPieces(),
                     "createMapper: a corner condition is set on patch "<<it->patch
                     <<", but there are "<<bases.nPieces()<<" patches.");
        const index_t d = bases.basis(it->patch).domainDim();
        GISMO_ENSURE(1 <= it->corner.m_index && it->corner.m_index <= (index_t(1) << d),
                     "createMapper: a corner condition is set on corner "<<it->corner.m_index
                     <<" of patch "<<it->patch<<", which has "<<(index_t(1) << d)<<" corners.");
        if (allComponents)
            checkComponentSelection(it->component, nComp, "a corner condition");
    }
}

// Patch and side references of the interfaces the conforming loop visits.
template<class T>
void checkInterfaces(const gsBoxTopology & topology, const gsFunctionSet<T> & bases)
{
    for ( gsBoxTopology::const_iiterator it = topology.iBegin();
          it != topology.iEnd(); ++it )
    {
        if (it->type() == interaction::contact) continue;
        checkPatchSide(bases, it->first() , "an interface of the topology");
        checkPatchSide(bases, it->second(), "an interface of the topology");
    }
}

// Component c of one patch is matched with component c of its neighbour,
// which presumes that the component's parametric direction is the same on
// both sides.  An interface whose direction map is not the identity (a
// rotated neighbour) would pair, e.g., a component normal to the interface
// on one side with one tangential to it on the other.  For distinct
// per-component bases that pairing is rejected rather than guessed; it is
// only detected otherwise when the two traces happen to differ in size.
//
// An interface that joins the same side of both patches (e.g. east to
// east) keeps every direction, but the two normal parametric directions
// point opposite ways in physical space.  A plain vector field is matched
// correctly across it; a Piola-mapped component space such as
// Raviart-Thomas needs the normal coefficients identified with a sign flip,
// which a dof mapper cannot express.  Distinct per-component bases are what
// such spaces are built from, so this pairing is rejected as well.
inline void checkAlignedInterfaces(const gsBoxTopology & topology)
{
    for ( gsBoxTopology::const_iiterator it = topology.iBegin();
          it != topology.iEnd(); ++it )
    {
        if (it->type() == interaction::contact) continue;
        const gsVector<index_t> & dirs = it->dirMap();
        for (index_t d = 0; d != dirs.size(); ++d)
            GISMO_ENSURE(dirs[d] == d,
                         "createMapper: the interface between "<<it->first()<<" and "
                         <<it->second()<<" maps parametric direction "<<d<<" onto direction "
                         <<dirs[d]<<".  Distinct per-component bases are matched component by "
                         "component, which requires every interface to keep each direction.");
        GISMO_ENSURE(it->first().side().index() != it->second().side().index(),
                     "createMapper: the interface between "<<it->first()<<" and "
                     <<it->second()<<" joins the same side of both patches.  Distinct "
                     "per-component bases are matched component by component, without the "
                     "sign flip of the normal component that a Piola-mapped space needs there.");
    }
}

// Glues every conforming interface of \a topology, each component against
// its own basis: a component's interface trace is determined by that
// component's knot vectors alone.  The traces are computed once per run of
// components sharing one function set, so a single-basis mapper costs one
// matchWith() per interface whatever its number of components.
template<class T>
void matchInterfaces(gsDofMapper & mapper,
                     const std::vector<const gsFunctionSet<T>*> & bases,
                     const gsBoxTopology & topology)
{
    gsMatrix<index_t> b1, b2;
    for ( gsBoxTopology::const_iiterator it = topology.iBegin();
          it != topology.iEnd(); ++it )
    {
        // Dofs must not be matched across a contact interface: the two
        // sides are physically distinct (cf. gsMultiBasis::repairInterfaces)
        if (it->type() == interaction::contact) continue;

        for (size_t c = 0; c != bases.size(); ++c)
        {
            if (0 == c || bases[c] != bases[c-1])
            {
                const gsBasis<T> & basis1 = bases[c]->basis(it->first().patch);
                const gsBasis<T> & basis2 = bases[c]->basis(it->second().patch);
                basis1.matchWith(*it, basis2, b1, b2);
                GISMO_ENSURE(b1.rows() == b2.rows(),
                             "createMapper: component "<<c<<" has "<<b1.rows()<<" dofs on side "
                             <<it->first()<<" but "<<b2.rows()<<" on side "<<it->second()
                             <<" of the same interface; they cannot be matched.");
            }
            mapper.matchDofs(it->first().patch, b1, it->second().patch, b2, c);
        }
    }
}

// Applies every boundary-condition kind the mapper knows of to the
// components each condition selects, taking boundary dofs and corner
// functions from the selected component's own basis.  Expects
// checkConditions() to have passed.  Within one condition the boundary dofs
// are extracted again only when the function set changes (\a prev), so the
// components of a single-basis mapper share one extraction.
template<class T>
void applyConditions(gsDofMapper & mapper,
                     const std::vector<const gsFunctionSet<T>*> & bases,
                     const gsBoundaryConditions<T> & bc, index_t unk)
{
    const index_t nComp = static_cast<index_t>(bases.size());
    gsMatrix<index_t> bnd, bnd1;
    const gsFunctionSet<T> * prev = nullptr;

    // Strong Dirichlet conditions
    for (typename gsBoundaryConditions<T>::const_iterator
         it = bc.begin("Dirichlet") ; it != bc.end("Dirichlet"); ++it )
    {
        if (unk!=-1 && it->unknown() != unk) continue;
        const index_t cc = it->unkComponent();
        prev = nullptr;
        for (index_t c = 0; c!=nComp; c++)
        {
            if (c!=cc && cc!=-1) continue;
            if (bases[c] != prev)
            {
                prev = bases[c];
                bnd = prev->basis(it->ps.patch).boundary(it->ps.side());
            }
            mapper.markBoundary(it->ps.patch, bnd, c);
        }
    }

    // Clamped boundary condition (per DoF)
    for (typename gsBoundaryConditions<T>::const_iterator
         it = bc.begin("Clamped") ; it != bc.end("Clamped"); ++it )
    {
        if (unk!=-1 && it->unknown() != unk) continue;
        const index_t cc = it->unkComponent();
        prev = nullptr;
        for (index_t c = 0; c!=nComp; c++)
        {
            if (c!=cc && cc!=-1) continue;
            if (bases[c] != prev)
            {
                prev = bases[c];
                bnd = prev->basis(it->ps.patch).boundary(it->ps.side());
                bnd1= prev->basis(it->ps.patch).boundaryOffset(it->ps.side(), 1);
                if (!it->ps.parameter())
                    bnd.swap(bnd1);
            }
            for (index_t k = 0; k < bnd.size(); ++k)
                mapper.matchDof(it->ps.patch, (bnd)(k, 0),
                                it->ps.patch, (bnd1)(k, 0), c);
        }
    }

    // Collapsed
    for (typename gsBoundaryConditions<T>::const_iterator
         it = bc.begin("Collapsed") ; it != bc.end("Collapsed"); ++it )
    {
        if (unk!=-1 && it->unknown() != unk) continue;
        const index_t cc = it->unkComponent();
        prev = nullptr;
        for (index_t c = 0; c!=nComp; c++)
        {
            if (c!=cc && cc!=-1) continue;
            if (bases[c] != prev)
            {
                prev = bases[c];
                bnd = prev->basis(it->ps.patch).boundary(it->ps.side());
            }
            // match all DoFs to the first one of the side
            for (index_t k = 0; k < bnd.size() - 1; ++k)
                mapper.matchDof(it->ps.patch, (bnd)(0, 0),
                                it->ps.patch, (bnd)(k + 1, 0), c);
        }
    }

    // Coupled boundary condition (per DoF)
    for (typename gsBoundaryConditions<T>::const_cpliterator
         it = bc.coupledBegin(); it != bc.coupledEnd(); ++it )
    {
        if (unk!=-1 && it->unknown!=-1 && it->unknown != unk) continue;
        const index_t cc = it->component;
        prev = nullptr;
        for (index_t c = 0; c!=nComp; c++)
        {
            if (c!=cc && cc!=-1) continue;
            if (bases[c] != prev)
            {
                prev = bases[c];
                bnd = prev->basis(it->ifc.first().patch).boundary(it->ifc.first().side());
                bnd1= prev->basis(it->ifc.second().patch).boundary(it->ifc.second().side());
                GISMO_ENSURE(bnd.rows() == bnd1.rows(),
                             "createMapper: a coupled condition couples "<<bnd.rows()<<" dofs of side "
                             <<it->ifc.first()<<" with "<<bnd1.rows()<<" dofs of side "
                             <<it->ifc.second()<<" in component "<<c<<".");
            }

            // match all DoFs to the first one of the side
            for (index_t k = 0; k < bnd.size() -1; ++k)
                mapper.matchDof(it->ifc.first() .patch, (bnd)(0, 0),
                                it->ifc.first() .patch, (bnd)(k + 1, 0), c);
            for (index_t k = 0; k < bnd1.size(); ++k)
                mapper.matchDof(it->ifc.second().patch, (bnd1)(k, 0),
                                it->ifc.first().patch,  (bnd)(k, 0), c);
        }
    }

    // Corners
    for (typename gsBoundaryConditions<T>::const_citerator
         it = bc.cornerBegin() ; it != bc.cornerEnd(); ++it )
    {
        if (unk!=-1 && it->unknown!=-1 && it->unknown != unk) continue;
        for (index_t c = 0; c!=nComp; ++c)
        {
            if (it->component!=-1 && it->component!=c) continue;
            mapper.eliminateDof(bases[c]->basis(it->patch).functionAtCorner(it->corner),
                                it->patch, c);
        }
    }
}

// The interfaces of \a t, each with its smaller patch side first, sorted.
inline std::vector<boundaryInterface> orientedInterfaces(const gsBoxTopology & t)
{
    std::vector<boundaryInterface> r;
    r.reserve(t.nInterfaces());
    for (gsBoxTopology::const_iiterator it = t.iBegin(); it != t.iEnd(); ++it)
        r.push_back(it->second() < it->first() ? it->getInverse() : *it);
    std::sort(r.begin(), r.end());
    return r;
}

// Topologies are compared as sets of records: patch count, dimension,
// boundary sides and interfaces, the latter oriented first side first and
// compared with their direction maps, orientations and interaction type.
// Labels are ignored; they play no part in mapper creation.
inline bool sameTopology(const gsBoxTopology & a, const gsBoxTopology & b)
{
    if (a.nBoxes() != b.nBoxes() || a.dim() != b.dim() ||
        a.nInterfaces() != b.nInterfaces() || a.nBoundary() != b.nBoundary())
        return false;

    std::vector<patchSide> ba = a.boundaries(), bb = b.boundaries();
    std::sort(ba.begin(), ba.end());
    std::sort(bb.begin(), bb.end());
    if (ba != bb) return false;

    const std::vector<boundaryInterface> ia = orientedInterfaces(a), ib = orientedInterfaces(b);
    for (size_t i = 0; i != ia.size(); ++i)
        if ( !(ia[i] == ib[i]) || ia[i].type() != ib[i].type() )
            return false;
    return true;
}

} // namespace internal

template<class T>
gsDofMapper createMapper(const gsFunctionSet<T>        & bases,
                         const gsBoxTopology           & topology,
                         const gsBoundaryConditions<T> & bc,
                         index_t nComp,
                         index_t unk,
                         bool    conforming,
                         bool    finalize)
{
    const bool hasBCs = bc.size() != 0;

    gsDofMapper mapper;

    // Only 2D and 3D mapped bases are numbered as a global identity; a 1D
    // one has always taken the patch-concatenated branch below.
    const short_t mappedDim = internal::mappedBasisDim(bases);
    if (2 == mappedDim || 3 == mappedDim)
    {
        mapper.setIdentity(bases.nPieces(), bases.size(), nComp);
    }
    else
    {
        GISMO_ASSERT(nComp>0,"Zero components");
        // Initialize offsets and dof holder: identical to the prologue of the
        // former gsDofMapper::init(), realized through the public
        // gsDofMapper(gsVector<index_t>,nComp) constructor.
        const index_t nPatches = bases.nPieces();
        gsVector<index_t> sz(nPatches);
        for (index_t k = 0; k != nPatches; ++k)
            sz[k] = bases.basis(k).size();

        mapper = gsDofMapper(sz, nComp);

        if (conforming)
        {
            internal::checkInterfaces(topology, bases);
            internal::matchInterfaces(mapper,
                std::vector<const gsFunctionSet<T>*>(mapper.numComponents(), &bases),
                topology);
        }
    }

    if (hasBCs)
    {
        internal::checkConditions(bc, unk, bases, mapper.numComponents(), false);
        internal::applyConditions(mapper,
            std::vector<const gsFunctionSet<T>*>(mapper.numComponents(), &bases),
            bc, unk);
    }

    if (finalize)
        mapper.finalize();
    return mapper;
}

template<class T>
gsDofMapper createMapper(const std::vector<const gsFunctionSet<T>*> & basesPerComp,
                         const gsBoxTopology           & topology,
                         const gsBoundaryConditions<T> & bc,
                         index_t unk, bool conforming, bool finalize)
{
    GISMO_ENSURE(!basesPerComp.empty(),
                 "createMapper: expecting one function set per component, got none.");
    const index_t nComp = static_cast<index_t>(basesPerComp.size());

    for (index_t c = 0; c != nComp; ++c)
    {
        GISMO_ENSURE(basesPerComp[c] != nullptr,
                     "createMapper: the function set of component "<<c<<" is a null pointer.");
        // A mapped basis reports mapped (global) indices from active() and
        // boundary(), while the conditions are applied here in the
        // patch-local indices of the underlying bases, so there is no
        // per-component meaning to give it.  The single-basis overload keeps
        // handling it.
        GISMO_ENSURE(0 == internal::mappedBasisDim(*basesPerComp[c]),
                     "createMapper: the function set of component "<<c<<" is a gsMappedBasis, "
                     "which the per-component overload does not accept; use the single-basis "
                     "createMapper(bases, topology, bc, nComp, ...) instead.");
    }

    const gsFunctionSet<T> & first = *basesPerComp.front();
    const index_t nPatches = first.nPieces();
    GISMO_ENSURE(nPatches > 0, "createMapper: the function set of component 0 has no patches.");
    for (index_t c = 1; c != nComp; ++c)
    {
        GISMO_ENSURE(basesPerComp[c]->nPieces() == nPatches,
                     "createMapper: component "<<c<<" has "<<basesPerComp[c]->nPieces()
                     <<" patches, component 0 has "<<nPatches<<"; all components must share "
                     "the same patches.");
        for (index_t k = 0; k != nPatches; ++k)
            GISMO_ENSURE(basesPerComp[c]->basis(k).domainDim() == first.basis(k).domainDim(),
                         "createMapper: patch "<<k<<" of component "<<c<<" has domain dimension "
                         <<basesPerComp[c]->basis(k).domainDim()<<", component 0 has "
                         <<first.basis(k).domainDim()<<".");
    }

    // Validated here, before any delegation, so that the same input is
    // rejected whether or not its components share one function set.
    if (conforming)
        internal::checkInterfaces(topology, first);
    internal::checkConditions(bc, unk, first, nComp, true);

    // One function set for every component is the single-basis case.  It is
    // delegated, so both overloads produce the same mapper for it -- one
    // not declared to have distinct component spaces.
    bool shared = true;
    for (index_t c = 1; c != nComp && shared; ++c)
        shared = (basesPerComp[c] == basesPerComp.front());
    if (shared)
        return createMapper(first, topology, bc, nComp, unk, conforming, finalize);

    // Distinct function sets: declared distinct even where every size
    // coincides (a Raviart-Thomas pair on a square mesh), because the
    // mapper cannot tell after the fact.
    std::vector<gsVector<index_t> > sz(nComp);
    for (index_t c = 0; c != nComp; ++c)
    {
        sz[c].resize(nPatches);
        for (index_t k = 0; k != nPatches; ++k)
            sz[c][k] = basesPerComp[c]->basis(k).size();
    }
    gsDofMapper mapper(sz, /*hasDistinctComponentSpaces=*/true);

    if (conforming)
    {
        internal::checkAlignedInterfaces(topology);
        internal::matchInterfaces(mapper, basesPerComp, topology);
    }

    if (bc.size() != 0)
        internal::applyConditions(mapper, basesPerComp, bc, unk);

    if (finalize)
        mapper.finalize();
    return mapper;
}

template<class T>
gsDofMapper createMapper(const std::vector<gsMultiBasis<T> > & basesPerComp,
                         const gsBoundaryConditions<T> & bc,
                         index_t unk, bool conforming, bool finalize)
{
    GISMO_ENSURE(!basesPerComp.empty(),
                 "createMapper: expecting one gsMultiBasis per component, got none.");

    std::vector<const gsFunctionSet<T>*> ptrs(basesPerComp.size());
    for (size_t c = 0; c != basesPerComp.size(); ++c)
    {
        GISMO_ENSURE(internal::sameTopology(basesPerComp[c].topology(),
                                            basesPerComp.front().topology()),
                     "createMapper: the topology of component "<<c
                     <<" differs from that of component 0.");
        ptrs[c] = &basesPerComp[c];
    }

    return createMapper(ptrs, basesPerComp.front().topology(), bc, unk, conforming, finalize);
}

template<class T>
gsDofMapper createMapper(const gsFunctionSet<T> & bases,
                         index_t nComp, bool conforming,
                         bool finalize)
{
    if (const gsMultiBasis<T> * mb = dynamic_cast<const gsMultiBasis<T>*>(&bases))
        return createMapper(bases, mb->topology(), gsBoundaryConditions<T>(),
                            nComp, 0, conforming, finalize);
    else
        return createMapper(bases, gsBoxTopology(), gsBoundaryConditions<T>(),
                            nComp, 0, conforming, finalize);
}

template<class T>
gsDofMapper createMapper(const gsFunctionSet<T> & bases, const gsBoxTopology & topology,
                         index_t nComp, bool conforming,
                         bool finalize)
{
    return createMapper(bases, topology, gsBoundaryConditions<T>(),
                        nComp, 0, conforming, finalize);
}

template<class T>
gsDofMapper createMapper(const gsFunctionSet<T> & bases, const gsBoundaryConditions<T> & bc,
                         index_t nComp, index_t unk, bool conforming,
                         bool finalize)
{
    if (const gsMultiBasis<T> * mb = dynamic_cast<const gsMultiBasis<T>*>(&bases))
        return createMapper(bases, mb->topology(), bc, nComp, unk, conforming, finalize);
    else
        return createMapper(bases, gsBoxTopology(), bc, nComp, unk, conforming, finalize);
}

template<class T>
gsDofMapper createMapper(const gsFunctionSet<T> & bases, const gsBoundaryConditions<T> & bc,
                         dirichlet::strategy ds, iFace::strategy is,
                         index_t nComp, index_t unk, bool finalize)
{
    const bool conforming = (is == iFace::glue);
    if (dirichlet::elimination == ds)
        return createMapper(bases, bc, nComp, unk, conforming, finalize);
    else
        return createMapper(bases, gsBoundaryConditions<T>(), nComp, unk, conforming, finalize);
}

}//namespace gismo
