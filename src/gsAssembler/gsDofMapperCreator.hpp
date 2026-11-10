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

template<class T>
gsDofMapper createMapper(const std::vector<const gsFunctionSet<T>*> & bases,
                         const gsBoxTopology           & topology,
                         const gsBoundaryConditions<T> & bc,
                         index_t unk,
                         bool    conforming,
                         bool    finalize)
{
    GISMO_ENSURE(!bases.empty(), "createMapper: the list of component bases is empty.");
    const index_t nComp = static_cast<index_t>(bases.size());
    for (index_t c = 0; c != nComp; ++c)
        GISMO_ENSURE(bases[c] != nullptr,
                     "createMapper: the basis of component " << c << " is null.");
    const index_t nPatches = bases[0]->nPieces();
    for (index_t c = 1; c != nComp; ++c)
        GISMO_ENSURE(bases[c]->nPieces() == nPatches,
                     "createMapper: component " << c << " has " << bases[c]->nPieces()
                     << " pieces, but component 0 has " << nPatches << ".");

    const bool hasBCs = bc.size() != 0;

    gsDofMapper mapper;

    if (dynamic_cast<const gsMappedBasis<2,T>*>(bases[0]) != nullptr ||
        dynamic_cast<const gsMappedBasis<3,T>*>(bases[0]) != nullptr)
    {
        for (index_t c = 1; c != nComp; ++c)
            GISMO_ENSURE(bases[c] == bases[0],
                         "createMapper: a mapped basis requires the same basis in every component (component "
                         << c << " differs).");
        mapper.setIdentity(nPatches, bases[0]->size(), nComp);
    }
    else
    {
        gsMatrix<index_t> sz(nPatches, nComp);
        for (index_t c = 0; c != nComp; ++c)
            for (index_t k = 0; k != nPatches; ++k)
                sz(k,c) = bases[c]->basis(k).size();

        mapper = gsDofMapper(sz);

        if (conforming)
        {
            gsMatrix<index_t> b1, b2;
            for ( gsBoxTopology::const_iiterator it = topology.iBegin();
                    it != topology.iEnd(); ++it )
            {
                // Dofs must not be matched across a contact interface: the two
                // sides are physically distinct (cf. gsMultiBasis::repairInterfaces)
                if (it->type() == interaction::contact) continue;

                for (index_t c = 0; c != nComp; ++c)
                {
                    // b1,b2 are reused for consecutive components sharing a basis
                    if (c == 0 || bases[c] != bases[c-1])
                        bases[c]->basis(it->first().patch).matchWith(
                            *it, bases[c]->basis(it->second().patch), b1, b2);
                    mapper.matchDofs(it->first().patch, b1, it->second().patch, b2, c);
                }
            }
        }
    }

    if (hasBCs)
    {
        // Strong Dirichlet conditions
        gsMatrix<index_t> bnd, bnd1;
        for (typename gsBoundaryConditions<T>::const_iterator
             it = bc.begin("Dirichlet") ; it != bc.end("Dirichlet"); ++it )
        {
            if (unk!=-1 && it->unknown() != unk) continue;

            GISMO_ASSERT(static_cast<size_t>(it->ps.patch) < mapper.numPatches(),
                        "Problem: a boundary condition is set on a patch id which does not exist.");

            const index_t cc = it->unkComponent();
            GISMO_ASSERT(cc >= -1 && cc < nComp,
                         "Problem: a boundary condition is set on a component which does not exist.");
            for (index_t c = 0; c != nComp; ++c)
            {
                if (cc == -1 || c == cc)
                {
                    bnd = bases[c]->basis(it->ps.patch).boundary(it->ps.side());
                    mapper.markBoundary(it->ps.patch, bnd, c);
                }
            }
        }

        // Clamped boundary condition (per DoF)
        for (typename gsBoundaryConditions<T>::const_iterator
             it = bc.begin("Clamped") ; it != bc.end("Clamped"); ++it )
        {
            if (unk!=-1 && it->unknown() != unk) continue;

            GISMO_ASSERT(static_cast<size_t>(it->ps.patch) < mapper.numPatches(),
                            "Problem: a boundary condition is set on a patch id which does not exist.");

            const index_t cc = it->unkComponent();
            for (index_t c = 0; c!=nComp; c++) // for all components
            {
                if (c==cc || cc==-1 )
                {
                    bnd = bases[c]->basis(it->ps.patch).boundary(it->ps.side());
                    bnd1= bases[c]->basis(it->ps.patch).boundaryOffset(it->ps.side(), 1);
                    if (!it->ps.parameter())
                        bnd.swap(bnd1);
                    for (index_t k = 0; k < bnd.size(); ++k)
                        mapper.matchDof(it->ps.patch, (bnd)(k, 0),
                                        it->ps.patch, (bnd1)(k, 0), c);
                }
            }
        }

        // Collapsed
        for (typename gsBoundaryConditions<T>::const_iterator
             it = bc.begin("Collapsed") ; it != bc.end("Collapsed"); ++it )
        {
            if (unk!=-1 && it->unknown() != unk) continue;

            GISMO_ASSERT(static_cast<size_t>(it->ps.patch) < mapper.numPatches(),
                            "Problem: a boundary condition is set on a patch id which does not exist.");
            const index_t cc = it->unkComponent();
            for (index_t c = 0; c!=nComp; c++) // for all components
            {
                if (c==cc || cc==-1)
                {
                    bnd = bases[c]->basis(it->ps.patch).boundary(it->ps.side());
                    // match all DoFs to the first one of the side
                    for (index_t k = 0; k < bnd.size() - 1; ++k)
                        mapper.matchDof(it->ps.patch, (bnd)(0, 0),
                                        it->ps.patch, (bnd)(k + 1, 0), c);
                }
            }
        }

        // Coupled boundary condition (per DoF)
        for (typename gsBoundaryConditions<T>::const_cpliterator
             it = bc.coupledBegin(); it != bc.coupledEnd(); ++it )
        {
            if (unk!=-1 && it->unknown!=-1 && it->unknown != unk) continue;

            GISMO_ASSERT(static_cast<size_t>(it->ifc.first().patch) < mapper.numPatches(),
                            "Problem: a boundary condition is set on a patch id which does not exist.");
            GISMO_ASSERT(static_cast<size_t>(it->ifc.second().patch) < mapper.numPatches(),
                            "Problem: a boundary condition is set on a patch id which does not exist.");

            const index_t cc = it->component;
            for (index_t c = 0; c!=nComp; c++) // for all components
            {
                if (c==cc || cc==-1)
                {
                    bnd = bases[c]->basis(it->ifc.first().patch).boundary(it->ifc.first().side());
                    bnd1= bases[c]->basis(it->ifc.second().patch).boundary(it->ifc.second().side());
                    GISMO_ASSERT(bnd.rows() == bnd1.rows(),
                                    "Problem: trying to couple boundaries of different size.");

                    // match all DoFs to the first one of the side
                    for (index_t k = 0; k < bnd.size() -1; ++k)
                        mapper.matchDof(it->ifc.first() .patch, (bnd)(0, 0),
                                        it->ifc.first() .patch, (bnd)(k + 1, 0), c);
                    for (index_t k = 0; k < bnd1.size(); ++k)
                        mapper.matchDof(it->ifc.second().patch, (bnd1)(k, 0),
                                        it->ifc.first().patch,  (bnd)(k, 0), c);
                }
            }
        }

        // Corners
        for (typename gsBoundaryConditions<T>::const_citerator
                it = bc.cornerBegin() ; it != bc.cornerEnd(); ++it )
        {
            if (unk!=-1 && it->unknown!=-1 && it->unknown != unk) continue;
            for (index_t r = 0; r!=nComp; ++r)
            {
                if (it->component!=-1 && it->component!=r) continue;
                GISMO_ASSERT(static_cast<size_t>(it->patch) < mapper.numPatches(),
                                "Problem: a corner boundary condition is set on a patch id which does not exist.");
                mapper.eliminateDof(bases[r]->basis(it->patch).functionAtCorner(it->corner),
                                    it->patch, r);
            }
        }
    }

    if (finalize)
        mapper.finalize();
    return mapper;
}

template<class T>
gsDofMapper createMapper(const gsFunctionSet<T>        & bases,
                         const gsBoxTopology           & topology,
                         const gsBoundaryConditions<T> & bc,
                         index_t nComp,
                         index_t unk,
                         bool    conforming,
                         bool    finalize)
{
    return createMapper<T>(std::vector<const gsFunctionSet<T>*>(static_cast<size_t>(nComp), &bases),
                           topology, bc, unk, conforming, finalize);
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
