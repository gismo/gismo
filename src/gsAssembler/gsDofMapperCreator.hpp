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

    if (hasBCs &&
        (dynamic_cast<const gsMappedBasis<2,T>*>(&bases) != nullptr ||
         dynamic_cast<const gsMappedBasis<3,T>*>(&bases) != nullptr))
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
            gsMatrix<index_t> b1, b2;
            for ( gsBoxTopology::const_iiterator it = topology.iBegin();
                    it != topology.iEnd(); ++it )
            {
                // Dofs must not be matched across a contact interface: the two
                // sides are physically distinct (cf. gsMultiBasis::repairInterfaces)
                if (it->type() == interaction::contact) continue;

                const gsBasis<T> & basis1 = bases.basis(it->first().patch);
                const gsBasis<T> & basis2 = bases.basis(it->second().patch);
                basis1.matchWith(*it, basis2, b1, b2);
                for (size_t i=0; i!=mapper.componentsSize(); ++i)
                    mapper.matchDofs(it->first().patch, b1,it->second().patch, b2,i);
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

            bnd = bases.basis(it->ps.patch).boundary(it->ps.side());
            mapper.markBoundary(it->ps.patch, bnd, it->unkComponent());
        }

        // Clamped boundary condition (per DoF)
        for (typename gsBoundaryConditions<T>::const_iterator
             it = bc.begin("Clamped") ; it != bc.end("Clamped"); ++it )
        {
            if (unk!=-1 && it->unknown() != unk) continue;

            GISMO_ASSERT(static_cast<size_t>(it->ps.patch) < mapper.numPatches(),
                            "Problem: a boundary condition is set on a patch id which does not exist.");

            const index_t cc = it->unkComponent();
            bnd = bases.basis(it->ps.patch).boundary(it->ps.side());
            bnd1= bases.basis(it->ps.patch).boundaryOffset(it->ps.side(), 1);
            if (!it->ps.parameter())
                bnd.swap(bnd1);
            for (index_t c = 0; c!=nComp; c++) // for all components
            {
                if (c==cc || cc==-1 )
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

            GISMO_ASSERT(static_cast<size_t>(it->ps.patch) < mapper.numPatches(),
                            "Problem: a boundary condition is set on a patch id which does not exist.");
            const index_t cc = it->unkComponent();
            bnd = bases.basis(it->ps.patch).boundary(it->ps.side());
            // match all DoFs to the first one of the side
            for (index_t c = 0; c!=nComp; c++) // for all components
            {
                if (c==cc || cc==-1)
                    for (index_t k = 0; k < bnd.size() - 1; ++k)
                        mapper.matchDof(it->ps.patch, (bnd)(0, 0),
                                        it->ps.patch, (bnd)(k + 1, 0), c);
            }
        }

        // Coupled boundary condition (per DoF)
        for (typename gsBoundaryConditions<T>::const_cpliterator
             it = bc.coupledBegin(); it != bc.coupledEnd(); ++it )
        {
            if (unk!=-1 && it->unknown != unk) continue;

            GISMO_ASSERT(static_cast<size_t>(it->ifc.first().patch) < mapper.numPatches(),
                            "Problem: a boundary condition is set on a patch id which does not exist.");
            GISMO_ASSERT(static_cast<size_t>(it->ifc.second().patch) < mapper.numPatches(),
                            "Problem: a boundary condition is set on a patch id which does not exist.");

            const index_t cc = it->component;
            bnd = bases.basis(it->ifc.first().patch).boundary(it->ifc.first().side());
            bnd1= bases.basis(it->ifc.second().patch).boundary(it->ifc.second().side());
            GISMO_ASSERT(bnd.rows() == bnd1.rows(),
                            "Problem: trying to couple boundaries of different size.");

            // match all DoFs to the first one of the side
            for (index_t c = 0; c!=nComp; c++) // for all components
            {
                if (c==cc || cc==-1)
                {
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
            if (unk!=-1 && it->unknown != unk) continue;
            for (index_t r = 0; r!=nComp; ++r)
            {
                if (it->component!=-1 && it->component!=r) continue;
                GISMO_ASSERT(static_cast<size_t>(it->patch) < mapper.numPatches(),
                                "Problem: a corner boundary condition is set on a patch id which does not exist.");
                mapper.eliminateDof(bases.basis(it->patch).functionAtCorner(it->corner),
                                    it->patch, it->component);
            }
        }
    }

    if (finalize)
        mapper.finalize();
    return mapper;
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
