/** @file gsIetiMapper.h

    @brief Algorithms for the assembling of matrices required for IETI-Solvers

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#pragma once

#include <gsCore/gsDofMapper.h>
#include <gsNurbs/gsKnotVector.h>
#include <gsIO/gsOptionList.h>
#include <gsUtils/gsCombinatorics.h>

namespace gismo
{
#define DEBUGVAR(a) gsInfo << "  " << #a << ": " << a << std::endl
#define DEBUGMATRIX(a) gsInfo << "  " << #a << ": " << a.rows() << " x " << a.cols() << std::endl

/** @brief
    Ieti Mapper

    This class contains algorithms for the assembling of matrices required for IETI-Solvrs.

    \ingroup Solver
*/
template< typename T >
class gsIetiMapper
{
    typedef gsSparseMatrix<T,RowMajor> Transfer;
    typedef memory::shared_ptr<Transfer> TransferPtr;

private:

    typedef std::vector< std::vector< std::pair<index_t,index_t> > > CouplingInfo;

    CouplingInfo getCoupling() const
    {
        CouplingInfo result;
        result.reserve(dm_global.coupledSize());
        gsVector<index_t> tmp;
        tmp.setZero(dm_global.freeSize(),1);
        const index_t numPatches = dm_global.numPatches();
        for (index_t k=0; k<numPatches; ++k)
        {
            const index_t patchSize = dm_global.patchSize(k);
            for (index_t i=0; i<patchSize; ++i)
            {
                const index_t j=dm_global.index(i,k);
                if (dm_global.is_coupled_index(j))
                {
                    if (tmp[j]==0)
                    {
                        std::vector< std::pair<index_t,index_t> > vec;
                        vec.push_back(std::pair<index_t,index_t>(k,i));
                        result.push_back(give(vec));
                        tmp[j] = result.size();
                    }
                    else
                        result[tmp[j]-1].push_back(std::pair<index_t,index_t>(k,i));
                }
            }
        }
        return result;
    }

public:


    /// Get jump matrices for local subproblems
    void computeJumpMatrices()
    {
        m_jumpMatrices.resize(0);
        CouplingInfo coupling = getCoupling();
        const index_t couplingSize = coupling.size();

        const index_t numPatches = dm_global.numPatches();
        GISMO_ASSERT( numPatches == dm_local.size(), "gsIetiMapper::jumpMatrices: "
            "The number of patches and the number of local dof mappers must match." );

        index_t numLagrangeMult = 0;
        for (index_t i=0; i<couplingSize; ++i)
        {
            const index_t n = coupling[i].size();
            numLagrangeMult += (n * (n-1))/2;
        }

        for (index_t i=0; i<numPatches; ++i)
            m_jumpMatrices.push_back(Transfer(numLagrangeMult, dm_local[i].freeSize()));

        index_t multiplier = 0;
        for (index_t i=0; i<couplingSize; ++i)
        {
            const index_t sz = coupling[i].size();
            GISMO_ASSERT( sz>1, "Found coupled dof that is not coupled to any other dof." );
            // We implement fully redundant. TODO: Allow alternatives.
            for (index_t j1=0; j1<sz-1; ++j1)
            {
                const index_t patch1 = coupling[i][j1].first;
                const index_t localIndex1 = coupling[i][j1].second;
                const index_t localMappedIndex1 = dm_local[patch1].index(localIndex1,0);
                for (index_t j2=j1+1; j2<sz; ++j2)
                {
                    const index_t patch2 = coupling[i][j2].first;
                    const index_t localIndex2 = coupling[i][j2].second;
                    const index_t localMappedIndex2 = dm_local[patch2].index(localIndex2,0);
                    GISMO_ASSERT(multiplier<numLagrangeMult, "bug." );
                    m_jumpMatrices[patch1](multiplier,localMappedIndex1) = 1.;
                    m_jumpMatrices[patch2](multiplier,localMappedIndex2) = -1.;
                    ++multiplier;
                }
            }
        }
        GISMO_ASSERT( multiplier == numLagrangeMult, "Have:"<<multiplier<<"!="<<numLagrangeMult);

    }

    Transfer jumpMatrix(index_t i) { return m_jumpMatrices[i]; }

    gsMatrix<T> constructGlobalSolutionFromLocalSolutions( const std::vector< gsMatrix<T> >& localContribs )
    {
        const index_t numPatches = dm_global.numPatches();
        GISMO_ASSERT( numPatches == dm_local.size(), "");
        GISMO_ASSERT( numPatches == localContribs.size() || numPatches == localContribs.size()-1, "");

        // TODO: incorporate values at primal dofs

        gsMatrix<T> result;
        result.setZero( dm_global.freeSize(), 1/*==localContribs[0].cols()*/ );


        for (index_t k=0; k<numPatches; ++k)
        {
            const index_t sz=dm_local[k].size();
            for (index_t i=0; i<sz; ++i)
            {
                if (dm_local[k].is_free(i,0) && dm_global.is_free(i,k))
                    result(dm_global.index(i,k),0) = localContribs[k](dm_local[k].index(i,0),0);
            }
        }

        return result;
    }

    void cornersAsPrimalConstraints(
        const std::vector< std::vector<index_t> >&       corners
    )
    {
        const index_t nPatches = dm_local.size();

        if (primalConstraints.empty())
            primalConstraints.resize(nPatches);

        if (primalConstraintsMapper.empty())
            primalConstraintsMapper.resize(nPatches);

        GISMO_ASSERT(primalConstraints.size() == nPatches, "The number of primal constraints does not fit.");
        GISMO_ASSERT(primalConstraints.size() == nPatches, "The number of primal constraints mappers does not fit.");

        std::vector<index_t> corner_list;
        for (index_t i=0; i<nPatches; ++i)
        {
            std::vector< gsSparseVector<real_t> > localPrimalConstraints;
            std::vector< index_t > localPrimalConstraintsMapper;

            const index_t nCorners = corners[i].size();
            for (index_t j=0; j<nCorners; ++j)
            {
                const index_t localIndex = dm_local[i].index(corners[i][j],0);
                const index_t globalIndex = dm_global.index(corners[i][j],i);
                const bool is_free = dm_global.is_free_index(globalIndex);

                if (is_free) {

                    const index_t sz = corner_list.size();
                    index_t index; //std::find or better?
                    for (index=0; index<sz && corner_list[index]!=globalIndex; ++index)
                    {}

                    if (index==sz)
                        corner_list.push_back(globalIndex);

                    gsSparseVector<> constr(dm_local[i].freeSize()); //!
                    constr[localIndex] = 1;

                    primalConstraints[i].push_back(give(constr));
                    primalConstraintsMapper[i].push_back(index);
                }
            }
        }
        nrPrimalConstraints = corner_list.size();

    }



public:
    gsDofMapper dm_global;
    std::vector<gsDofMapper> dm_local;
    std::vector< gsSparseMatrix<real_t,RowMajor> > m_jumpMatrices;

    index_t nrPrimalConstraints;
    std::vector< std::vector< gsSparseVector<real_t> > > primalConstraints;
    std::vector< std::vector< index_t > >                primalConstraintsMapper;
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsIetiMapper.hpp)
#endif
