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
#include <gsCore/gsMultiBasis.h>
#include <gsNurbs/gsKnotVector.h>
#include <gsIO/gsOptionList.h>
#include <gsUtils/gsCombinatorics.h>
#include <gsAssembler/gsExprAssembler.h>

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

public:
    gsIetiMapper() : m_nrPrimalConstraints(0) {}

    void init(gsDofMapper dofMapperGlobal, const gsMatrix<T>& fixedPart)
    {
        const index_t nPatches = dofMapperGlobal.numPatches();

        m_dofMapperGlobal = give(dofMapperGlobal);
        m_dofMapperLocal.resize(nPatches);
        m_fixedPart.resize(nPatches);

        for (index_t k=0; k<nPatches; ++k)
        {
            const index_t nDofs = m_dofMapperGlobal.patchSize(k);
            m_dofMapperLocal[k].setIdentity(1,nDofs);
            for (index_t i=0; i<nDofs; ++i)
            {
                const index_t idx = m_dofMapperGlobal.index(i,k);
                if (m_dofMapperGlobal.is_boundary_index(idx))
                    m_dofMapperLocal[k].eliminateDof(i,0);
            }
            m_dofMapperLocal[k].finalize();

            const index_t szFixedPart = m_dofMapperLocal[k].boundarySize();
            m_fixedPart[k].setZero(szFixedPart,1);
            for (index_t i=0; i<nDofs; ++i)
            {
                const index_t idx = m_dofMapperGlobal.index(i,k);
                if (m_dofMapperGlobal.is_boundary_index(idx))
                {
                    const index_t globalBoundaryIdx = m_dofMapperGlobal.bindex(i,k);
                    const index_t localBoundaryIdx = m_dofMapperLocal[k].bindex(i,0);
                    m_fixedPart[k](localBoundaryIdx,0) = fixedPart(globalBoundaryIdx,0);
                }
            }
        }
        computeJumpMatrices();
    }

    void initFeSpace(const typename gsExprAssembler<T>::space& u, index_t k)
    {
        GISMO_ASSERT( u.mapper().size() ==  m_dofMapperLocal[k].size(),
            "gsIetiMapper::initFeSpace: The sizes do not agree." );
        /*if (m_dofMapperLocal[k].freeSize() == u.mapper().freeSize())
        {
            gsInfo << "Patch " << k << ": The number of Dirichlet values does match; the error of their values is "
                << (u.fixedPart()-m_fixedPart[k]).norm() << "\n";
        }
        else
        {
            gsInfo << "Patch " << k << ": The number of Dirichlet values does not match.\n";
            gsInfo << u.fixedPart().transpose() << "\n";
            gsInfo << m_fixedPart[k].transpose() << "\n";
        }*/
        const_cast<expr::gsFeSpace<real_t>&>(u).mapper() = m_dofMapperLocal[k];
        const_cast<expr::gsFeSpace<real_t>&>(u).fixedPart() = m_fixedPart[k];
    }


private:

    typedef std::vector< std::vector< std::pair<index_t,index_t> > > CouplingInfo;

    CouplingInfo getCoupling() const
    {
        CouplingInfo result;
        result.reserve(m_dofMapperGlobal.coupledSize());
        gsVector<index_t> tmp;
        tmp.setZero(m_dofMapperGlobal.freeSize(),1);
        const index_t numPatches = m_dofMapperGlobal.numPatches();
        for (index_t k=0; k<numPatches; ++k)
        {
            const index_t patchSize = m_dofMapperGlobal.patchSize(k);
            for (index_t i=0; i<patchSize; ++i)
            {
                const index_t j=m_dofMapperGlobal.index(i,k);
                if (m_dofMapperGlobal.is_coupled_index(j))
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

    void computeJumpMatrices()
    {
        m_jumpMatrices.resize(0);
        CouplingInfo coupling = getCoupling();
        const index_t couplingSize = coupling.size();

        const index_t numPatches = m_dofMapperGlobal.numPatches();
        GISMO_ASSERT( numPatches == m_dofMapperLocal.size(), "gsIetiMapper::jumpMatrices: "
            "The number of patches and the number of local dof mappers must match." );

        index_t numLagrangeMult = 0;
        for (index_t i=0; i<couplingSize; ++i)
        {
            const index_t n = coupling[i].size();
            numLagrangeMult += (n * (n-1))/2;
        }

        for (index_t i=0; i<numPatches; ++i)
            m_jumpMatrices.push_back(Transfer(numLagrangeMult, m_dofMapperLocal[i].freeSize()));

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
                const index_t localMappedIndex1 = m_dofMapperLocal[patch1].index(localIndex1,0);
                for (index_t j2=j1+1; j2<sz; ++j2)
                {
                    const index_t patch2 = coupling[i][j2].first;
                    const index_t localIndex2 = coupling[i][j2].second;
                    const index_t localMappedIndex2 = m_dofMapperLocal[patch2].index(localIndex2,0);
                    GISMO_ASSERT(multiplier<numLagrangeMult, "bug." );
                    m_jumpMatrices[patch1](multiplier,localMappedIndex1) = 1.;
                    m_jumpMatrices[patch2](multiplier,localMappedIndex2) = -1.;
                    ++multiplier;
                }
            }
        }
        GISMO_ASSERT( multiplier == numLagrangeMult, "Have:"<<multiplier<<"!="<<numLagrangeMult);

    }

public:

    Transfer jumpMatrix(index_t i) { return m_jumpMatrices[i]; }

    gsMatrix<T> constructGlobalSolutionFromLocalSolutions( const std::vector< gsMatrix<T> >& localContribs )
    {
        const index_t numPatches = m_dofMapperGlobal.numPatches();
        GISMO_ASSERT( numPatches == m_dofMapperLocal.size(), "");
        GISMO_ASSERT( numPatches == localContribs.size(), "");

        gsMatrix<T> result;
        result.setZero( m_dofMapperGlobal.freeSize(), 1/*==localContribs[0].cols()*/ );


        for (index_t k=0; k<numPatches; ++k)
        {
            const index_t sz=m_dofMapperLocal[k].size();
            for (index_t i=0; i<sz; ++i)
            {
                if (m_dofMapperLocal[k].is_free(i,0) && m_dofMapperGlobal.is_free(i,k))
                    result(m_dofMapperGlobal.index(i,k),0) = localContribs[k](m_dofMapperLocal[k].index(i,0),0);
            }
        }

        return result;
    }

    void cornersAsPrimalConstraints( const gsMultiBasis<T>& mb )
    {
        const index_t nPatches = m_dofMapperLocal.size();

        std::vector< std::vector<index_t> > corners;
        corners.reserve(nPatches);
        for (index_t i=0; i<nPatches; ++i)
        {
            std::vector<index_t> local_corners;
            for (boxCorner it = boxCorner::getFirst(mb.dim()); it!=boxCorner::getEnd(mb.dim()); ++it)
                local_corners.push_back( mb[i].functionAtCorner(it) );
            corners.push_back(give(local_corners));
        }


        if (m_primalConstraints.empty())
            m_primalConstraints.resize(nPatches);

        if (m_primalConstraintsMapper.empty())
            m_primalConstraintsMapper.resize(nPatches);

        GISMO_ASSERT(m_primalConstraints.size() == nPatches, "The number of primal constraints does not fit.");
        GISMO_ASSERT(m_primalConstraintsMapper.size() == nPatches, "The number of primal constraints mappers does not fit.");

        std::vector<index_t> corner_list;
        for (index_t i=0; i<nPatches; ++i)
        {
            std::vector< gsSparseVector<real_t> > localPrimalConstraints;
            std::vector< index_t > localPrimalConstraintsMapper;

            const index_t nCorners = corners[i].size();
            for (index_t j=0; j<nCorners; ++j)
            {
                const index_t localIndex = m_dofMapperLocal[i].index(corners[i][j],0);
                const index_t globalIndex = m_dofMapperGlobal.index(corners[i][j],i);
                const bool is_free = m_dofMapperGlobal.is_free_index(globalIndex);

                if (is_free) {

                    const index_t sz = corner_list.size();
                    index_t index; //std::find or better?
                    for (index=0; index<sz && corner_list[index]!=globalIndex; ++index)
                    {}

                    if (index==sz)
                        corner_list.push_back(globalIndex);

                    gsSparseVector<> constr(m_dofMapperLocal[i].freeSize()); //!
                    constr[localIndex] = 1;

                    m_primalConstraints[i].push_back(give(constr));
                    m_primalConstraintsMapper[i].push_back(index);
                }
            }
        }
        m_nrPrimalConstraints = corner_list.size();

    }

    index_t nrLagrangeMultipliers() {
        GISMO_ASSERT(! m_jumpMatrices.empty(), "Not yet known.");
        return m_jumpMatrices[0].rows();
    }

    index_t nrPrimalConstraints() {
        return m_nrPrimalConstraints;
    }

    const std::vector< gsSparseVector<real_t> > & primalConstraints(index_t i) const
    { return m_primalConstraints[i]; }

    const std::vector<index_t> & primalConstraintsMapper(index_t i) const
    { return m_primalConstraintsMapper[i]; }

private:
    gsDofMapper                                          m_dofMapperGlobal;
    std::vector<gsDofMapper>                             m_dofMapperLocal;
    std::vector< gsMatrix<T> >                           m_fixedPart;

    std::vector< gsSparseMatrix<real_t,RowMajor> >       m_jumpMatrices;

    index_t                                              m_nrPrimalConstraints;
    std::vector< std::vector< gsSparseVector<real_t> > > m_primalConstraints;
    std::vector< std::vector< index_t > >                m_primalConstraintsMapper;
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsIetiMapper.hpp)
#endif
