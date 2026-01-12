/** @file gsHDomain.hpp

    @brief Provides implementation of the HDomain class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): G. Kiss, A. Mantzaflaris
*/

#pragma once

#include <gsDomain/gsHDomain.h>
#include <gsHSplines/gsHTensorBasis.h>
#include <memory>

namespace gismo
{

template<short_t d, class T, class Z>
gsHDomain<d,T,Z>::gsHDomain(const gsHTree<d,Z>& tree,
                       const gsHTensorBasis<d,T>& basis)
    :
    m_tree(tree),
    m_basis(basis)
{
}

template<short_t d, class T, class Z>
typename gsHDomain<d,T,Z>::domainIter gsHDomain<d,T,Z>::beginAll() const
{
    return domainIter(new gsHDomainIterator<T,d,Z>(m_tree,m_basis, this->patchId()));
}

template<short_t d, class T, class Z>
typename gsHDomain<d,T,Z>::domainIter gsHDomain<d,T,Z>::beginBdr(const boxSide bs) const
{
    return domainIter(new gsHDomainBoundaryIterator<T,d,Z>(m_tree,m_basis, bs, this->patchId()));
}

template<short_t d, class T, class Z>
size_t gsHDomain<d,T,Z>::numElements() const
{
    leafIterator it = m_tree.beginLeafIterator();
    size_t nel(0);
    while (it.good())
    {
        if (m_basis.manualLevels() )
        {
            index_t ll, uu;
            size_t nel_local = 1;
            for (short_t i = 0; i < d; ++i)
            {
                ll = it.lowerCorner()[i];
                uu = it.upperCorner()[i];
                m_basis._diadicIndexToKnotIndex(it.level(),i,ll);
                m_basis._diadicIndexToKnotIndex(it.level(),i,uu);
                nel_local *= uu - ll;
            }
            nel += nel_local;
        }
        else
            nel += ( it.upperCorner() - it.lowerCorner() ).prod();
        it.next();
    }
    return nel;
}

template<short_t d, class T, class Z>
size_t gsHDomain<d,T,Z>::numElementsBdr(boxSide const & s) const
{
    GISMO_ASSERT(s != boundary::none, "Not implemented");
    leafIterator it = m_tree.beginLeafIterator();
    size_t nel(0);
    size_t nel_local;
    while (it.good())
    {
        if  (leafOnBoundary(s,it))
        {
            nel_local = 1;
            for (short_t i = 0; i < d; ++i)
                if (i != s.direction())
                {
                    if (m_basis.manualLevels() )
                    {
                        index_t ll = it.lowerCorner()[i];
                        index_t uu = it.upperCorner()[i];
                        m_basis._diadicIndexToKnotIndex(it.level(),s.direction(),ll);
                        m_basis._diadicIndexToKnotIndex(it.level(),s.direction(),uu);
                        nel_local *= uu - ll;
                    }
                    else
                        nel_local *= it.upperCorner()[i] - it.lowerCorner()[i];
                }
            nel +=  nel_local;
        }
        it.next();
    }
    return nel;
}

template<short_t d, class T, class Z>
short_t gsHDomain<d,T,Z>::degree(short_t i) const
{
    return m_basis.degree(i);
}

template<short_t d, class T, class Z>
short_t gsHDomain<d,T,Z>::dim() const { return d; }

template<short_t d, class T, class Z>
gsMatrix<T> gsHDomain<d,T,Z>::boundingBox() const
{
    return m_basis.support();
}

template<short_t d, class T, class Z>
const gsHTree<d,Z> & gsHDomain<d,T,Z>::tree() const { return m_tree; }

template<short_t d, class T, class Z>
bool gsHDomain<d,T,Z>::leafOnBoundary(const boxSide & s, const leafIterator leaf) const
{
    if ( s.parameter() )
    {
        // AM: a little ugly for now, to be improved
        size_t diadicSize;
        if (m_basis.manualLevels() )
        {
            const gsKnotVector<T> & kv = m_basis.tensorLevel(leaf.level()).knots(s.direction());
            index_t start = 0;
            index_t end  = kv.uSize()-1;
            m_basis._diadicIndexToKnotIndex(leaf.level(),s.direction(),start);
            m_basis._diadicIndexToKnotIndex(leaf.level(),s.direction(),end);
            diadicSize = end - start;
        }
        else
            diadicSize = m_basis.tensorLevel(leaf.level()).knots(s.direction()).uSize() - 1;
        return static_cast<size_t>(leaf.upperCorner().at(s.direction()) ) == diadicSize;// todo: more efficient
    }
    else
        return leaf.lowerCorner().at(s.direction()) == 0;
}


template<short_t d, class T, class Z>
typename gsDomain<T>::Ptr gsHDomain<d,T,Z>::decompose(size_t npieces) const
{
    auto result = memory::make_shared(new gsCompositeDomain<T>());

    // Get the number of leaves
    size_t nleaves = 0;
    leafIterator it_count = m_tree.beginLeafIterator();
    while (it_count.good())
    {
        nleaves++;
        it_count.next();
    }

    // Get the number of elements per leaf and other info
    std::vector<LeafInfo<T>> all_leaf_infos;
    all_leaf_infos.reserve(nleaves); 
    size_t total_elements = 0;
    index_t leafIdx = 0;
    leafIterator it = m_tree.beginLeafIterator();
    while (it.good())
    {
        size_t nel = 1;
        gsVector<index_t> lower = it.lowerCorner();
        gsVector<index_t> upper = it.upperCorner();

        if (m_basis.manualLevels() )
        {
            for (short_t i = 0; i < d; ++i)
            {
                m_basis._diadicIndexToKnotIndex(it.level(),i,lower[i]);
                m_basis._diadicIndexToKnotIndex(it.level(),i,upper[i]);
                nel *= upper[i] - lower[i];
            }
        }
        else
        {
            nel *= (upper - lower).prod();
        }
        
        total_elements += nel;
        all_leaf_infos.push_back({lower, upper, leafIdx, it.level(), nel});
        
        it.next();
        leafIdx++;
    }
    GISMO_ENSURE(total_elements > npieces || npieces == 0, "Number of pieces exceeds number of elements.");

    // if the number of pieces is smaller than the number of leaves, we merge leaves
    if (nleaves > npieces)
    {
        // Karmarkar Karp merging of patches until we reach npieces
        std::vector<std::vector<index_t>> leafIndices_per_piece(npieces);
        std::vector<size_t> piece_sizes(npieces, 0);
        for (const auto& leaf_info : all_leaf_infos)
        {
            // Find piece with minimum size
            size_t minIdx = 0;
            for (size_t j = 1; j < npieces; ++j) {
                if (piece_sizes[j] < piece_sizes[minIdx]) {
                    minIdx = j;
                }
            }
            leafIndices_per_piece[minIdx].push_back(leaf_info.global_id);
            piece_sizes[minIdx] += leaf_info.num_elements;
        }

        // Now create subdomains from the pieces
        for (const auto& piece_indices : leafIndices_per_piece)
        {
            gsCompositeDomain<T> piece_domain;
            for (const auto& leaf_id : piece_indices)
            {
                const auto& leaf_info = all_leaf_infos[leaf_id];
                // Get the tensor domain for the specific level of this leaf
                // We need to create it manually from the knot vectors at this level
                const auto& tensorBasis = m_basis.tensorLevel(leaf_info.level);
                // Create ranges for the subdomain
                gsVector<index_t> ranges_start(d);
                gsVector<index_t> ranges_end(d);
                for(short_t i = 0; i < d; ++i)
                {
                    // For manual levels, convert diadic indices to knot indices
                    // For non-manual levels, diadic indices are already knot indices
                    index_t start = leaf_info.lower[i];
                    index_t end = leaf_info.upper[i];
                    if (m_basis.manualLevels())
                    {
                        m_basis._diadicIndexToKnotIndex(leaf_info.level, i, start);
                        m_basis._diadicIndexToKnotIndex(leaf_info.level, i, end);
                    }
                    // For non-manual levels, use diadic indices directly as knot indices
                    ranges_start[i] = start;
                    ranges_end[i] = end;
                }
                // Create the tensor domain
                auto parentDomain = std::dynamic_pointer_cast<gsTensorDomain<d,T>>(tensorBasis.domain());
                GISMO_ASSERT(parentDomain, "Leaf domain is not a tensor domain.");
                result->keepAlive(parentDomain);
                auto leaf_domain = memory::make_shared(new gsTensorSubDomain<d,T>(*parentDomain, ranges_start, ranges_end, this->patchId(), parentDomain));
                piece_domain.addDomain(leaf_domain);
            }
            result->addDomain(memory::make_shared(new gsCompositeDomain<T>(piece_domain)));
        }
    }
    else // npieces >= num_leaves => assign leaves to pieces
    {
        // Strategy: Decompose each leaf into a number of pieces such that every piece has roughly `elementsPerPiece` elements.
        // We use the Karmarkar-Karp (KK) algorithm to distribute the pieces among the leaves.
        std::vector<size_t> pieces_per_leave(nleaves, 1);
        size_t current_total_pieces = nleaves;
        while (current_total_pieces < npieces)        
        {
            // Find leaf with maximum number of elements per piece
            size_t maxIdx = 0;
            T max_ratio = -1.0;
            for (size_t i = 0; i < nleaves; ++i)
            {
                T ratio = (T)all_leaf_infos[i].num_elements / pieces_per_leave[i];
                if (ratio > max_ratio)
                {
                    max_ratio = ratio;
                    maxIdx = i;
                }
            }
            pieces_per_leave[maxIdx]++;
            current_total_pieces++;
        }

        // Sum the pieces_per_leave to verify it matches npieces
        size_t verify_total = std::accumulate(pieces_per_leave.begin(), pieces_per_leave.end(), 0);
        // Now decompose each leaf accordingly and add to result
        for (size_t i = 0; i < nleaves; ++i)
        {
            const auto& leaf_info = all_leaf_infos[i];
            const auto& tensorBasis = m_basis.tensorLevel(leaf_info.level);
            // Create ranges for the subdomain
            gsVector<index_t> ranges_start(d);
            gsVector<index_t> ranges_end(d);
            for(short_t j = 0; j < d; ++j)
            {
                // For manual levels, convert diadic indices to knot indices
                // For non-manual levels, diadic indices are already knot indices
                index_t start = leaf_info.lower[j];
                index_t end = leaf_info.upper[j];
                if (m_basis.manualLevels())
                {
                    m_basis._diadicIndexToKnotIndex(leaf_info.level, j, start);
                    m_basis._diadicIndexToKnotIndex(leaf_info.level, j, end);
                }
                // For non-manual levels, use diadic indices directly as knot indices
                ranges_start[j] = start;
                ranges_end[j] = end;
            }
            // Create the tensor domain and decompose it directly
            auto parentDomain = std::dynamic_pointer_cast<gsTensorDomain<d,T>>(tensorBasis.domain());
            GISMO_ASSERT(parentDomain, "Leaf domain is not a tensor domain.");
            result->keepAlive(parentDomain);
            auto decomposed_leaf = gsTensorSubDomain<d,T>::decompose(*parentDomain, ranges_start, ranges_end, pieces_per_leave[i], this->patchId(), parentDomain);
            // Add to result
            if (auto* decomposed_composite = dynamic_cast<gsCompositeDomain<T>*>(decomposed_leaf.get()))
            {
                for (size_t j = 0; j < decomposed_composite->nPieces(); ++j)
                {
                    result->addDomain(decomposed_composite->subdomain(j));
                }
            }
            else
            {
                GISMO_ERROR("Decomposition of leaf did not return a composite domain.");
            }
        }
    }
    return result;
}

} // namespace gismo
