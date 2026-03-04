/** @file gsTrimmedDomainIterator.h

    @brief Iterator over the elements of a trimmed domain

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst, A. Mantzaflaris
*/

#pragma once

#include <gsDomain/gsDomainIterator.h>
#include <gsDomain/gsTrimmedDomain.h>

namespace gismo
{
/**
 * @brief Iterator over the elements of a trimmed domain.
 *
 * Iterates over leaves of the KdTree whose sign satisfies SignOp.
 * For each leaf at level L spanning [lc,uc) in level-0 element coords,
 * iterates over all level-L sub-elements within that range using
 * per-leaf break vectors built from m_breaks[L][dim].
 *
 * This follows the same pattern as gsHDomainIterator: each leaf may
 * contain multiple sub-elements, and the iterator steps through them
 * lexicographically.
 *
 * \ingroup Domain
 */

template<typename SignOp, short_t d, class T, class Z /*= size_t*/>
class gsTrimmedDomainIterator : public gsDomainIterator<T>
{
public:

    typedef gsTrimmedDomain<d,T,Z> Domain_t;
    typedef typename Domain_t::TData_t TData_t;
    typedef typename Domain_t::Tree_t Tree_t;
    typedef typename Tree_t::const_literator leafIterator;
    typedef typename Tree_t::point_t point_t;

    typedef typename std::vector<T>::const_iterator uiter;

    typedef typename gsDomainIterator<T>::uPtr domainIter;

private:
    const Domain_t & m_tdomain;

    leafIterator m_leaf;

    // Per-leaf break vectors (built in updateLeaf())
    std::vector< std::vector<T> > m_leafBreaks;

    // Iterators into m_leafBreaks for current element position
    gsVector<uiter, d> m_curElement;
    gsVector<uiter, d> m_meshStart;
    gsVector<uiter, d> m_meshEnd;

    // Parametric bounds of the current element
    gsVector<T> m_lower;
    gsVector<T> m_upper;

public:

    explicit gsTrimmedDomainIterator(const Domain_t & _tdomain)
    :
    gsDomainIterator<T>(),
    m_tdomain(_tdomain),
    m_leafBreaks(d),
    m_lower(d),
    m_upper(d)
    {
        reset();
    }

    gsTrimmedDomainIterator(const gsTrimmedDomainIterator & other) = default;
    domainIter clone() const override { return domainIter(new gsTrimmedDomainIterator(*this)); }

    void next() override
    {
        bool isGood = nextLexicographic(m_curElement, m_meshStart, m_meshEnd);
        if (isGood)
            updateCorners();
        else
            nextLeaf();
    }

    void next(index_t increment) override
    {
        bool isGood(m_leaf.good());
        for (index_t i = 0; i < increment && isGood; ++i)
        {
            isGood = nextLexicographic(m_curElement, m_meshStart, m_meshEnd);
            if (!isGood)
                isGood = nextLeaf();
        }
        if (isGood)
            updateCorners();
    }

    void reset() override
    {
        m_leaf = m_tdomain.tree().beginLeafIterator();
        if (m_leaf.good() && !SignOp()(m_leaf.data().sign()))
            nextLeaf();
        else
            updateLeaf();
    }

    gsVector<T> lowerCorner() const override { return m_lower; }

    gsVector<T> upperCorner() const override { return m_upper; }

    short_t sign() const { return m_leaf.data().sign(); }

private:

    /// Advances to the next leaf satisfying SignOp; calls updateLeaf() on success.
    bool nextLeaf()
    {
        bool isGood = m_leaf.next();
        while (isGood)
        {
            if (SignOp()(m_leaf.data().sign()))
            {
                updateLeaf();
                return true;
            }
            isGood = m_leaf.next();
        }
        return false;
    }

    /// Builds per-leaf break vectors from the domain's m_breaks at the leaf's level,
    /// covering the parametric range of the leaf. Resets the element cursor.
    void updateLeaf()
    {
        if (!m_leaf.good()) return;

        const point_t & lc = m_leaf.data().lowerCorner();
        const point_t & uc = m_leaf.data().upperCorner();
        const index_t L = m_leaf.data().level();

        const unsigned nLevels = m_tdomain.numLevels();
        const unsigned usedLevel = (static_cast<unsigned>(L) < nLevels)
                                   ? static_cast<unsigned>(L) : nLevels - 1;

        // With level-L leaf coordinates, they directly index the level-L break vector.
        // No need to multiply by subdiv.
        for (short_t j = 0; j < d; ++j)
        {
            const std::vector<T> & bk = m_tdomain.breaks(usedLevel, j);
            const index_t start = static_cast<index_t>(lc[j]);
            const index_t end   = static_cast<index_t>(uc[j]);

            m_leafBreaks[j].clear();
            for (index_t idx = start; idx <= end; ++idx)
                m_leafBreaks[j].push_back(bk[idx]);

            m_curElement[j] = m_meshStart[j] = m_leafBreaks[j].begin();
            m_meshEnd[j]    = m_leafBreaks[j].end() - 1;
        }

        updateCorners();
    }

    /// Updates m_lower and m_upper from the current element iterators.
    void updateCorners()
    {
        for (short_t j = 0; j < d; ++j)
        {
            m_lower[j] = *m_curElement[j];
            m_upper[j] = *(m_curElement[j] + 1);
        }
    }
};

} // namespace gismo
