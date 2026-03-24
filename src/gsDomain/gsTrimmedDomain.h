/** @file gsTrimmedDomain.h

    @brief TODO

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst, A. Mantzaflaris
*/

#pragma once

#include <gsDomain/gsDomain.h>
#include <gsDomain/gsCutTreeData.h>
#include <gsAssembler/gsQuadrature.h>
#include <gsHSplines/gsHTensorBasis.h>

namespace gismo
{

// Cells with sign 0
struct BoundarySign { bool operator()(short_t s) const { return s==0; } };
// Cells with sign -1
struct InteriorSign { bool operator()(short_t s) const { return s< 0; } };
// Cells with sign  1
struct ExteriorSign { bool operator()(short_t s) const { return s> 0; } };
// Cells with signs 0 or -1
struct AllSign      { bool operator()(short_t s) const { return s<=0; } };
// Cells with any sign
struct AnySign      { bool operator()(short_t s) const { return true; } };

/**
 * @brief TODO
 *
 * \ingroup Domain
 */

template<short_t d, class T, class Z=size_t>
class gsTrimmedDomain : public gsDomain<T>
{
public:
    typedef gsCutTreeData<d,Z> TData_t;
    typedef gsKdTree<d,Z,TData_t > Tree_t;
    typedef gsDomainIteratorWrapper<T> domainIter;
    typedef typename Tree_t::const_literator leafIterator;
    typedef typename Tree_t::point_t point_t;

    template <typename SignOp, short_t _d, class _T, class _Z>
    friend class gsTrimmedDomainIterator;

protected:
    Tree_t m_tree; ///< The KdTree partition characterizing in/out/cut elements
    /// m_breaks[level][dim] = break vector at that level for that direction.
    /// Level 0 is the coarsest (background) grid; higher levels are dyadic refinements.
    /// Leaf coordinates (lowerCorner, upperCorner) are stored as indices into m_breaks[L][j]
    /// for their respective level L, making the domain self-contained and independent of any basis.
    std::vector< std::vector< std::vector<T> > > m_breaks;

    short_t m_deg; // For now: the degree is here: we need somehow to get the quadrature order needed for integrals...

public: // virtual interface

    virtual gsMatrix<T> boundingBox() const override
    {
        GISMO_NO_IMPLEMENTATION
    }

    /// Computes the sign at given points \a u in the domain
    virtual gsVector<short_t> sign(const gsMatrix<T> & u)
    {
        GISMO_NO_IMPLEMENTATION
    }

public: // non-virtual interface

    inline std::vector<bool> inDomain(const gsMatrix<T> & u)
    {
        std::vector<bool> res(u.cols());
        gsVector<short_t> val = this->sign(u);
        bool * r = res.data();
        for (short_t * a = val.data(); a != val.data()+val.size(); ++a)
            *(r++) = (*a<0 ? false : true);
        return res;
    }

    inline std::vector<bool> onBoundary(const gsMatrix<T> & u)
    {
        std::vector<bool> res(u.cols());
        gsVector<short_t> val = this->sign(u);
        bool * r = res.data();
        for (short_t * a = val.data(); a != val.data()+val.size(); ++a)
            *(r++) = (*a!=0 ? false : true);
        return res;
    }

    domainIter beginAll() const override
    {
        return begin<AllSign>();
    }

    domainIter beginBdr(const boxSide bs) const override
    {
        return begin<BoundarySign>();
    }

    domainIter beginActive() const
    {
        return beginAll();
    }

    domainIter beginExterior() const
    {
        return begin<ExteriorSign>();
    }

    domainIter beginInterior() const
    {
        return begin<InteriorSign>();
    }

    domainIter beginAny() const
    {
        return begin<AnySign>();
    }

    template<typename SignOp>
    domainIter begin() const
    {
        return domainIter(new gsTrimmedDomainIterator<SignOp,d,T,Z>(*this));
    }

    template<typename SignOp>
    domainIter end() const
    {
        return domainIter(new gsDomainIteratorEnd<T>(this->template numElements<SignOp>()));
    }

    /// Returns the number of elements whose sign satisfies SignOp.
    /// Each leaf at level L spanning [lc,uc) in level-0 coords contains
    /// (uc-lc).prod() * 2^(L*d) sub-elements (one per level-L break cell).
    template<typename SignOp>
    size_t numElements() const
    {
        leafIterator it = m_tree.beginLeafIterator();
        size_t nel(0);
        while (it.good())
        {
            if (SignOp()(it.data().sign()))
            {
                const point_t & uc = it.data().upperCorner();
                const point_t & lc = it.data().lowerCorner();
                // With level-L leaf coordinates, element count is simply the span.
                // No need to multiply by 2^L since coordinates index the level-L break vector.
                size_t leafNel = 1;
                for (short_t j = 0; j < d; ++j)
                    leafNel *= static_cast<size_t>(uc[j] - lc[j]);
                nel += leafNel;
            }
            it.next();
        }
        return nel;
    }

    size_t numElements() const override
    {
        return numElements<AllSign>();
    }

    size_t numElementsBdr(boxSide const & s = boundary::none) const override
    {
        return numElements<BoundarySign>();
    }

    short_t dim() const override { return d; }

    short_t degree(short_t i = 0) const
    {GISMO_UNUSED(i); return m_deg;}

    const Tree_t & tree() const { return m_tree; }

    /// Returns break vector at given \a level for direction \a j.
    const std::vector<T> & breaks(unsigned level, short_t j) const
    {
        GISMO_ASSERT(level < m_breaks.size(), "Level " << level << " out of range");
        return m_breaks[level][j];
    }

    unsigned numLevels() const { return static_cast<unsigned>(m_breaks.size()); }

protected:

    /// Dyadically subdivides breaks up to and including \a level.
    /// m_breaks[0] must already be populated; each subsequent level
    /// inserts midpoints between every consecutive pair of parent breaks.
    void _ensureLevel(unsigned level)
    {
        while (m_breaks.size() <= static_cast<size_t>(level))
        {
            const unsigned parentL = static_cast<unsigned>(m_breaks.size()) - 1;
            m_breaks.push_back(std::vector<std::vector<T>>(d));
            for (short_t j = 0; j < d; ++j)
            {
                const std::vector<T> & parentBk = m_breaks[parentL][j];
                std::vector<T> & childBk = m_breaks.back()[j];
                childBk.reserve(parentBk.size() * 2 - 1);
                for (size_t i = 0; i + 1 < parentBk.size(); ++i)
                {
                    childBk.push_back(parentBk[i]);
                    childBk.push_back((parentBk[i] + parentBk[i+1]) / T(2));
                }
                childBk.push_back(parentBk.back());
            }
        }
    }

    /// Returns the parametric lower and upper corners of \a node from m_breaks.
    void _getParamCorners(const Tree_t * node,
                          gsVector<T,d> & u1, gsVector<T,d> & u2) const
    {
        const point_t & lc = node->nodeData().lowerCorner();
        const point_t & uc = node->nodeData().upperCorner();
        const index_t   L  = node->nodeData().level();
        for (short_t j = 0; j < d; ++j)
        {
            u1[j] = m_breaks[L][j][static_cast<index_t>(lc[j])];
            u2[j] = m_breaks[L][j][static_cast<index_t>(uc[j])];
        }
    }

    /// Evaluates sign at quadrature points within \a node's parametric box and
    /// assigns the result to the node's sign field. Uses \a QuRule for quadrature.
    void _classifyLeaf(Tree_t * node, gsLobattoRule<T> & QuRule)
    {
        gsVector<T,d>  u1, u2;
        gsVector<T>    wts;
        gsMatrix<T>    pts;
        _getParamCorners(node, u1, u2);
        QuRule.mapTo(u1, u2, pts, wts);
        gsVector<short_t> sgn = this->sign(pts);
        node->nodeData().sign() =
            (sgn.array() ==  1).all() ?  1 :
            (sgn.array() == -1).all() ? -1 : 0;
    }

    /// Classifies all leaves of m_tree using m_breaks[0] for index->parametric
    /// conversion. Mixed leaves spanning more than one level-0 element are split
    /// recursively; single-element mixed leaves are marked as boundary (sign=0).
    void _classifyTree(index_t samples)
    {
        gsVector<index_t> numNodes(d);
        numNodes.setConstant(samples);
        gsLobattoRule<T> QuRule(numNodes);

        std::vector<Tree_t*> stack;
        stack.reserve(std::pow(4, d));
        stack.push_back(&m_tree);

        while (!stack.empty())
        {
            Tree_t * cur = stack.back();
            stack.pop_back();

            if (!cur->isLeaf())
            {
                stack.push_back(cur->left);
                stack.push_back(cur->right);
                continue;
            }

            _classifyLeaf(cur, QuRule);

            if (cur->nodeData().sign() == 0)
            {
                const point_t & lc = cur->nodeData().lowerCorner();
                const point_t & uc = cur->nodeData().upperCorner();
                if ((uc - lc).prod() > 1)
                {
                    cur->anyMidSplit(1);
                    stack.push_back(cur->left);
                    stack.push_back(cur->right);
                }
            }
        }
    }

    /// Adaptive refinement of the cut-tree driven by a stopping predicate.
    ///
    /// Traverses all leaves. For each leaf, evaluates the sign at \a samples
    /// quadrature points per direction. Assigns sign (-1/0/1). If the predicate
    /// \a shouldStop returns false, splits the leaf and promotes it to the next
    /// dyadic level (doubling the index coords and incrementing level by 1).
    ///
    /// \a shouldStop is called as shouldStop(sign, u1, u2) and should return
    /// true when no further refinement is needed for this leaf.
    template<typename StopPred>
    void _classifyTreeAdaptive(index_t samples, StopPred shouldStop)
    {
        gsVector<index_t> numNodes(d);
        numNodes.setConstant(samples);
        gsLobattoRule<T> QuRule(numNodes);

        std::vector<Tree_t*> stack;
        stack.reserve(std::pow(4, d));
        stack.push_back(&m_tree);

        while (!stack.empty())
        {
            Tree_t * cur = stack.back();
            stack.pop_back();

            if (!cur->isLeaf())
            {
                stack.push_back(cur->left);
                stack.push_back(cur->right);
                continue;
            }

            const index_t L = cur->nodeData().level();
            _ensureLevel(static_cast<unsigned>(L));

            gsVector<T,d> u1, u2;
            _getParamCorners(cur, u1, u2);
            _classifyLeaf(cur, QuRule);
            const short_t leafSign = cur->nodeData().sign();

            if (!shouldStop(leafSign, u1, u2))
            {
                const index_t newL = L + 1;
                _ensureLevel(static_cast<unsigned>(newL));

                point_t & lc = cur->nodeData().lowerCorner();
                point_t & uc = cur->nodeData().upperCorner();
                for (short_t j = 0; j < d; ++j)
                {
                    lc[j] = static_cast<Z>(lc[j]) * 2;
                    uc[j] = static_cast<Z>(uc[j]) * 2;
                }
                cur->nodeData().level() = newL;

                short_t splitAxis = 0;
                Z maxSpan = 0;
                for (short_t j = 0; j < d; ++j)
                {
                    Z span = uc[j] - lc[j];
                    if (span > maxSpan) { maxSpan = span; splitAxis = j; }
                }
                Z splitPos = lc[splitAxis] + (uc[splitAxis] - lc[splitAxis]) / 2;
                cur->split(splitAxis, splitPos);
                stack.push_back(cur->left);
                stack.push_back(cur->right);
            }
        }
    }

    /// Size-based adaptive init. Refines until each leaf has volume <= maxElementSize,
    /// or is a cut cell with volume <= minElementSize. Leaves are tracked at their
    /// current dyadic level; m_breaks is extended on demand via _ensureLevel().
    void init(T maxElementSize = T(1), T minElementSize = T(0.1), index_t samples = 10)
    {
        gsMatrix<T> bb = this->boundingBox();

        m_breaks.clear();
        m_breaks.resize(1);
        m_breaks[0].resize(d);
        for (short_t j = 0; j < d; ++j)
        {
            m_breaks[0][j].push_back(static_cast<T>(bb(0, j)));
            m_breaks[0][j].push_back(static_cast<T>(bb(1, j)));
        }

        point_t rootUpper;
        for (short_t j = 0; j < d; ++j)
            rootUpper[j] = static_cast<Z>(1);

        m_tree = Tree_t(TData_t(point_t::Zero(), rootUpper, 0));

        _classifyTreeAdaptive(samples, [maxElementSize, minElementSize]
            (short_t leafSign, const gsVector<T,d> & u1, const gsVector<T,d> & u2) -> bool
        {
            const T volume = (u2 - u1).prod();
            if (leafSign != 0)
                return volume <= maxElementSize;
            else
                return volume <= minElementSize;
        });
    }

    /// Initializes from a bounding box + uniform cell count per direction.
    void init(const gsMatrix<T>         & bbox,
              const gsVector<index_t,d> & numCells,
              index_t samples = 5)
    {
        m_breaks.clear();
        m_breaks.resize(1);
        m_breaks[0].resize(d);
        for (short_t j = 0; j < d; ++j)
        {
            const T lo = bbox(0, j);
            const T hi = bbox(1, j);
            gsKnotVector<T> kv(lo, hi, numCells[j] - 1, 2);
            m_breaks[0][j] = kv.breaks();
        }

        point_t upperIndex;
        for (short_t j = 0; j < d; ++j)
            upperIndex[j] = static_cast<Z>(numCells[j]);

        m_tree = Tree_t(TData_t(point_t::Zero(), upperIndex, 0));
        _classifyTree(samples);
    }

    /// Initializes from a tensor B-spline basis (uses its knot vectors as the level-0 grid).
    void init(const gsTensorBSplineBasis<d,T> & tbasis, index_t samples = 5)
    {
        m_deg = tbasis.maxDegree();
        m_breaks.clear();
        m_breaks.resize(1);
        m_breaks[0].resize(d);
        for (short_t j = 0; j < d; ++j)
            m_breaks[0][j] = tbasis.knots(j).breaks();

        point_t upperIndex;
        for (short_t j = 0; j < d; ++j)
            upperIndex[j] = static_cast<Z>(tbasis.knots(j).numElements());

        m_tree = Tree_t(TData_t(point_t::Zero(), upperIndex, 0));
        _classifyTree(samples);
    }

    /// Mirrors the HTensorBasis tree structure into m_tree.
    /// Populates m_breaks with per-level break vectors from the HTB's knot vectors.
    /// Sets leaf geometry (lowerCorner, upperCorner, level) to level-L element indices.
    /// Does NOT classify leaves (no sign assignment).
    ///
    /// Leaf coordinates are stored as level-L element indices (indexing into m_breaks[L][j]).
    /// This makes the trimmed domain self-contained and independent of the basis.
    /// Internal node split positions are stored in level-0 element indices for consistency
    /// with the paired traversal logic, though these are not used after construction.
    ///
    /// \param htbasis The hierarchical tensor basis to mirror.
    void _mirrorHTBTree(const gsHTensorBasis<d,T> & htbasis)
    {
        typedef gsHTree<d,index_t>                         HTree_t;
        typedef gsKdTree<d,index_t,gsHTreeData<d,index_t>> HNode_t;
        typedef gsVector<index_t,d>                        HPoint_t;

        const HTree_t & htree  = htbasis.tree();
        const unsigned  idxLv  = htree.getIndexLevel();
        const unsigned  maxLvl = htbasis.maxLevel();

        // Extract break vectors for all levels from the HTB knot vectors
        m_breaks.clear();
        m_breaks.resize(maxLvl + 1);
        for (unsigned L = 0; L <= maxLvl; ++L)
        {
            const gsTensorBSplineBasis<d,T> * tbL =
                dynamic_cast<const gsTensorBSplineBasis<d,T>*>(htbasis.getBases()[L]);
            GISMO_ASSERT(tbL != nullptr,
                         "gsHTensorBasis level-" << L << " basis is not a gsTensorBSplineBasis");
            m_breaks[L].resize(d);
            for (short_t j = 0; j < d; ++j)
                m_breaks[L][j] = tbL->knots(j).breaks();
        }

        // Root of cut-tree covers same level-0 element range as the HTB.
        // (Temporarily stored; will be overwritten when leaves are populated.)
        const point_t upp0 = htree.upperCornerIndex().template cast<Z>();
        m_tree = Tree_t(TData_t(point_t::Zero(), upp0, 0));

        // Paired traversal: mirror HTB tree structure into cut-tree
        typedef std::pair<const HNode_t*, Tree_t*> NodePair;
        std::vector<NodePair> stack;
        stack.reserve(128);
        stack.emplace_back(&htree, &m_tree);

        while (!stack.empty())
        {
            auto [htNode, cutNode] = stack.back();
            stack.pop_back();

            if (htNode->isLeaf())
            {
                const int L = htNode->nodeData().level();
                // Stride: how many global-index units per one level-L unit
                const index_t h = static_cast<index_t>(1) << (idxLv - L);

                const HPoint_t & gk1 = htNode->nodeData().lowerCorner();
                const HPoint_t & gk2 = htNode->nodeData().upperCorner();

                // Convert global-index-level coords to level-L local indices.
                // Leaf coordinates are stored as level-L element indices.
                point_t lcL, ucL;
                for (short_t j = 0; j < d; ++j)
                {
                    const index_t k1L = static_cast<index_t>(gk1[j]) / h;
                    const index_t k2L = static_cast<index_t>(gk2[j]) / h;
                    // Clamp to valid range within level-L break vector
                    const index_t nElemL = static_cast<index_t>(m_breaks[L][j].size()) - 1;
                    lcL[j] = static_cast<Z>(k1L);
                    ucL[j] = static_cast<Z>(std::min(k2L, nElemL));
                }

                cutNode->nodeData().lowerCorner() = lcL;
                cutNode->nodeData().upperCorner() = ucL;
                cutNode->nodeData().level() = static_cast<index_t>(L);
            }
            else
            {
                // Convert the HTB split position (global-index-level) to level-0 coords
                const int splitAxis  = htNode->axis;
                const index_t gPos   = static_cast<index_t>(htNode->pos);

                // We need to figure out the level of the children to compute the
                // level-0 split position. Use the left child's level as reference.
                // For internal nodes, we simply split the parent's level-0 range
                // at the mid-point corresponding to the HTB split.
                const point_t & parentLc = cutNode->nodeData().lowerCorner();
                const point_t & parentUc = cutNode->nodeData().upperCorner();

                // Determine the level this split operates at.
                // The HTB internal node splits at a position expressed in global-index coords.
                // The left child spans [parentGk1, gPos], right spans [gPos, parentGk2].
                // We compute the level-0 element index at gPos using the coarsest applicable level.
                // Use level of the cut-node (inherited from its construction) as a guide.
                const index_t parentL = cutNode->nodeData().level();
                const index_t h = static_cast<index_t>(1) << (idxLv - static_cast<unsigned>(parentL));
                const index_t kL  = gPos / h;              // level-parentL index
                const index_t c   = kL >> parentL;          // level-0 element index

                Z cutPos = static_cast<Z>(c);
                // Clamp to a valid split (must be strictly between parent corners)
                if (cutPos <= parentLc[splitAxis]) cutPos = parentLc[splitAxis] + static_cast<Z>(1);
                if (cutPos >= parentUc[splitAxis]) cutPos = parentUc[splitAxis] - static_cast<Z>(1);

                cutNode->split(static_cast<int>(splitAxis), cutPos);

                stack.emplace_back(htNode->left,  cutNode->left);
                stack.emplace_back(htNode->right, cutNode->right);
            }
        }
    }

    /// Initializes from a hierarchical tensor basis.
    ///
    /// First mirrors the HTensorBasis tree topology via _mirrorHTBTree(), then
    /// classifies all leaves via _classifyTreeAdaptive() with no refinement
    /// (since HTB structure is already fully specified).
    ///
    /// Leaf coordinates are stored as level-L element indices (indexing into m_breaks[L][j]),
    /// making the domain self-contained and independent of the original basis.
    ///
    /// \param htbasis The hierarchical tensor basis to initialize from.
    /// \param samples Number of quadrature points per direction for sign classification (default: 5).
    void init(const gsHTensorBasis<d,T> & htbasis,
              index_t samples = 5)
    {
        _mirrorHTBTree(htbasis);
        _classifyTreeAdaptive(samples, [](short_t sign, const gsVector<T,d>& u1, const gsVector<T,d>& u2) { return true; });
    }
};

} // namespace gismo
