/** @file gsHTree.h

    @brief Provides the class gsHTree, the hierarchical box tree of
    hierarchical spline spaces, built on gsKdTree.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): G. Kiss, A. Mantzaflaris
*/

#pragma once

#include <algorithm>
#include <limits>
#include <list>
#include <queue>
#include <stack>
#include <vector>

#include <gsCore/gsBoundary.h>
#include <gsDomain/gsKdTree.h>
#include <gsHSplines/gsAAPolyline.h>
#include <gsHSplines/gsHTreeData.h>
#include <gsHSplines/gsVSegment.h>

namespace gismo {

template<short_t d, class Z = index_t> class gsHTree;

namespace internal {
    /// True iff all leaves hit by the query box are of level \a level
    struct query1_visitor
    {
        typedef bool return_type;

        // initialize result as true
        static return_type init() {return true;}

        template<short_t d, class Z>
        static void visitLeaf(const gismo::gsKdTree<d,Z,gismo::gsHTreeData<d,Z> > * leafNode, int level, return_type & res)
        {
            if ( leafNode->nodeData().level() != level )
                res = false;
        }
    };

    /// True iff all leaves hit by the query box are of level above \a level
    struct query2_visitor
    {
        typedef bool return_type;

        // initialize result as true
        static return_type init() {return true;}

        template<short_t d, class Z>
        static void visitLeaf(const gismo::gsKdTree<d,Z,gismo::gsHTreeData<d,Z> > * leafNode, int level, return_type & res)
        {
            if ( leafNode->nodeData().level() <= level )
                res = false;
        }
    };

    /// Minimum level of the leaves hit by the query box
    struct query3_visitor
    {
        typedef int return_type;

        // initialize result as a max possible value, since we are looking
        // for a minimum
        static return_type init() {return 1000000;}

        template<short_t d, class Z>
        static void visitLeaf(const gismo::gsKdTree<d,Z,gismo::gsHTreeData<d,Z> > * leafNode, int , return_type & res)
        {
            if ( leafNode->nodeData().level() < res )
                res = leafNode->nodeData().level();
        }
    };

    /// Maximum level of the leaves hit by the query box
    struct query4_visitor
    {
        typedef int return_type;

        // initialize result as a minimum possible value, since we are
        // looking for a maximum
        static return_type init() {return -1;}

        template<short_t d, class Z>
        static void visitLeaf(const gismo::gsKdTree<d,Z,gismo::gsHTreeData<d,Z> > * leafNode , int , return_type & res)
        {
            if ( leafNode->nodeData().level() > res )
                res = leafNode->nodeData().level();
        }
    };

    /// Decreases the level by 1 for all leaves
    struct levelDown_visitor
    {
        typedef int return_type;
        static return_type init() {return 0;}

        template<short_t d, class Z>
        static void visitLeaf(gismo::gsKdTree<d,Z,gismo::gsHTreeData<d,Z> > * leafNode, return_type &)
        {
            leafNode->nodeData().level()--;
        }
    };

    struct maxLevel_visitor
    {
        typedef int return_type;
        static return_type init() {return 0;}

        template<short_t d, class Z>
        static void visitLeaf(gismo::gsKdTree<d,Z,gismo::gsHTreeData<d,Z> > * leafNode, return_type &i)
        {
            if (leafNode->nodeData().level()>i) i=leafNode->nodeData().level();
        }
    };

    /// Multiplies by 2 the split positions and the leaf data of every node
    struct liftCoordsOneLevel_visitor
    {
        typedef int return_type;
        static return_type init() {return 0;}

        template<short_t d, class Z>
        static void visitNode(gismo::gsKdTree<d,Z,gismo::gsHTreeData<d,Z> > * n, return_type &)
        {
            n->multiplyByTwo();
        }
    };

    /// Divides by 2 the split positions and the leaf data of every node
    struct reduceCoordsOneLevel_visitor
    {
        typedef int return_type;
        static return_type init() {return 0;}

        template<short_t d, class Z>
        static void visitNode(gismo::gsKdTree<d,Z,gismo::gsHTreeData<d,Z> > * n, return_type &)
        {
            n->divideByTwo();
        }
    };

} // namespace internal

/// Leaf iterator of gsHTree: corners are returned in the level-local
/// index space of the leaf's level; data() exposes the raw leaf data,
/// whose corners are global indices at the tree's index level.
template<short_t d, class Z, bool isconst = false>
class gsHTreeLeafIter : public gsKdTreeLeafIter<gsKdTree<d,Z,gsHTreeData<d,Z> >, isconst>
{
    typedef gsKdTreeLeafIter<gsKdTree<d,Z,gsHTreeData<d,Z> >, isconst> Base;
public:
    typedef gsVector<Z,d> point;

    gsHTreeLeafIter() : Base(), m_index_level(0) { }

    gsHTreeLeafIter(const gsKdTree<d,Z,gsHTreeData<d,Z> > * root, unsigned index_level)
    : Base(root->axis == -2 ? Base() : Base(root)), m_index_level(index_level) { }

    /// Level of the current leaf
    int level() const { return this->data().level(); }

    /// Lower corner of the current leaf, in the indices of its own level
    point lowerCorner() const
    {
        point p = this->data().lowerCorner();
        const unsigned s = m_index_level - this->data().level();
        for (short_t i = 0; i != d; ++i)
            p[i] >>= s;
        return p;
    }

    /// Upper corner of the current leaf, in the indices of its own level
    point upperCorner() const
    {
        point p = this->data().upperCorner();
        const unsigned s = m_index_level - this->data().level();
        for (short_t i = 0; i != d; ++i)
            p[i] >>= s;
        return p;
    }

    /// The index level of the tree
    unsigned indexLevel() const { return m_index_level; }

    /// True iff the current leaf is aligned to the grid of its level
    bool isAligned() const { return this->data().isAligned(m_index_level); }

private:
    unsigned m_index_level;
};

/**
    @brief The hierarchical box tree: a gsKdTree whose leaves carry
    gsHTreeData, i.e., a box and a level.

    This is the root node of the tree. Its leaves partition the box
    [0, upperCorner()] of the parameter domain into the cells of the
    hierarchical domain, each with its level. The boxes inside the tree
    are stored in global indices on the finest grid (the index level, see
    getIndexLevel()). The input of insertBox(), sinkBox(), clearBox() and the
    query functions is given in level-local indices, and the leaf
    iterator returns level-local indices.

    \tparam d the dimension
    \tparam Z the coordinate index type

    \ingroup HSplines
*/
template<short_t d, class Z>
class gsHTree : public gsKdTree<d,Z,gsHTreeData<d,Z> >
{
public:
    typedef gsKdTree<d,Z,gsHTreeData<d,Z> > node;
    typedef typename node::data_t data_t;

    typedef typename node::point_t point; // it's a gsVector<index_t,d>

    typedef typename data_t::kdBox box; // it's a gsAABB<d,Z>

    typedef gsHTreeLeafIter<d,Z,false> literator;
    typedef gsHTreeLeafIter<d,Z,true>  const_literator;

private:

    /// Keeps the highest upper indices (at level gsHTree::m_indexLevel)
    point m_upperIndex;

    /// The level of the box representation (global indices)
    unsigned m_indexLevel;

    /// Maximum level present in the tree
    unsigned m_maxInsLevel;

    unsigned m_maxPath;

public:
    gsHTree() : m_indexLevel(0)
    {
        m_maxInsLevel = 0;
        m_maxPath = 0;
    }

    gsHTree(point const & upp)
    {
        m_maxInsLevel = 0;
        m_maxPath = 0;
        init(upp);
    }

    /// Copy constructor (makes a deep copy)
    gsHTree(const gsHTree & o) : node(o),
        m_upperIndex(o.m_upperIndex),
        m_indexLevel(o.m_indexLevel),
        m_maxInsLevel(o.m_maxInsLevel),
        m_maxPath(o.m_maxPath)
    { }

    /// Assignment operator (makes a deep copy)
    gsHTree& operator=( const gsHTree & o)
    {
        if ( this != &o )
        {
            node::operator=(o);
            m_upperIndex  = o.m_upperIndex;
            m_indexLevel  = o.m_indexLevel;
            m_maxInsLevel = o.m_maxInsLevel;
            m_maxPath    = o.m_maxPath;
        }
        return *this;
    }

    /// Move constructor. The moved-from tree is left as an empty node:
    /// it has no leaves, so it may be traversed (leaf iteration,
    /// leafSize(), size(), getBoxes*()), copied, destroyed, assigned to or
    /// re-initialised with init(); it must not be queried (query1-4,
    /// queryLevelCell, levelOf) or modified (insertBox, sinkBox, clearBox)
    /// before init().
    gsHTree(gsHTree&& o) noexcept
    : node(give(static_cast<node&>(o))),
      m_upperIndex(o.m_upperIndex),
      m_indexLevel(o.m_indexLevel),
      m_maxInsLevel(o.m_maxInsLevel),
      m_maxPath(o.m_maxPath)
    { }

    /// Move assignment. The moved-from tree is left as an empty node;
    /// see the move constructor for what it may be used for.
    gsHTree& operator=( gsHTree&& o) noexcept
    {
        if ( this != &o )
        {
            node::operator=(give(static_cast<node&>(o)));
            m_upperIndex  = o.m_upperIndex;
            m_indexLevel  = o.m_indexLevel;
            m_maxInsLevel = o.m_maxInsLevel;
            m_maxPath     = o.m_maxPath;
        }
        return *this;
    }

    /// Initialize the tree
    void init(point const & upp, unsigned index_level)
    {
        m_indexLevel = index_level;
        m_maxInsLevel = 0;

        for (short_t i=0; i<d; ++i)
            m_upperIndex[i] = (upp[i]<< m_indexLevel);

        m_maxPath = 1;

        // Reset this root to a single leaf covering [0, m_upperIndex]
        delete this->data;
        delete this->left;
        delete this->right;
        this->data = new gsHTreeData<d,Z>(m_upperIndex);
        this->axis = -1;
        this->parent = this->left = this->right = nullptr;
    }

    /// Initialize the tree with computing the index_level.
    void init(point const & upp)
    {
        // Idea: find index_level so that for each i, upp[i] * (2 ^ index_level) does not overflow in Z.
        // See issue #672 for more details.

        // backwards compatibility
        Z oldMax = 13;

        Z numMax = std::numeric_limits<Z>::max();

        std::vector<Z> logUpps(d);
        for(short_t i=0; i<d; i++)
        {
            // prevent division by zero
            if(upp[i] == 1)
                logUpps[i] = oldMax;
            else
            {
                // floor of log_2 (numMax / upp[i]):
                logUpps[i] = math::floor( (math::log(numMax) - math::log(upp[i])) / math::log(2) );
            }
        }

        // If the computed number would be too big we take 13 as we used to.
        init(upp, std::min( *std::min_element(logUpps.begin(), logUpps.end()), oldMax) );
    }

    /// Destructor deletes the whole tree
    ~gsHTree() { }

    /// Clones the object
    gsHTree * clone() const { return new gsHTree(*this); }

    /// Leaf iterator with level-local corners, see gsHTreeLeafIter
    literator beginLeafIterator()
    { return literator(this, m_indexLevel); }

    const_literator beginLeafIterator() const
    { return const_literator(this, m_indexLevel); }

    /// Returns the upper corner of the tree, in global indices on the index level
    const point & upperCorner() const { return m_upperIndex; }

    /// Multiplies all indices of the tree by 2
    void multiplyByTwo()
    {
        m_upperIndex *= 2;
        this->template nodeSearch< internal::liftCoordsOneLevel_visitor >();
    }

    /// Divides all indices of the tree by 2
    void divideByTwo()
    {
        m_upperIndex /= 2;
        this->template nodeSearch< internal::reduceCoordsOneLevel_visitor >();
    }

    /// Returns the number of distinct knots in direction \a k of level \a lvl
    int numBreaks(int lvl, int k) const
    {
        return m_upperIndex[k] >> (m_indexLevel - lvl);
    }

    inline unsigned getIndexLevel() const
    {
        return m_indexLevel;
    }

    inline unsigned getMaxInsLevel() const
    {
        return m_maxInsLevel;
    }

    void computeMaxInsLevel()
    {
        m_maxInsLevel = this->template leafSearch< internal::maxLevel_visitor >();
    }

    inline void decrementLevel()
    {
        m_maxInsLevel--;
        this->template leafSearch< internal::levelDown_visitor >();
    }

    /** \brief Merge two trees into their common refinement.
    *
    * Both trees must describe the same domain and use the same index level.
    * The returned tree contains the union of the leaf boxes from both inputs,
    * compressed to remove siblings that end up at the same level.
    */
    static gsHTree merge(const gsHTree & tree1, const gsHTree & tree2);

    /** \brief Returns the boxes which make up the hierarchical domain
    * and the respective levels.
    *
    * Returns a list of boxes defined by left-bottom (b1) and
    * right-top (b2) corners for the splitting in B-spline patches
    * together with the corresponding levelOf()
    *
    * The boxes \em b1 and \em b2 are given as matrices
    * of size <em>n</em> x <em>d</em>, where \em d is the dimension
    * of the domain, and where \em n is the number of boxes of
    * the gsHTree.\n
    *
    * The numbers in \em b1 and \em b2 are given as
    * unique knot indices of gsHTree::m_maxInsLevel
    *
    * \param[out] b1 <em>n</em> x <em>d</em>-matrix, left bottom corners of boxes
    * \param[out] b2 <em>n</em> x <em>d</em>-matrix, right upper corners of boxes
    * \param[out] level vector of length \em n, corresponding levels
    */
    // getBoxes-functions might get removed at some point of time.
    // Use iterators instead whenever possible.
    void getBoxes(gsMatrix<Z>& b1,
                  gsMatrix<Z>& b2,
                  gsVector<Z>& level) const;

    /** \brief Returns the boxes which make up the hierarchical domain
    * and the respective levels touching side \em s.
    *
    * Similar to getBoxes() with additional input of the side \em s.
    * Returns only those boxes which touch the side \em s of the domain.
    *
    * \param[in] s Side (as in boundary::side) specifying the side of the domain
    * \param[out] b1 <em>n</em> x <em>d</em>-matrix, left bottom corners of boxes
    * \param[out] b2 <em>n</em> x <em>d</em>-matrix, right upper corners of boxes
    * \param[out] level vector of length \em n, corresponding levels
    */
    // getBoxes-functions might get removed at some point of time.
    // Use iterators instead whenever possible.
    void getBoxesOnSide(boundary::side s,
                        gsMatrix<Z>& b1,
                        gsMatrix<Z>& b2,
                        gsVector<Z>& level) const;

    /// Returns a list of boxes defined by left-bottom (b1) and
    /// right-top (b2) corners for the splitting in B-spline patches
    /// together with the corresponding levelOf
    /// b1 and b2 are indexing in the corresponding level indices
    /// \param b1 left bottom corners of boxes
    /// \param b2 right upper corners of boxs
    /// \param level corresponding levels
    // getBoxes-functions might get removed at some point of time.
    // Use iterators instead whenever possible.
    void getBoxesInLevelIndex(gsMatrix<Z>& b1,
                              gsMatrix<Z>& b2,
                              gsVector<Z>& level) const;

    ///return a list of polylines- boundaries of each connected
    ///component for all levels in the parameter space
    std::vector<std::vector<std::vector<std::vector<Z> > > > getPolylines() const;

    std::vector<std::vector<std::vector<Z> > > getPolylinesSingleLevel(std::vector<gsVSegment<Z> >& seg) const;

    /// Returns the level of the point \a p
    int levelOf(point const & p, int level) const
    { return this->knotSearch(p,level)->nodeData().level();}

    /// Return the upper corner of the tree in level 0
    const point upperCornerIndex() const
    {
        point ind = m_upperIndex;
        for (short_t i=0; i<d; ++i)
            ind[i] = (ind[i] >> m_indexLevel);
        return ind;
    }
    /// Inserts the box [\a k1, \a k2] (level-local indices of level \a lvl) at level \a lvl
    void insertBox(point const & k1, point const & k2, int lvl)
    {
        GISMO_ENSURE( lvl <= static_cast<int>(m_indexLevel), "Max index level reached..");

        // Make a box
        data_t iData(k1,k2,lvl);
        if( iData.aabb().isDegenerate() )
            return;

        // Represent box in the index level
        local2globalIndex( iData.lowerCorner(), static_cast<unsigned>(lvl), iData.lowerCorner() );
        local2globalIndex( iData.upperCorner(), static_cast<unsigned>(lvl), iData.upperCorner());

        // Ensure that the box is within the valid limits
        if ( ( iData.lowerCorner().array() >= m_upperIndex.array() ).any() )
        {
            gsWarn<<" Invalid box coordinate "<<  k1.transpose() <<" at level" <<lvl<<".\n";
            return;
        }

        // Initialize stack
        std::vector<node*> stack;
        stack.reserve( 2 * (m_maxPath + d) );
        stack.push_back(this);

        node * curNode;
        while ( ! stack.empty() )
        {
            curNode = stack.back();
            stack.pop_back();

            if ( curNode->isLeaf() ) // reached a leaf
            {
                // Since we reached a leaf, it should overlap with iBox

                // If this leaf is already in level lvl, then we have nothing to do
                if ( lvl <= curNode->nodeData().level() )
                    continue;

                // Split the leaf (if possible)
                const unsigned h = 1 << (m_indexLevel - curNode->nodeData().level()) ;
                node * newLeaf = curNode->adaptiveAlignedSplit(iData, h);

                // If curNode is still a leaf, its domain is almost
                // contained in iBox
                if ( !newLeaf ) //  curNode->isLeaf()
                {
                    // Increase level and recurse
                    if ( ++curNode->nodeData().level() != lvl)
                        stack.push_back(curNode);
                }
                else // treat new child
                {
                    stack.push_back(newLeaf);
                }
            }
            else // roll down the tree
            {
                if ( iData.upperCorner()[curNode->axis] <= curNode->pos)
                    // iBox overlaps only left child of this split-node
                    stack.push_back(curNode->left);
                else if  ( iData.lowerCorner()[curNode->axis] >= curNode->pos)
                    // iBox overlaps only right child of this split-node
                    stack.push_back(curNode->right);
                else
                {
                    // iBox overlaps both children of this split-node
                    stack.push_back(curNode->left );
                    stack.push_back(curNode->right);
                }
            }
        }

        // Update maximum inserted level
        if ( static_cast<unsigned>(lvl) > m_maxInsLevel)
            m_maxInsLevel = lvl;
    }

    /// Raises by one level the leaves inside the box [\a k1, \a k2] (level-local indices of level \a lvl)
    void sinkBox (point const & k1,
                  point const & k2, int lvl)
    {
        GISMO_ENSURE( m_maxInsLevel+1 <= m_indexLevel,
                      "Max index level might be reached..");

        // Make a box
        data_t iData(k1,k2,lvl);
        if( iData.aabb().isDegenerate() )
            return;

        // Represent box in the index level
        local2globalIndex( iData.lowerCorner() , static_cast<unsigned>(lvl), iData.lowerCorner() );
        local2globalIndex( iData.upperCorner(), static_cast<unsigned>(lvl), iData.upperCorner());

        // Ensure that the box is within the valid limits
        if ( ( iData.lowerCorner().array() >= m_upperIndex.array() ).any() )
        {
            return;
        }

        // Initialize stack
        std::stack<node*, std::vector<node*> > stack;
        stack.push(this);

        node * curNode;
        while ( ! stack.empty() )
        {
            curNode = stack.top();
            stack.pop();

            if ( curNode->isLeaf() ) // reached a leaf
            {
                // Since we reached a leaf, it should overlap with iBox.
                // Split the leaf (if possible)
                const unsigned h = 1 << (m_indexLevel - curNode->nodeData().level()) ;
                node * newLeaf = curNode->adaptiveAlignedSplit(iData, h);

                // If curNode is still a leaf, its domain is almost
                // contained in iBox
                if ( !newLeaf ) //  implies curNode was a leaf
                {
                    // Increase level
                    if ( ++curNode->nodeData().level() > static_cast<int>(m_maxInsLevel) )
                        m_maxInsLevel = curNode->nodeData().level();
                }
                else // treat new child
                {
                    stack.push(newLeaf);
                }
            }
            else // walk down the tree
            {
                if ( iData.upperCorner()[curNode->axis] <= curNode->pos)
                    // iBox overlaps only left child of this split-node
                    stack.push(curNode->left);
                else if  ( iData.lowerCorner()[curNode->axis] >= curNode->pos)
                    // iBox overlaps only right child of this split-node
                    stack.push(curNode->right);
                else
                {
                    // iBox overlaps both children of this split-node
                    stack.push(curNode->left );
                    stack.push(curNode->right);
                }
            }
        }
    }

    /// Lowers the leaves inside the box [\a k1, \a k2] (level-local indices of level \a lvl) to level \a lvl, at least
    void clearBox ( point const & k1, point const & k2,
                    int lvl) // CONSTRAINT: lvl is "minimum level"
    {
        GISMO_ENSURE( lvl <= static_cast<int>(m_indexLevel), "Max index level reached..");

        // Make a box
        data_t iData(k1,k2,lvl);
        if( iData.aabb().isDegenerate() )
            return;

        // Represent box in the index level
        local2globalIndex( iData.lowerCorner() , static_cast<unsigned>(lvl), iData.lowerCorner() );
        local2globalIndex( iData.upperCorner(), static_cast<unsigned>(lvl), iData.upperCorner());

        // Ensure that the box is within the valid limits
        if ( ( iData.lowerCorner().array() >= m_upperIndex.array() ).any() )
        {
            gsWarn<<" Invalid box coordinate "<<  k1.transpose() <<" at level" <<lvl<<".\n";
            return;
        }

        // Initialize stack
        std::vector<node*> stack;
        stack.reserve( 2 * (m_maxPath + d) );
        stack.push_back(this);

        node * curNode;
        while ( ! stack.empty() )
        {
            curNode = stack.back();
            stack.pop_back();

            if ( curNode->isLeaf() ) // reached a leaf
            {
                // Since we reached a leaf, it should overlap with iBox

                // If this leaf is already in level lvl, then we have nothing to do
                if ( curNode->nodeData().level() <= lvl )
                    continue;

                // Split the leaf (if possible)
                const unsigned h = 1 << (m_indexLevel - curNode->nodeData().level()) ;
                node * newLeaf = curNode->adaptiveAlignedSplit(iData, h);
                // If curNode is still a leaf, its domain is almost
                // contained in iBox
                if ( !newLeaf ) //  curNode->isLeaf()
                {
                    // Decrease level and reccurse
                    if ( --curNode->nodeData().level() != lvl)
                        stack.push_back(curNode);
                }
                else // treat new child
                {
                    stack.push_back(newLeaf);
                }
            }
            else // roll down the tree
            {
                if ( iData.upperCorner()[curNode->axis] <= curNode->pos)
                    // iBox overlaps only left child of this split-node
                    stack.push_back(curNode->left);
                else if  ( iData.lowerCorner()[curNode->axis] >= curNode->pos)
                    // iBox overlaps only right child of this split-node
                    stack.push_back(curNode->right);
                else
                {
                    // iBox overlaps both children of this split-node
                    stack.push_back(curNode->left );
                    stack.push_back(curNode->right);
                }
            }
        }

        m_maxInsLevel = this->template leafSearch< internal::maxLevel_visitor >(); // compute again globally
    }

    bool query1(point const & k1, point const & k2,
                int level) const
    {
        return boxSearch< internal::query1_visitor >(k1,k2,level);
    }

    bool query2(point const & k1, point const & k2,
                int level) const
    {
        return boxSearch< internal::query2_visitor >(k1,k2,level);
    }

    int query3(point const & k1, point const & k2,
               int level) const
    {
        return boxSearch< internal::query3_visitor >(k1,k2,level);
    }

    int query4(point const & k1, point const & k2,
               int level) const
    {
        return boxSearch< internal::query4_visitor >(k1,k2,level);
    }

    std::pair<point, point>
    queryLevelCell(point const & lower, point const & upper,
                   int level) const
    {
        std::pair<point,point> tmp = boxSearch< get_cell_visitor >(lower,upper,level);
        global2localIndex(tmp.first,level,tmp.first);
        global2localIndex(tmp.second,level,tmp.second);
        return tmp;
    }

    inline void local2globalIndex(gsVector<Z, d> const & index,
                                  unsigned lvl,
                                  gsVector<Z, d> & result) const
    {
        for(short_t i = 0; i!=d; ++i)
            result[i] = index[i] << (m_indexLevel-lvl) ;
    }

    inline void global2localIndex(gsVector<Z, d> const & index,
                                  unsigned lvl,
                                  gsVector<Z, d> & result) const
    {
        for(short_t i = 0; i!=d; ++i)
            result[i] = index[i] >> (this->m_indexLevel-lvl) ;
    }

    void computeFinestIndex(gsVector<Z, d> const & index,
                            unsigned lvl,
                            gsVector<Z, d> & result) const
    {
        for(short_t i = 0; i!=d; ++i)
            result[i] = index[i] << (m_maxInsLevel-lvl) ;
    }

    void computeLevelIndex(gsVector<Z, d> const & index,
                           unsigned lvl,
                           gsVector<Z, d> & result) const
    {
        for(short_t i = 0; i!=d; ++i)
            result[i] = index[i] >> (m_maxInsLevel-lvl) ;
    }

private:

    /// \brief Returns the boxes of all leaves, as vectors
    /// [ [lower corner],[upper corner], Level ], with corners given as
    /// knot indices on level gsHTree::m_maxInsLevel.
    void getBoxes_vec(std::vector<std::vector<Z> >& boxes) const;

    /// \brief connect the boxes returned from quadtree getBoxes_vec()
    ///
    /// If two neighbouring boxes could be represented by a single
    /// box (i.e., if the union of two axis-aligned boxes is again
    /// an axis-aligned box), then these two are merged into a single box.
    ///
    /// \param[in,out] boxes Format as of getBoxes_vec(), i.e., each box
    /// is represented as vector of size <em>2*d + 1</em> containing
    /// [ [lower corner],[upper corner], Level ], where the corners
    /// are defined by the knot indices on level gsHTree::m_maxInsLevel.
    void connect_Boxes(std::vector<std::vector<Z> > &boxes) const;

    /// For each x-coordinate delete repeated parts of vertical segments.
    /// \a vert_seg_lists list, where for each x coordinate (sorted increasingly) a list of vertical segments
    /// with that coordinate is saved.
    void getRidOfOverlaps( std::list< std::list< gsVSegment<Z> > >& vert_seg_lists ) const;

    /// Sweepline algorithm.
    void sweeplineConnectAndMerge( std::vector< std::vector< std::vector<Z> > >& result,
                                   std::list< std::list< gsVSegment<Z> > >& vert_seg_lists ) const;

    /// Applies \a visitor to the leaves hit by the box [\a k1, \a k2],
    /// given in the indices of level \a level.  The visitor controls the
    /// return type and the update of the result type at every leaf node
    template<typename visitor>
    typename visitor::return_type
    boxSearch(point k1, point k2, int level) const
    {
        local2globalIndex( k1, static_cast<unsigned>(level), k1 );
        local2globalIndex( k2, static_cast<unsigned>(level), k2 );
        return this->template rangeSearch< visitor>(k1,k2,level, m_maxPath);
    }

    /// Returns the leaf containing the point \a p, given in the indices of level \a level
    const node * knotSearch(const point & p, int level) const
    {
        point pp;
        local2globalIndex(p, static_cast<unsigned>(level), pp);

        GISMO_ASSERT( ( pp.array() <= m_upperIndex.array() ).all(),
                      "pointSearch: Wrong input: "<< p.transpose()<<", level "<<level
                      << ", i.e. "<< pp.transpose()
                      << " should not exceed "<< m_upperIndex.transpose()<< ".\n" );
        return node::pointSearch(pp);
    }

    /// Returns an cell/element box of a requested level
    // TODO: stop traverse as soon as it is found for the first time..
    struct get_cell_visitor
    {
        typedef std::pair<point,point> return_type;

        // initialize result
        static return_type init()
        {
            return std::make_pair(point::Zero(),point::Zero());
        }

        static void visitLeaf(const gismo::gsKdTree<d,Z,gsHTreeData<d,Z> > * leafNode , int level, return_type & res)
        {
            if ( leafNode->nodeData().level() == level )
            {
                res.first  = leafNode->nodeData().lowerCorner();
                res.second = leafNode->nodeData().upperCorner();
            }
        }
    };

public:
    // see http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsHTree.hpp)
#endif
