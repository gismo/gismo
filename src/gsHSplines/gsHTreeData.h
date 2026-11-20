/** @file gsHTreeData.h

    @brief Leaf data policy of gsHTree: an axis-aligned box in global
    indices together with its hierarchical level.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <ostream>

#include <gsDomain/gsKdTree.h>
#include <gsHSplines/gsAABB.h>

namespace gismo {

/**
    @brief Leaf data of gsHTree.

    Stores the box covered by a leaf and the hierarchical level of the
    leaf. The corners of the box are global indices on the index level of
    the owning gsHTree (the finest grid), not level-local indices.

    \tparam d the dimension
    \tparam Z the coordinate index type

    \ingroup HSplines
*/
template<short_t d, class Z>
class gsHTreeData : public gsKdTreeData<d,Z>
{
public:
    typedef          gsAABB<d, Z> kdBox;
    typedef typename kdBox::point point;
private:

    int m_level;
    kdBox box;

public:
    gsHTreeData(point const & k1, point const & k2, int lvl = 0)
    : m_level(lvl), box(k1,k2)
    { }

    explicit gsHTreeData(point const & k2)
    : m_level(0), box(k2)
    { }

    friend std::ostream & operator<<(std::ostream & os, const gsHTreeData & td)
    {
        os << "Box("<<td.m_level<<"): ["<< td.box.first.transpose()
           <<"]~["<<td.box.second.transpose() <<"]\n";
        return os;
    }

    bool check() const
    { return (box.second.array() >= box.first.array()).all(); }

    /// Makes the children of the freshly split \a node meet at \a node.pos
    static void split(const gsKdTree<d,Z,gsHTreeData> & node)
    {
        node.left->nodeData().upperCorner()[node.axis] =
            node.right->nodeData().lowerCorner()[node.axis] = node.pos;
    }

    /// Merges the data of right child into left child of the node
    static void mergeToLeft(const gsKdTree<d,Z,gsHTreeData> & node)
    {
        node.left ->nodeData().upperCorner()[node.axis] =
            node.right->nodeData().upperCorner()[node.axis];
    }

    void multiplyByTwo()
    {
        box.first .array() *= 2;
        box.second.array() *= 2;
    }

    void divideByTwo()
    {
        box.first .array() /= 2;
        box.second.array() /= 2;
    }

    /// Chooses a split of \a node aligned to the grid of spacing \a h
    /// that separates the box of \a node from the box of \a insData.
    /// Sets \a doSplit to 1 (resp. -1) if the right (resp. left) child
    /// intersects the box \a insData, and to 0 if no split is needed.
    static void adaptiveAlignedSplit(const gsHTreeData & insData, unsigned h,
                                     gsKdTree<d,Z,gsHTreeData> &  node,
                                     int & doSplit)
    {
        for (short_t i = 0; i < d; ++i)
        {
            const Z c1 = insData.lowerCorner()[i] - insData.lowerCorner()[i] % h; //floor
            const Z cc = insData.upperCorner()[i] % h;
            const Z c2 = insData.upperCorner()[i] + (cc ? h-cc : 0 ); // ceil

            if ( c1 > node.nodeData().lowerCorner()[i] )
            {
                // right child intersects insBox
                node.axis = i;
                node.pos = c1;
                doSplit  = 1;
                return;
            }
            else if ( c2 <node.nodeData().upperCorner()[i]  )
            {
                // left child intersects insBox
                node.axis = i;
                node.pos = c2;
                doSplit  = -1;
                return;
            }
        }
        doSplit  = 0;
    }

    int level() const {return m_level;}
    int & level() {return m_level;}

    const point & lowerCorner() const { return box.first; }
          point & lowerCorner()       { return box.first; }

    const point & upperCorner() const { return box.second; }
          point & upperCorner()       { return box.second; }

    const kdBox & aabb() const { return box; }
    kdBox & aabb() { return box; }

    /// True iff the box is aligned to the grid of its own level, i.e.
    /// its global corners are multiples of 2^(index_level - level)
    bool isAligned(unsigned index_level) const
    {
        const unsigned h = 1 << (index_level - m_level);

        for ( index_t i = 0; i!=box.first.size(); ++i )
        {
            if (box.second[i] % h != 0 ||
                box.first[i]  % h != 0 )
                return false;
        }
        return true;
    }

    // see http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen
};

} // namespace gismo
