
#pragma once

#include <gsDomain/gsKdTree.h>
#include <gsHSplines/gsAABB.h>

namespace gismo 
{

template<short_t d, class Z>
class gsCutTreeData
{
    typedef gsCutTreeData Self_t;
    typedef gsKdTree<d,Z,Self_t> Tree_t;
public:
    typedef          gsAABB<d, Z> kdBox;
    typedef typename kdBox::point point;
private:

    index_t m_level; ///< Uniform refinement level (same for all directions)
    kdBox box;
    short_t m_sign;

public:
    gsCutTreeData(point const & k1, point const & k2,
                  index_t lvl = 0)
    : m_level(lvl), box(k1,k2)
    { }

    explicit gsCutTreeData(point const & k2)
    : m_level(0), box(k2)
    { }

    friend std::ostream & operator<<(std::ostream & os, const Self_t & td)
    {
        os << "Box(level="<< td.m_level <<"): ["
           << td.box.first.transpose()
           <<"]~["<<td.box.second.transpose() <<"]\n";
        return os;
    }

    short_t sign() const { return m_sign; }

    short_t & sign() { return m_sign; }



    bool check() const
    { return (box.second.array() >= box.first.array()).any(); }

    static void split(const Tree_t & node)
    {
        node.left->nodeData().upperCorner()[node.axis] =
            node.right->nodeData().lowerCorner()[node.axis] = node.pos;
    }

    static void nextMidSplit(Tree_t &  node)
    {
        node.axis = ( node.parent == 0 ? 0 : (node.parent->axis+1)%d );
        node.pos  = node.box->first[node.axis] +
            (node.box->second[node.axis] - node.box->first[node.axis])/2 ;
    }

    static void anyMidSplit(unsigned h, Tree_t &  node,
                            int & doSplit)
    {
        const unsigned mask = ~(h - 1);
        for ( unsigned i = 0; i < d; ++i )
        {
            const unsigned c =
                (node.nodeData().box.first [i] + (node.nodeData().box.second[i] - node.nodeData().box.first[i])/2) & mask ;
            if ( c != (unsigned)node.nodeData().box.first [i] ) // avoid degenerate split
            {
                node.axis = i;
                node.pos  = c;
                doSplit   = 1;
                return;
            }
        }
        doSplit = 0;
    }

    // Merges the data of right child into left child of the node
    static void mergeToLeft(const gsKdTree<d,Z,Self_t> & node)
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

    index_t   level() const { return m_level; }
    index_t & level()       { return m_level; }

    const point & lowerCorner() const { return box.first; }
          point & lowerCorner()       { return box.first; }

    const point & upperCorner() const { return box.second; }
          point & upperCorner()       { return box.second; }

    const kdBox & aabb() const { return box; }
    kdBox & aabb() { return box; }

};

}
