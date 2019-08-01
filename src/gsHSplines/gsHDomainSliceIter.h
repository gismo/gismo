/** @file gsHDomainSliceIter.h

    @brief Provides declaration of HDomainSliceIter class.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger, A. Mantzaflaris
*/

#pragma once

#include <gsHSplines/gsKdNode.h>
#include <gsCore/gsTemplateTools.h>

namespace gismo
{

/** 
    @brief Iterates over the leaves of an gsHDomain (tree) that
    intersect with a slice position.
    
    A slice in this context is defined by: 
    - The direction index \a dir which is perpendicular to the slice position
    - A span index \a pos along direction \a dir which indicates the position 
      of the "cut"

    \ingroup HSplines
*/

template<typename node, bool isconst = false>
class gsHDomainSliceIter
{
public:
    //typedef kdnode<d, index_t> node;
    typedef typename choose<isconst, const node&, node&>::type reference;
    typedef typename choose<isconst, const node*, node*>::type pointer;

    typedef typename node::point::Projection_t point;

    static const index_t d = point::_Rows;// d is the slice dimension

public:
    reference operator*() const { return *curNode; }
    pointer  operator->() const { return  curNode; }

public:

    gsHDomainSliceIter() 
    : m_dir(0), m_pos(0), m_last(0), curNode(0), m_index_level(0)
    { }

    gsHDomainSliceIter( node * const root_node, 
                        const unsigned _dir, 
                        const unsigned _pos,
                        const unsigned _last,
                        const unsigned index_level)
    : m_dir(m_dir), m_pos(_pos), m_last(_last), m_index_level(index_level)
    { 
        m_stack.push(root_node);

        // Go to the first leaf
        next();
    }

    // Next leaf
    bool next()
    {
        while ( ! m_stack.empty() )
        {
            curNode = m_stack.top();
            m_stack.pop();
            
            if ( curNode->isLeaf() )
            {
                // does this box intersect the slice ?  
                // note: boxes are considered half-open, eg. products
                // of intervals [a,b), except from the rightmost
                // interval which is closed
                // if ( (curNode->box->lowCorner()[m_dir] <= m_pos) && 
                //      (curNode->box->uppCorner()[m_dir] >  m_pos  ||
                //      (m_pos == m_last && curNode->box->uppCorner()[m_dir] == m_pos) )
                //    )
                     return true;
            }
            else // this is a split-node
            {
                if ( curNode->axis == m_dir )
                {
                    if (m_pos < curNode->pos)
                        m_stack.push(curNode->left);
                    else 
                        m_stack.push(curNode->right);
                }
                else
                {
                    m_stack.push(curNode->left);
                    m_stack.push(curNode->right);
                }
            }
        }

        // Leaves exhausted
        curNode = NULL;
        return false;
    }

    /// Returns true iff we are still pointing at a valid leaf
    bool good() const   { return curNode != 0; }

    /// The iteration is done in the sub-tree hanging from node \em root_node
    void startFrom( node * const root_node)
    {
        m_stack.clear();
        m_stack.push(root_node);
    }

    int level() const { return curNode->level; }

    /// \brief The lower corner of the sliced box.
    /// Note that \a m_dir is skipped
    point lowerCorner() const
    { 
        point result;
        result.topRows   (m_dir  ) = curNode->box->first.topRows(m_dir     );
        result.bottomRows(d-m_dir) = curNode->box->first.bottomRows(d-m_dir);

        const int lvl = curNode->level;

        for ( index_t i = 0; i!= result.size(); ++i )
            result[i] = result[i] >> (m_index_level-lvl) ;

        return result; 
    }

    /// \brief The upper corner of the sliced box.
    /// Note that \a m_dir is skipped
    point upperCorner() const
    { 
        point result = curNode->box->second;
        result.topRows   (m_dir  ) = curNode->box->second.topRows(m_dir     );
        result.bottomRows(d-m_dir) = curNode->box->second.bottomRows(d-m_dir);

        const int lvl = curNode->level;

        for ( index_t i = 0; i!=result.size(); ++i )
            result[i] = result[i] >> (m_index_level-lvl) ;

        return result; 
    }

    unsigned indexLevel() const {return m_index_level;}

private:

    // Slice information
    unsigned m_dir; // direction normal to the slice
    unsigned m_pos; // slice position index (span-index)
    unsigned m_last; // last (span-index) in direction \a m_dir --> most likely not needed

    // current node
    node * curNode;

    /// The level of the box representation
    unsigned m_index_level;

    // stack of pointers to tree nodes, used in next()
    std::stack<node*> m_stack;// to do: change type to std::vector
};



} // end namespace gismo
