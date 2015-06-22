/** @file gsHDomainSliceIter.h

    @brief Provides declaration of HDomainLeafIter class.

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
class gsHDomainLeafIter
{
public:
    //typedef kdnode<d, unsigned> node;
    typedef typename choose<isconst, const node&, node&>::type reference;
    typedef typename choose<isconst, const node*, node*>::type pointer;

    typedef typename node::point point;

    typedef typename node::point::Projection_t slicedPoint;// to do: this must have dimension = d-1

public:
    reference operator*() const { return *m_curNode; }
    pointer  operator->() const { return  m_curNode; }

public:

    gsHDomainSliceIter() 
    : m_dir(0), m_pos(0), m_curNode(0), m_index_level(0)
    { }

    // to do: constructor by gsHDomain instead
    gsHDomainSliceIter( node * const root_node, 
                        const unsigned _dir, 
                        const unsigned _pos,
                        const unsigned index_level)
    : m_dir(m_dir), m_pos(_pos), m_index_level(index_level)
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
            m_curNode = m_stack.top();
            m_stack.pop();
            
            if ( m_curNode->isLeaf() )
            {
                // does this box intersect the slice ?  
                // note: boxes are considered half-open, eg. products
                // of intervals [a,b), except from the rightmost
                // interval which is closed
                if ( (m_curNode->box->lowCorner()[dir] <= m_pos) && 
                     (m_curNode->box->uppCorner()[dir] >  m_pos  ||
                     (m_pos == m_last && m_curNode->box->uppCorner()[dir] == m_pos) )
                   )
                     return true;
            }
            else // this is a split-node
            {
                m_stack.push(m_curNode->left );
                m_stack.push(m_curNode->right);
            }
        }

        // Leaves exhausted
        m_curNode = NULL;
        return false;
    }

    /// Returns true iff we are still pointing at a valid leaf
    bool good() const   { return m_curNode != 0; }

    /// The iteration is done in the sub-tree hanging from node \start
    void startFrom( node * const root_node)
    {
        m_stack.clear();
        m_stack.push(root_node);
    }

    int level() const { return m_curNode->level; }

    /// \brief The lower corner of the sliced box
    /// TODO: Note that m_dir is skipped
    slicedPoint lowerCorner() const
    { 
        point result = m_curNode->box->first;
        const int lvl = m_curNode->level;

        for ( index_t i = 0; i!= result.size(); ++i )
            result[i] = result[i] >> (m_index_level-lvl) ;

        return result; 
    }

    /// \brief The upper corner of the sliced box
    /// TODO: Note that m_dir is skipped
    slicedPoint upperCorner() const
    { 
        point result = m_curNode->box->second;
        const int lvl = m_curNode->level;

        for ( index_t i = 0; i!=result.size(); ++i )
            result[i] = result[i] >> (m_index_level-lvl) ;

        return result; 
    }

    unsigned indexLevel() {return m_index_level;}

private:

    // Slice information
    unsigned m_dir; // direction normal to the slice
    unsigned m_pos; // slice position index (span-index)
    unsigned m_last; // last (span-index) in direction \a m_dir

    // current node
    node * m_curNode;

    /// The level of the box representation
    unsigned m_index_level;

    // stack of pointers to tree nodes, used in next()
    std::stack<node*> m_stack;// to do: change type to std::vector
};



} // end namespace gismo
