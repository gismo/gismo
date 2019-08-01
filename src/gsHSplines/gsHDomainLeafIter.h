/** @file gsHDomainLeafIter.h

    @brief Provides declaration of HDomainLeafIter class.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsHSplines/gsKdNode.h>
#include <gsCore/gsTemplateTools.h>

namespace gismo
{

/** 
    @brief Iterates over the leaves of an gsHDomain (tree).
    
    \ingroup HSplines
*/

template<typename node, bool isconst = false>
class gsHDomainLeafIter
{
public:
    //typedef kdnode<d, index_t> node;
    typedef typename util::conditional<isconst, const node&, node&>::type reference;
    typedef typename util::conditional<isconst, const node*, node*>::type pointer;

    typedef typename node::point point;

public:
    reference operator*() const { return *curNode; }
    pointer  operator->() const { return  curNode; }

public:

    gsHDomainLeafIter() : curNode(0)
    { }

    explicit gsHDomainLeafIter( node * const root_node, index_t index_level)
        : m_index_level(index_level)
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
                return true;
            }
            else // this is a split-node
            {
                m_stack.push(curNode->left );
                m_stack.push(curNode->right);
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

    point lowerCorner() const
    { 
        point result = curNode->box->first;
        const int lvl = curNode->level;

        //result = result.array() / (1>> (m_index_level-lvl)) ;
        for ( index_t i = 0; i!= result.size(); ++i )
            result[i] = result[i] >> (m_index_level-lvl) ;

        return result; 
    }

    point upperCorner() const
    { 
        point result = curNode->box->second;
        const int lvl = curNode->level;

        for ( index_t i = 0; i!=result.size(); ++i )
            result[i] = result[i] >> (m_index_level-lvl) ;

        return result; 
    }

    index_t indexLevel() const {return m_index_level;}

private:

    // current node
    node * curNode;

    /// The level of the box representation
    index_t m_index_level;

    // stack of pointers to tree nodes, used in next()
    std::stack<node*> m_stack; // to do: change type to std::vector
};



} // end namespace gismo
