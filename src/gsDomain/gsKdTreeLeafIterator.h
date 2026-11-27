/** @file gsKdTreeLeafIterator.h

    @brief Provides a leaf iterator over a gsKdTree (sub)tree.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <stack>

#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsTemplateTools.h>

namespace gismo
{

/**
    @brief Iterates over the leaves of a gsKdTree (sub)tree.

    \tparam node the tree node type (an instance of gsKdTree), which must
    provide the member types \a data_t and \a point_t

    \ingroup HSplines
*/

template<typename node, bool isconst = false>
class gsKdTreeLeafIter
{
public:
    typedef typename util::conditional<isconst, const node&, node&>::type reference;
    typedef typename util::conditional<isconst, const node*, node*>::type pointer;
    
    typedef typename node::data_t data_t;
    typedef typename node::point_t point_t;

public:
    reference operator*() const { return *curNode; }
    pointer  operator->() const { return  curNode; }

public:

    gsKdTreeLeafIter() : curNode(0)
    { }

    explicit gsKdTreeLeafIter( const node * root_node) : curNode(nullptr) // const node *: treat node as const. node * const: treat pointer const
    {
        m_stack.push(const_cast<node*>(root_node));

        next();
    }

    bool next()
    {
        while ( ! m_stack.empty() )
        {
            curNode = m_stack.top();
            m_stack.pop();

            if ( curNode->axis == -2 ) // empty node: no leaves
                continue;

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

    /// The iteration is done in the sub-tree hanging from node \em root_node.
    /// The caller must call next() to reach the first leaf.
    void startFrom( node * const root_node)
    {
        m_stack = std::stack<node*>();
        m_stack.push(root_node);
    }

    data_t & data() const { return *curNode->data; }

private:

    // current node
    node * curNode;

    // stack of pointers to tree nodes, used in next()
    std::stack<node*> m_stack; // TODO: change type to std::vector
};



} // end namespace gismo
