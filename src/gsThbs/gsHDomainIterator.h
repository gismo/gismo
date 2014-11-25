/** @file gsHDomainIterator.h

    @brief Provides declaration of iterator of hierarchical basis.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): J. Speh
*/

#pragma once

#include <gsThbs/gsHDomain.h>

#include <gsNurbs/gsCompactKnotVector.h>
#include <gsNurbs/gsTensorBSplineBasis.h>

#include <gsCore/gsDomainIterator.h>

namespace gismo
{

/** 
 @brief Re-implements gsDomainIterator for iteration over all elements
 of a <b>hierarchical</b> parameter domain.\n

  <em>See gsDomainIterator for more detailed documentation and an
  example of the typical use!!!</em>\n Used, e.g., for basis of
  classes gsHTensorBasis or gsTHBSplineBasis.
 */

template<typename T, unsigned d>
class gsHDomainIterator: public gsDomainIterator<T>
{
public:
    // here I have unsigned, maybe put another template parameter for this
    // unsigned...
    typedef kdnode<d, unsigned> node;

    typedef typename node::point point; 

    typedef typename std::vector<T>::const_iterator  uiter;


public:
    gsHDomainIterator(const gsHTensorBasis<d, T>& hbs)
        : gsDomainIterator<T>(hbs),
          m_tree(hbs.tree()),
          m_bases(hbs.getBases())
    {
        init();
    }

    gsHDomainIterator(const gsTHBSplineBasis<d, T>& thbs)
        : gsDomainIterator<T>(thbs),
          m_tree(thbs.tree()),
          m_bases(thbs.getBases())
    {
        init();
    }


    // ---> Documentation in gsDomainIterator.h
    bool next()
    {
        this->m_isGood = nextLexicographic(m_curElement, m_meshStart, m_meshEnd);

        if (this->m_isGood) // new element in m_cnode
            update();
        else // went through all elements in m_cnode
            this->m_isGood = selectFirstElementInNextBox();

        return this->m_isGood;
    }

    // if we want to iterate again throuh all elements, we must reset iterator
    /// Resets the iterator so that it can be used for another iteration through all elements.
    void reset()
    {
        this->m_isGood = selectFirstElementInNextBox();
    }



    // ---> Documentation in gsDomainIterator.h
    // Compute a suitable quadrature rule of the given order for the current element
    void computeQuadratureRule(const gsVector<int>& numIntNodes)
    {
        m_quadrature.setNodes(numIntNodes);
        m_quadrature.mapTo(m_lower, m_upper, this->quNodes, this->quWeights);
    }


    // get the basis function indices which are active in the current element
    void getActiveFunctions(gsMatrix<unsigned>& act)
    {
        this->m_basis.active_into(center, act);
    }


    const gsMatrix<unsigned>& computeActiveFunctions()
    {
        this->m_basis.active_into(center, this->activeFuncs);
        return this->activeFuncs;
    }

    const gsVector<T>& lowerCorner() const { return m_lower; }

    const gsVector<T>& upperCorner() const { return m_upper; }

    int getLevel() const
    {
        return m_cnode->level;
    }

private:
    /// proceed to the next box; returns true if the end is not reached yet
    ///
    // typical usage
    // while (gsHDomainIterator::nextBox())
    // {
    //      use current box via m_cnode
    // }
    //
    // you can reach current box via pointer m_cnode that points to a node which
    // defines a box
    bool nextBox()
    {
        if (m_stopNextBox)
        {
            // set next node for new iteration over whole domain
            m_nextNode = m_tree.m_root;
            m_stopNextBox = false;
            return false;
        }

        // current node
        node* cnode = m_nextNode;

        while (true)
        {
            if (cnode != NULL)
            {
                m_stack.push(cnode);
                cnode = cnode->left;
            }
            else
            {
                if (m_stack.empty())
                {
                    // this should happen only when the tree is empty
                    gsWarn << "We have an empty tree\n";
                    return false;
                }

                cnode = m_stack.top();
                m_stack.pop();

                // current node should be a leaf (we updated implementation)
                if (cnode->right != NULL)
                {
                    gsWarn << "This should not happen in this version of "
                        "KD-tree: cnode is not a leaf.\n";

                    return false;
                }

                // cnode is a leaf // TO DO: use m_cnode->isLeaf()
                this->m_cnode = cnode;

                // we pop out "a parent" of cnode and continue to the right
                if (m_stack.empty())
                {
                    // this is the last iteration, we return false next call
                    // of the nextBox function
                    m_stopNextBox = true;
                    return true;
                }

                cnode = m_stack.top();
                m_stack.pop();
                m_nextNode = cnode->right;

                return true;

            }
        }

        gsWarn << "Something horribly wrong happend.\n";
        return false;
    }


    /// returns true if there is a box with an element
    bool selectFirstElementInNextBox()
    {
        m_boxSelected = false;
        while (!m_boxSelected)
        {
            if (!nextBox())
                return false;
            else
            {
                m_boxSelected = true;
                isAligned();
                if (m_isAligned)
                    setAlignedMesh();
                else
                    setNotAlignedMesh();
            }
        }

        update();
        return true;
    }

    /// Computes lower, upper and center point of the current element, maps the reference
    /// quadrature nodes and weights to the current element, and computes the
    /// active functions.
    void update()
    {
        for (unsigned dim = 0; dim < d; ++dim)
        {
            m_lower(dim) = *(m_curElement(dim));
            m_upper(dim) = *(m_curElement(dim) + 1);
            center(dim) = T(0.5) * (m_lower(dim) + m_upper(dim));
        }

        // Update quadrature rule
        m_quadrature.mapTo(m_lower, m_upper, this->quNodes, this->quWeights);
        computeActiveFunctions();
    }

    void setNotAlignedMesh()
    {
        const point & lower = m_cnode->lowCorner();
        const point & upper = m_cnode->uppCorner();

        const int mult = this->mult();

        for (unsigned dim = 0; dim < d; dim++)
        {
            m_breaks[dim].clear();

            unsigned start = lower(dim);
            const unsigned end = upper(dim);

            if (start == end)
            {
                m_boxSelected = false;
                return;
            }

            const gsCompactKnotVector<T> & kv =
                    m_bases[m_tree.m_indexLevel]->component(dim).knots();

            if (start % mult != 0)
            {
                m_breaks[dim].push_back(kv.uValue(start));

                start = start / mult;
                start = (start + 1) * mult;
            }


            for (unsigned mid = start; mid < end; mid += mult)
            {
                m_breaks[dim].push_back(kv.uValue(mid));
            }

            m_breaks[dim].push_back(kv.uValue(end));

            updateMeshVariables(dim);
        }
    }


    void setAlignedMesh()
    {
        const point & lower = m_cnode->lowCorner();
        const point & upper = m_cnode->uppCorner();

        // 2^(maxLevel - cnodeLevel)
        const int mult = this->mult();
        const int level = m_cnode->level;

        for (unsigned dim = 0; dim < d; ++dim)
        {
            const unsigned start = lower(dim) / mult;
            const unsigned end = upper(dim) / mult;
            const gsCompactKnotVector<T> & kv =
                    m_bases[level]->component(dim).knots();

            m_breaks[dim].clear();
            for (unsigned index = start; index <= end; ++index)
                m_breaks[dim].push_back(kv.uValue(index));

            updateMeshVariables(dim);
        }
    }


    void updateMeshVariables(unsigned dim)
    {
        m_curElement(dim) = m_breaks[dim].begin();
        m_meshStart(dim) = m_breaks[dim].begin();

        // for n breaks, we have n - 1 elements (spans)
        m_meshEnd(dim) = m_breaks[dim].end() - 1;

        // selected box is degenerated
        if (m_meshEnd(dim) == m_meshStart(dim))
            m_boxSelected = false;
    }

    /// Initializes some variables.
    void init()
    {
        m_stopNextBox = false;
        m_nextNode = m_tree.m_root;

        m_meshEnd.resize(d);
        m_meshStart.resize(d);
        m_curElement.resize(d);

        m_lower.resize(d);
        m_upper.resize(d);

        m_breaks = std::vector<std::vector<T> >(d, std::vector<T>());

        // Set to one quadrature point by default
        m_quadrature.setNodes( gsVector<int>::Ones(d) );

        reset();
    }


    /// Checks is m_cnode is aligned in level m_cnode.level
    bool isAligned()
    {
        const point & lower = m_cnode->lowCorner();
        const point & upper = m_cnode->uppCorner();
        
        const int mult = this->mult();

        for (unsigned dim = 0; dim < d; ++dim)
        {
            if (lower(dim) % mult != 0 ||
                upper(dim) % mult != 0)
            {
                m_isAligned = false;
                return false;
            }
        }
        m_isAligned = true;
        return true;
    }

    /// compute 2^(maxLevel - cnodeLevel)
    int mult()
    {
        const int level = m_cnode->level;
        const int maxLevel = static_cast<int>(m_tree.m_indexLevel);

        return 1 << (maxLevel - level);
    }


// =============================================================================
// members
// =============================================================================

public:

    using gsDomainIterator<T>::center;

    // Alignment operator Eigen-derived members of this class which
    // have static size
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
    // current node, points to current box in nextBox iterator
    node* m_cnode;

    // kd tree
    const gsHDomain< d >& m_tree;

    // the list of nestes spaces
    const std::vector<gsTensorBSplineBasis<d,T,gsCompactKnotVector<T> >* >& m_bases;

    // '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    // variables for nextBox() iteration
    // .........................................................................

    // next node, we will use this node in next iteration of nextBox function
    node* m_nextNode;

    // stack of pointers to nodes, used in nextBox function
    std::stack<node*> m_stack;

    // stop the next iteration of nextBox function
    bool m_stopNextBox;

    // if the box (pointed by m_cnode) is aligned with level - 1, m_isAligned
    // is True
    bool m_isAligned;

    // true is box is selected, false otherwise
    bool m_boxSelected;

    // '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    // variables for next() iteration
    // .........................................................................

    // extent of the tensor grid
    gsVector<uiter, d> m_meshStart, m_meshEnd;

    // coordinates of the grid cell boundaries
    std::vector< std::vector<T> > m_breaks;

    // current element as pointers to it's supporting mesh-lines
    gsVector<uiter, d> m_curElement;

    // '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    // variables for element
    // .........................................................................

    // Quadrature rule
    gsGaussRule<T> m_quadrature;

    // parameter coordinates of current grid cell
    gsVector<T> m_lower, m_upper;


};

} // end namespace gismo
