/** @file gsHDomainIterator.h

    @brief Provides declaration of iterator of hierarchical domain.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): J. Speh
*/

#pragma once

#include <gsHSplines/gsHDomain.h>

#include <gsNurbs/gsTensorBSplineBasis.h>

#include <gsCore/gsDomainIterator.h>

namespace gismo
{

// Documentation in gsDomainIterator
/** @brief Re-implements gsDomainIterator for iteration over all boundary
  * elements of a <b>hierarchical</b> parameter domain.
  * 
  * <em>See
  * gsDomainIterator for more detailed documentation and an example of
  * the typical use!!!</em>\n Used, e.g., for basis of classes
  * gsHTensorBasis or gsTHBSplineBasis.
  *
  * \ingroup HSplines
  */

template<typename T, unsigned d>
class gsHDomainIterator: public gsDomainIterator<T>
{
public:

    typedef kdnode<d, unsigned> node;

    typedef typename node::point point; 

    typedef typename std::vector<T>::const_iterator  uiter;

    typedef gsHDomain<d,unsigned> hDomain;

    typedef typename hDomain::const_literator leafIterator;

public:

    gsHDomainIterator(const gsHTensorBasis<d, T> & hbs)
    : gsDomainIterator<T>(hbs)
    {
        // Initialize mesh data
        m_meshStart.resize(d);
        m_meshEnd  .resize(d);
        
        // Initialize cell data
        m_curElement.resize(d);
        m_lower     .resize(d);
        m_upper     .resize(d);

        // Allocate breaks
        m_breaks = std::vector<std::vector<T> >(d, std::vector<T>());

        // Set to one quadrature point by default
        m_quadrature.setNodes( gsVector<int>::Ones(d) );

        m_leaf = hbs.tree().beginLeafIterator();
        updateLeaf();
    }

    // ---> Documentation in gsDomainIterator.h
    bool next()
    {
        this->m_isGood = nextLexicographic(m_curElement, m_meshStart, m_meshEnd);

        if (this->m_isGood) // new element in m_leaf
            updateElement();
        else // went through all elements in m_leaf
            this->m_isGood = nextLeaf();

        return this->m_isGood;
    }

    // ---> Documentation in gsDomainIterator.h
    bool next(index_t increment)
    {
      for (index_t i = 0; i < increment; i++)
            this->m_isGood = nextLexicographic(m_curElement, m_meshStart, m_meshEnd);

        if (this->m_isGood) // new element in m_leaf
            updateElement();
        else // went through all elements in m_leaf
            this->m_isGood = nextLeaf();

        return this->m_isGood;
    }
    
    /// Resets the iterator so that it can be used for another
    /// iteration through all boundary elements.
    void reset()
    {
        const gsHTensorBasis<d, T>* hbs =  dynamic_cast<const gsHTensorBasis<d, T> *>(m_basis);
        m_leaf = hbs->tree().beginLeafIterator();
        updateLeaf();
    }
    
    // ---> Documentation in gsDomainIterator.h Compute a suitable
    // quadrature rule of the given order for the current element
    void computeQuadratureRule(const gsVector<int>& numIntNodes)
    {
        m_quadrature.setNodes(numIntNodes);
        m_quadrature.mapTo(m_lower, m_upper, this->quNodes, this->quWeights);
    }

    void computeQuadratureRuleDefault()
    {
        // uses same formula as gsGaussAssembler::getNumIntNodesFor( gsBasis )
        gsVector<int> numIntNodes( m_basis->dim() );
        for (int i = 0; i < m_basis->dim(); ++i)
            numIntNodes[i] = m_basis->degree(i) + 1;

        computeQuadratureRule( numIntNodes );
    }
    
    // get the basis function indices which are active in the current
    // element
    void getActiveFunctions(gsMatrix<unsigned>& act)
    {
        this->m_basis->active_into(center, act);
    }
    
    const gsMatrix<unsigned>& computeActiveFunctions()
    {
        this->m_basis->active_into(center, this->activeFuncs);
        return this->activeFuncs;
    }

    const gsVector<T>& lowerCorner() const { return m_lower; }

    const gsVector<T>& upperCorner() const { return m_upper; }

    int getLevel() const
    {
        return m_leaf.level();
    }

private:

    gsHDomainIterator();

    /// returns true if there is a another leaf with a boundary element
    bool nextLeaf()
    {
        this->m_isGood = m_leaf.next();

        if ( m_leaf.good() )
            updateLeaf();

        return this->m_isGood;
    }

    /// Computes lower, upper and center point of the current element, maps the reference
    /// quadrature nodes and weights to the current element, and computes the
    /// active functions.
    void updateLeaf()
    {
        const point & lower = m_leaf.lowerCorner();
        const point & upper = m_leaf.upperCorner();
        // gsDebug<<"leaf "<<  lower.transpose() <<", " 
        //        << upper.transpose() <<"\n";

        const int level2 = m_leaf.level();

        // Update leaf box
        for (unsigned dim = 0; dim < d; ++dim)
        {
            const unsigned start = lower(dim);
            const unsigned end  = upper(dim) ;

            const gsKnotVector<T> & kv =
                static_cast<const gsHTensorBasis<d,T>*>(m_basis)
                ->tensorLevel(level2).component(dim).knots();

            // knotVals = kv.unique()

            m_breaks[dim].clear();
            for (unsigned index = start; index <= end; ++index)
                m_breaks[dim].push_back(kv.uValue(index));

            m_curElement(dim) = 
            m_meshStart(dim)  = m_breaks[dim].begin();

            // for n breaks, we have n - 1 elements (spans)
            m_meshEnd(dim) =  m_breaks[dim].end() - 1;
        }

        // We are at a new element, so update cell data
        updateElement();
    }

    /// Computes lower, upper and center point of the current element, maps the reference
    /// quadrature nodes and weights to the current element, and computes the
    /// active functions.
    void updateElement()
    {
        // Update cell data
        for (unsigned i = 0; i < d ; ++i)
        {
            m_lower[i]  = *m_curElement[i];
            m_upper[i]  = *(m_curElement[i]+1);
            center[i] = T(0.5) * (m_lower[i] + m_upper[i]);
        }

        // Update quadrature rule
        m_quadrature.mapTo(m_lower, m_upper, this->quNodes, this->quWeights);

        // Update Active basis functions
        computeActiveFunctions();
    }

// =============================================================================
// members
// =============================================================================

public:

    using gsDomainIterator<T>::center;
    using gsDomainIterator<T>::m_basis;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:

    // The current leaf node of the tree
    leafIterator m_leaf;

    // Coordinates of the grid cell boundaries
    // \todo remove this member
    std::vector< std::vector<T> > m_breaks;

    // Extent of the tensor grid
    gsVector<uiter, d> m_meshStart, m_meshEnd;

    // Current element as pointers to it's supporting mesh-lines
    gsVector<uiter, d> m_curElement;

    // parameter coordinates of current grid cell
    gsVector<T> m_lower, m_upper;

    // Quadrature rule
    gsGaussRule<T> m_quadrature;

};

} // end namespace gismo
