/** @file gsTensorDomainBoundaryIterator.h

    @brief Iterator over the boundary elements of a tensor-structured grid

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsDomainIterator.h>

#include <gsAssembler/gsGaussRule.h>

#include <gsUtils/gsCombinatorics.h>

namespace gismo
{

/** 
 * @brief Re-implements gsDomainIterator for iteration over all elements of the boundary of a <b>tensor product</b> parameter domain.\n
 * <em>See gsDomainIterator for more detailed documentation and an example of the typical use!!!</em>
 *
 * \ingroup Tensor
 */

// Class which enables iteration over all elements of a tensor product parameter domain
// Documentation in gsDomainIterator.h

template<class T, int D, typename uiter>
//template<class T, int D=Dynamic, typename uiter = typename gsKnotVector<T>::const_uiterator>
class gsTensorDomainBoundaryIterator : public gsDomainIterator<T>
{
public:

    gsTensorDomainBoundaryIterator( const gsBasis<T>& b, const boxSide & s )
    : gsDomainIterator<T>(b, s),
      d( m_basis->dim() ),
      lower ( gsVector<T, D>::Zero(d) ),
      upper ( gsVector<T, D>::Zero(d) )
    {
        center =  gsVector<T, D>::Zero(d);
        par = s.parameter();
        dir = s.direction();
        meshBegin.resize(d);
        meshEnd.resize(d);
        curElement.resize(d);
        breaks.reserve(d);

        for (int i=0; i < dir; ++i) 
        {
            breaks.push_back( m_basis->component(i).domain()->breaks() );
            meshEnd[i]   = breaks[i].end() - 1;
            meshBegin[i] = curElement[i] = breaks[i].begin();
            //meshEnd[i]    = m_basis->component(i).domain()->uend() - 1;
            //meshBegin[i]  = 
            //curElement[i] = m_basis->component(i).knots().ubegin();

            if (meshEnd[i] == curElement[i])
                m_isGood = false;
        }

        // Fixed direction
        breaks.push_back( m_basis->component(dir).domain()->breaks() );

        meshEnd[dir]    = ( par ? breaks[dir].end() - 1 : breaks[dir].begin() + 1 );
        curElement[dir] =
        meshBegin[dir]  = ( par ? breaks[dir].end() - 2 : breaks[dir].begin()     );
        tindex = curElement[dir] - breaks[dir].begin();
        //meshEnd[dir]    = ( par ? m_basis->component(dir).knots().uend() - 1   : 
        //                          m_basis->component(dir).knots().ubegin() + 1 );
        //curElement[dir] =
        //meshBegin[dir]  = ( par ? 
        //                    m_basis->component(dir).knots().uend() - 2 : 
        //                    m_basis->component(dir).knots().ubegin()   );
        //tindex = curElement[dir] - m_basis->component(i).knots().ubegin();

        for (int i=dir+1; i < d; ++i) 
        {
            breaks.push_back( m_basis->component(i).domain()->breaks() );
            meshEnd[i]   = breaks[i].end() - 1;
            meshBegin[i] = curElement[i] = breaks[i].begin();
            //meshEnd[i]    = m_basis->component(i).knots().uend() - 1;
            //meshBegin[i]  = 
            //curElement[i] = m_basis->component(i).knots().ubegin();

            if (meshEnd[i] == curElement[i])
                m_isGood = false;
        }

        // Set to one quadrature point by default
        m_quadrature.setNodes( gsVector<index_t>::Ones(d) );

        if (m_isGood)
            update();
    }

    // ---> Documentation in gsDomainIterator.h
    // proceed to the next element; returns true if end not reached yet
    bool next()
    {
        m_isGood = m_isGood && nextLexicographic(curElement, meshBegin, meshEnd);
        if (m_isGood)
            update();
        return m_isGood;
    }

    // ---> Documentation in gsDomainIterator.h
    // proceed to the next element (skipping #increment elements);
    // returns true if end not reached yet
    bool next(index_t increment)
    {
        for (index_t i = 0; i < increment; i++)
            m_isGood = m_isGood && nextLexicographic(curElement, meshBegin, meshEnd);
        if (m_isGood)
            update();
        return m_isGood;
    }

    /// Resets the iterator implementation copied from the constructor
    /// note that it fails for sides containing 1 element in any direction
    /// do not know the rationale for it
    void reset()
    {
        curElement=meshBegin;
        m_isGood = true;
        for(int i=0; i < d; ++i)
        {
            if (i!=dir && curElement[i]==meshEnd[i])
                m_isGood=false;
        }
        if (m_isGood)
            update();
    }

    /// Return the tensor index of the current element
    gsVector<unsigned, D> index() const
    {
        gsVector<unsigned, D> curr_index(d);  
        for (int i = 0; i < dir; ++i)
            curr_index[i]  = curElement[i] - meshBegin[i];
        for (int i = dir+1; i < d; ++i)
            curr_index[i]  = curElement[i] - meshBegin[i];
        curr_index[dir]  = tindex;
        return curr_index; 
    }

    const gsVector<T> & lowerCorner() const
    { return lower; }

    const gsVector<T> & upperCorner() const
    { return upper; }

    const T getPerpendicularCellSize() const
    {
        return *(curElement[dir]+1) - *curElement[dir];
    }

    /// Returns the number of elements.
    size_t numElements() const
    {
        size_t result = 1;
        for (short_t i = 0; i < dir; ++i)
            result *= breaks[i].size() - 1;
        for (short_t i = dir+1; i < d; ++i)
            result *= breaks[i].size() - 1;
        
        return result;
    }

    void adjacent( const gsVector<bool> & orient, 
                   gsDomainIterator<T>  & other )
    {
        // 2D only for now

        gsTensorDomainBoundaryIterator & other_ = 
            static_cast< gsTensorDomainBoundaryIterator &>(other);

        int a1 = !dir;
        int a2 = !other_.dir;

        other_.curElement[a2] = std::lower_bound( 
            other_.breaks[a2].begin(), other_.breaks[a2].end(), 
            orient[0] ? *curElement[a1] : *(curElement[a1]+1) );
        other_.update();
    }

    /// Function to set the breakpoints in direction \a i manually
    void setBreaks(std::vector<T> newBreaks, index_t i) // i: direction
    {
        breaks[i].swap(newBreaks);
        meshEnd[i]   = breaks[i].end() - 1;
        meshBegin[i] = curElement[i] = breaks[i].begin();
        reset();
    }

private:

    /// Computes lower, upper and center point of the current element, maps the reference
    /// quadrature nodes and weights to the current element, and computes the
    /// active functions. Plus some additional, boundary-iterator-specific things.
    void update()
    {
        for (int i = 0; i < dir ; ++i)
        {
            lower[i]  = *curElement[i];
            upper[i]  = *(curElement[i]+1);
            center[i] = T(0.5) * (lower[i] + upper[i]);
        }
        lower[dir]  = 
        upper[dir]  =
        center[dir] = (par ? *(curElement[dir]+1) : *curElement[dir] );
        for (int i = dir+1; i < d; ++i)
        {
            lower[i]  = *curElement[i];
            upper[i]  = *(curElement[i]+1);
            center[i] = T(0.5) * (lower[i] + upper[i]);
        }

        //gsDebug<<"lower: "<< lower.transpose() <<", upper="<<upper.transpose() <<"\n";
    }

// Data members
protected:
    using gsDomainIterator<T>::m_basis;
    using gsDomainIterator<T>::m_isGood;
    using gsDomainIterator<T>::center;
    using gsDomainIterator<T>::m_side;

private:

    // the dimension of the parameter space
    short_t d;

    // Boundary parameters
    short_t  dir;
    bool par;
    unsigned tindex;

    // coordinates of the grid cell boundaries
    std::vector< std::vector<T> > breaks;

    // Quadrature rule
    gsGaussRule<T> m_quadrature;

    // First mesh-line on the tensor grid
    gsVector<uiter, D> meshBegin;

    // Last mesh-line on the tensor grid
    gsVector<uiter, D> meshEnd;

    // Current element as pointers to it's supporting mesh-lines
    gsVector<uiter, D> curElement;

    // parameter coordinates of current grid cell
    gsVector<T> lower, upper;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
}; // class gsTensorDomainBoundaryIterator


} // namespace gismo
