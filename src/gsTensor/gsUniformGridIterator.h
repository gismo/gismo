/** @file gsUniformGridIterator.h

    @brief Provides iteration over points in a uniform tensor grid

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsUtils/gsCombinatorics.h>


namespace gismo
{

/** 
    Iterator over the points of a tensor-product point grid.
*/
template<class T, unsigned d>
class gsUniformGridIterator
{
private:
    gsUniformGridIterator()  { }

public:
    
    gsUniformGridIterator(gsVector<T> const & a, 
                          gsVector<T> const & b, 
                          gsVector<unsigned, d> const & np)
    : m_low(a), m_upp(b), m_numPoints(np)
    {
        GISMO_ASSERT( (m_numPoints.array()>0).all(), "Number of points is zero.");

        m_step.array() = (b-a).array() / (m_numPoints.array() - 1)
            .matrix().cwiseMax(1).template cast<T>().array() ;
        m_size = m_numPoints.prod();
        reset();
    }

    gsUniformGridIterator(gsMatrix<T> const & ab, 
                          gsVector<unsigned, d> const & np)
    : m_low(ab.col(0)), m_upp(ab.col(1)), m_numPoints(np)
    {
        GISMO_ASSERT( (m_numPoints.array()>0).all(), "Number of points is zero.");

        m_step.array() = (ab.col(1)-ab.col(0)).array() / (m_numPoints.array() - 1)
            .matrix().cwiseMax(1).template cast<T>().array() ;
        m_size = m_numPoints.prod();
        reset();
    }

    gsUniformGridIterator(gsMatrix<T> const & ab, unsigned numPoints)
    : m_low(ab.col(0)), m_upp(ab.col(1))
    {
        // deduce number of points per direction
        const gsVector<T,d> span = m_upp - m_low;
        const gsVector<T,d> wght = span / span.sum();
        const T h = math::pow( span.prod() / numPoints, 1.0 / d);
        for (index_t i = 0; i != m_low.size(); ++i)
            m_numPoints[i] = cast<T,unsigned>(math::ceil( span[i] / (h*wght[i]) ) );
        GISMO_ASSERT( (m_numPoints.array()>0).all(), "Number of points is zero.");

        m_step.array() = (ab.col(1)-ab.col(0)).array() / (m_numPoints.array() - 1)
            .matrix().cwiseMax(1).template cast<T>().array() ;
        m_size = m_numPoints.prod();
        reset();
    }
        

public:

    operator bool() const {return m_count < m_size;}
    
    gsUniformGridIterator & operator++()
    {
        if ( nextLexicographic(currIndex, m_numPoints) )
            update();
        return *this;        
    }
    
    const gsMatrix<T> & operator*() const {return currPoint;}

    const gsMatrix<T> * operator->() const {return &currPoint;}

    const gsMatrix<T> & point() const {return *this;}

    const gsVector<unsigned> & tensorIndex() const { return currIndex;}
    
    unsigned flatIndex() const { return m_count;  }

    const gsVector<T,d> & step() const {return m_step;}
    
    unsigned size() const {return m_size;}

    const gsVector<unsigned,d> & numPoints() const {return m_numPoints;}

    void reset()
    {
        currIndex.setZero(d);
        currPoint = m_low;
        m_count = 0;
    }
    
private:

    void update()
    {
        ++m_count;
        for (unsigned i = 0; i < d; ++i)
        {
            if ( currIndex[i] == 0 )
                currPoint.at(i) = m_low[i]; // avoid numerical error at first val
            else if ( currIndex[i] == m_numPoints[i]-1 )
                currPoint.at(i) = m_upp[i]; // avoid numerical error at last val
            else
                currPoint.at(i) = m_low[i] + currIndex[i] * m_step[i];
        }
    }

private:

    gsVector<T,d>  m_low;
    gsVector<T,d>  m_upp;
    gsVector<T,d> m_step;
    
    gsVector<unsigned,d> m_numPoints;
    unsigned m_size;

    gsVector<unsigned,d> currIndex;
    gsMatrix<T>          currPoint;
    unsigned m_count;
};



}; // namespace gismo
