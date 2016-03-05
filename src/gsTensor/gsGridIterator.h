/** @file gsIntegerGridIterator.h

    @brief Provides iteration over integer or numeric points in a (hyper-)cube

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once


namespace gismo
{

/** 
    Iterator over the integer points of a tensor-product point-grid.

    The iteration is done in lexicographic order.

    closed == true   :  iteration over [a, b]
    closed == false  :  iteration over [a, b)
*/
template<class Z, int d, bool closed>
class gsGridIterator<Z,d,closed,true>
{
public:
    typedef gsVector<Z,d> point;

public:

    gsGridIterator() { }
            
    gsGridIterator(point const & a, point const & b)
    { reset(a,b); }

    gsGridIterator(point const & b)
    { reset(point::Zero(b.size()), b); }

    gsGridIterator(gsMatrix<Z,d,2> const & ab)
    { reset(ab.col(0), ab.col(1) ); }

    inline void reset(point const & a, point const & b)
    {
        GISMO_ASSERT(a.rows() == b.rows(), "Invalid endpoint dimensions");
        m_low = m_cur = a;
        m_upp = b;
        m_dim = closed ?
            ( (m_low.array() <= m_upp.array()).all() ? a.rows() : 0 ) : 
            ( (m_low.array() <  m_upp.array()).all() ? a.rows() : 0 ) ;
    }

    void reset() { reset(m_low,m_upp); }

    // See http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF( (sizeof(point)%16)==0 );

public:

    operator bool() const {return 0!=m_dim;}
        
    const point & operator*() const {return m_cur;}

    const point * operator->() const {return &m_cur;}
    
    inline gsGridIterator & operator++()
    {
        if (closed)  // Iteration over [m_low, m_upp]
        {
            for (int i = 0; i != (d==-1?m_dim:d); ++i)
                if ( m_cur[i] != m_upp[i] )
                {
                    ++m_cur[i];
                    return *this;
                }
                else
                    m_cur[i] = m_low[i];
            m_dim = 0;//done
            return *this;
        }
        else // Iteration over [m_low, m_upp)
        {
            for (int i = 0; i != (d==-1?m_dim:d); ++i)
                if (++m_cur[i] == m_upp[i])
                {
                    if (i + 1 == m_dim)
                    {
                        m_dim = 0;//done
                        return *this;
                    }
                    else
                        m_cur[i] = m_low[i];
                }
                else
                    return *this;
            GISMO_ERROR("Fatal error in gsGridIterator.");
        }
        // else mode==2 (todo: boundary iterator)
    }

    inline bool isFloor(int i) const { return m_cur[i] == m_low[i];}

    inline bool isCeil (int i) const
    { return closed ? m_cur[i] == m_upp[i] : m_cur[i] + 1 == m_upp[i];}

    point numPoints() const {return m_upp - m_low;}

private:
    
    point m_low, m_upp;
    
    point  m_cur;

    int m_dim;
};


/** 
    Iterator over uniformly distributed numeric points inside a
    (hyper-)cube.
*/
template<class T, int d>
class gsGridIterator<T,d,1,false>
{   // note: closed = true
public:

    typedef gsVector<T,d> point;
    
    typedef gsGridIterator<index_t, d, false> integer_iterator;

    typedef typename integer_iterator::point point_index;
public:
    
    gsGridIterator(point const & a, 
                   point const & b, 
                   point_index const & np)
    : m_iter(np)
    {
        reset(a, b);
    }

    gsGridIterator(gsMatrix<T,2> const & ab, 
                   point_index const & np)
    : m_iter(np)
    {
        reset(ab.col(0), ab.col(1));
    }

    gsGridIterator(gsMatrix<T,2> const & ab, unsigned numPoints)
    : m_low(ab.col(0)), m_upp(ab.col(1))
    {
        // deduce the number of points per direction
        const point span = ab.col(1) - ab.col(0);
        const point wght = span / span.sum();
        const T h = math::pow( span.prod() / numPoints, 1.0 / (d!=-1?d:ab.rows()) );
        point_index npts(ab.rows());
        for (index_t i = 0; i != (d!=-1?d:ab.rows()); ++i)
            npts[i] = cast<T,index_t>(math::ceil( span[i] / (h*wght[i]) ) );
        m_iter = integer_iterator(npts);

        reset(ab.col(0), ab.col(1));
    }

    void reset() { m_cur = m_low;}

    void reset(point const & a, 
               point const & b)
    {
        GISMO_ASSERT( (m_iter.numPoints().array()>0).all(), "Number of points is zero.");
        m_cur = m_low = a;
        m_upp = b;
        m_step = (b-a).array() / (m_iter.numPoints().array() - 1)
            .matrix().cwiseMax(1).template cast<T>().array() ;

    }

    // See http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF( (sizeof(point)%16)==0 );

public:

    operator bool() const {return m_iter;}
    
    inline gsGridIterator & operator++()
    {
        if ( ++m_iter )
            for (int i = 0; i != (d==-1?m_low.size():d); ++i)
                if ( m_iter.isFloor(i) )
                    m_cur.at(i) = m_low[i]; // avoid numerical error at first val
                else if ( m_iter.isCeil(i) )
                    m_cur.at(i) = m_upp[i]; // avoid numerical error at last val
                else
                    m_cur.at(i) = m_low[i] + m_iter->at(i) * m_step[i];
        return *this;        
    }
    
    const gsMatrix<T> & operator*() const {return m_cur;}

    const gsMatrix<T> * operator->() const {return &m_cur;}

    const point_index & tensorIndex() const { return *m_iter;}
    
    const point & step() const {return m_step;}
    
    index_t size() const {return numPoints().prod();}

    point_index numPoints() const {return m_iter.numPoints();}

    inline bool isFloor(int i) const { return m_iter.isFloor(i);}
    
    inline bool isCeil (int i) const { return m_iter.isCeil(i);}

    const integer_iterator & index_iterator() const { return m_iter;}
    
private:

    point  m_low, m_upp, m_step;

    integer_iterator m_iter;    
    gsMatrix<T>      m_cur;
};


} // namespace gismo
