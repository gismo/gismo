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

// note: default arguments are found in gsForwardDeclarations.h
template<typename Z, int d, int mode, bool> class gsGridIterator { };


/** 
    \brief Iterator over the integer points of a tensor-product point-grid.
    
    The iteration is done in lexicographic order.

    - mode = 0 : iteration over [a, b) or [a, b]
    - mode = 1 : iteration over the boundry points of [a, b) or [a, b]
    - mode = 2 : iteration over the vertices of [a, b) or [a, b]

    The open or closed case is determined by a constructor flag

    \note Iteration over the boundary including offsets is possible
    using the free functions in gsUtils/gsCombinatorics.h

    \tparam Z type of the integer coordinates of the index vector

    \tparam d statically known dimension, or dynamic dimension if d =
    -1 (default value)

    \tparam mode 0: all points in [a,b], 1: all points in [a,b), 2:
    vertices of cube [a,b], 0: boundary of cube [a,b],

    \ingroup Tensor
*/
template<class Z, int d, int mode>
class gsGridIterator<Z,d,mode,true>
{
public:
    typedef gsVector<Z,d> point;

public:

    gsGridIterator()
    { GISMO_STATIC_ASSERT(std::numeric_limits<Z>::is_integer,INCONSISTENT_INSTANTIZATION); }
            
    gsGridIterator(point const & a, point const & b, bool open = true)
    { reset(a, b, open); }

    gsGridIterator(point const & b, bool open = true)
    { reset(point::Zero(b.size()), b, open); }

    gsGridIterator(gsMatrix<Z,d,2> const & ab, bool open = true)
    { reset(ab.col(0), ab.col(1), open); }

    inline void reset(point const & a, point const & b, bool open = true)
    {
        GISMO_ASSERT(a.rows() == b.rows(), "Invalid endpoint dimensions");
        m_low = m_cur = a;
        if (open) m_upp = b.array() - 1; else m_upp = b;
        m_dim = ( (m_low.array() <= m_upp.array()).all() ? a.rows() : 0 );
    }

    void reset() { reset(m_low,m_upp, false); }

    // See http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF( (sizeof(point)%16)==0 );

public:

    operator bool() const {return 0!=m_dim;}
        
    const point & operator*() const {return m_cur;}

    const point * operator->() const {return &m_cur;}
    
    inline gsGridIterator & operator++()
    {
        // Multiple implementation according to mode
        switch (mode)
        {
        case 0: // ----------- Iteration over [m_low, m_upp)]
        case 2: // iteration over vertices of the cube [m_low, m_upp)]
            for (int i = 0; i != (d==-1?m_dim:d); ++i)
                if ( m_cur[i] != m_upp[i] )
                {
                    if (0==mode) ++m_cur[i]; else m_cur[i] = m_upp[i];
                    return *this;
                }
                else
                    m_cur[i] = m_low[i];
            m_dim = 0;//done
            return *this;

        case 1: // ----------- Iteration over boundary of [m_low, m_upp)]
            for (int i = 0; i != (d==-1?m_dim:d); ++i)
            {        
                if ( m_cur[i] != m_upp[i] )
                {
                    if ( m_cur[i] == m_low[i] && ( i+1!=m_dim || m_dim==1) )
                    {
                        int c = i+1;
                        for (int j = c; j!=(d==-1?m_dim:d); ++j)
                            if ( (m_cur[j] == m_low[j]) || 
                                 (m_cur[j] == m_upp[j]) )
                                ++c;
                        
                        if ( c==1 )
                            m_cur[i] = m_upp[i];
                        else
                            ++m_cur[i];
                    }
                    else
                        m_cur[i]++;
                    
                    for (int k = i-1; k!=-1; --k)
                        m_cur[k] = m_low[k];
                    return *this;
                }
            }
            m_dim = 0;//done
            return *this;

        default:
            GISMO_STATIC_ASSERT(mode > -1 && mode<3, INCONSISTENT_INSTANTIZATION);
        }
    }

    inline bool isFloor(int i) const { return m_cur[i] == m_low[i];}

    inline bool isCeil (int i) const
    { return mode ? m_cur[i] == m_upp[i] : m_cur[i] + 1 == m_upp[i];}

    bool isBoundary() const
    {
        if ( mode )
            return (m_cur.array() == m_low.array()).any() ||
                   (m_cur.array() == m_upp.array()).any() ;
        else
            return (m_cur.array()     == m_low.array()).any() ||
                   (m_cur.array() + 1 == m_upp.array()).any() ;
    }
    
    point numPoints() const {return m_upp - m_low;}

private:
    
    point m_low, m_upp;
    
    point  m_cur;

    int m_dim;
};


/** 
    \brief Iterator over uniformly distributed numeric points inside a
    (hyper-)cube.
*/
template<class T, int d>
class gsGridIterator<T,d,0,false>
{   // note: the iteration is over the volume of [a,b]
public:

    typedef gsVector<T,d> point;

    typedef gsGridIterator<index_t, d, 0> integer_iterator;

    typedef typename integer_iterator::point point_index;
public:
    
    gsGridIterator(point const & a, 
                   point const & b, 
                   point_index const & np)
    : m_iter(np, 1)
    {
        reset(a, b);
    }

    gsGridIterator(gsMatrix<T,2> const & ab, 
                   point_index const & np)
    : m_iter(np, 1)
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
        m_iter = integer_iterator(npts, 1);

        reset(ab.col(0), ab.col(1));
    }

    void reset() { m_cur = m_low; m_iter.reset(); }

    void reset(point const & a, 
               point const & b)
    {
        GISMO_ASSERT( (m_iter.numPoints().array()>0).all(), "Number of points is zero.");
        m_cur = m_low = a;
        m_upp = b;
        m_step = (b-a).array() / (m_iter.numPoints().array() - 1)
            .matrix().cwiseMax(1).template cast<T>().array() ;
        m_iter.reset();
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
                else
                {
                    if ( m_iter.isCeil(i) )
                        m_cur.at(i) = m_upp[i]; // avoid numerical error at last val
                    else
                        m_cur.at(i) = m_low[i] + m_iter->at(i) * m_step[i];
                }
        
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

    inline bool isBoundary() const { return m_iter.onBoundary();}
    
    const integer_iterator & index_iterator() const { return m_iter;}
    
private:

    point  m_low, m_upp, m_step;

    integer_iterator m_iter;    
    gsMatrix<T>      m_cur;
};


} // namespace gismo
