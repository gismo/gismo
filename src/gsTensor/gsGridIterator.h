/** @file gsGridIterator.h

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
   \brief Specifies aliases describing the modes for gismo::gsGridIterator

   \ingroup enums
*/
enum gsGridIteratorMode
{
    CUBE   = 0, ///< Cube mode iterates over all lattice points inside a cube
    BDR    = 1, ///< Boundary mode iterates over boundary lattice points only
    VERTEX = 2, ///< Vertex mode iterates over cube vertices only
    CWISE  = 3  ///< Coordinate-wise mode iterates over a grid given by coordinate vectors
};
    
// note: default arguments are found in gsForwardDeclarations.h
template<typename Z, int mode, short_t d, bool> class gsGridIterator
{using Z::GISMO_ERROR_gsGridIterator_has_invalid_template_arguments;};

/** 
    \brief Iterator over the Cartesian product of integer points in a
    tensor-product grid.
    
    The iteration is done in lexicographic order.

    - mode = 0 : iteration over [a, b) or [a, b]
    gsGridIterator<T,CUBE> grid(a,b);

    - mode = 1 : iteration over the boundry points of [a, b) or [a, b]
    gsGridIterator<T,BDR> grid(a,b);
        
    - mode = 2 : iteration over the vertices of [a, b) or [a, b]
    gsGridIterator<T,VERTEX> grid(a,b);

    
    The open or closed case is determined by a constructor flag.
    A third argument, eg.  gsGridIterator<T,CUBE,d> grid(a,b);
    specifies fixed size dimension \a d
    
    \note Iteration over the boundary including offsets is possible
    using the free functions in gsUtils/gsCombinatorics.h

    \tparam Z type of the integer coordinates of the index vector

    \tparam d statically known dimension, or dynamic dimension if d =
    -1 (default value)

    \tparam mode 0: all integer points in [a,b], 1: all points in
    [a,b), 2: vertices of cube [a,b]

    \ingroup Tensor
*/
template<class Z, int mode, short_t d>
class gsGridIterator<Z,mode,d,true>
{
public:
    typedef gsVector<Z,d> point;

public:

    /**
       \brief Empty constructor
    */
    gsGridIterator()
    {
        GISMO_STATIC_ASSERT(std::numeric_limits<Z>::is_integer,"The template parameter needs to be an integer type.");
        GISMO_STATIC_ASSERT(mode > -1 && mode < 3, "The mode of gsGridIterator needs to be 0, 1 or 2.");
    }

    /**
       \brief Constructor using lower and upper limits

       \param a lower limit (vertex of a cube)
       \param b upper limit (vertex of a cube)
       \param open if true, the iteration is performed for the points in $[a_i,b_i)$
    */
    gsGridIterator(point const & a, point const & b, bool open = true)
    { reset(a, b, open); }

    /**
       \brief Constructor using upper limit. The iteration starts from zero.

       \param b upper limit (vertex of a cube)
       \param open if true, the iteration is performed for the points in $[0,b_i)$
    */
    explicit gsGridIterator(point const & b, bool open = true)
    { reset(point::Zero(b.size()), b, open); }

    /**
       \brief Constructor using lower and upper limits

       \param ab Matrix with two columns, corresponding to lower and upper limits
       \param open if true, the iteration is performed for the points in $[a_i,b_i)$
    */
    explicit gsGridIterator(gsMatrix<Z,d,2> const & ab, bool open = true)
    { reset(ab.col(0), ab.col(1), open); }

    /**
       \brief Resets the iterator using new lower and upper limits

       \param a lower limit (vertex of a cube)
       \param b upper limit (vertex of a cube)
       \param open if true, the iteration is performed for the points in $[a_i,b_i)$
    */    
    inline void reset(point const & a, point const & b, bool open = true)
    {
        GISMO_ASSERT(a.rows() == b.rows(), "Invalid endpoint dimensions");
        m_low = m_cur = a;
        if (open) m_upp = b.array() - 1; else m_upp = b;
        m_valid = ( (m_low.array() <= m_upp.array()).all() ? true : false );
    }

    /**
       \brief Resets the iterator, so that a new iteration over the
       points may start
    */
    void reset() { reset(m_low,m_upp, false); }

    // See http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
    //EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF( (sizeof(point)%16)==0 );

public:

    operator bool() const {return m_valid;}
    
    const point & operator*() const {return m_cur;}

    const point * operator->() const {return &m_cur;}
    
    inline gsGridIterator & operator++()
    {
        // Multiple implementations according to mode
        switch (mode)
        {
        case 0: // ----------- Iteration over [m_low, m_upp]
        case 2: // iteration over vertices of the cube [m_low, m_upp]
            for (index_t i = 0; i != m_cur.size(); ++i)
                if ( m_cur[i] != m_upp[i] )
                {
                    if (0==mode) ++m_cur[i]; else m_cur[i] = m_upp[i];
                    return *this;
                }
                else
                    m_cur[i] = m_low[i];
            m_valid = 0;//done
            return *this;

        case 1: // ----------- Iteration over boundary of [m_low, m_upp]
            for (index_t i = 0; i != m_cur.size(); ++i)
            {        
                if ( m_cur[i] != m_upp[i] )
                {
                    if ( m_cur[i] == m_low[i] && (i+1!=m_cur.size() || 1==m_cur.size()) )
                    {
                        index_t c = i+1;
                        for (index_t j = c; j!=m_cur.size(); ++j)
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

                    m_cur.head(i) = m_low.head(i);
                    return *this;
                }
            }
            /*fall through to default*/
        default:
            m_valid = 0;//done
            return *this;
        }
    }

    /**
       \brief Returns true if the \a i-th coordinate has minimal value
    */
    inline bool isFloor(int i) const { return m_cur[i] == m_low[i];}

    /**
       \brief Returns true if the \a i-th coordinate has maximal value
    */
    inline bool isCeil (int i) const
    { return mode ? m_cur[i] == m_upp[i] : m_cur[i] + 1 == m_upp[i];}

    /**
       \brief Returns true if the current point lies on a boundary
    */
    bool isBoundary() const
    {
        return (m_cur.array() == m_low.array()).any() ||
            (m_cur.array() == m_upp.array()).any() ;
    }

    /**
       \brief Returns the total number of points that are iterated
    */
    index_t numPoints() const
    {
        switch (mode)
        {
        case 0: return ((m_upp-m_low).array() + 1).prod();
        case 1: return ((m_upp-m_low).array() + 1).prod() - ((m_upp-m_low).array()-1).prod();
        case 2: return ( 1 << m_low.rows() );
        default: GISMO_ERROR("gsGridIterator::numPoints");
        }
    }

    /**
       \brief Returns the total number of points per coordinate which
       are iterated
    */
    point numPointsCwise() const { return (m_upp-m_low).array() + 1;}

    /**
       \brief Returns the first point in the iteration
    */
    const point & lower() const {return m_low;}

    /**
       \brief Returns the last point in the iteration
    */
    const point & upper() const {return m_upp;}

    /**
       \brief Utility function which returns the vector of strides

       \return An integer vector (stride vector). the entry \f$s_i\f$
       implies that the current \f$i\f$-th coordinate of the iterated
       point will appear again after \f$s_i\f$ increments.  Moreover,
       the dot product of (*this - this->lower()) with the stride
       vector results in the "flat index", i.e. the cardinal index of
       the point in the iteration sequence.

       If the iteration starts from the point zero
       (this->lower()==zero), then the stride vector has the property
       that that when we add it to *this we obtain the next point in
       the iteration sequence. In this case the flat index is obtained
       by taking the dot product with *this.
    */
    point strides() const
    {
        point res(m_low.rows());
        res[0] = 1;
        for (index_t i = 1; i != res.size(); ++i )
            res[i] = res[i-1] * ( m_upp[i-1] - m_low[i-1] + 1 );
        return res;
    }

private:
    
    point m_low, m_upp; ///< Iteration lower and upper limits
    point m_cur;        ///< Current point pointed at by the iterator
    bool  m_valid;      ///< Indicates the state of the iterator
};


/** 
    \brief Iterator over a Cartesian product of uniformly distributed
    numeric points inside a (hyper-)cube.

    The iteration is done in the natural lexicographic order.

    - mode = 0 : iteration over uniform samples in [a, b]
    gsGridIterator<T,CUBE> grid(a,b);

    - mode = 1 : iteration over uniform samples on the boundary points [a, b]
    gsGridIterator<T,BDR> grid(a,b);

    - mode = 2 : iteration over the vertices [a, b]
    gsGridIterator<T,VERTEX> grid(a,b);

    A third argument, eg.  gsGridIterator<T,CUBE,d> grid(a,b);
    specifies fixed size dimension \a d.

    \tparam T type of the numeric coordinates of the index vector

    \tparam d statically known dimension, or dynamic dimension if d =
    -1 (default value)

    \tparam mode 0: all sample points in [a,b], 1: all points in
    [a,b), 2: vertices of cube [a,b]

    \ingroup Tensor
*/
template<class T, int mode, short_t d>
class gsGridIterator<T,mode,d,false>
{
public:

    typedef gsVector<T,d> point;

    typedef gsGridIterator<index_t, mode, d> integer_iterator;

    typedef typename integer_iterator::point point_index;
public:

    /**
       \brief Constructor using lower and upper limits

       Uniformly sampled points will be generated

       \param a lower limit (vertex of a cube)
       \param b upper limit (vertex of a cube)
       \param np number of sample points per coordinate
    */
    gsGridIterator(point const & a, 
                   point const & b, 
                   point_index const & np)
    : m_iter(np, 1)
    {
        reset(a, b);
        GISMO_STATIC_ASSERT(mode > -1 && mode < 3, "The mode of gsGridIterator needs to be 0, 1 or 2.");
    }

    /**
       \brief Constructor using lower and upper limits
     
       Uniformly sampled points will be generated

       \param ab Matrix with two columns, corresponding to lower and upper limits  
       \param np number of sample points per coordinate
    */
    gsGridIterator(gsMatrix<T,d,2> const & ab, 
                   point_index const & np)
    : m_iter(np, 1)
    {
        reset(ab.col(0), ab.col(1));
    }

    /**
       \brief Constructor using lower and upper limits
       
       Uniformly sampled points will be generated. The number of the
       points is approximately \a numPoints

       \param ab Matrix with two columns, corresponding to lower and upper limits  
       \param numPoints the number sample points to be used (in total, approximately)
    */
    gsGridIterator(gsMatrix<T,d,2> const & ab, unsigned numPoints)
    : m_low(ab.col(0)), m_upp(ab.col(1))
    {
        // deduce the number of points per direction for an approximately uniform grid
        const gsVector<double,d> L = (ab.col(1) - ab.col(0)).template cast<double>();
        const double h = std::pow(L.prod()/numPoints, 1.0 / m_low.rows());
        const point_index npts = (L/h).unaryExpr((double(*)(double))std::ceil).template cast<index_t>();
        m_iter = integer_iterator(npts, 1);
        reset(ab.col(0), ab.col(1));
    }

    /**
       \brief Resets the iterator, so that a new iteration over the
       points may start
    */
    void reset() { m_cur = m_low; m_iter.reset(); }

    /**
       \brief Resets the iterator using new lower and upper limits
       
       \param a lower limit (vertex of a cube)
       \param b upper limit (vertex of a cube)
    */    
    void reset(point const & a, 
               point const & b)
    {
        m_cur = m_low = a;
        m_upp = b;
        m_step = (b-a).array() / (m_iter.numPointsCwise().array() - 1)
            .matrix().cwiseMax(1).template cast<T>().array() ;
        m_iter.reset();
    }

    gsMatrix<T> toMatrix() const
    {
        gsMatrix<T> res(m_cur.rows(), numPoints());
        integer_iterator iter = m_iter;
        iter.reset();

        for(index_t c = 0; iter; ++iter, ++c)
            update(*iter,res.col(c).data());
        
        return res;
    }
    
    // See http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
    //EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF( (sizeof(point)%16)==0 );

public:

    operator bool() const {return m_iter;}
        
    const gsMatrix<T> & operator*() const {return m_cur;}

    const gsMatrix<T> * operator->() const {return &m_cur;}

    inline gsGridIterator & operator++()
    {
        if ( ++m_iter ) update(*m_iter, m_cur.data());
        return *this;        
    }

    /**
       \brief Returns true if the \a i-th coordinate has minimal value
    */
    inline bool isFloor(int i) const { return m_iter.isFloor(i);}

    /**
       \brief Returns true if the \a i-th coordinate has maximal value
    */
    inline bool isCeil (int i) const { return m_iter.isCeil(i);}

    /**
       \brief Returns true if the current point lies on a boundary
    */
    inline bool isBoundary() const { return m_iter.isBoundary();}

    /**
       \brief Returns the total number of points that are iterated
    */
    index_t numPoints() const {return m_iter.numPoints();}
    
    /**
       \brief Returns the total number of points per coordinate which
       are iterated
    */
    point_index numPointsCwise() const {return m_iter.numPointsCwise();}

    /**
       \brief Returns the first point in the iteration
    */
    const point & lower() const {return m_low;}

    /**
       \brief Returns the last point in the iteration
    */
    const point & upper() const {return m_upp;}

    /**
       \brief Returns the tensor index of the current point
    */
    const point_index & tensorIndex() const { return *m_iter;}

    /**
       \brief Returns the \a i-th index of the current point
    */
    index_t index(const index_t i) const { return m_iter->at(i);}

    /**
       \brief Returns the coordinate-wise stepping between the samples
    */
    const point & step() const {return m_step;}

    /**
       \brief Returns a reference to the underlying integer lattice
       iterator
    */    
    const integer_iterator & index_iterator() const { return m_iter;}

    /**
       \brief Returns the point corresponding to tensor index \a ti
       
       \note Unlikely to std iterators, the position is not
       counted from the current position of \a this, but from the
       starting point of the iteration
       
       \param ti a valid tensor index of a point in the iteration sequence
    */
    inline const gsMatrix<T> operator [] (const point_index & ti) const
    {
        gsMatrix<T> res(m_low.rows(),1);
        update(ti, res.data());
        return res;
    }

private:
    
    // Update the point \a pt to the position \a ti
    inline void update(const point_index & ti, T * pt) const
    {
        for (index_t i = 0; i != m_low.size(); ++i)
            if ( ti[i] == m_iter.lower()[i] )
                *(pt++) = m_low[i]; // avoid numerical error at first val
            else if ( ti[i] == m_iter.upper()[i] )
                *(pt++) = m_upp[i]; // avoid numerical error at last val
            else
                *(pt++) = m_low[i] + ti[i] * m_step[i];
    }

private:

    point  m_low, m_upp, m_step; ///< Iteration lower and upper limits and stepsize
    
    integer_iterator m_iter;     ///< Underlying integer lattice iterator
    gsMatrix<T>      m_cur;      ///< Current point pointed at by the iterator
};


/** 
    \brief Iterator over a Cartesian product of points, which is
    given by coordinate-wise point sets

    The iteration is done in lexicographic order. Usage:

    gsGridIterator<T,CWISE> grid(a,b);

    or 

    gsGridIterator<T,CWISE,d> grid(a,b);

    \tparam T type of the numeric coordinates of the index vector

    \tparam d statically known dimension, or dynamic dimension if d =
    -1 (default value)

    \ingroup Tensor
*/
template<class T, short_t d> // mode == 3 == CWISE
class gsGridIterator<T,CWISE,d,false> 
{
public:

    typedef gsGridIterator<index_t, 0, d> integer_iterator;

    typedef typename integer_iterator::point point_index;

    typedef gsVector<const T *,d, 2> CwiseData; //non-aligned array
    
public:

    /**
       \brief Constructor using references to the coordinate vectors
       
       \param cwise A container of matrices or vectors, each
       containing the sample points in the respective coordinate
    */
    template<class CwiseContainer>
    explicit gsGridIterator(const CwiseContainer & cwise)
    : m_cwise(cwise.size()), m_cur(m_cwise.size(),1)
    {
        point_index npts(m_cwise.size());
        for (index_t i = 0; i != npts.size(); ++i)
        {
            m_cwise[i] = cwise[i].data();
            npts[i]    = cwise[i].size() - 1;
            GISMO_ASSERT(cwise[i].cols()==1 || cwise[i].rows()==1, "Invalid input");
        }
        m_iter = integer_iterator(npts, 0);
        //if (1==cwise.front().cols())
        //    m_cur.derived().resize(1, npts.size());
        //else
        //m_cur.derived().resize(npts.size(), 1);
        update(*m_iter, m_cur.data());
    }

    template<class Matrix_t>
    explicit gsGridIterator(const std::vector<Matrix_t*> & cwise)
    : m_cwise(cwise.size()), m_cur(m_cwise.size(),1)
    {
        point_index npts(m_cwise.size());
        for (index_t i = 0; i != npts.size(); ++i)
        {
            m_cwise[i] = cwise[i]->data();
            npts[i]    = cwise[i]->size() - 1;
            GISMO_ASSERT(cwise[i]->cols()==1 || cwise[i]->rows()==1, "Invalid input");
        }
        m_iter = integer_iterator(npts, 0);
        update(*m_iter, m_cur.data());
    }

    /**
       \brief Resets the iterator, so that a new iteration over the
       points may start
    */
    void reset() { m_iter.reset(); update(*m_iter, m_cur.data());}

    //void restart() { m_iter.reset(); update(*m_iter, m_cur.data());}

public:

    operator bool() const {return m_iter;}
        
    const gsMatrix<T> & operator*() const {return m_cur;}

    const gsMatrix<T> * operator->() const {return &m_cur;}
    
    inline gsGridIterator & operator++()
    {
        if (++m_iter) update(*m_iter, m_cur.data());
        return *this;
    }

    /**
       \brief Returns true if the \a i-th coordinate has minimal value
    */
    inline bool isFloor(int i) const { return m_iter.isFloor(i);}

    /**
       \brief Returns true if the \a i-th coordinate has maximal value
    */
    inline bool isCeil (int i) const { return m_iter.isCeil(i);}

    /**
       \brief Returns true if the current point lies on a boundary
    */
    inline bool isBoundary() const { return m_iter.isBoundary();}

    /**
       \brief Returns the total number of points that are iterated
    */
    index_t numPoints() const {return m_iter.numPoints();}
    
    /**
       \brief Returns the total number of points per coordinate which
       are iterated
    */
    point_index numPointsCwise() const {return m_iter.numPointsCwise();}

    /**
       \brief Returns the tensor index of the current point
    */
    const point_index & tensorIndex() const { return *m_iter;}

    /**
       \brief Returns the \a i-th index of the current point
    */
    index_t index(const index_t i) const { return m_iter->at(i);}

    /**
       \brief Returns a reference to the underlying integer lattice
       iterator
    */    
    const integer_iterator & index_iterator() const { return m_iter;}

    /**
       \brief Returns the point corresponding to tensor index \a ti
       
       \note Unlikely to std iterators, the position is not
       counted from the current position of \a this, but from the
       starting point of the iteration
       
       \param ti a valid tensor index of a point in the iteration sequence
    */
    inline const gsMatrix<T> operator [] (const point_index & ti) const
    {
        gsMatrix<T> res(m_cur.rows(),1);
        update(ti, res.data());
        return res;
    }

private:

    // Update the point \a pt to the position \a ti
    inline void update(const point_index & ti, T * pt)
    {
        for (index_t i = 0; i != m_cur.rows(); ++i)
            *(pt++) = m_cwise[i][ti[i]];
    }

private:

    CwiseData m_cwise; ///< List of coordinate-wise values
    
    integer_iterator m_iter; ///< Underlying integer lattice iterator
    gsMatrix<T>      m_cur;  ///< Current point pointed at by the iterator
};


} // namespace gismo
