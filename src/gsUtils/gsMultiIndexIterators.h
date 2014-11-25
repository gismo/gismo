/*
* gsMultiIndexIterarors.h created on 09.07.2014
*
* This file contains iterators over discrete set of points in Z^n.
* In particular it contains a tensorGrid and a simplexGrid iterators.
*
* The content was in gsCombinatorics
*
* Author: Andrea Bressan
*
* This file is part of the G+SMO library
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>

namespace gismo {



/**
    \brief interface for random access iterators on multi-indexes

    they iterate on structured subset of f$\mathbb{Z}^nf$.
    the standard loop structure for an iterator it is

    for (it.first() ; it.good(); it.next())
    {
        do something ...
    }

    or backward

    for (it.last() ; it.good(); it.previous())
    {
        do something ...
    }

    or with predetermined endpoints

    for (it.setTo(startingIndex) ;
            it.good() && it.multiIndex()!=lastIndex;
            it.next())
    {
        do something ...
    }

    \tparam Flat is the type of the scalar index
    \tparam dim  is the dimension of the multi-index

    \ingroup combinatorics
**/

template<typename Flat=index_t, int dim=-1>
class gsMultiIndexIterator
{
public:
    typedef gsVector<Flat,dim> MIndexT;
    typedef Flat               FIndexT;
    enum {
        Tdim=dim,
        d=dim
        };
public:
    gsMultiIndexIterator(int _dim)
        : m_good(false)
    {}
    virtual ~gsMultiIndexIterator(){}
    /**
        \brief tells if the iterator is pointing to a valid multi-index
    **/
    bool               good() {return m_good;}

    /**
        \brief reset the iterator to the first position
    **/
    bool               reset() {this->first(); return good();}
    /**
        \brief returns the coordinate in Z^n of the current position
    **/
    const   MIndexT & multiIndex() const
    {
        GISMO_ASSERT(m_good,"do not ask for indexes without asking for validity");
        return m_vpos;
    }
    /**
        \brief advance the iterator by 1 position
    **/
    virtual bool      next() = 0;
    /**
        \brief decrease the iterator by 1 position
    **/
    virtual bool      previous() = 0;

    /**
        \brief reposition the iterator to the first element
        \return true if the first element exists, false if the set is empty
    **/
    virtual bool       first()=0;
    /**
        \brief reposition the iterator to the last element
        \return true if the last element exists, false if the set is empty
    **/
    virtual bool       last  ()=0;

    /**
        \brief returns the flat index of the current position
    **/
    virtual Flat       flatIndex() const = 0;
    /**
        \brief move the current position to the supplied point in Z^n

        \return true if the supplied position is lecit and the operation
                was performend
    **/
    virtual bool       setTo (const MIndexT & here) {m_good=false; return m_good;}
    /**
        \brief move the current position to the supplied flat index

        \return true if the supplied index is lecit and the operation
                was performend
    **/
    virtual bool       setTo (Flat here)   {m_good=false; return m_good;}

    /**
        \brief convert a multi-index to a flat index
      **/
    virtual Flat    multiToFlat (MIndexT const & index) const = 0;

    /**
        \brief convert a flat index to a multi-index
    **/
    virtual MIndexT flatToMulti (Flat index) const = 0;

    /**
        \brief check validity of a multi-index
    **/
    virtual bool   multiIsValid (MIndexT const & index) const = 0;
    /**
        \brief check validity of a multi-index
    **/
    virtual bool   flatIsValid (Flat index) const = 0;


    // See http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF( (sizeof(MIndexT)%16)==0 )
    //EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:
    /// iterator validity, set to false when passing the last element
    /// and restored only on setTo, first, last or reset functions if
    /// the iterator is not empty
    bool    m_good;
    /// current multiIndex
    MIndexT m_vpos;
};



// partial specialization for runtime dimension
template<typename Flat>
class gsMultiIndexIterator<Flat,-1>
{
public:
    typedef gsVector<Flat,-1>  MIndexT;
    typedef Flat               FIndexT;
    enum {
        Tdim=-1
        };
public:
    gsMultiIndexIterator(int dim)
        : d(dim), m_good(false)
    {
        m_vpos.resize(d);
    }

    virtual ~gsMultiIndexIterator(){}

    bool               good() {return m_good;}
    bool               reset() {this->first(); return good();}
    const   MIndexT & multiIndex() const
    {
        GISMO_ASSERT(m_good,"do not ask for indexes without asking for validity");
        return m_vpos;
    }
    virtual Flat       flatIndex() const = 0;

    virtual bool       next() = 0;
    virtual bool       previous() = 0;

    virtual bool       first()=0;
    virtual bool       last  ()=0;
    virtual bool       setTo (const MIndexT & here) {m_good=false; return m_good;}
    virtual bool       setTo (Flat here)   {m_good=false; return m_good;}

    virtual Flat       multiToFlat (MIndexT const & index) const = 0;
    virtual MIndexT    flatToMulti (Flat index) const = 0;

    virtual bool       multiIsValid (MIndexT const & index) const = 0;
    virtual bool       flatIsValid (Flat index) const = 0;


    // See http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF( (sizeof(MIndexT)%16)==0 )
    //EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:
    /// dimension
    const index_t d;
    /// bool to represent iterator validity, true when the iterator
    /// points to a new value
    bool    m_good;
    /// current multiIndex
    MIndexT m_vpos;
};


template<typename Flat, int d> class gsTensorGridVertexIterator;
template<typename Flat, int d> class gsTensorGridBoundaryIterator;


/**
    \brief iterates on tensor grids

    With a tensor grid we mean a subset of Z^n containing the
    points f$(a_0,...,a_n)f$ such that f$m_i\leq a_i< M_if$.

    The iteration is in lexicographic order.
    The multi-index is the vector of point coordinates and the
    flat index is the position during iteration.

  \ingroup combinatorics
**/
template<typename Flat=index_t, int dim=-1>
class gsTensorGridIterator : public gsMultiIndexIterator<Flat,dim>
{
public:
    typedef typename gsMultiIndexIterator<Flat,dim>::MIndexT MIndexT;
    typedef typename gsMultiIndexIterator<Flat,dim>::FIndexT FIndexT;
    using gsMultiIndexIterator<Flat,dim>::good;

public:
    /**
     * \brief constructor form the maximum of the coodinates
     *
     * iterates over f$(a_0,...,a_n)f$ such that f$0\leq a_i\leq M_if$.
     * \param M
     */
    gsTensorGridIterator(MIndexT const & M);
    /**
     * \brief construct an iterator by setting the minimum and
     *
     * iterates over f$(a_0,...,a_n)f$ such that f$m_i\leq a_i\leq M_if$.
     * \param m
     * \param M
     */
    gsTensorGridIterator(MIndexT const & m, MIndexT const & M);

public:


    /**
     * \brief construct an iterator by setting lower and upper corners
     *
     * Iterates over f$(a_0,...,a_n)f$ such that f$m_i\leq a_i< M_if$.
     * the returned flat indexes of the SubGridIterator coincide with those
     * of the original iterators on the same coordinate.
     * \param m
     * \param M
     */
    virtual gsTensorGridIterator<Flat,dim>
    * makeSubGridIterator (MIndexT const & m, MIndexT const & M);
    /**
     * \brief construct a boundary iterator by setting the thickness
     *
     * Iterates over f$(a_0,...,a_n)f$ such that f$m_i\leq a_i<m_i+mOffset_i f$
     * or f$M_i-mOffset_i\leq a_i < M_i f$.
     * the returned flat indexes of the boundary iterator coincide with those
     * of the original iterators on the same coordinate.
     * \param m
     * \param M
     */
    virtual gsTensorGridIterator<Flat,dim>
    * makeBoundaryIterator (MIndexT const & mOffset, MIndexT const & MOffset);

    /**
     * \brief construct an iterator
     *
     * Iterates over the vertex.
    **/
    virtual gsTensorGridIterator<Flat,dim>
    * makeVertexIterator ();
public: // gsMultiIndexIterator interface
    virtual bool       previous     ();
    virtual bool       next         ();
    virtual bool       first        ();
    virtual bool       last         ();
    virtual bool       setTo        (MIndexT const & here);
    virtual bool       setTo        (Flat    here);
    virtual Flat       multiToFlat  (MIndexT const & index)     const;
    virtual MIndexT    flatToMulti  (Flat    index)     const;
    virtual bool       multiIsValid (MIndexT const & index)     const;
    virtual bool       flatIsValid  (Flat    index)     const;
    virtual Flat       flatIndex    ()                  const;

    // See http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
    //EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF( (sizeof(MIndexT)%16)==0 )
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    //EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:

    using gsMultiIndexIterator<Flat,dim>::d;
    using gsMultiIndexIterator<Flat,dim>::m_good;
    using gsMultiIndexIterator<Flat,dim>::m_vpos;

    /// lower corner, included in the range
    MIndexT       m_vmin;
    /// upper corner, out of the range m_vmax-1 is the last index
    MIndexT       m_vmax;
    /// dimension of the grid, used to convert from scalar to multi-index
    MIndexT       m_vinc;
    /// minimum scalar index
    Flat          m_fmin;
    /// maximum scalar index
    Flat          m_fmax;
};


/**
    \brief provides iterators on the vertexes of a tensor grid

    Note that make subgridIterator and makeBoundaryIterator
    returns iterators to the original tensorGrid

  \ingroup combinatorics
**/
template<typename Flat=index_t, int dim=-1>
class gsTensorGridVertexIterator
        : public gsTensorGridIterator<Flat,dim>
{
public:
    using typename gsMultiIndexIterator<Flat,dim>::MIndexT;
    using typename gsMultiIndexIterator<Flat,dim>::FIndexT;
    using gsMultiIndexIterator<Flat,dim>::good;
public:
    gsTensorGridVertexIterator( gsTensorGridIterator<Flat,dim> &area);
    gsTensorGridVertexIterator( MIndexT const & lower, MIndexT const & upper);
    virtual bool next ();
    virtual bool previous();
    virtual bool multiIsValid (MIndexT const & index)  const;
protected:
    using  gsMultiIndexIterator<Flat,dim>::d;
    using  gsMultiIndexIterator<Flat,dim>::m_good;

    using  gsTensorGridIterator<Flat,dim>::m_vmin;
    using  gsTensorGridIterator<Flat,dim>::m_vmax;
    using  gsTensorGridIterator<Flat,dim>::m_vpos;
};

/**
    \brief provides iterators on the surface of a tensor grid

    Note that makeSubgridIterator and makeBoundaryIterator
    returns iterators to the original tensorGrid

  \ingroup combinatorics
**/
template<typename Flat=index_t, int dim=-1>
class gsTensorGridBoundaryIterator
        : public gsTensorGridIterator<Flat,dim>
{
public:
    typedef typename gsMultiIndexIterator<Flat,dim>::MIndexT MIndexT;
    typedef typename gsMultiIndexIterator<Flat,dim>::FIndexT FIndexT;
    using gsMultiIndexIterator<Flat,dim>::good;
public:
    /**
     * \brief constructor form the maximum of the coodinates
     *
     * iterates over f$(a_0,...,a_n)f$ such that f$0\leq a_i\leq M_if$.
     * \param M
     */
    gsTensorGridBoundaryIterator(
            gsTensorGridIterator<Flat,dim> &area,
            MIndexT                    const & mOffset,
            MIndexT                    const & MOffset);

    virtual bool       first();
    virtual bool       last();
    virtual bool       previous();
    virtual bool       next();
    virtual bool       multiIsValid (MIndexT const & index)  const;

    // See http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
    // EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF( (sizeof(MIndexT)%16)==0 )
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    //EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:

    using  gsMultiIndexIterator<Flat,dim>::d;
    using  gsMultiIndexIterator<Flat,dim>::m_good;

    using  gsTensorGridIterator<Flat,dim>::m_vmin;
    using  gsTensorGridIterator<Flat,dim>::m_vmax;
    using  gsTensorGridIterator<Flat,dim>::m_vinc;
    using  gsTensorGridIterator<Flat,dim>::m_vpos;
    using  gsTensorGridIterator<Flat,dim>::m_fmin;
    using  gsTensorGridIterator<Flat,dim>::m_fmax;

    /// thickness of the lower faces orthogonal to all coordinates
    MIndexT        m_lowOff;
    /// thickness of the upper faces orthogonal to all coordinates
    MIndexT        m_uppOff;
    /// index of the first coordinate in which either m_lowOff or m_uppOff
    /// is not 0
    int            m_mne; // minimum not empty
};

/**
    \brief iterates over the compositions of \a sum as the sum of dim
    different numbers.

    If the dimension is specified in the template parameter then the
    dimension argument in the constructor is ignored.

    \ingroup combinatorics
**/
template<typename Flat=index_t, int dim=-1>
class gsCompositionIterator
        : public gsMultiIndexIterator<Flat,dim>
{
public:
    typedef typename gsMultiIndexIterator<Flat,dim>::MIndexT MIndexT;
    typedef typename gsMultiIndexIterator<Flat,dim>::FIndexT FIndexT;
    using gsMultiIndexIterator<Flat,dim>::good;
public:
    gsCompositionIterator(Flat sum, int _dim=dim);


    virtual bool       previous();
    virtual bool       next();
    virtual bool       first();
    virtual bool       last();
    virtual Flat       flatIndex()                     const;
    virtual bool       setTo (MIndexT const & here);
    virtual bool       setTo (Flat       here);
    virtual Flat       multiToFlat (MIndexT const & index)     const;
    virtual MIndexT    flatToMulti (Flat index)        const;
    virtual bool       multiIsValid (MIndexT const & index)    const;
    virtual bool       flatIsValid (Flat index)        const;

    // See http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
    //EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF( (sizeof(MIndexT)%16)==0 )
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    //EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:
    /// first multi-index
    MIndexT m_vmin;
    /// last multi-index
    MIndexT m_vmax;
    /// minimum scalar-index
    Flat    m_fmin;
    /// maximum scalar-index
    Flat    m_fmax;
    /// current scalar position, saved because it is hard to compute
    Flat    m_fpos;
    /// sum of the coordinates
    Flat    m_fsum;
    using  gsMultiIndexIterator<Flat,dim>::m_good;
    using  gsMultiIndexIterator<Flat,dim>::d;
    using  gsMultiIndexIterator<Flat,dim>::m_vpos;
};


/**
    \brief iterates over the simplex

    returns all points f$(a_0,\dots,a_{n-1})f$ in f$Z^nf$ such that
    f$\sum_{i=0}^{n-1}a_i \leq k f$

  \ingroup combinatorics
**/

template<typename Flat=index_t, int dim=-1>
class gsSimplexIterator
        : public gsMultiIndexIterator<Flat,dim>
{
public:
    typedef typename gsMultiIndexIterator<Flat,dim>::MIndexT MIndexT;
    typedef typename gsMultiIndexIterator<Flat,dim>::FIndexT FIndexT;
    using gsMultiIndexIterator<Flat,dim>::good;
public:
    gsSimplexIterator(Flat k, int _dim=dim);


    virtual bool       previous();
    virtual bool       next();
    virtual bool       first();
    virtual bool       last();
    virtual Flat       flatIndex()                     const;
    virtual bool       setTo (const MIndexT & here);
    virtual bool       setTo (Flat       here);
    virtual Flat       multiToFlat (const MIndexT &index)     const;
    virtual MIndexT    flatToMulti (Flat index)        const;
    virtual bool       multiIsValid (const MIndexT &index)    const;
    virtual bool       flatIsValid (Flat index)        const;

    // See http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
    // EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF( (sizeof(MIndexT)%16)==0 )
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    //EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:

    /// minimum valid multiIndex
    MIndexT m_vmin;
    /// maximum valid multiIndex
    MIndexT m_vmax;
    /// scalar that keeps track of current sum: m_comp+m_vpos.sum()=m_fsum
    Flat    m_comp;
    /// minimum scalar index
    Flat    m_fmin;
    /// maximum scalar index
    Flat    m_fmax;
    /// current scalar position
    Flat    m_fpos;
    /// required sum
    Flat    m_fsum;
    /// members fron the multiIndexIterator, see there for description
    using  gsMultiIndexIterator<Flat,dim>::m_good;
    using  gsMultiIndexIterator<Flat,dim>::d;
    using  gsMultiIndexIterator<Flat,dim>::m_vpos;
};


} // namespace gismo


#ifndef GISMO_HEADERS_ONLY
#include GISMO_HPP_HEADER(gsMultiIndexIterators.hpp)
#endif
