/** @file gsPointIterator.h

    @brief Provides declaration of point iterators.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Bressan
*/

#pragma once

#include <gsUtils/gsMultiIndexIterators.h>

namespace gismo
{


/**
    \brief a class that wraps a multiIndex iterator and produces
    a grid of linearly spaced points

    The interFace is the same of a multiIndexIterator except that the
    current point coordinates can be queried with currPoint;

    Calling functions not in the multiIndexIteratorRandomAccess interface
    is undefined behaviour except for when the underlyingItegerIteratorType
    is gsTensorGridIterator and the dimension of the iterator and the
    target space are the same.
    In this case makeSubGridIterator, makeBoundaryIterator and makeVertexIterator
    are supported.
**/
template <typename T, typename IntegerIterator, typename PointProducer>
class gsPointIterator
    : public  IntegerIterator
{
public:
    using   typename IntegerIterator::MIndexT;
    using   typename IntegerIterator::FIndexT;
    typedef typename PointProducer::PointT PointT;
    using IntegerIterator::good;
public:
    gsPointIterator (
                const IntegerIterator & indexes,
                const PointProducer   & points)
        : IntegerIterator(indexes), m_conv(points)
    {
        update();
    }

    virtual bool       previous     ();
    virtual bool       next         ();
    virtual bool       first        ();
    virtual bool       last         ();
    virtual bool       setTo        (MIndexT const & here);
    virtual bool       setTo        (FIndexT here);

    const PointT & currPoint() const {return m_conv.m_curr;}
protected:
    PointProducer          m_conv;

    using IntegerIterator::d;
    using IntegerIterator::m_good;
    using IntegerIterator::m_vpos;

    void update ()
    {
        if(m_good)
        {
            m_conv.update(m_vpos);
        }
    }
};

// ?
template <typename T, typename IntegerIterator, int Tdim=-1>
class PPSimplexAveraging
{
public:
    friend class gsPointIterator<T,IntegerIterator,PPSimplexAveraging>;
    typedef gsVector<T,Tdim> PointT;
    typedef typename IntegerIterator::MIndexT MIndexT;
    enum {
        Ddim=IntegerIterator::Tdim
    };
public:
    PPSimplexAveraging (const gsMatrix<T> &vertex, T sum )
        : m_vert(vertex), m_sum(sum)
    {}
    void update(MIndexT const & newPos)
    {
        m_curr=m_vert*(newPos.template cast<T>().array()/m_sum).matrix();
    }
protected:
    gsMatrix<T,Tdim,Ddim>    m_vert;
    PointT                   m_curr;
    T                        m_sum;
};


// ?
template <typename T,typename IntegerIterator, int Tdim=-1>
class PPTensorAveraging
{
public:
    friend class gsPointIterator<T,IntegerIterator,PPTensorAveraging>;

    template <typename T2,typename IntegerIterator2, int Tdim2>
    friend class PPTensorAveraging;

    typedef gsVector<T,Tdim> PointT;
    typedef typename IntegerIterator::MIndexT MIndexT;
public:
    PPTensorAveraging (const MIndexT & max, const PointT & low, const PointT & upp)
       : m_low(low), m_upp(upp), m_max(max.template cast<T>().array()-1)
    {}

    template <typename IntegerIterator2>
    PPTensorAveraging (const PPTensorAveraging<T,IntegerIterator2,Tdim> & other)
       : m_low(other.m_low), m_upp(other.m_upp), m_max(other.m_max)
    {}

    void update(MIndexT newPos)
    {
        PointT tp=newPos.template cast<T>();
        m_curr= m_low.array()    * ( (m_max-tp).array() / m_max.array() )
                  +m_upp.array() * (        tp.array()  / m_max.array() ) ;
    }
protected:
    const PointT          m_low;
    const PointT          m_upp;
    const PointT          m_max;
    PointT                m_curr;
};


// documentation ?
template <typename T, int dim=-1, typename Flat=index_t>
class gsTensorPointGridIterator
    : public gsPointIterator<T,gsTensorGridIterator<Flat,dim>, PPTensorAveraging<T,gsTensorGridIterator<Flat,dim>,dim> >
{
    protected:
        using gsPointIterator<T,gsTensorGridIterator<Flat,dim>, PPTensorAveraging<T,gsTensorGridIterator<Flat,dim>,dim> >::m_conv;
        typedef gsTensorGridBoundaryIterator<Flat,dim>  boundaryItT;
        typedef PPTensorAveraging<T,boundaryItT,dim>    boundaryPMT;
        typedef gsTensorGridIterator<Flat,dim>          subGridItT;
        typedef PPTensorAveraging<T,subGridItT,dim>     subGridPMT;
    public:
        using   typename gsTensorGridIterator<Flat,dim>::MIndexT;
        typedef typename PPTensorAveraging<T,gsTensorGridIterator<Flat,dim>,dim>::PointT PointT;
        typedef gsPointIterator<T,boundaryItT, boundaryPMT>                          BoundaryPointIt;
        typedef gsPointIterator<T,subGridItT, subGridPMT>                            SubGridPointIt;
    public:
        gsTensorPointGridIterator(const MIndexT & M, const PointT & low,const PointT & upp)
            : gsPointIterator<T,gsTensorGridIterator<Flat,dim>, PPTensorAveraging<T,gsTensorGridIterator<Flat,dim>,dim> >
            ( gsTensorGridIterator<Flat,dim>(M),
              PPTensorAveraging<T,gsTensorGridIterator<Flat,dim>, dim> (M,low,upp) )

        {}

        virtual SubGridPointIt * makeSubGridIterator (MIndexT const & m, MIndexT const & M)
        {

            subGridItT * newB =  gsTensorGridIterator<Flat,dim>::makeSubGridIterator(m,M);
            gsPointIterator<T,subGridItT, subGridPMT> *result= new gsPointIterator<T,subGridItT, subGridPMT>(*newB,m_conv);
            delete newB;
            return result;
        }

        virtual BoundaryPointIt * makeBoundaryIterator (MIndexT const &mOffset, 
                                                        MIndexT const &MOffset)
        {
            boundaryItT * newB =  static_cast<boundaryItT*>( gsTensorGridIterator<Flat,dim>::makeBoundaryIterator(mOffset,MOffset));
            BoundaryPointIt *result= new gsPointIterator<T,boundaryItT, boundaryPMT>(*newB,m_conv);
            delete newB;
            return result;
        }
};


template <typename T, int dim=-1, typename Flat=index_t>
class gsSimplexPointGridIterator
    : public gsPointIterator<T, gsCompositionIterator<Flat,dim+1>, PPSimplexAveraging<T,gsCompositionIterator<Flat,dim+1>,dim > >
{
    public:
        typedef typename gsCompositionIterator<Flat,dim+1>::MIndexT MIndexT;
        typedef typename PPSimplexAveraging<T,gsCompositionIterator<Flat,dim+1>,dim >::PointT PointT;
    public:
        gsSimplexPointGridIterator( Flat pointPerSide, const gsMatrix<T,dim,dim+1> & vertex)
            : gsPointIterator< T,gsCompositionIterator<Flat,dim+1>, PPSimplexAveraging<T,gsCompositionIterator<Flat,dim+1>, dim > >
            ( gsCompositionIterator<Flat,dim+1>(pointPerSide, vertex.cols()),
              PPSimplexAveraging<T,gsCompositionIterator<Flat,dim+1>,dim >(vertex,pointPerSide) )

        {}
    protected:
        using gsPointIterator<T, gsCompositionIterator<Flat,dim+1>, PPSimplexAveraging<T,gsCompositionIterator<Flat,dim+1>,dim > >::m_conv;
};

// ?
template <typename T, typename Flat>
class gsSimplexPointGridIterator<T,-1,Flat>
    : public gsPointIterator<T, gsCompositionIterator<Flat>, PPSimplexAveraging<T,gsCompositionIterator<Flat> > >
{
    public:
        typedef typename gsCompositionIterator<Flat,-1>::MIndexT MIndexT;
        typedef typename PPSimplexAveraging<T,gsCompositionIterator<Flat> >::PointT PointT;
    public:
        gsSimplexPointGridIterator( Flat pointPerSide, const gsMatrix<T> & vertex)
            : gsPointIterator< T,gsCompositionIterator<Flat>, PPSimplexAveraging<T,gsCompositionIterator<Flat> > >
            ( gsCompositionIterator<Flat>(pointPerSide, vertex.cols()),
              PPSimplexAveraging<T,gsCompositionIterator<Flat> >(vertex,pointPerSide) )

        {}
    protected:
        using gsPointIterator<T, gsCompositionIterator<Flat>, PPSimplexAveraging<T,gsCompositionIterator<Flat> > >::m_conv;
};


///////////////////////////////////////////////////
///////////////////////////////////////////////////


template <typename T, typename IntegerIterator, typename PointProducer>
bool gsPointIterator<T,IntegerIterator,PointProducer>::previous()
{
    IntegerIterator::previous();
    update();
    return m_good;
}

template <typename T, typename IntegerIterator, typename PointProducer>
bool gsPointIterator<T,IntegerIterator,PointProducer>::next()
{
    IntegerIterator::next();
    update();
    return m_good;
}

template <typename T, typename IntegerIterator, typename PointProducer>
bool gsPointIterator<T,IntegerIterator,PointProducer>::first()
{
    IntegerIterator::first();
    update();
    return m_good;
}


template <typename T, typename IntegerIterator, typename PointProducer>
bool gsPointIterator<T,IntegerIterator,PointProducer>::last()
{
    IntegerIterator::last();
    update();
    return m_good;
}

template <typename T, typename IntegerIterator, typename PointProducer>
bool gsPointIterator<T,IntegerIterator,PointProducer>::setTo(MIndexT const & here)
{
    IntegerIterator::setTo(here);
    update();
    return m_good;
}

template <typename T, typename IntegerIterator, typename PointProducer>
bool gsPointIterator<T,IntegerIterator,PointProducer>::setTo(FIndexT here)
{
    IntegerIterator::setTo(here);
    update();
    return m_good;
}

} // namespace gismo
