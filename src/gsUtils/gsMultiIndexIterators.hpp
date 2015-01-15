/*
* gsMultiIndexIterarors.h created on 09.07.2014
*
* Author: Andrea Bressan
*
* This file is part of the G+SMO library
*/

#pragma once

#include<gsUtils/gsMultiIndexIterators.h>
#include<gsUtils/gsCombinatorics.h>

namespace gismo {


///////////////////////////////////////////
/// TENSOR GRID ITERATOR FUNCTIONS
/// 
/// \ingroup Utils

template<typename Flat, int dim>
gsTensorGridIterator<Flat,dim>::gsTensorGridIterator (MIndexT const & M)
    : gsMultiIndexIterator<Flat,dim>(M.size()),
      m_vmax(M),
      m_vinc(M),
      m_fmin(0),
      m_fmax(M.prod())
{
    m_vmin.setConstant(d,0);
    m_vpos=m_vmin;
    m_vinc[0]=1;
    for (int i=1;i<M.size();i++)
    {
        m_vinc[i]=m_vinc[i-1]*M[i-1];
    }
    GISMO_ASSERT(m_fmax==m_vinc[M.size()-1]*M[M.size()-1], "absurd error");
    if ( (M.array()>0).all() )
        m_good=true;
}

template<typename Flat, int dim>
gsTensorGridIterator<Flat,dim>::gsTensorGridIterator (MIndexT const & m, MIndexT const & M)
    : gsMultiIndexIterator<Flat,dim>(M.size()),
      m_vmin(m),
      m_vmax(M),
      m_vinc(M),
      m_fmin(0),
      m_fmax((M-m).prod())
{
    m_vpos=m;
    MIndexT size = M-m;
    m_vinc(0)=1;
    for (int i=1;i<M.size();i++)
    {
        m_vinc(i)=m_vinc(i-1)*size(i-1);
    }
    GISMO_ASSERT(m_fmax==m_vinc(M.size()-1)*size(M.size()-1), "absurd error");
    if ( (M.array()>m.array()).all() )
        m_good=true;
}

template<typename Flat, int dim>
bool  gsTensorGridIterator<Flat,dim>::next()
{
    if (!m_good)
        return false;
    int i;
    for (i=0; i<d && m_vpos(i)==m_vmax(i)-1; ++i)
    {
        m_vpos(i)=m_vmin(i);
    }
    if(i!=d)
        m_vpos(i)+=1;
    else
        m_good=false;
    return  m_good;
}

template<typename Flat, int dim>
bool  gsTensorGridIterator<Flat,dim>::previous()
{
    if (!m_good)
        return false;
    int i;
    for (i=0; i<d && m_vpos(i)==m_vmin(i); ++i)
    {
        m_vpos(i)=m_vmax(i)-1;
    }
    if(i!=d)
        m_vpos(i)-=1;
    else
        m_good=false;
    return  m_good;
}

template<typename Flat, int dim>
bool  gsTensorGridIterator<Flat,dim>::first()
{
    m_vpos=m_vmin;
    if ( (m_vmax.array()>m_vmin.array()).all() )
        m_good=true;
    return m_good;
}

template<typename Flat, int dim>
bool  gsTensorGridIterator<Flat,dim>::last()
{
    m_vpos=m_vmax.array()-1;
    if ( (m_vmax.array()>m_vmin.array()).all() )
        m_good=true;
    return m_good;
}

template<typename Flat, int dim>
Flat  gsTensorGridIterator<Flat,dim>::flatIndex()   const
{
    return m_fmin+m_vinc.dot(m_vpos-m_vmin);
}

template<typename Flat, int dim>
bool  gsTensorGridIterator<Flat,dim>::setTo (MIndexT const & here)
{
    if ( this->multiIsValid(here) )
    {
        m_vpos=here;
        m_good=true;
    }
    else
    {
        m_good=false;
    }
    return m_good;
}

template<typename Flat, int dim>
bool  gsTensorGridIterator<Flat,dim>::setTo (Flat here)
{
    return setTo(flatToMulti(here));
}

template<typename Flat, int dim>
gsTensorGridIterator<Flat,dim>
        * gsTensorGridIterator<Flat,dim>::makeSubGridIterator (MIndexT const & m, MIndexT const & M)
{
    MIndexT mm = m.cwiseMax( m_vmin ),
            MM = M.cwiseMin( m_vmax );

    gsTensorGridIterator<Flat,dim> *result;
    result = new gsTensorGridIterator<Flat,dim>(mm,MM);
    result->m_vinc = m_vinc;
    result->m_fmin = multiToFlat(m);
    return result;
}

template<typename Flat, int dim>
gsTensorGridIterator<Flat,dim>
        * gsTensorGridIterator<Flat,dim>::makeBoundaryIterator (MIndexT const & mOffset, MIndexT const & MOffset)
{
    return new gsTensorGridBoundaryIterator<Flat,dim>(*this, mOffset,MOffset);
}


template<typename Flat, int dim>
gsTensorGridIterator<Flat,dim>
        * gsTensorGridIterator<Flat,dim>::makeVertexIterator ()
{
    return new gsTensorGridVertexIterator<Flat,dim>(*this);
}


template<typename Flat, int dim>
Flat gsTensorGridIterator<Flat,dim>::multiToFlat(MIndexT const & index)  const
{
    return m_fmin+m_vinc.dot(index-m_vmin);
}

template<typename Flat, int dim>
typename gsMultiIndexIterator<Flat,dim>::MIndexT gsTensorGridIterator<Flat,dim>::flatToMulti(Flat index)  const
{
    MIndexT result(d);

    index-=m_fmin;
    for (int i=d-1; i>=0;i--)
    {
        result(i)=index/m_vinc(i);
        index%=m_vinc(i);
    }
    return result+m_vmin;
}

template<typename Flat, int dim>
bool gsTensorGridIterator<Flat,dim>::multiIsValid(MIndexT const & index)  const
{
    return (index.array()>=m_vmin.array()).all() && (index.array()<m_vmax.array()).all();
}

template<typename Flat, int dim>
bool gsTensorGridIterator<Flat,dim>::flatIsValid(Flat index)  const
{
    MIndexT multi = flatToMulti(index);
    return  this->multiIsValid(multi) && ( index==multiToFlat(multi) );
}


///////////////////////////////////////////
/// TENSOR GRID BOUNDARY ITERATOR FUNCTIONS
///

template<typename Flat, int dim>
gsTensorGridVertexIterator<Flat,dim>::gsTensorGridVertexIterator(gsTensorGridIterator<Flat,dim> &area)
    : gsTensorGridIterator<Flat,dim>(area)
{
    this->reset();
    if (multiIsValid(m_vpos))
        m_good=true;
}

template<typename Flat, int dim>
gsTensorGridVertexIterator<Flat,dim>::gsTensorGridVertexIterator(MIndexT const & lower, MIndexT const & upper)
    : gsTensorGridIterator<Flat,dim>(lower,upper)
{
    this->reset();
    if (multiIsValid(m_vpos))
        m_good=true;
}

template<typename Flat, int dim>
bool gsTensorGridVertexIterator<Flat,dim>::next()
{
    if (!m_good)
        return false;
    int i;
    for (i=0; i<d && m_vpos(i)==m_vmax(i)-1; ++i)
    {
        m_vpos(i)=m_vmin(i);
    }
    if(i!=d)
        m_vpos(i)=m_vmax(i)-1;
    else
        m_good=false;
    return  m_good;
}

template<typename Flat, int dim>
bool gsTensorGridVertexIterator<Flat,dim>::previous()
{
    if (!m_good)
        return false;
    int i;
    for (i=0; i<d && m_vpos(i)==m_vmin(i); ++i)
    {
        m_vpos(i)=m_vmax(i)-1;
    }
    if(i!=d)
        m_vpos(i)=m_vmin(i);
    else
        m_good=false;
    return  m_good;
}

template<typename Flat, int dim>
bool gsTensorGridVertexIterator<Flat,dim>::multiIsValid (MIndexT const & index)  const
{
    return ((index.array()==(m_vmax.array()-1))||(index.array()==m_vmin.array())).all();
}


///////////////////////////////////////////
/// TENSOR GRID BOUNDARY ITERATOR FUNCTIONS
///

template<typename Flat, int dim>
gsTensorGridBoundaryIterator<Flat,dim>::gsTensorGridBoundaryIterator(
            gsTensorGridIterator<Flat,dim> &area,
            MIndexT                    const & mOffset,
            MIndexT                    const & MOffset)
            : gsTensorGridIterator<Flat,dim>(area), m_lowOff(mOffset), m_uppOff(MOffset),m_mne(d)
{
    for (int i=d-1;i>=0;--i)
    {
        m_lowOff(i) = m_lowOff(i)<0 ? 0:m_lowOff(i);
        m_uppOff(i) = m_uppOff(i)<0 ? 0:m_uppOff(i);
        if (m_lowOff(i)+m_uppOff(i)>m_vmax(i)-m_vmin(i))
        {
            m_lowOff(i) = m_vmax(i)-m_vmin(i);
            m_uppOff(i) = 0;
        }
        if ( m_lowOff(i)+m_uppOff(i))
            m_mne=i;
    }
    if (m_mne==0 && m_lowOff(0)+m_uppOff(0)==0)
    {
        m_mne=d;
    }
    first();
}

template<typename Flat, int dim>
bool  gsTensorGridBoundaryIterator<Flat,dim>::first()
{
    if(m_mne==d)
        m_good=false;
    else
    {
        m_vpos=m_vmin;
        int i;
        for (i=m_mne;i<d && !m_lowOff(i);++i)
            ;
        if (i==d)
            m_vpos(m_mne)=m_vmax(m_mne)-m_uppOff(m_mne);
        m_good=true;
    }
    return m_good;
}

template<typename Flat, int dim>
bool  gsTensorGridBoundaryIterator<Flat,dim>::last()
{
    if(m_mne==d)
        m_good=false;
    else
    {
        m_vpos=m_vmax.array()-1;
        int i;
        for (i=m_mne;i<d && !m_uppOff(i);++i)
            ;
        if (i==d)
            m_vpos(m_mne)=m_vmin(m_mne)+m_lowOff(m_mne)-1;
        m_good=true;
    }
    return m_good;
}

template<typename Flat, int dim>
bool  gsTensorGridBoundaryIterator<Flat,dim>::previous()
{
    if (!m_good)
        return false;
    int  tc;  // to change
    // mne = minimum non empty: the index of the first direction in which the offset is not 0

    for (tc=0; tc < d && m_vpos(tc)==m_vmin(tc); ++tc)
    {
            m_vpos(tc)=m_vmax(tc)-1;
    }
    if (tc == d)
    {
        m_good=false;
    }
    else if (tc < m_mne)
    {
        m_vpos(tc)-=1;
    }
    else if (tc > m_mne)
    {
        m_vpos(tc)-=1;
        for ( tc=m_mne; tc<d; ++tc )
        {
            if ( m_vpos(tc)-m_vmin(tc)<m_lowOff(tc) || m_vmax(tc)-m_vpos(tc)<=m_uppOff(tc) )
                break;
        }
        if (tc==d && m_uppOff(m_mne)==0 )
            m_vpos(m_mne)=m_vmin(m_mne)+m_lowOff(m_mne)-1;
    }
    else
    {
        m_vpos(tc)-=1;
        int i;
        for (i=tc; i<d; ++i)
        {
            if ( m_vpos(i)-m_vmin(i)<m_lowOff(i) || m_vmax(i)-m_vpos(i)<=m_uppOff(i) )
                break;
        }
        if ( i==d ) // not in boundary
        {
            if( m_lowOff(tc) && m_vpos(tc)>= m_vmin(tc)+m_lowOff(tc) )
                m_vpos(tc)=m_vmin(tc)+m_lowOff(tc)-1;
            else
            {
                if (m_uppOff(tc))
                    m_vpos(tc)=m_vmax(tc)-1;
                else
                    m_vpos(tc)=m_vmin(tc)+m_lowOff(tc)-1;
                for (tc+=1; tc < d && m_vpos(tc)==m_vmin(tc); ++tc)
                {
                    m_vpos(tc)=m_vmax(tc)-1;
                }
                if (tc == d)
                {
                    m_good=false;
                }
                else
                    m_vpos(tc)-=1;
            }
        }
    }
    return m_good;
}

template<typename Flat, int dim>
bool  gsTensorGridBoundaryIterator<Flat,dim>::next()
{
    if (!m_good)
        return false;
    int  tc;  // to change
    // mne = minimum non empty: the index of the first direction in which the offset is not 0

    for (tc=0; tc < d && m_vpos(tc)==m_vmax(tc)-1; ++tc)
    {
            m_vpos(tc)=m_vmin(tc);
    }
    if (tc == d)
    {
        m_good=false;
    }
    else if (tc < m_mne)
    {
        // we are below the mne direction -> some direction with higher index is in the boundary
        // remember m_good=true -> not empty
        m_vpos(tc)+=1;
    }
    else if (tc > m_mne)
    {
        // we are above the mne -> increment by one and if not in other boundary set mne to the
        //                         minimum inside a boundary
        m_vpos(tc)+=1;
        for ( tc=m_mne; tc<d; ++tc )
        {
            if ( m_vpos(tc)-m_vmin(tc)<m_lowOff(tc) || m_vmax(tc)-m_vpos(tc)<=m_uppOff(tc) )
                break;
        }
        if (tc==d && m_lowOff(m_mne)==0 )
            m_vpos(m_mne)=m_vmax(m_mne)-m_uppOff(m_mne);
    }
    else // if (tc==mne)
    {
        // we are at mne -> we increment by one, if we are in the boundary -> ok
        //                                       else we check if we can increment vpos(tc) to get in the boundary
        //                                       else we decrease vpos(tc) to the boundary and increase vpos(tc+1)
        m_vpos(tc)+=1;
        int i;
        for (i=tc; i<d; ++i)
        {
            if ( m_vpos(i)-m_vmin(i)<m_lowOff(i) || m_vmax(i)-m_vpos(i)<=m_uppOff(i) )
                break;
        }
        if ( i==d ) // not in boundary
        {
            if( m_uppOff(tc) && m_vpos(tc)<m_vmax(tc)-m_uppOff(tc) )
                m_vpos(tc)=m_vmax(tc)-m_uppOff(tc);
            else
            {
                if (m_lowOff(tc))
                    m_vpos(tc)=m_vmin(tc);
                else
                    m_vpos(tc) =m_vmax(tc)-m_uppOff(tc);
                for (tc+=1; tc < d && m_vpos(tc)==m_vmax(tc)-1; ++tc)
                {
                    m_vpos(tc)=m_vmin(tc);
                }
                if (tc == d)
                {
                    m_good=false;
                }
                else
                    m_vpos(tc)+=1;
            }
        }
    }
    return m_good;
}

template<typename Flat, int dim>
bool gsTensorGridBoundaryIterator<Flat,dim>::multiIsValid(MIndexT const & index)  const
{
    return
       ( index.array()>=m_vmin.array() && index.array()<(m_vmin+m_lowOff).array() ).any()
       ||
       ( index.array()>=(m_vmax-m_uppOff).array() && index.array()<m_vmax.array()).any() ;
}

///////////////////////////////////////////
/// COMPOSITION ITERATOR FUNCTIONS
///
///

template<typename Flat, int dim>
gsCompositionIterator<Flat,dim>::gsCompositionIterator(Flat sum, int _dim)
    : gsMultiIndexIterator<Flat,dim>(_dim),
      m_fmin(0),
      m_fpos(m_fmin),
      m_fsum(sum)
{
    m_vmin.setConstant(d,0);
    m_vmax=m_vmin;
    m_vmin(0)=m_fsum;
    m_vmax(d-1)=m_fsum;
    m_vpos=m_vmin;
    m_fmax=m_fmin+binomial<Flat>(m_fsum+d-1, m_fsum)-1;
    m_good=true;
}


template<typename Flat, int dim>
bool  gsCompositionIterator<Flat,dim>::next()
{
    if (!m_good || m_vpos==m_vmax)
    {
        m_good=false;
        return m_good;
    }
    int i=0;
    Flat t;
    while (m_vpos(i)==0)
        i++;
    t=m_vpos(i); // if i=0 save the value
    m_vpos(i)=0;
    m_vpos(0)=t-1;
    m_vpos(i+1)++;
    m_fpos++;
    return m_good;
}

template<typename Flat, int dim>
bool  gsCompositionIterator<Flat,dim>::previous()
{
        if ( !m_good || m_vpos==m_vmin)
        {
            m_good=false;
            return m_good;
        }
        int i=1;
        while (m_vpos(i)==0)
            i++;
        m_vpos(i-1)=m_vpos(0)+1;
        m_vpos(i)--;
        if (i!=1)
            m_vpos(0)=0;
        m_fpos--;
        return m_good;
}

template<typename Flat, int dim>
bool  gsCompositionIterator<Flat,dim>::first() {return setTo(m_vmin);}

template<typename Flat, int dim>
bool  gsCompositionIterator<Flat,dim>::last()  {return setTo(m_vmax);}

template<typename Flat, int dim>
Flat  gsCompositionIterator<Flat,dim>::flatIndex() const {return m_fpos;}

template<typename Flat, int dim>
bool  gsCompositionIterator<Flat,dim>::setTo (MIndexT const & here)
{
    if (multiIsValid(here))
    {
        m_vpos=here;
        m_fpos=multiToFlat(here);
        m_good=true;
    }
    else
        m_good=false;
    return m_good;
}

template<typename Flat, int dim>
bool  gsCompositionIterator<Flat,dim>::setTo (Flat here)
{
    if (flatIsValid(here))
    {
        m_vpos=flatToMulti(here);
        m_fpos=here;
        m_good=true;
    }
    else
        m_good=false;
    return m_good;
}


template<typename Flat, int dim>
Flat  gsCompositionIterator<Flat,dim>::multiToFlat (MIndexT const & index) const
{
    GISMO_ASSERT(multiIsValid(index), "Invalid index");

    Flat f_index= m_fmax;
    Flat to_sum = m_fsum-index(d-1);

    for (int i=d-1; i>0 && to_sum>0;)
    {
         f_index -= binomial<Flat>( to_sum+i-1, i);
         to_sum  -= index(--i);
    }
    return f_index;
}

template<typename Flat, int dim>
typename gsMultiIndexIterator<Flat,dim>::MIndexT gsCompositionIterator<Flat,dim>::flatToMulti (Flat index) const
{
    GISMO_ASSERT(flatIsValid(index), "Invalid index");

    int cur_digit, temp_digit, i;
    Flat cur_limit, cur_increment, missing(m_fsum);
    MIndexT result;
    result.setConstant(d,0);
    temp_digit=d;
    while(missing>0)
    {
        // find the position of the mosts significant non-null digit
        cur_digit=0;
        cur_limit=1;
        // cur_limit=binomial(missing+cur_digit,cur_digit);
        while (cur_digit<temp_digit && index>=cur_limit)
        {
            cur_digit++;
            cur_limit= cur_limit *(missing+cur_digit)/cur_digit;
            // we are computing binomial (missing+cur_digit,cur_digit)
            // from binomial (missing+cur_digit-1,cur_digit-1)
        }

        // find the most significant digit
        i=missing;
        cur_increment = 1;
        // cur_increment = binomial(missing+cur_digit-i,cur_digit)=1;
        while (i>0 && index+cur_increment<cur_limit)
        {
            i--;
            cur_increment = cur_increment*(missing+cur_digit-i)/(missing-i);
            // cur_increment = binomial(missing+cur_digit-i,cur_digit);
        }

        result[cur_digit]=i;
        missing-=i;
        index+=cur_increment;
        index-=cur_limit;
        temp_digit=cur_digit;
    }
    return result;
}


template<typename Flat, int dim>
bool  gsCompositionIterator<Flat,dim>::multiIsValid (MIndexT const & index) const {return index.sum()==m_fsum;}

template<typename Flat, int dim>
bool  gsCompositionIterator<Flat,dim>::flatIsValid (Flat index) const {return m_fmin<=index && index<=m_fmax;}


///////////////////////////////////////////
/// SIMPLEX ITERATOR FUNCTIONS
///
///

template<typename Flat, int dim>
gsSimplexIterator<Flat,dim>::gsSimplexIterator(Flat sum, int _dim)
    : gsMultiIndexIterator<Flat,dim>(_dim),
      m_fmin(0),
      m_fsum(sum)
{
    m_vmin.setConstant(d,0);
    m_vmax=m_vmin;
    m_vmax(d-1)=m_fsum;
    m_fmax=m_fmin+binomial<Flat>(m_fsum+d, m_fsum)-1;
    first();
}


template<typename Flat, int dim>
bool  gsSimplexIterator<Flat,dim>::next()
{
    if (!m_good || m_fpos==m_fmax)
    {
        m_good=false;
        return m_good;
    }
    if (m_comp!=0)
    {
        m_comp-=1;
        m_vpos(0)+=1;
    }
    else
    {
        int i=0;
        while (m_vpos(i)==0)
            i++;
        m_comp=m_vpos(i)-1;
        m_vpos(i)=0;
        m_vpos(i+1)+=1;
    }
    m_fpos++;
    return m_good;
}

template<typename Flat, int dim>
bool  gsSimplexIterator<Flat,dim>::previous()
{
        if ( !m_good || m_fpos==m_fmin)
        {
            m_good=false;
            return m_good;
        }
        if (m_vpos(0)!=0)
        {
            m_comp+=1;
            m_vpos(0)-=1;
        }
        else
        {
            int i=1;
            while (m_vpos(i)==0)
                i++;
            m_vpos(i-1)=m_comp+1;
            m_vpos(i)-=1;
            m_comp=0;
        }
        m_fpos--;
        return m_good;
}

template<typename Flat, int dim>
bool  gsSimplexIterator<Flat,dim>::first() {return setTo(m_vmin);}

template<typename Flat, int dim>
bool  gsSimplexIterator<Flat,dim>::last()  {return setTo(m_vmax);}

template<typename Flat, int dim>
Flat  gsSimplexIterator<Flat,dim>::flatIndex() const {return m_fpos;}

template<typename Flat, int dim>
bool  gsSimplexIterator<Flat,dim>::setTo (const MIndexT & here)
{
    if (multiIsValid(here))
    {
        m_vpos=here;
        m_fpos=multiToFlat(here);
        m_comp=m_fsum-here.sum();
        m_good=true;
    }
    else
        m_good=false;
    return m_good;
}

template<typename Flat, int dim>
bool  gsSimplexIterator<Flat,dim>::setTo (Flat here)
{
    if (flatIsValid(here))
    {
        m_vpos=flatToMulti(here);
        m_fpos=here;
        m_comp=m_fsum-m_vpos.sum();
        m_good=true;
    }
    else
        m_good=false;
    return m_good;
}


template<typename Flat, int dim>
Flat  gsSimplexIterator<Flat,dim>::multiToFlat (const MIndexT &index) const
{
    GISMO_ASSERT(multiIsValid(index), "Invalid index");

    Flat f_index= m_fmax;
    Flat to_sum = m_fsum-index(d-1);
    for (int i=d-1; i>0 && to_sum>0;)
    {
         f_index -= binomial<Flat>( to_sum+i, i+1);
         to_sum  -= index(--i);
    }
    f_index-=m_comp;
    return f_index;
}

template<typename Flat, int dim>
typename gsMultiIndexIterator<Flat,dim>::MIndexT gsSimplexIterator<Flat,dim>::flatToMulti (Flat index) const
{
        int cur_digit, temp_digit, i;
        Flat cur_limit, cur_increment, missing(m_fsum);
        gsVector<Flat> result(d+1);
        result.setConstant(d+1,0);
        temp_digit=d;
        while(missing>0)
        {
            // find the position of the mosts significant non-null digit
            cur_digit=0;
            cur_limit=1;
            // cur_limit=binomial(missing+cur_digit,cur_digit);
            while (cur_digit<temp_digit && index>=cur_limit)
            {
                cur_digit++;
                cur_limit= cur_limit *(missing+cur_digit)/cur_digit;
                // we are computing binomial (missing+cur_digit,cur_digit)
                // from binomial (missing+cur_digit-1,cur_digit-1)
            }
            // find the most significant digit
            i=missing;
            cur_increment = 1;
            // cur_increment = binomial(missing+cur_digit-i,cur_digit)=1;
            while (i>0 && index+cur_increment<cur_limit)
            {
                i--;
                cur_increment = cur_increment*(missing+cur_digit-i)/(missing-i);
                // cur_increment = binomial(missing+cur_digit-i,cur_digit);
            }

            result(cur_digit)=i;
            missing-=i;
            index+=cur_increment;
            index-=cur_limit;
            temp_digit=cur_digit;
        }
        return result.block(1,0,d,1);
}


template<typename Flat, int dim>
bool  gsSimplexIterator<Flat,dim>::multiIsValid (const MIndexT &index) const
{
    Flat sum=index.sum();
    return 0<=sum && sum<=m_fsum && (index.array()>=0).all();
}

template<typename Flat, int dim>
bool  gsSimplexIterator<Flat,dim>::flatIsValid (Flat index) const {return m_fmin<=index && index<m_fmax;}


} // namespace gismo
