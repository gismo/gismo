/** @file gsMultiBasis.hpp

    @brief Provides declaration of MultiBasis class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <iterator>

#include <gsCore/gsMultiPatch.h>
#include <gsUtils/gsCombinatorics.h>

namespace gismo
{

template<class T>
gsMultiBasis<T>::gsMultiBasis( const gsBasis<T> & bb )
: m_topology( bb.dim() )
{
    m_bases.push_back( bb.clone() );
    m_topology.addBox();
    m_topology.addAutoBoundaries();
}

template<class T>
gsMultiBasis<T>::gsMultiBasis( const gsMultiPatch<T> & mpatch )
: m_topology( mpatch )
{
    m_bases = mpatch.basesCopy();
}
  
template<class T>
gsMultiBasis<T>::gsMultiBasis( const gsMultiBasis& other )
: m_bases         ( other.m_bases        ),
  m_topology      ( other.m_topology     )
{
    cloneAll( other.m_bases.begin(), other.m_bases.end(),
              this->m_bases.begin() );
}

template<class T>
gsMultiBasis<T>::~gsMultiBasis()
{
    freeAll(m_bases);
}

template<class T>
std::ostream& gsMultiBasis<T>::print( std::ostream& os ) const
{
    gsInfo<<"Topology: "<< m_topology <<"\n";
    
    return os;
}

 
template<class T>
void gsMultiBasis<T>::addBasis( gsBasis<T>* g ) 
{
    gsDebug<< "TO DO\n";
    if ( m_topology.dim() == -1 ) 
    {
        m_topology.setDim( g->dim() );
    } else {
        assert( g->dim() == m_topology.dim() );
    }
    m_bases.push_back( g ) ;
    m_topology.addBox();
}
  
template<class T>
int gsMultiBasis<T>::findBasisIndex( gsBasis<T>* g ) const {
    typename BasisContainer::const_iterator it
        = std::find( m_bases.begin(), m_bases.end(), g );
    assert( it != m_bases.end() );
    return it - m_bases.begin();
}
  
template<class T>
void gsMultiBasis<T>::addInterface( gsBasis<T>* g1, boxSide s1,
                                    gsBasis<T>* g2, boxSide s2 ) 
{
    int p1 = findBasisIndex( g1 );
    int p2 = findBasisIndex( g2 );
    m_topology.addInterface( p1, s1, p2, s2);
}
  

template<class T>
int gsMultiBasis<T>::maxDegree(int k) const
{
    GISMO_ASSERT(m_bases.size(), "Empty multibasis.");
    int result = m_bases[0]->degree(k);
    for (size_t i = 0; i < m_bases.size(); ++i)
        if (m_bases[i]->degree(k) > result )
            result = m_bases[i]->degree(k);
    return result;
}

template<class T>
int gsMultiBasis<T>::maxCwiseDegree() const
{
    GISMO_ASSERT(m_bases.size(), "Empty multibasis.");
    int result = m_bases[0]->maxDegree();
    for (size_t i = 0; i < m_bases.size(); ++i)
        result = math::max(m_bases[i]->maxDegree(), result);
    return result;
}

template<class T>
int gsMultiBasis<T>::minCwiseDegree() const
{
    GISMO_ASSERT(m_bases.size(), "Empty multibasis.");
    int result = m_bases[0]->minDegree();
    for (size_t i = 0; i < m_bases.size(); ++i)
        result = math::min(m_bases[i]->minDegree(), result);
    return result;
}

template<class T>
int gsMultiBasis<T>::minDegree(int k) const
{
    GISMO_ASSERT(m_bases.size(), "Empty multibasis.");
    int result = m_bases[0]->degree(k);
    for (size_t i = 0; i < m_bases.size(); ++i)
        if (m_bases[i]->degree(k) < result )
            result = m_bases[i]->degree(k);
    return result;
}

template<class T>
void gsMultiBasis<T>::getMapper(bool conforming, gsDofMapper & mapper) const
{
    mapper = gsDofMapper(*this);//.init(*this);
    
    if ( conforming )  // Conforming boundaries ?
    {
        for ( gsBoxTopology::const_iiterator it = m_topology.iBegin();
              it != m_topology.iEnd(); ++it )
        {
            matchInterface(*it,mapper);
        }
    }
    
    mapper.finalize();
}


template<class T>
void gsMultiBasis<T>::getMapper(bool conforming, 
                                const gsBoundaryConditions<T> & bc, 
                                int unk,
                                gsDofMapper & mapper) const
{
    mapper = gsDofMapper(*this, bc, unk); //.init(*this, bc, unk);
    
    if ( conforming ) // Conforming boundaries ?
    {
        for ( gsBoxTopology::const_iiterator it = m_topology.iBegin();
              it != m_topology.iEnd(); ++it )
        {
            matchInterface(*it,mapper);
        }
    }

    mapper.finalize();
}
    
template<class T>
void gsMultiBasis<T>::matchInterface(const boundaryInterface & bi, gsDofMapper & mapper) const
{
    // Grab the indices to be matched
    gsMatrix<unsigned>
        * b1= m_bases[bi.first() .patch]->boundary( bi.first() .side() ),
        * b2= m_bases[bi.second().patch]->boundary( bi.second().side() );

    GISMO_ASSERT( b1->rows() == b2->rows(), 
                  "Input error, sizes do not match: "<<b1->rows()<<"!="<<b2->rows() );


    // Compute tensor structure of b1 -- to do move to tensor basis
    const index_t d = dim();
    const index_t p1 = bi.first().patch;
    const index_t s1 = bi.first().direction();
    gsVector<int>  bSize(d-1);
    index_t c = 0;
    for (index_t k = 0; k<d; ++k )
    {
        if ( k == s1 ) 
            continue;
        bSize[c] = m_bases[p1]->component(k).size();
        c++;
    }

    // Reorder the indices so that they match on the interface
    bi.matchDofs(bSize, *b1, *b2);

    // All set, match interface dofs
    for (c = 0; c<b1->size(); ++c)
        mapper.matchDof(bi.first().patch, (*b1)(c,0), bi.second().patch, (*b2)(c,0) );

    delete b1;
    delete b2;
}



} // namespace gismo
