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

#include <gsUtils/gsCombinatorics.h>
#include <gsCore/gsMultiPatch.h>

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
    gsVector<bool> orient( m_topology.dim() - 1 );
    orient.setConstant( true );
    m_topology.addInterface( p1, s1, p2, s2, orient );
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
            gsMatrix<unsigned>
                * b1= m_bases[it->first().patch]->boundary( it->first().side() ),
                * b2= m_bases[it->second().patch]->boundary( it->second().side() );
            
            mapper.matchInterface( it->first().patch, it->second().patch, *b1, *b2, it->orient());
            
            delete b1;
            delete b2;
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
    mapper.init(*this, bc, unk);
    
    if ( conforming ) // Conforming boundaries ?
    {
        for ( gsBoxTopology::const_iiterator it = m_topology.iBegin();
              it != m_topology.iEnd(); ++it )
        {
            gsMatrix<unsigned>
                * b1= m_bases[it->first().patch]->boundary( it->first().side() ),
                * b2= m_bases[it->second().patch]->boundary( it->second().side() );
            
            mapper.matchInterface( it->first().patch, it->second().patch, *b1, *b2, it->orient());
            
            delete b1;
            delete b2;
        }
    }

    mapper.finalize();
}
    
} // namespace gismo
