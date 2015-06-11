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
//#include <gsHSplines/gsTHBSpline.h>
#include <gsHSplines/gsHTensorBasis.h>
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
void gsMultiBasis<T>::getMapper(bool conforming, 
                                gsDofMapper & mapper, 
                                bool finalize) const
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
    
    if (finalize)
        mapper.finalize();
}


template<class T>
void gsMultiBasis<T>::getMapper(bool conforming, 
                                const gsBoundaryConditions<T> & bc, 
                                int unk,
                                gsDofMapper & mapper, 
                                bool finalize) const
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

    if (finalize)
        mapper.finalize();
}
    
template<class T>
void gsMultiBasis<T>::matchInterface(const boundaryInterface & bi, gsDofMapper & mapper) const
{
    const gsHTensorBasis<2,T> * bas0 = dynamic_cast< const gsHTensorBasis<2,T> * >( m_bases[ bi.first().patch ] );
    const gsHTensorBasis<2,T> * bas1 = dynamic_cast< const gsHTensorBasis<2,T> * >( m_bases[ bi.second().patch ] );

    // Check, if an gsHTensorBasis is involved, and if so,
    // call matchInterfaceHTensor
    if( bas0 != 0 && bas1 != 0 )
        matchInterfaceHTensor( bi, mapper );
    else if( bas0 != 0 || bas1 != 0 )
        GISMO_ASSERT(false, "One Basis is HTensor, the other is not. Or dimension is not 2. Cannot handle this. You should implement that.");
    else
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
    } // if-else
}

template<class T>
void gsMultiBasis<T>::matchInterfaceHTensor(const boundaryInterface & bi, gsDofMapper & mapper) const
{
    const gsHTensorBasis<2,T> * bas0 = dynamic_cast< const gsHTensorBasis<2,T> * >( m_bases[ bi.first().patch ] );
    const gsHTensorBasis<2,T> * bas1 = dynamic_cast< const gsHTensorBasis<2,T> * >( m_bases[ bi.second().patch ] );

    // see if the orientation is preserved on side second()
    const gsVector<bool> dirOrient = bi.dirOrientation();
    bool orientPreserv = dirOrient[ ( bi.first().direction() + 1 ) % 2 ];

    // get the global indices of the basis functions which are
    // active on the interface
    gsMatrix<unsigned> b0b = * bas0->boundary( bi.first().side() );
    //gsMatrix<unsigned> b1b = * bas1->boundary( bi.second().side() );

    for( unsigned i=0; i < b0b.rows(); i++)
    {
        // get the level of the basis function on side first()
        unsigned L = bas0->levelOf( b0b(i,0) );
        // get the flat tensor index
        // (i.e., the single-number-index on level L)...
        unsigned flat0 = bas0->flatTensorIndexOf( b0b(i,0) );
        // ... and change it to the tensor-index.
        gsVector<unsigned> tens0 = bas0->tensorLevel(L).tensorIndex( flat0 );

        // tens1 will store the tensor-index on side second()...
        gsVector<unsigned> tens1(2);
        // ...and flat1 the corresponding flat index.
        unsigned flat1;

        // idxUse stores the "relevant" one of the tensor indices
        // on side first(). The other index is constant for the
        // whole domain side, right?
        unsigned idxUse;
        if( bi.first().direction() == 0 ) // west or east
            idxUse = tens0[1];
        else
            idxUse = tens0[0];

        // Size of the tensor-product basis of level L on
        // interface side second().
        unsigned N0 = bas1->tensorLevel(L).component(0).size();
        unsigned N1 = bas1->tensorLevel(L).component(1).size();

        switch( bi.second().side().index() )
        {
        case 1: // west
            tens1[0] = 0;
            if( orientPreserv )
                tens1[1] = idxUse;
            else
                tens1[1] = N1-1 - idxUse;
            flat1 = bas1->tensorLevel(L).index( tens1 );
            break;
        case 2: // east
            tens1[0] = N0-1;
            if( orientPreserv )
                tens1[1] = idxUse;
            else
                tens1[1] = N1-1 - idxUse;
            flat1 = bas0->tensorLevel(L).index( tens1 );
            break;
        case 3: // south
            tens1[1] = 0;
            if( orientPreserv )
                tens1[0] = idxUse;
            else
                tens1[0] = N0-1 - idxUse;
            flat1 = bas1->tensorLevel(L).index( tens1 );
            break;
        case 4: // north
            tens1[1] = N1-1;
            if( orientPreserv )
                tens1[0] = idxUse;
            else
                tens1[0] = N0-1 - idxUse;
            flat1 = bas0->tensorLevel(L).index( tens1 );
            break;
        default:
            GISMO_ASSERT(false,"3D, huh? Not implemented yet. You can do it!");
        }

        // compute the "continuous" index on second(), i.e., the index
        // in the numbering which is global over all levels.
        unsigned cont1 = bas1->flatTensorIndexToHierachicalIndex( flat1, L );

        // finally, match these two degrees of freedom
        mapper.matchDof( bi.first().patch, b0b(i,0), bi.second().patch, cont1 );
    }
}


template<class T>
bool gsMultiBasis<T>::repairInterface( const boundaryInterface & bi )
{
    // get direction and orientation maps
    const gsVector<bool> dirOrient = bi.dirOrientation();

    // get the bases of both sides as gsHTensorBasis
    const gsHTensorBasis<2,T> * bas0 = dynamic_cast< const gsHTensorBasis<2,T> * >( m_bases[ bi.first().patch ] );
    const gsHTensorBasis<2,T> * bas1 = dynamic_cast< const gsHTensorBasis<2,T> * >( m_bases[ bi.second().patch ] );

    gsMatrix<unsigned> lo;
    gsMatrix<unsigned> up;
    gsVector<unsigned> level;

    int dir;
    unsigned idxExponent;

    // get the higher one of both indexLevels
    unsigned indexLevelUse = ( bas0->tree().getIndexLevel() > bas1->tree().getIndexLevel() ? bas0->tree().getIndexLevel() : bas1->tree().getIndexLevel() );
    unsigned indexLevelDiff0 = indexLevelUse - bas0->tree().getIndexLevel();
    unsigned indexLevelDiff1 = indexLevelUse - bas1->tree().getIndexLevel();

    // get the box-representation of the gsHDomain on the interface
    bas0->tree().getBoxesOnSide( bi.first().side(), lo, up, level);

    dir = ( bi.first().direction() + 1 ) % 2;
    bool orientPreserv = dirOrient[ dir ];
    // for mapping the indices to the same
    idxExponent = ( indexLevelUse - bas0->tree().getMaxInsLevel());
    gsMatrix<unsigned> intfc0( lo.rows(), 3 );
    for( index_t i=0; i < lo.rows(); i++)
    {
        intfc0(i,0) = lo(i,dir) << idxExponent;
        intfc0(i,1) = up(i,dir) << idxExponent;
        intfc0(i,2) = level[i];
    }
    intfc0.sortByColumn(0);

    // get the box-representation of the gsHDomain on the interface
    bas1->tree().getBoxesOnSide( bi.second().side(), lo, up, level);
    dir = ( bi.second().direction() + 1 ) % 2;
    idxExponent = ( indexLevelUse - bas1->tree().getMaxInsLevel());
    gsMatrix<unsigned> intfc1( lo.rows(), 3 );
    for( index_t i=0; i < lo.rows(); i++)
    {
        intfc1(i,0) = lo(i,dir) << idxExponent;
        intfc1(i,1) = up(i,dir) << idxExponent;
        intfc1(i,2) = level[i];
    }

    // now the knot indices in intfc0 and intfc1 both correspond to
    // numbering on level "indexLevelUse"

    // get upper corners, but w.r.t. level "indexLevelUse"
    gsVector<unsigned,2> upperCorn0 = bas0->tree().upperCorner();
    upperCorn0[0] = upperCorn0[0] << indexLevelDiff0;
    upperCorn0[1] = upperCorn0[1] << indexLevelDiff0;

    gsVector<unsigned,2> upperCorn1 = bas1->tree().upperCorner();
    upperCorn1[0] = upperCorn1[0] << indexLevelDiff1;
    upperCorn1[1] = upperCorn1[1] << indexLevelDiff1;

    if( !orientPreserv )
    {
        // flip the knot indices
        for( index_t i=0; i < lo.rows(); i++)
        {
            unsigned tmp = upperCorn1[dir] - intfc1(i, 1);
            intfc1(i,1)  = upperCorn1[dir] - intfc1(i, 0);
            intfc1(i,0)  = tmp;
        }
    }
    intfc1.sortByColumn(0);

    GISMO_ASSERT(intfc0( intfc0.rows()-1, 1) == intfc1( intfc1.rows()-1, 1)," Something wrong with interfaces! Mark 264");

    // Merge the knot spans from both sides into intfcU
    // intfcU[i][0]: end-knot-index
    // intfcU[i][1]: level on first()
    // intfcU[i][2]: level on second()
    int i0 = 0; int i1 = 0;
    std::vector< std::vector< unsigned > > intfcU;
    while( i0 < intfc0.rows() && i1 < intfc1.rows() )
    {
        std::vector<unsigned> tmp(3);

        if( intfc0( i0, 1 ) == intfc1( i1, 1 ) )
        {
            tmp[0] = intfc0(i0,1);
            tmp[1] = intfc0(i0,2);
            tmp[2] = intfc1(i1,2);
            intfcU.push_back( tmp );
            i0++;
            i1++;
        }
        else if( intfc0( i0, 1 ) > intfc1( i1, 1 ) )
        {
            tmp[0] = intfc1(i1,1);
            tmp[1] = intfc0(i0,2);
            tmp[2] = intfc1(i1,2);
            intfcU.push_back( tmp );
            i1++;
        }
        else
        {
            tmp[0] = intfc0(i0,1);
            tmp[1] = intfc0(i0,2);
            tmp[2] = intfc1(i1,2);
            intfcU.push_back( tmp );
            i0++;
        }
    }

    // create the refineboxes needed for
    // reparing the interface
    unsigned knot0;
    unsigned knot1 = 0;
    std::vector<unsigned> refElts0;
    std::vector<unsigned> refElts1;

    for( unsigned i=0; i < intfcU.size(); i++)
    {
        knot0 = knot1;
        knot1 = intfcU[i][0];
        unsigned L0 = intfcU[i][1];
        unsigned L1 = intfcU[i][2];

        if( L0 < L1 ) // refine first()
        {
            refElts0.push_back( L1 );

            // knot indices on level L1:
            unsigned knot0L = knot0 >> ( indexLevelUse - L1 );
            unsigned knot1L = knot1 >> ( indexLevelUse - L1 );

            if( bi.first().side() % 2 == 1 ) // west or south
            {
                if( bi.first().direction() == 0 ) // west
                {
                    refElts0.push_back( 0 );
                    refElts0.push_back( knot0L );
                    refElts0.push_back( 1 );
                    refElts0.push_back( knot1L );
                }
                else // south
                {
                    refElts0.push_back( knot0L );
                    refElts0.push_back( 0 );
                    refElts0.push_back( knot1L );
                    refElts0.push_back( 1 );
                }
            }
            else // east or north
            {
                gsVector<unsigned> upperCornOnLevel(2);
                upperCornOnLevel[0] = ( upperCorn0[0] >> ( indexLevelUse - L1 ) );
                upperCornOnLevel[1] = ( upperCorn0[1] >> ( indexLevelUse - L1 ) );

                if( bi.first().direction() == 0 ) //east
                {
                    refElts0.push_back( upperCornOnLevel[0]-1 );
                    refElts0.push_back( knot0L );
                    refElts0.push_back( upperCornOnLevel[0] );
                    refElts0.push_back( knot1L );
                }
                else
                {
                    refElts0.push_back( knot0L );
                    refElts0.push_back( upperCornOnLevel[1]-1 );
                    refElts0.push_back( knot1L );
                    refElts0.push_back( upperCornOnLevel[1] );
                }
            }
        }
        else if( L0 > L1 ) // refine second()
        {
            refElts1.push_back( L0 );

            // knot indices on level "indexLevelUse":
            unsigned knot0L = knot0;
            unsigned knot1L = knot1;
            // flip, if necessary
            if( !orientPreserv )
            {
                unsigned tmp = knot0L;
                knot0L = ( upperCorn1[dir] - knot1L );
                knot1L = ( upperCorn1[dir] - tmp );
            }
            // push to level L0
            knot0L = knot0L >> ( indexLevelUse - L0 );
            knot1L = knot1L >> ( indexLevelUse - L0 );

            if( bi.second().side() % 2 == 1 ) // west or south
            {
                if( bi.second().direction() == 0 ) // west
                {
                    refElts1.push_back( 0 );
                    refElts1.push_back( knot0L );
                    refElts1.push_back( 1 );
                    refElts1.push_back( knot1L );
                }
                else // south
                {
                    refElts1.push_back( knot0L );
                    refElts1.push_back( 0 );
                    refElts1.push_back( knot1L );
                    refElts1.push_back( 1 );
                }
            }
            else // east or north
            {
                gsVector<unsigned> upperCornOnLevel(2);
                upperCornOnLevel[0] = ( upperCorn1[0] >> ( indexLevelUse - L0 ) );
                upperCornOnLevel[1] = ( upperCorn1[1] >> ( indexLevelUse - L0 ) );

                if( bi.second().direction() == 0 ) //east
                {
                    refElts1.push_back( upperCornOnLevel[0]-1 );
                    refElts1.push_back( knot0L );
                    refElts1.push_back( upperCornOnLevel[0] );
                    refElts1.push_back( knot1L );
                }
                else // north
                {
                    refElts1.push_back( knot0L );
                    refElts1.push_back( upperCornOnLevel[1]-1 );
                    refElts1.push_back( knot1L );
                    refElts1.push_back( upperCornOnLevel[1] );
                }
            }
        }
    }

    if( refElts0.size() > 0 )
        m_bases[ bi.first().patch ]->refineElements( refElts0 );
    if( refElts1.size() > 0 )
        m_bases[ bi.second().patch ]->refineElements( refElts1 );


    //delete bas0;
    //delete bas1;

    return ( ( refElts0.size() > 0 ) || ( refElts1.size() > 0 ) );

}

} // namespace gismo
