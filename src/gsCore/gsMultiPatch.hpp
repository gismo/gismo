#pragma once

#include <iterator>

#include <gsCore/gsBasis.h>
#include <gsCore/gsGeometry.h>

#include <gsUtils/gsCombinatorics.h>

namespace gismo
{

template<class T>
gsMultiPatch<T>::gsMultiPatch(const gsGeometry<T> & geo )
: gsBoxTopology( geo.parDim() )
{
    m_patches.push_back( geo.clone() );
    addBox();
    this->addAutoBoundaries();
}
  
template<class T>
gsMultiPatch<T>::gsMultiPatch( const gsMultiPatch& other )
: gsBoxTopology( other ), m_patches( other.m_patches.size() )
{
    // clone all geometries
    cloneAll( other.m_patches.begin(), other.m_patches.end(),
              this->m_patches.begin() );
}

template<class T>
gsMultiPatch<T>::gsMultiPatch( const std::vector<gsGeometry<T> *>& patches )
        : gsBoxTopology( patches[0]->parDim(), patches.size() ) , m_patches( patches )
{
    this->addAutoBoundaries();
}

template<class T>
gsMultiPatch<T>::gsMultiPatch( const PatchContainer& patches,
            const std::vector<patch_side>& boundary,
            const std::vector<boundaryInterface>& interfaces )
        : gsBoxTopology( patches[0]->parDim(), patches.size(), boundary, interfaces ),
          m_patches( patches ) 
{ }

template<class T>
gsMultiPatch<T>::~gsMultiPatch()
{
    freeAll(m_patches);
}

template<class T>
int gsMultiPatch<T>::geoDim() const 
{
    GISMO_ASSERT( m_patches.size() > 0 , "Empty multipatch object.");
    return m_patches[0]->geoDim();
}

template<class T>
int gsMultiPatch<T>::coDim() const 
{
    GISMO_ASSERT( m_patches.size() > 0 , "Empty multipatch object.");
    return m_patches[0]->geoDim() - m_dim;
}

template<class T>
gsMatrix<T> 
gsMultiPatch<T>::parameterRange(int i) const
{
    return m_patches[i]->basis().support();
}

template<class T>
gsBasis<T> &
gsMultiPatch<T>::basis( std::size_t i ) const
{
    GISMO_ASSERT( i < m_patches.size(), "Invalid patch index requested from gsMultiPatch" );
    return m_patches[i]->basis();
}

template<class T>
std::vector<gsBasis<T> *> gsMultiPatch<T>::basesCopy() const
{
    std::vector<gsBasis<T> *> bb;
    for ( typename Base::const_iterator it = m_patches.begin();
          it != m_patches.end(); ++it ) {
        bb.push_back( ( *it )->basis().clone() );
    }
    return bb ;
}
  
template<class T>
void gsMultiPatch<T>::addPatch( gsGeometry<T>* g ) 
{
    if ( m_dim == -1 ) {
        m_dim = g->parDim();
    } else {
        assert( m_dim == g->parDim() );
    }
    m_patches.push_back( g ) ;
    addBox();
}
  
template<class T>
int gsMultiPatch<T>::findPatchIndex( gsGeometry<T>* g ) const {
    typename PatchContainer::const_iterator it
        = std::find( m_patches.begin(), m_patches.end(), g );
    assert( it != m_patches.end() );
    return it - m_patches.begin();
}
  
template<class T>
void gsMultiPatch<T>::addInterface( gsGeometry<T>* g1, boundary::side s1,
                                    gsGeometry<T>* g2, boundary::side s2 ) {
    int p1 = findPatchIndex( g1 );
    int p2 = findPatchIndex( g2 );
    gsVector<bool> orient( m_dim - 1 );
    orient.setConstant( true );
    gsBoxTopology::addInterface( p1, s1, p2, s2, orient );
}

template<class T>
void gsMultiPatch<T>::uniformRefine(int numKnots)
{
    for ( typename Base::const_iterator it = m_patches.begin();
          it != m_patches.end(); ++it ) 
    {
        ( *it )->uniformRefine(numKnots);
    }
}

  
template<class T>
bool gsMultiPatch<T>::computeTopology( T tol )
{
    // Clear boundary and interface containers
    gsBoxTopology::clear();
    const std::size_t np   = m_patches.size();
    const int         n    = this->geoDim();
    const index_t     sz   = 1 << m_dim;
    // number of vertices on one patch side
    const std::size_t nver = 1 << ( m_dim - 1 );
    bool match(0);

    std::vector<patch_side> all_sides;    // All sides
    std::vector<gsMatrix<T> > all_verts;  // All vertices in the multipatch
    std::vector<std::vector<int> > cverts;// indices of vertices of a hypercube

    gsVector<bool> orient(m_dim-1);     
    gsMatrix<T>    vm1( n, nver ), vm2( n, nver );

    // assumes all patches have the same parameter domain
    gsBasis<T> & basis = m_patches[0]->basis();
    gsMatrix<T>  para  = basis.support();
    gsMatrix<T>  pts( m_dim, sz );

    // Produce all corners of the 01-cube in lex-order
    gsVector<unsigned> v( m_dim ), ones( m_dim );
    ones.setConstant( 2 );
    v.setZero();
    cverts.resize( 2 * m_dim + 1 ); // ignore cverts[0]
    unsigned r = 0;
    do 
    { // Generate vertices of hypercube
        for ( int i = 0; i != m_dim; ++i ) 
        {
            pts( i, r ) = para( i, v[i] );
            cverts[ sideOf( i, v[i] ) ].push_back( r );
        }
        r++;
    } while ( nextLexicographic( v, ones ) );

    // Compute all sides and vertices
    all_verts.reserve( np );
    all_sides.reserve( np * 2 * m_dim );
    for ( std::size_t patch = 0; patch != np; ++patch ) 
    {
        all_verts.push_back( m_patches[patch]->eval( pts ) );
        for ( int i = 1; i <= 2 * m_dim; ++i ) 
        {
            all_sides.push_back( patch_side( patch, i ) );
        }
    }

    while ( all_sides.size() ) // For all remaining sides
    {
        const patch_side side1 = all_sides.back();
        all_sides.pop_back();

        // Get vertices of side1
        for ( std::size_t r = 0; r != nver; ++r )
            vm1.col( r ) = all_verts[side1.patch].col( cverts[side1.side][r] );

        // For all remaining non-matched sides
        for ( typename std::vector<patch_side>::iterator side2 = all_sides.begin();
              side2 != all_sides.end(); ++side2 )
        {
            // Initialize data for of side2
            std::vector<int> verts2(cverts[side2->side]);
            for ( std::size_t r = 0; r != nver; ++r )
                vm2.col( r ) = all_verts[side2->patch].col( cverts[side2->side][r] );

            // Try to match side1 and side2 vertex by vertex
            std::size_t r1(0);
            do
            {
                match = false;
                for ( std::size_t r2 = r1; r2 != nver; ++r2 )
                    if ( ( vm1.col(r1) - vm2.col(r2) ).norm() < tol )
                    {
                        vm2.col(r1).swap(vm2.col(r2));
                        std::swap(verts2[r1], verts2[r2]);
                        match = true ;
                        break;
                    }
            }
            while ( (++r1 < nver) && match);

            if (match) // Did we get an interface ?
            {
                // Compute orientation
                int c(0);
                for ( int i = 0; i<m_dim-1; ++i )
                {
                    orient[i] = ( verts2[0] < verts2[(1<<i)] );

                    if ( ( vm1.col( 0 ) - vm1.col( (1<<i) ) ).norm() < tol ) 
                        c++;

                    if ( ( vm1.col(nver-1) - vm1.col(nver-1-(1<<i) ) ).norm() < tol ) 
                        c++;                      
                }

                if (c) // Check for collapsing edge/face
                {
                    gsWarn << "Detected "<< (c==2*m_dim-2 ? 1 : c ) <<" collapsing "
                           << (c==2*m_dim-2 ? "face" : "edges ") 
                           <<" on the side, topology might be incorrect: "
                           << side1 << " "<< *side2<<  "\n";
                }

                // Add interface
                gsBoxTopology::addInterface( boundaryInterface( side1, *side2, orient ) );
                  
                // remove this side and  break loop
                std::swap( *side2, all_sides.back() );
                all_sides.pop_back();
                break;
            }
        }

        // Side not matched, so it is a boundary side
        if ( !match ) 
        {
            gsBoxTopology::addBoundary( side1 );
            //gsDebug<<"Added boundary "<< m_boundary.back() <<"\n";
        }
    }
    return true;
}
  
}
