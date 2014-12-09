/** @file gsHTensorBasis.h

    @brief Provides implementation of HTensorBasis common operations.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): G. Kiss, A. Mantzaflaris, J. Speh
*/

#pragma once 

#include <gsUtils/gsMesh/gsMesh.h>

namespace gismo
{

template<unsigned d, class T>
gsMatrix<T> gsHTensorBasis<d,T>::support() const
{
    return m_bases[0]->support();
}

//i in continuous numbering
template<unsigned d, class T>
gsMatrix<T> gsHTensorBasis<d,T>::support(const unsigned & i) const
{
    // Get the level
    int lvl = levelOf(i);
    // Return the the support
    return m_bases[lvl]->support( m_xmatrix[lvl][ i - m_xmatrix_offset[lvl] ] );
}

template<unsigned d, class T> inline
int gsHTensorBasis<d,T>::get_max_inserted_level() const
{
    return m_tree.getMaxInsLevel();
}

// S.K.
template<unsigned d, class T> inline
gsMatrix<int> gsHTensorBasis<d,T>::getLevelAtPoint( gsMatrix<T> Pts ) const
{
    gsMatrix<int> ret( 1, Pts.cols() );
    gsVector<unsigned int> loIdx( Dim );
    gsVector<unsigned int> upIdx( Dim );

    int maxLevel = m_tree.getMaxInsLevel();

    for( int j = 0; j < Pts.cols(); j++ )
    {
        loIdx.setZero();
        upIdx.setZero();
        for( int i =0; i < Dim; i++)
        {
            loIdx[i] = m_bases[maxLevel]->component(i).knots().Uniquefindspan( Pts(i,j) );
            upIdx[i] = loIdx[i] + 1;
        }
        ret(0,j) = m_tree.query3( loIdx, upIdx, maxLevel);
    }
    return ret;
}

template<unsigned d, class T> inline
void gsHTensorBasis<d,T>::numActive(const gsMatrix<T> & u, gsVector<unsigned>& result) const 
{
    result.resize( u.cols() );
    result.setZero();

    gsVector<unsigned,d> low, upp;
    gsCombinat<gsVector<unsigned,d> > c;

    for(index_t p = 0; p < u.cols(); p++ ) //for all input points
    {
        for(int i = 0; i != d; ++i)
        {
            low[i] = m_bases.back()->component(i).knots().Uniquefindspan(u(i,p));
            upp[i] = low[i]+1;
        }        
        //identify the maximum necessary level
        int lvl = m_tree.query3(low, upp,m_bases.size()-1);

        for(int i = 0; i <= lvl; i++)
        {
            m_bases[i]->active_cwise(u.col(p), low, upp);
            c.first_lattice_point(low,upp,low);
            do
            {
                CMatrix::const_iterator it = 
                    m_xmatrix[i].find_it_or_fail( m_bases[i]->index(low) );
                
                if( it != m_xmatrix[i].end() )// if index is found
                    result[p]++;
            }
            while( c.next_lattice_point(low) );
        }
    }
}


template<unsigned d, class T>
void gsHTensorBasis<d,T>::connectivity(const gsMatrix<T> & nodes, gsMesh<T> & mesh) const 
{
    const index_t sz  = size();
    GISMO_ASSERT( nodes.rows() == sz, "Invalid input.");

    // Add vertices
    for(index_t i = 0; i< sz; ++i )
        mesh.addVertex( nodes.row(i).transpose() );

    const gsVector<unsigned,d> & low = gsVector<unsigned,d>::Zero();    

    // For all levels
    for(int lvl = 0; lvl <= get_max_inserted_level(); lvl++)
    {
        const gsTensorBSplineBasis<d,T,gsCompactKnotVector<T> > & bb = *m_bases[lvl];
        const CMatrix & cmat = m_xmatrix[lvl];
        
        // Last tensor-index in level lvl
        gsVector<unsigned, d> end(d);
        for (unsigned i = 0; i < d; ++i)
            end(i) = bb.size(i) - 1;

        unsigned k, s;
        gsVector<unsigned,d> v, upp;
        for (unsigned i = 0; i < d; ++i) // For all axes
        {
            s      = bb.stride(i);
            v      = low;
            upp    = end;
            upp[i] = 0; // suppress to face v[i]==0
            
            do // Insert all edges normal to axis i 
            {
                k = bb.index(v);
                for (unsigned j = 0; j != end[i]; ++j)
                {
                    if ( cmat.bContains( k )  && cmat.bContains( k+s ) )
                    {
                        // inefficient for now
                        const index_t kInd =  m_xmatrix_offset[lvl] + 
                            (std::lower_bound(cmat.begin(), cmat.end(), k ) 
                             - cmat.begin() );
                        
                        // inefficient for now
                        const index_t kNextInd =  m_xmatrix_offset[lvl] + 
                            (std::lower_bound(cmat.begin(), cmat.end(), k+s ) 
                             - cmat.begin() );
                        
                        mesh.addEdge(kInd, kNextInd );
                    }
                    k += s ;
                }
            }
            while ( nextCubePoint(v, low, upp) );
        }   
    }
}
    
template<unsigned d, class T>
int gsHTensorBasis<d,T>::size() const 
{
    return m_xmatrix_offset.back();
}
template<unsigned d, class T>
void gsHTensorBasis<d,T>::refine_withCoefs(gsMatrix<T> & coefs, gsMatrix<T> const & boxes)
{
    std::vector<gsSortedVector<unsigned> > OX = m_xmatrix;
    refine(boxes);
    gsMatrix<> transf;
    this->transfer(OX, transf);
    coefs = transf*coefs;
}

template<unsigned d, class T>
void gsHTensorBasis<d,T>::refineElements_withCoefs(gsMatrix<T> & coefs,std::vector<unsigned> const & boxes){
    std::vector<gsSortedVector<unsigned> > OX = m_xmatrix;
    refineElements(boxes);
    gsMatrix<> transf;
    this->transfer(OX, transf);
    coefs = transf*coefs;
}

template<unsigned d, class T>
void gsHTensorBasis<d,T>::uniformRefine_withCoefs(gsMatrix<T>& coefs, int numKnots){
    std::vector<gsSortedVector<unsigned> > OX = m_xmatrix;
    //uniformRefine(numKnots);
    //gsMatrix<> transf;
    //CMatrix temp;
    //this->m_xmatrix.insert(this->m_xmatrix.begin(), temp);
    //gsVector<unsigned> level;
    //gsMatrix<unsigned> p1, p2;
    //this->m_tree.getBoxes(p1,p2,level);
    std::vector<unsigned> boxes;
    unsigned lvl;
    
    for ( typename hdomain_type::literator it = m_tree.beginLeafIterator(); it.good(); it.next() )
    {
//        gsInfo <<" level : "<< it.level() <<"\n";
//        gsInfo <<" lower : "<< it.lowerCorner() <<"\n";
//        gsInfo <<" upper : "<< it.upperCorner() <<"\n";

        lvl = it.level() + 1; 
        const point l = it.lowerCorner();
        const point u = it.upperCorner();

        boxes.push_back(lvl);
        boxes.push_back( l(0) * 2);
        boxes.push_back( l(1) * 2);
        boxes.push_back( u(0) * 2);
        boxes.push_back( u(1) * 2);
    }

    safe(this->clone())->refineElements_withCoefs(coefs, boxes);
    this->uniformRefine(numKnots);
    //this->m_xmatrix.erase(this->m_xmatrix.begin(),this->m_xmatrix.begin()+1);
    //coefs = transf*coefs;
}

template<unsigned d, class T>
void gsHTensorBasis<d,T>::refine(gsMatrix<T> const & boxes) 
{
    GISMO_ASSERT(boxes.rows() == d, "refine() needs d rows of boxes.");
    GISMO_ASSERT(boxes.cols()%2 == 0, "Each box needs two corners but you don't provied refine() with them.");

    gsMatrix<T> para = support();
    for(int i = 0; i < boxes.cols()/2; i++)
    {
        for( unsigned j = 0; j < d; j++ )
        {
            GISMO_ASSERT( para(j,0) <= boxes(j, 2*i) , 
                          "In refine() the first corner is outside the computational domain.");
            GISMO_ASSERT( para(j,1) >= boxes(j, 2*i+1), 
                          "In refine() the second corner is outside the computational domain." );
        }
    }

    //enlarge for boxes of maxlevel
    gsVector<unsigned,d> k1;
    gsVector<unsigned,d> k2;

    for(int i = 0; i < boxes.cols()/2; i++)
    {
        for(int j = 0; j < k1.size();j++){
            k1[j] = m_bases[m_bases.size()-1]->component(j).knots().Uniquefindspan(boxes(j,2*i));
        }
        for(int j = 0; j < k2.size();j++){
            k2[j] = m_bases[m_bases.size()-1]->component(j).knots().Uniquefindspan(boxes(j,2*i+1))+1;
        }

        int level = m_tree.query3(k1,k2,m_bases.size()-1);
        for(int j = 0; j < k1.size();j++){
            k1[j] = m_bases[level+1]->component(j).knots().Uniquefindspan(boxes(j,2*i));
        }
        for(int j = 0; j < k2.size();j++){
            k2[j] = m_bases[level+1]->component(j).knots().Uniquefindspan(boxes(j,2*i+1))+1;
        }

        insert_box(k1, k2, level+1);
    }

/*
// to do
for(index_t i = 0; i < boxes.cols()/2; i++)
{
// Locate box in finest level
for(unsigned j = 0; j < d ;j++)
{
k1[j] = m_bases.back()->component(j).knots().Uniquefindspan(boxes(j,2*i));
k2[j] = m_bases.back()->component(j).knots().Uniquefindspan(boxes(j,2*i+1))+1;
}

// Get the actual level of the box
const int level = m_tree.query3(k1, k2, m_bases.size()-1 );

// Sink coordinates one level higher
for(unsigned j = 0; j < d ;j++)
{
k1[j] = ( k1[j] << 1 );
// equiv: = m_bases[level+1]->component(j).knots().Uniquefindspan(boxes(j,2*i));
k2[j] = ( k2[j] << 1 );
// equiv: = m_bases[level+1]->component(j).knots().Uniquefindspan(boxes(j,2*i+1))+1;
}

// Insert the box to the next level
insert_box(k1, k2, level+1);
}
*/

    update_structure();
}

template<unsigned d, class T>
void gsHTensorBasis<d,T>::refineElements(std::vector<unsigned> const & boxes)
{
    gsVector<unsigned int, d> i1;
    gsVector<unsigned int, d> i2;

    GISMO_ASSERT( (boxes.size()%(2*d + 1))==0,
                  "The points did not define boxes properly. The boxes were not added to the basis.");
    for( unsigned int i = 0; i < (boxes.size())/(2*d+1); i++)
    {
        for( unsigned j = 0; j < d; j++ )
        {
            i1[j] = boxes[(i*(2*d+1))+j+1];
            i2[j] = boxes[(i*(2*d+1))+d+j+1];
        }
        insert_box(i1,i2,boxes[i*(2*d+1)]);
    }
    
    update_structure();
}

//protected functions



/// Construct the characteristic matrix of level \a level ; i.e., set all the matrix entries corresponding to active functions to one and the rest to zero.
template<unsigned d, class T>
void gsHTensorBasis<d,T>::set_activ1(int level)
{
    //gsDebug<<" Setting level "<< level <<"\n";
    gsVector<unsigned,d> low, upp;

    CMatrix & cmat = m_xmatrix[level];
    
    // Clear previous entries
    cmat.clear();

    gsCombinat<gsVector<typename gsCompactKnotVector<T>::const_iterator,d> > it; // vector of iterators. one iterator for each dimension
    gsVector<typename gsCompactKnotVector<T>::const_iterator,d> ends, curr;
    gsVector<unsigned,d> ind;

    for(unsigned i = 0; i != d; ++i)
    {
        curr[i] = m_bases[level]->component(i).knots().begin();         // beginning of the iteration in i-th direction
        ends[i] = m_bases[level]->component(i).knots().end()-(m_deg[i]+2);// end of the iteration in i-th direction
    }

    it.first_lattice_point(curr, ends, curr); // This is crucial, since it sets the ends to be the ends of iteration.
    do
    {
        for(unsigned i = 0; i != d; ++i)
        {
            low[i]  = curr[i].span;   // lower left corner of the span of the function
            upp[i]  = (curr[i]+m_deg[i]+1).span; // upper right corner of the span of the function
            ind[i]  = curr[i].index; // index of the function in the matrix
        }

        if ( m_tree.query3(low, upp,level) == level) //if active
            cmat.push_unsorted( m_bases[level]->index( ind ) );
    }
    while (  it.next_lattice_point(curr) ); // while there are some functions (i.e., some combinations of iterators) left

    cmat.sort();
}

template<unsigned d, class T>
void gsHTensorBasis<d,T>::setActiveToLvl(int level, std::vector<gsSortedVector<unsigned> >& x_matrix_lvl){
    x_matrix_lvl.resize(level+1);

    gsVector<unsigned,d> low, upp;
    for(int j =0; j < level+1; j++){
        int counterA = 0;
        int counterB = 0;
        //CMatrix & cmat = x_matrix_lvl[j];

        // Clear previous entries
        x_matrix_lvl[j].clear();

        gsCombinat<gsVector<typename gsCompactKnotVector<T>::const_iterator,d> > it; // vector of iterators. one iterator for each dimension
        gsVector<typename gsCompactKnotVector<T>::const_iterator,d> ends, curr;
        gsVector<unsigned,d> ind;

        for(unsigned i = 0; i != d; ++i)
        {
            curr[i] = m_bases[j]->component(i).knots().begin();         // beginning of the iteration in i-th direction
            ends[i] = m_bases[j]->component(i).knots().end()-(m_deg[i]+2);// end of the iteration in i-th direction
        }

        it.first_lattice_point(curr, ends, curr); // This is crucial, since it sets the ends to be the ends of iteration.
        do
        {
            for(unsigned i = 0; i != d; ++i)
            {
                low[i]  = curr[i].span;   // lower left corner of the span of the function
                upp[i]  = (curr[i]+m_deg[i]+1).span; // upper right corner of the span of the function
                ind[i]  = curr[i].index; // index of the function in the matrix
            }
            if(j < level){
                if ( m_tree.query3(low, upp,j) == j){ //if active
                    counterA ++;
                    x_matrix_lvl[j].push_unsorted( m_bases[j]->index( ind ) );
                }
            }else{
                if ( m_tree.query3(low, upp,j) >= j){ //if active
                    counterB ++;
                    x_matrix_lvl[j].push_unsorted( m_bases[j]->index( ind ) );
                }
            }
        }
        while (  it.next_lattice_point(curr) ); // while there are some functions (i.e., some combinations of iterators) left

        x_matrix_lvl[j].sort();
    }

}

///private functions
template<unsigned d, class T> inline
void gsHTensorBasis<d,T>::insert_box(gsVector<unsigned,d> const & k1,
                                     gsVector<unsigned,d> const & k2, int lvl)
{
    // Remember box in History (for debugging)
    // m_boxHistory.push_back( box(k1,k2,lvl) );

    m_tree.insertBox(k1,k2, lvl);

    // updateTensorLevels();

/*
// Note: We allow insertion boxes not to be fully in the domain
for (unsigned dim = 0; dim < d; dim++)
{
const gsCompactKnotVector<T>& kv =
this->m_bases[lvl]->component(dim).knots();

if (kv.uSize() <= k2(dim) || kv.uSize() <= k1(dim))
{
gsWarn << "You want insert a box outside the domain. "
<< "Box [" << k1.transpose()
<< "] x [" << k2.transpose() << "] is not inserted."
<< std::endl;
return;
}
}
*/

}

template<unsigned d, class T>
void gsHTensorBasis<d,T>::makeCompressed()
{
    // Compress the tree
    // m_tree.makeCompressed();

    while ( ! m_xmatrix_offset[1] )
    {
        delete m_bases.front();
        m_bases.erase( m_bases.begin() );
        m_tree.decrementLevel();
        m_xmatrix.erase( m_xmatrix.begin() );
        m_xmatrix_offset.erase( m_xmatrix_offset.begin() );
    }
    
}

template<unsigned d, class T>
void gsHTensorBasis<d,T>::flatTensorIndexesToHierachicalIndexes(gsSortedVector< int > & indexes,const int level) const
{
    CMatrix::const_iterator xmat_pointer = m_xmatrix[level].begin();
    CMatrix::const_iterator xmat_end = m_xmatrix[level].end();
    gsSortedVector< int >::iterator ind_pointer = indexes.begin();
    gsSortedVector< int >::iterator ind_end = indexes.end();
    unsigned index = 0;
    while(ind_pointer!=ind_end&&xmat_pointer!=xmat_end)
    {
        if(*ind_pointer<static_cast<int>(*xmat_pointer))
        {
            (*ind_pointer)=-1;
            ++ind_pointer;
        }
        else if(*ind_pointer==static_cast<int>(*xmat_pointer))
        {
            (*ind_pointer)=m_xmatrix_offset[level]+index;
            ++xmat_pointer;
            ++index;
            ++ind_pointer;
        }
        else
        {
            ++xmat_pointer;
            ++index;
        }
    }
    while(ind_pointer!=ind_end)
    {
        (*ind_pointer)=-1;
        ++ind_pointer;
    }
}

template<unsigned d, class T>
int gsHTensorBasis<d,T>::flatTensorIndexToHierachicalIndex(unsigned index,const int level) const
{
    if( m_xmatrix.size()<=static_cast<unsigned>(level) )
        return -1;
    CMatrix::const_iterator it = std::lower_bound(m_xmatrix[level].begin(), m_xmatrix[level].end(), index );
    if(it == m_xmatrix[level].end() || *it != index)
        return -1;
    else
        return m_xmatrix_offset[level] + ( it - m_xmatrix[level].begin() );
}

template<unsigned d, class T>
void gsHTensorBasis<d,T>::activeBoundaryFunctionsOfLevel(const unsigned level,const boxSide & s,std::vector<bool>& actives) const
{
    const gsMatrix<unsigned> * bound = m_bases[level]->boundary(s);
    const index_t sz = bound->rows();
    //gsSortedVector< int > indexes(bound->data(),bound->data()+sz);
    gsSortedVector< int > indexes;
    indexes.resize(sz,-1);
    if(static_cast<int>(level)<=get_max_inserted_level())
    {
        for(index_t i = 0;i<sz;++i)
            indexes[i]=(*bound)(i,0);
        flatTensorIndexesToHierachicalIndexes(indexes,level);
    }
    actives.resize(indexes.size(),false);
    std::fill (actives.begin(),actives.end(),false);
    for(unsigned i = 0;i<indexes.size();i++)
        if(indexes[i]!=-1)
            actives[i]=true;
    delete bound;
}

template<unsigned d, class T>
void gsHTensorBasis<d,T>::update_structure() // to do: rename as updateCharMatrices
{
    // Setup the characteristic matrices
    m_xmatrix.clear();
    m_xmatrix.resize(m_tree.getMaxInsLevel()+1);

    // Compress the tree
    m_tree.makeCompressed();

    for(std::size_t i = 0; i != m_xmatrix.size(); i ++)
        set_activ1(i);
    
    // Compute offsets
    m_xmatrix_offset.clear();
    m_xmatrix_offset.reserve(m_xmatrix.size()+1);
    m_xmatrix_offset.push_back(0);
    for (std::size_t i = 0; i != m_xmatrix.size(); i++)
    {
        m_xmatrix_offset.push_back( 
            m_xmatrix_offset.back() + m_xmatrix[i].size() );
    }
}

template<unsigned d, class T>
void gsHTensorBasis<d,T>::createMoreLevels(int numLevels) const
{
    // to do: check for overflow

    for ( int i = 0; i < numLevels; ++i )
    {
        gsTensorBSplineBasis<d,T,gsCompactKnotVector<T> > 
            * next_basis = m_bases.back()->clone();
        next_basis->uniformRefine(1);
        m_bases.push_back (next_basis);
    }   
}

template<unsigned d, class T>
void gsHTensorBasis<d,T>::updateTensorLevels()
{
    const int extraLevels = get_max_inserted_level() + 1 - m_bases.size();

    createMoreLevels(extraLevels);
}

template<unsigned d, class T>
void gsHTensorBasis<d,T>::initialize_class(gsBasis<T> const&  tbasis, int nlevels)
{
    // Degrees
    m_deg.resize(d);

    for( unsigned i = 0; i < d; i++)
        m_deg[i] = tbasis.degree(i);

    // Construct the initial basis
    if ( const gsTensorBSplineBasis<d,T,gsKnotVector<T> > * tb = 
         dynamic_cast<const gsTensorBSplineBasis<d,T,gsKnotVector<T> >*>(&tbasis) )
    {
        std::vector<gsBSplineBasis<T, gsCompactKnotVector<T> > * > cw_bases(d);
        
        for ( unsigned i = 0; i!=d; ++i )
        {
            cw_bases[i]=
                new gsBSplineBasis<T, gsCompactKnotVector<T> >( 
                    gsCompactKnotVector<T>( tb->component(i).knots()) );
        }

        m_bases.push_back(
            new gsTensorBSplineBasis<d,T,gsCompactKnotVector<T> >(cw_bases) );
    }
    else if ( const gsTensorBSplineBasis<d,T,gsCompactKnotVector<T> > * tb = 
              dynamic_cast<const gsTensorBSplineBasis<d,T,gsCompactKnotVector<T> >*>(&tbasis) )
    {
        m_bases.push_back( tb->clone() );
    }
    else
    {
        GISMO_ERROR("Cannot construct a Hierarchical basis from "<< tbasis );
    }
   
    // Initialize the binary tree
    point upp;
    for ( unsigned i = 0; i!=d; ++i )
        upp[i] = m_bases[0]->knots(i).uSize()-1;

    m_tree.init(upp, nlevels-1);

// /*
    // Produce tensor-product spaces by dyadic refinement for all levels
    m_bases.reserve(nlevels);
    for(unsigned int i = 1; i <=m_tree.getIndexLevel(); i++)
    {
        gsTensorBSplineBasis<d,T,gsCompactKnotVector<T> > 
            * next_basis = m_bases[i-1]->clone();
        next_basis->uniformRefine(1);
        m_bases.push_back( next_basis );
    }
// */

    undefined_value =unsigned (std::numeric_limits<int>::max());

}


template<unsigned d, class T>
void gsHTensorBasis<d,T>::active_into(const gsMatrix<T> & u, gsMatrix<unsigned>& result) const
{
    gsVector<unsigned,d> low, upp;
    gsCombinat<gsVector<unsigned,d> > c;

    std::vector<std::vector<unsigned> > temp_output;//collects the outputs
    temp_output.resize( u.cols() );

    for(index_t p = 0; p < u.cols(); p++) //for all input points
    {
        const gsMatrix<T> & curr = u.col(p);
        for(unsigned i = 0; i != d; ++i)
        {
            low[i] = m_bases.back()->component(i).knots().Uniquefindspan( curr(i,0) );
            upp[i] = low[i]+1;
        }
        // Identify the maximum necessary level
        int lvl = m_tree.query3(low, upp, m_bases.size()-1);

        for(int i = 0; i <= lvl; i++)
        {
            m_bases[i]->active_cwise(curr, low, upp);
            c.first_lattice_point(low,upp,low);
            do
            {
                CMatrix::const_iterator it = 
                    m_xmatrix[i].find_it_or_fail( m_bases[i]->index(low) );
                        
                if( it != m_xmatrix[i].end() )// if index is found
                {
                    temp_output[p].push_back(
                        this->m_xmatrix_offset[i] + (it - m_xmatrix[i].begin() )
                        );
                }
            }
            while( c.next_lattice_point(low) );
        }
    }

    int index =0;
    //finding the longest vector to set the matrix size
    for(std::size_t i = 0; i < temp_output.size();i++)
        if(temp_output[i].size()> temp_output[index].size())
            index = i;

    result.resize(temp_output[index].size(),u.cols());
    for(index_t i = 0; i < result.cols(); i++)
        for (index_t j = 0; j < result.rows();j++)
            if (std::size_t(j) < temp_output[i].size())
                result(j,i) = temp_output[i][j];
            else
                result(j,i) = 0 ;
}

template<unsigned d, class T>
gsMatrix<unsigned> *  gsHTensorBasis<d,T>::boundary( ) const
{
    std::vector<unsigned> temp;
    gsVector<unsigned,d>  ind;
    for(int i = 0; i <= this->get_max_inserted_level(); i++)
        for (CMatrix::const_iterator it = m_xmatrix[i].begin(); 
             it != m_xmatrix[i].end(); it++)
        {
            ind = this->m_bases[i]->tensorIndex(*it);
            for (unsigned j=0; j!=d; ++j )
                if ( (ind[j]==0) || (ind[j]==unsigned(this->m_bases[i]->size(j)-1)) )
                {
                    temp.push_back(m_xmatrix_offset[i] + (it-m_xmatrix[i].begin()) );
                    break;
                }
        }
    return makeMatrix<unsigned>(temp.begin(),temp.size(),1 ).release();
}

template<unsigned d, class T>
gsMatrix<unsigned> *  gsHTensorBasis<d,T>::boundary(boxSide const & s ) const
{ 
    //get information on the side
    int k   = direction(s);
    int par = parameter(s);
    
    std::vector<unsigned> temp;
    gsVector<unsigned,d>  ind;
    // i goes through all levels of the hierarchical basis
    for(int i = 0; i <= this->get_max_inserted_level(); i++)
    {
        unsigned r = ( par ? this->m_bases[i]->size(k) - 1 : 0);
        for (CMatrix::const_iterator it = m_xmatrix[i].begin(); 
             it != m_xmatrix[i].end(); it++)
        {
            ind = this->m_bases[i]->tensorIndex(*it);
            if ( ind[k]==r )
                temp.push_back( 
                    m_xmatrix_offset[i] +  (it - m_xmatrix[i].begin() )
                    );
        }
    }
    return makeMatrix<unsigned>(temp.begin(),temp.size(),1 ).release();
}

template<unsigned d, class T>
void gsHTensorBasis<d,T>::evalAllDers_into(const gsMatrix<T> & u, int n, 
                                           gsMatrix<T>& result) const
{
    GISMO_ASSERT( n<=2, "gsHTensorBasis::evalAllDers() not implemented for n > 2." );
    
    gsMatrix<T> values;
    gsVector<unsigned> tmp;
    numActive(u, tmp);
    unsigned nb = tmp.maxCoeff();
    
    if( n == 0 )
        result.resize( nb, u.cols());
    else if( n == 1 )
        result.resize((1 + d)*nb, u.cols());
    else if (n == 2 )
        result.resize((1 + d + d*(d+1)/2 )*nb, u.cols());
    result.setZero();

    this->eval_into(u, values);
    result.topRows(nb) = values;
      
    if ( n >= 1 )
    {
        this->deriv_into(u, values);
        result.middleRows(nb,d*nb) = values;
    }

    if ( n >= 2 )
    {
        this->deriv2_into(u, values);
        result.bottomRows( nb*d*(d+1)/2) = values;
    }

}

template<unsigned d, class T>
void gsHTensorBasis<d,T>::uniformRefine(int numKnots)
{
    GISMO_ASSERT(numKnots == 1, "Only implemented for numKnots = 1"); 

    GISMO_ASSERT( m_tree.getMaxInsLevel() < static_cast<unsigned>(m_bases.size()-1),
                  "Problem with max level, increase number of levels (to do).");

    // Delete the first level
    delete m_bases.front();
    m_bases.erase( m_bases.begin() );

    // Keep consistency of finest level
    gsTensorBSplineBasis<d,T,gsCompactKnotVector<T> > * last_basis 
        = m_bases.back()->clone();
    last_basis->uniformRefine(1);
    m_bases.push_back( last_basis );

    // Lift all indices in the tree by one level
    m_tree.multiplyByTwo();

    update_structure();
}

template<unsigned d, class T>
std::vector<std::vector< std::vector<unsigned int> > > gsHTensorBasis<d,T>::domain_boundariesInKnotIndices(std::vector< std::vector<std::vector< std::vector<unsigned int> > > >& result)const{
    result.clear();
    std::vector<std::vector< std::vector<int> > > res_aabb;
    std::vector<std::vector< std::vector<unsigned int> > > res_aabb_unsigned;


    //std::vector< std::vector<std::vector< std::vector<unsigned int> > > > temp;
    result =  this->m_tree.getPolylines();
    res_aabb.resize( result.size() );

    for(unsigned int i0 = 0; i0 < result.size(); i0++)
    {

        res_aabb[i0].resize( result[i0].size() );
        for(unsigned int i1 = 0; i1 < result[i0].size(); i1++)
        {

            res_aabb[i0][i1].resize( 4 );
            res_aabb[i0][i1][0] = 1000000;
            res_aabb[i0][i1][1] = 1000000;
            res_aabb[i0][i1][2] = -10000000;
            res_aabb[i0][i1][3] = -10000000;
            for(unsigned int i2 = 0; i2 < result[i0][i1].size(); i2++)
            {

                if(res_aabb[i0][i1][0]>signed(result[i0][i1][i2][0])){
                    //res_aabb[i0][i1][0] = result[i0][i1][i2][0];
                    res_aabb[i0][i1][0] = result[i0][i1][i2][0];
                }
                if(res_aabb[i0][i1][1]>signed(result[i0][i1][i2][1])){
                    res_aabb[i0][i1][1] = result[i0][i1][i2][1];
                }
                if(res_aabb[i0][i1][2]<signed(result[i0][i1][i2][2])){
                    res_aabb[i0][i1][2] = result[i0][i1][i2][2];
                }
                if(res_aabb[i0][i1][3]<signed(result[i0][i1][i2][3])){
                    res_aabb[i0][i1][3] = result[i0][i1][i2][3];
                }
            }
        }
    }

    res_aabb_unsigned.resize(res_aabb.size());
    for (unsigned int i = 0; i < res_aabb.size(); i++){
        res_aabb_unsigned[i].resize(res_aabb[i].size());
        for (unsigned int j = 0; j < res_aabb[i].size(); j++){
            res_aabb_unsigned[i][j].resize(res_aabb[i][j].size());
            for (unsigned int k = 0; k < res_aabb[i][j].size(); k++){
                if(res_aabb[i][j][k]<0){
                    gsWarn << "conversion form signed to unsigned"<< std::endl;
                }
                res_aabb_unsigned[i][j][k] = res_aabb[i][j][k];
            }
        }
    }
    return res_aabb_unsigned;
}



template<unsigned d, class T>
std::vector<std::vector< std::vector<unsigned int> > > gsHTensorBasis<d,T>::domain_boundaries(std::vector< std::vector<std::vector< std::vector<T> > > >& result)const{
    result.clear();
    std::vector<std::vector< std::vector<int> > > res_aabb;
    std::vector<std::vector< std::vector<unsigned int> > > res_aabb_unsigned;


    std::vector< std::vector<std::vector< std::vector<unsigned int> > > > temp;
    temp =  this->m_tree.getPolylines();
    res_aabb.resize( temp.size() );
    // We cannot simply say
    // result = this->m_tree.getPolylines();
    // because the return value of getPolylines() are vectors of ints and we have no implicit conversion of int to T.

    // We want result to be of the same size as temp. We achieve this here and in the for cycles.
    result.resize(temp.size());

    // Gabor understands these two lines, Dominik doesn't yet =).
    std::vector<T> x_dir(m_bases[m_bases.size()-1]->component(0).knots().unique());
    std::vector<T> y_dir(m_bases[m_bases.size()-1]->component(1).knots().unique());

    for(unsigned int i0 = 0; i0 < temp.size(); i0++)
    {
        result[i0].resize( temp[i0].size() );
        res_aabb[i0].resize( temp[i0].size() );
        for(unsigned int i1 = 0; i1 < temp[i0].size(); i1++)
        {
            result[i0][i1].resize( temp[i0][i1].size() );
            res_aabb[i0][i1].resize( 4 );
            res_aabb[i0][i1][0] = 1000000;
            res_aabb[i0][i1][1] = 1000000;
            res_aabb[i0][i1][2] = -10000000;
            res_aabb[i0][i1][3] = -10000000;
            for(unsigned int i2 = 0; i2 < temp[i0][i1].size(); i2++)
            {
                result[i0][i1][i2].resize(4);
                // We could as well successively push_back() them but this should be slightly more efficient.
                result[i0][i1][i2][0] = (x_dir[temp[i0][i1][i2][0]]);
                result[i0][i1][i2][1] = (y_dir[temp[i0][i1][i2][1]]);
                result[i0][i1][i2][2] = (x_dir[temp[i0][i1][i2][2]]);
                result[i0][i1][i2][3] = (y_dir[temp[i0][i1][i2][3]]);
                if(res_aabb[i0][i1][0]>signed(temp[i0][i1][i2][0])){
                    //res_aabb[i0][i1][0] = result[i0][i1][i2][0];
                    res_aabb[i0][i1][0] = temp[i0][i1][i2][0];
                }
                if(res_aabb[i0][i1][1]>signed(temp[i0][i1][i2][1])){
                    res_aabb[i0][i1][1] = temp[i0][i1][i2][1];
                }
                if(res_aabb[i0][i1][2]<signed(temp[i0][i1][i2][2])){
                    res_aabb[i0][i1][2] = temp[i0][i1][i2][2];
                }
                if(res_aabb[i0][i1][3]<signed(temp[i0][i1][i2][3])){
                    res_aabb[i0][i1][3] = temp[i0][i1][i2][3];
                }
            }
        }
    }

    res_aabb_unsigned.resize(res_aabb.size());
    for (unsigned int i = 0; i < res_aabb.size(); i++){
        res_aabb_unsigned[i].resize(res_aabb[i].size());
        for (unsigned int j = 0; j < res_aabb[i].size(); j++){
            res_aabb_unsigned[i][j].resize(res_aabb[i][j].size());
            for (unsigned int k = 0; k < res_aabb[i][j].size(); k++){
                if(res_aabb[i][j][k]<0){
                    gsWarn << "conversion form signed to unsigned"<< std::endl;
                }
                res_aabb_unsigned[i][j][k] = res_aabb[i][j][k];
            }
        }
    }
    return res_aabb_unsigned;
}


template<unsigned d, class T>
void  gsHTensorBasis<d,T>::transfer(const std::vector<gsSortedVector<unsigned> >& old, gsMatrix<T>& result){
    //gsMatrix<T> > result;
    gsTensorBSplineBasis<d,T, gsCompactKnotVector<T> > T_0_copy = this->tensorLevel(0);
    std::vector< gsSparseMatrix<T,RowMajor> > transfer;
    transfer.resize(this->get_max_inserted_level() );
    for(int i = 0; i < this->get_max_inserted_level();i++)
    {
        //T_0_copy.uniformRefine_withTransfer(transfer[i], 1);
        std::vector<std::vector<T> > knots;
        for(unsigned int dim = 0l; dim < d;dim++)
        {
            gsCompactKnotVector<T> & ckv = m_bases[i]->component(dim).knots();
            gsCompactKnotVector<T> & fkv = m_bases[i + 1]->component(dim).knots();

            std::vector<T> dirKnots;
            _differenceBetweenKnotVectors(ckv, 0, ckv.uSize() - 1,
                                          fkv, 0, fkv.uSize() - 1,
                                          dirKnots);
            knots.push_back(dirKnots);

            gsDebug << "level: " << i << "\n"
                    << "direction: " << dim << "\n";
            gsDebugVar(gsAsMatrix<T>(dirKnots));
        }
        T_0_copy.refine_withTransfer(transfer[i], knots);
    }
    result = this->coarsening_direct(old,this->m_xmatrix, transfer);
}

template<unsigned d, class T>
void gsHTensorBasis<d,T>::increaseMultiplicity(index_t lvl, int dir, T knotValue, int mult)
{
    if (m_bases[lvl]->component(dir).knots().has(knotValue))
    {
        for(unsigned int i =lvl;i < m_bases.size();i++)
            m_bases[i]->component(dir).insertKnot(knotValue,mult);
        update_structure();
    }
    else
    {
        gsWarn<<"Knot value not in the given knot vector."<<std::endl;
    }
}


} // namespace gismo


//template<unsigned d, class T> inline void
//gsHTensorBasis<d,T>::local2globalIndex( gsVector<unsigned,d> const & index,
//                    unsigned lvl,
//                    gsVector<unsigned,d> & result
//                    ) const
//{
//    for ( unsigned i = 0; i!=d; ++i )
//        result[i] = index[i] << (m_tree.m_index_level-lvl) ;
//}

//template<unsigned d, class T> inline void
//gsHTensorBasis<d,T>::global2localIndex( gsVector<unsigned,d> const & index,
//                                        unsigned lvl,
//                                        gsVector<unsigned,d> & result
//    ) const
//{
//    for ( unsigned i = 0; i!=d; ++i )
//        result[i] = index[i] >> (m_tree.m_index_level-lvl) ;
//}

//template<unsigned d, class T> inline
//int gsHTensorBasis<d,T>::get_level(unsigned function) const
//{
//    std::vector<unsigned>::const_iterator it =
//        std::upper_bound(max_size.begin(), max_size.end(), function);
//    GISMO_ASSERT( it < max_size.end(),
//                  "Something went wrong in get_level.\n");
//    return it - max_size.begin() -1;
//}

////get the x, y coordinates for a basis function from its number and level
//template<unsigned d, class T>
//std::vector<unsigned int> gsHTensorBasis<d,T>::get_x_y(int function, int level) const
//{
//    std::vector<unsigned int> result;
//    result.resize(2);
//    function = function - this->max_size[level];//get the numbering in the correct level
//    result[0] = function % (this->m_bases[level]->size(0));//(this->m_bases[level]->component(0).knots().size()-(this->m_deg[0]+1));
//    result[1] = (function - result[0]) / (this->m_bases[level]->size(0));
//    return result;
//}

////copy a 2D vector to a matrix
//template<unsigned d, class T>
//void gsHTensorBasis<d,T>::copy_to_matrix( std::vector<std::vector<T> >  temp, gsMatrix<T>& result) const
//{
//    int index =0;
//    for(unsigned int i = 0; i < temp.size();i++)//finding the longest vector to set the matrix size
//        if(temp[i].size()> temp[index].size())
//            index = i;
//    result.resize(temp[index].size(),temp.size());
//    for(int i = 0; i < result.cols(); i++){
//        for (int j = 0; j < result.rows();j++){
//            if (unsigned (j) < temp[i].size()){
//                (result)(j,i) = temp[i][j];
//            }else{
//                (result)(j,i) =0;
//            }
//        }
//    }
//}
