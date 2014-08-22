
#pragma once

#include <gsCore/gsBasis.h>

#include <gsThbs/gsHDomain.h>
#include <gsThbs/gsQuadIndex.h>
#include <gsThbs/gsHDomainIterator.h>
#include <gsThbs/gsHDomainBoundaryIterator.h>

#include <gsDataStructures/gsBoundary.h>

#include <gsNurbs/gsCompactKnotVector.h>
#include <gsNurbs/gsTensorBSplineBasis.h>

#include <limits>


namespace gismo
{


/// Small struct for a hierarchical box
/// Consists of two corners (lower and upper) and a level
template<class point>
struct gsHBox
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    gsHBox( const point & low, const point & upp, int lvl = 0)
    : lower(low), upper(upp), level(lvl)
    { }

    /// The lower left corner of the box
    point lower;

    /// The upper right corner of the box
    point upper;

    /// level in which the box lives
    int level ;

};

struct lvl_coef{
    int pos;
    int unsigned lvl;
    double coef;
};


/**
* Class representing a (scalar) hierarchical tensor
* basis of functions \f$ \mathbb R^d \to \mathbb R \f$.
*
*
* The principal idea for constructing the hierarchical basis is as follows
* (in simplified version):
*
* 1. Take a sequence of simple tensor-product bases \f$B^0,\ B^1,\ldots, B^L\f$.
* Each of these bases \f$ B^\ell \f$ defines a <em>level</em> \f$\ell\f$ of the hierarchy.
* Note that we assume that \f$ B^{k+1} \f$ is always a "finer" basis than \f$ B^k \f$.
* 2. From each of these basis \f$ B^\ell \f$, select a set of basis functions in a very smart way.
* This gives you a set of basis functions \f$S^\ell \subseteq B^\ell \f$ of level \f$\ell\f$.
* 3. Take the union of these sets \f$ H = \bigcup_{\ell = 0,\ldots,L} S^\ell \f$.
* This is your hierarchical basis \f$ H \f$ (assuming that you selected
* the sets of functions \f$ S^\ell \f$ in a smart and appropriate way).
*
* <em>Remark on the numbering of the basis functions of \f$ H \f$:</em>
*
* The functions in \f$ H \f$ have global indices \f$0, \ldots, N\f$.
* The numbering is sorted by levels in the following sense. Let \f$n^\ell\f$ be the number of basis functions
* selected from level \f$\ell\f$ (i.e., \f$ n^\ell = | S^\ell |\f$), then the global
* indices \f$0,\ldots,n^0-1\f$ correspond to functions which are taken from \f$ B^0 \f$,
* indices \f$ n^0,\ldots,n^0+n^1\f$ to functions from \f$B^1\f$ and so forth.
*
*
*    Template parameters
*    \param d is the domain dimension
*    \param T is the coefficient type
*
*    \ingroup basis
*/

template<unsigned d, class T>
class gsHTensorBasis: public gsBasis<T>
{

public:
    /// Shared pointer for gsHTensorBasis
    typedef memory::shared_ptr< gsHTensorBasis > Ptr;
    
    typedef T Scalar_t;

    typedef typename gsHDomain<d>::point point; 

    typedef gsHBox<point> box; 
    
    typedef std::vector< box > boxHistory;

    typedef gsSortedVector< unsigned > CMatrix; // charMatrix_

    typedef gsHDomain<d> hdomain_type;

    /// Dimension of the parameter domain
    static const int Dim = d;
    
public:
    /// Default empty constructor
    gsHTensorBasis() { }

    gsHTensorBasis( gsBasis<T> const&  tbasis, int nlevels)
    {      
        GISMO_ASSERT( nlevels > 0, "Invalid number of levels." );
        initialize_class(tbasis, nlevels);
        // Build the characteristic matrices
        update_structure();
    }

    gsHTensorBasis( gsTensorBSplineBasis<d,T> const&  tbasis, 
                    int nlevels, std::vector<unsigned> boxes)
    {
        GISMO_ASSERT( nlevels > 0, "Invalid number of levels." );
        initialize_class(tbasis, nlevels);
        gsVector<unsigned int> i1;
        gsVector<unsigned int> i2;
        i1.resize(d);
        i2.resize(d);
        // Set all functions to active
        GISMO_ASSERT( (boxes.size()%(2*d+1))==0,
        "The points did not define boxes properly. The basis was created without any domain structure.");

        for( unsigned i = 0; i < (boxes.size()/(2*d+1)); i++)
        {
            for( unsigned int j = 0; j < d; j++)
            {
                i1[j] = boxes[(2*d+1)*i+j+1];
                i2[j] = boxes[(2*d+1)*i+j+d+1];
            }
            insert_box(i1,i2,boxes[i*(2*d+1)]);
        }

        // Build the characteristic matrices
        update_structure();

    }

/**
 * @brief gsHTensorBasis
 * @param tbasis - tensor basis
 * @param nlevels - number of levels
 * @param boxes - matrix containing boxes - each 2x2 submatrix
 * contains the lover left and upper right corner of the box
 *
 * - the level where the box should be inserted is one higher as the
 *                 level where it is completely contained
 */
    gsHTensorBasis( gsTensorBSplineBasis<d,T> const&  tbasis, 
                    int nlevels, gsMatrix<T> const & boxes)
    {
        //assert(boxes.rows() == 2);    //can accept only 2D coordinates- Delete during the generalization of the lac to nD
        GISMO_ASSERT(boxes.rows() == d, "Points in boxes need to be of dimension d.");
        GISMO_ASSERT(boxes.cols()%2 == 0, "Each box needs two corners but you don't provied gsHTensorBasis constructor with them.");
        GISMO_ASSERT( nlevels > 0, "Invalid number of levels." );
        initialize_class(tbasis, nlevels);
        
        gsVector<unsigned int,d> k1;
        gsVector<unsigned int,d> k2;
        
        for(index_t i = 0; i < boxes.cols()/2; i++)
        {
            for(unsigned j = 0; j < d; j++)
            {
                k1[j] = this->m_bases.back()->component(j).knots().Uniquefindspan(boxes(j,2*i));
                k2[j] = this->m_bases.back()->component(j).knots().Uniquefindspan(boxes(j,2*i+1))+1;
            }
            int level = m_tree.query3(k1,k2,m_bases.size()-1);
            for(unsigned j = 0; j < d; j++)
            {
                k1[j] = this->m_bases[level+1]->component(j).knots().Uniquefindspan(boxes(j,2*i));
                k2[j] = this->m_bases[level+1]->component(j).knots().Uniquefindspan(boxes(j,2*i+1))+1;
            }

            insert_box(k1,k2,level+1);

            // Build the characteristic matrices
            update_structure();

        }
}

/**
 * @brief gsHTensorBasis
 * @param tbasis
 * @param nlevels - maximum number of levels
 * @param boxes - matrix containing boxes - eaxh 2x2 submatrix
 * contains the lover left and upper right corner of the box
 * @param levels
 */
    gsHTensorBasis( gsTensorBSplineBasis<d,T> const&  tbasis, int nlevels, 
                    gsMatrix<T> const & boxes, 
                    std::vector<unsigned int> levels)
    {
        GISMO_ASSERT(boxes.rows() == d, "Points in boxes need to be of dimension d.");
        GISMO_ASSERT(boxes.cols()%2 == 0, "Each box needs two corners but you don't provied gsHTensorBasis constructor with them.");
        for(std::size_t i = 0; i < levels.size(); i++){// no level should me bigger then maxlevel
            GISMO_ASSERT(levels[i]< unsigned(nlevels), "All the variables levels[i] needs to be smaller than nlevels.");
        }
        GISMO_ASSERT(unsigned (boxes.cols()/2) <= levels.size(), "We don't have enough levels for the boxes.");
        GISMO_ASSERT( nlevels > 0, "Invalid number of levels." );
        initialize_class(tbasis, nlevels);
        // Set all functions to active
        gsVector<unsigned,d> k1;
        gsVector<unsigned,d> k2;
        
        for(index_t i = 0; i < boxes.cols()/2; i++)
        {
            for(unsigned j = 0; j < d; j++)
            {
                k1[j] = m_bases[levels[i]]->component(j).knots().Uniquefindspan(boxes(j,2*i));
                k2[j] = m_bases[levels[i]]->component(j).knots().Uniquefindspan(boxes(j,2*i+1))+1;
            }

            /* m_boxHistory.push_back( box(k1,k2,levels[i]) );  */                      
            this->m_tree.insertBox(k1,k2, levels[i]);
            
            // Build the characteristic matrices
            update_structure();
        }
    }
    
    /// Copy constructor
    gsHTensorBasis( const gsHTensorBasis & o) : gsBasis<T>(o)
    {
        //max_size         = o.max_size;
        m_xmatrix_offset = o.m_xmatrix_offset;
        undefined_value  = o.undefined_value;
        m_deg            = o.m_deg;
        m_tree           = o.m_tree;
        m_xmatrix        = o.m_xmatrix;
        
        m_bases.reserve( m_bases.size() );
        for( typename std::vector<gsTensorBSplineBasis<d,T,gsCompactKnotVector<T> >* > 
                 ::const_iterator it = 
                 o.m_bases.begin(); it != o.m_bases.end(); ++it )
            m_bases.push_back( (*it)->clone() );
        }
    
    /// Destructor
    virtual ~gsHTensorBasis() 
    { 
        freeAll( m_bases );
    }
    

protected:

    // Members that will be changed/removed
    // TO DO: remove these members after they are not used anymore

    unsigned int undefined_value; // to do: make static..

    std::vector<int> m_deg;

//protected:
public:

    // The good members that will stay

    /// \brief The list of nestes spaces.
    ///
    /// See documentation for the class for details on the underlying
    /// structure.
    ///
    /// Recall that the hierarchical basis is built from
    /// a sequence of underlying bases \f$ B^0, B^1,\ldots, B^L\f$.
    /// These underlying bases are stored in gsHTensorBasis.m_bases,
    /// which is of type std::vector.\n
    /// <em>m_bases[k]</em> stores the pointer to the
    /// (global) tensor-product basis \f$ B^k\f$.
    mutable std::vector<gsTensorBSplineBasis<d,T,gsCompactKnotVector<T> >* > m_bases;

    /// \brief The characteristic matrices for each level.
    ///
    /// See documentation for the class for details on the underlying
    /// structure.
    ///
    /// Characteristic matrices provide information on the relation between\n
    /// the basis functions of this gsHTensorBasis \f$ H \f$ and\n
    /// the tensor-product basis functions of the underlying
    /// tensor-product bases \f$ B^\ell \f$.
    ///
    /// Let <em>vk = m_xmatrix[k]</em>. \em vk is a gsSortedVector. It contains
    /// a list of indices of the basis function of level \em k, i.e., of the
    /// basis functions which "are taken" from \f$B^k\f$.
    /// These indices are stored as the global indices in \f$B^k\f$.
    std::vector< CMatrix > m_xmatrix;

    /// The tree structure of the index space
    hdomain_type m_tree;

    // // Stores the coordinates of all inserted boxes
    // (for debugging purposes)
    // boxHistory m_boxHistory;

protected:

    /// \brief Stores the offsets of active functions for all levels
    ///
    /// See documentation for the class for details on the underlying
    /// structure. As mentioned there, the basis functions of the
    /// hierarchical basis \f$ H\f$ have a global numbering, where
    /// the functions from \f$B^0\f$ come first, then those from \f$B^1\f$,
    /// then \f$B^2\f$, and so forth.
    ///
    /// The entry <em>m_xmatrix_offset[k]</em>
    /// indicates the index from which the basis functions from
    /// level \em k (i.e., those taken from \f$ B^k \f$) start.
    std::vector<unsigned> m_xmatrix_offset;

    //------------------------------------
 
public:

    /// Returns the dimension of the parameter space
    virtual int dim() const 
    { return d; }
    
    /// Returns the number of breaks (distinct knot values) in direction \a k of level \a lvl
    int numBreaks(int lvl, int k) const
    {
        return m_tree.numBreaks(lvl,k);
    }

    /// Returns the number of knots in direction \a k of level \a lvl
    int numKnots(int lvl, int k) const
    {
        if ( lvl - m_bases.size() + 1 > 0 )
            createMoreLevels( lvl - m_bases.size() + 1);

        return m_bases[lvl]->component(k).knots().size();
    }

    /// Returns the \a i-th knot in direction \a k at level \a lvl
    T knot(int lvl, int k, int i)
    {
        if ( lvl - m_bases.size() + 1 > 0 )
            createMoreLevels( lvl - m_bases.size() + 1);

        return m_bases[lvl]->component(k).knot(i);
    }

    /// Returns the anchors points that represent the members of the basis
    virtual void anchors_into(gsMatrix<T> & result) const 
    {
        result.resize( d, this->size()) ;
        unsigned k(0);

        gsVector<unsigned, d> ind;
        for(std::size_t i = 0; i < m_xmatrix.size(); i++)
        {
            for( CMatrix::const_iterator it =
                 m_xmatrix[i].begin(); it != m_xmatrix[i].end(); it++)
            {
                ind = m_bases[i]->tensorIndex(*it);
                for ( unsigned r = 0; r!=d; ++r ) 
                    result(r,k) = m_bases[i]->component(r).knots().greville( ind[r] );
                k++;
            }
        }
    }

    virtual void connectivity(const gsMatrix<T> & nodes, gsMesh<T> & mesh) const;

    // Make gsBasis::connectivity visible
    virtual void connectivity(gsMesh<T> & mesh) const
    { gsBasis<T>::connectivity(mesh); }

  // Prints the characteristic matrices (ie. the indices of all basis
  // functions in the basis)
  void printCharMatrix(std::ostream &os = gsInfo) const
    {  
      os<<"Characteristic matrix:\n";
      for(int i = 0; i<= get_max_inserted_level(); i++)
      {
          if ( m_xmatrix[i].size() )
          {
              os<<"- level="<<i<< 
                  ", size="<<m_xmatrix[i].size() << ":\n";
              os << "("<< m_bases[i]->tensorIndex(*m_xmatrix[i].begin()).transpose() <<")";
              for( CMatrix::const_iterator  it = 
                   m_xmatrix[i].begin()+1; it != m_xmatrix[i].end(); it++)
              {
                os << ", ("<< m_bases[i]->tensorIndex(*it).transpose() <<")";
              }
              os <<"\n";
          }
          else
          {
              os<<"- level="<<i<<" is empty.\n";
          }
      }
  }

  /// Prints the spline-space hierarchy
  void printSpaces(std::ostream &os = gsInfo) const
    {  
      os<<"Spline-space hierarchy:\n";
      for(int i = 0; i<= get_max_inserted_level(); i++)
      {
          if ( m_xmatrix[i].size() )
          {
              os<<"- level="<<i<< 
                  ", size="<<m_xmatrix[i].size() << ":\n";
              os << "Space: "<< * m_bases[i] <<")";
          }
          else
          {
              os<<"- level="<<i<<" is empty.\n";
          }
      }
  }

  /// Prints the spline-space hierarchy
  void printBasic(std::ostream &os = gsInfo) const
  {
      os << "basis of dimension " <<this->dim()<<
          ",\nlevels="<< this->m_tree.max_ins_level+1 <<", size="<<
          this->size()<<", tree_nodes="<< this->m_tree.size() << ".\n";
      gsMatrix<T> supp  = this->support();
      os << "Domain: ["<< supp.col(0).transpose()<< "]..["<<
          supp.col(1).transpose()<< "].\n";
      os <<"Size per level: ";
      for(unsigned i = 0; i<= this->m_tree.max_ins_level; i++)
          os << this->m_xmatrix[i].size()<< " ";
      os<<"\n";
  }
    
    // Returns the active functions at points \a u
    virtual void active_into(const gsMatrix<T> & u, gsMatrix<unsigned>& result) const 
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


/// Returns the indices of the basis functions that touch the domain
/// boundary
virtual gsMatrix<unsigned> * boundary( ) const 
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

/// Returns the indices of the basis functions that touch the domain boundary
virtual gsMatrix<unsigned> * boundary(boundary::side const & s ) const 
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

  void evalAllDers_into(const gsMatrix<T> & u, int n, gsMatrix<T>& result) const
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


  const gsHDomain<d> & tree() const { return m_tree; }

  gsHDomain<d> &       tree()       { return m_tree; }

  /// Cleans the basis, removing any inactive levels 
  void makeCompressed();


  /// Returns the boundary basis for side s
  // virtual gsHTensorBasis<d,T> * boundaryBasis(boundary::side const & s ) const 

  /// Returns a bounding box for the basis' domain
  gsMatrix<T> support() const;

  gsMatrix<T> support(const unsigned & i) const;

  gsMatrix<T> supportInterval(unsigned dir) const;

  /// Clone function. Used to make a copy of a derived basis
  virtual gsHTensorBasis * clone() const = 0;

  /// The number of basis functions in this basis
  int size() const;

  /// The number of nodes in the tree representation
  int treeSize() const
  {
      return m_tree.size();
  }

  /// The number of active basis functions at points \a u
  void numActive(const gsMatrix<T> & u, gsVector<unsigned>& result) const;

  /// The 1-d basis for the i-th parameter component
  virtual gsBSplineBasis<T,gsCompactKnotVector<T> > & component(unsigned i) const
  {
      return m_bases[ this->get_max_inserted_level() ]->component(i);
  }

  /// Returns the tensor basis member of level i 
  gsTensorBSplineBasis<d,T,gsCompactKnotVector<T> > & tensorLevel(unsigned i) const
  {
      return *this->m_bases[i];
  }

  // Refine the basis uniformly by inserting \a numKnots new knots on each knot span
  //virtual void uniformRefine(int numKnots = 1)

  // Refine the basis uniformly and adjust the given matrix of coefficients accordingly
  //virtual void uniformRefine_withCoefs(gsMatrix<T>& coefs, int numKnots = 1)

  // Refine the basis uniformly and produce a sparse matrix which
  // maps coarse coefficient vectors to refined ones
  //virtual void uniformRefine_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, int numKnots = 1)

  /// @brief If the basis is of polynomial or piecewise polynomial
  /// type, then this function returns the polynomial degree.
  virtual inline int degree() const
    { return m_bases[0]->degree(); }

  /// @brief If the basis is a tensor product of (piecewise)
  /// polynomial bases, then this function returns the polynomial
  /// degree of the \a i-th component.
  virtual inline int degree(int i) const
    { return m_bases[0]->degree(i);}

  /// return maximum of inserted level. has to be smaller or equal to m_tree.m_index_level
  int get_max_inserted_level()const;

  /// Returns the level in which the indices are stored internally
  int maxAllowdLevel() const
  {
      return m_tree.m_index_level;
  }

  /// Returns the level of \a function, which is a hier. Id index
  int get_level(unsigned function) const; 

    /// Returns the level of the function indexed \a i (in continued indices)
    inline int levelOf(unsigned i) const
    {
        return std::upper_bound(m_xmatrix_offset.begin(), 
                                m_xmatrix_offset.end(), i)
            - m_xmatrix_offset.begin() - 1;
    }

/*
    const boxHistory & get_inserted_boxes() const
    {
        return m_boxHistory;
    }
*/

    virtual void refine(gsMatrix<T> const & boxes);

    virtual void refine(std::vector<unsigned> const & boxes);

    virtual void uniformRefine(int numKnots = 1)
    {
        GISMO_ASSERT(numKnots == 1, "Only implemented for numKnots = 1"); 

        GISMO_ASSERT( m_tree.max_ins_level < static_cast<unsigned>(m_bases.size()-1), 
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

    typename gsBasis<T>::domainIter makeDomainIterator() const 
    { 
        return typename gsBasis<T>::domainIter(new gsHDomainIterator<T, d>(*this));
    }

    typename gsBasis<T>::domainIter makeDomainIterator(const boundary::side & s) const
    {            
        // Next line iterates over boundary elements of max. inserted level
        // return typename gsBasis<T>::domainIter(
        //     new gsTensorDomainBoundaryIterator<T, d>(*this,s) );
        return typename gsBasis<T>::domainIter(
            new gsHDomainBoundaryIterator<T, d>(*this,s) );
    }

    /// @brief Returns the tensor index of the function indexed \a i
    /// (in continued indices).
    ///
    /// @param[in] i Global (continued) index of a basis function of
    /// the hierarchical basis.
    /// @return The tensor index of this basis function
    /// with respect to the tensor-product basis
    /// of the corresponding level.
    inline
    unsigned flatTensorIndexOf(const unsigned i) const
    {

        const int level = this->levelOf(i);

        const unsigned offset = this->m_xmatrix_offset[level];
        const unsigned ind_in_level = this->m_xmatrix[level][i - offset];

        return ind_in_level;
    }


    /// @brief Returns the tensor index of the function indexed \a i
    /// (in continued indices).
    ///
    /// @param[in] i Global (continued) index of a basis function of
    /// the hierarchical basis.
    /// @param level Level of the i-th basis function.
    /// @return The tensor index of this basis function
    /// with respect to the tensor-product basis
    /// of \em level.
    inline
    unsigned flatTensorIndexOf(const unsigned i, const unsigned level) const
    {

        const unsigned offset = this->m_xmatrix_offset[level];
        const unsigned ind_in_level = this->m_xmatrix[level][i - offset];

        return ind_in_level;
    }

    /// @brief Gives polylines on the boundaries between different levels of the mesh.
    /// @param result variable where to write the polylines in the form
    /// < levels < polylines_in_one_level < one_polyline < one_segment (x1, y1, x2, y2) > > > > ,
    /// where <x1, y1, x2, y2 > are so that (x1, y1) <=LEX  (x2, y2).
    /// @return bounding boxes of the polylines in the form
    /// < levels < polylines_in_one_level < x_ll, y_ll, x_ur, y_ur > > >, where "ur" stands for "upper right" and "ll" for "lower left".
    std::vector<std::vector< std::vector<unsigned int> > > domain_boundaries(std::vector< std::vector<std::vector< std::vector<T> > > >& result)const;

    /// @brief Gives polylines on the boundaries between different levels of the mesh.
    /// @param result variable where to write the polylines in the form
    /// < levels < polylines_in_one_level < one_polyline < one_segment (x1, y1, x2, y2) > > > > -kin knot vector indices!!!!!!!!
    /// where <x1, y1, x2, y2 > are so that (x1, y1) <=LEX  (x2, y2).
    /// @return bounding boxes of the polylines in the form
    /// < levels < polylines_in_one_level < x_ll, y_ll, x_ur, y_ur > > >, where "ur" stands for "upper right" and "ll" for "lower left".
    std::vector<std::vector< std::vector<unsigned int> > > domain_boundariesInKnotIndices(std::vector< std::vector<std::vector< std::vector<unsigned int> > > >& result)const;
    // this function probably needs better implementation...
    int numElements() const
    {
        gsHDomainIterator<T, d> domIter(*this);

        int numEl = 0;
        for (; domIter.good(); domIter.next())
        {
            numEl++;
        }

        return numEl;
    }

    /// @brief transformes a sortedVector \a indexes of flat tensor index
    /// of the bspline basis of \a level to hierachical indexes in place. If a flat
    /// tensor index is not found, it will transform to -1.
    ///
    /// @param[in] indexes flat tensor indexes of the function in level
    /// @param level Level of the basis.
    void flatTensorIndexesToHierachicalIndexes(gsSortedVector< int > & indexes,const int level) const;

    /// @brief takes a flat tensor \a index
    /// of the bspline basis of \a level and gives back the hierachical index. If a flat
    /// tensor index is not found, it will return -1.
    ///
    /// @param[in] index flat tensor index of the function in level
    /// @param level Level of the basis.
    /// @return hierachical index, or -1 if it was not found
    int flatTensorIndexToHierachicalIndex(unsigned index,const int level) const;

    /// Fills the vector actives with booleans, that determine if a function
    /// of the given level is active. The functions on the boundary are ordered
    /// in ascending patchindex order.
    ///
    /// \param[in] level : level of the boundary functions
    /// \param[in] s : boundary side
    /// \param[out] actives : the result, true if its active, false if not
    void activeBoundaryFunctionsOfLevel(const unsigned level,const boundary::side & s,std::vector<bool>& actives) const;

protected:

    /// @brief Updates the basis structure (eg. charact. matrices, etc), to
    /// be called after any modifications.
    void update_structure(); // to do: rename as updateCharMatrices
    
    /// @brief Makes sure that there are enough refinement levels to compute
    /// with the basis
    void updateTensorLevels();

    /// @brief Makes sure that there are enough refinement levels to compute
    /// with the basis
    void createMoreLevels(int numLevels) const;

    /// Computes difference between coarser knot vector (ckv) and finer knot
    /// vector (fkv). Difference is computed just between c_low, c_high
    /// indices f_low, f_high indices for ckv and fkv respectfully. Result
    /// is stored in vector knots.
    ///
    /// @param ckv coarse knot vector
    /// @param c_low low index of the interested area for the ckv
    /// @param c_high high index of the interested area for the ckv
    /// @param fkv finer knot vector
    /// @param f_low low index of the interested area for the fkv
    /// @param f_high high index of the interested area for the fkv
    /// @param knots {k | k \in fkv & k \notin ckv}
    static
    void _differenceBetweenKnotVectors(const gsCompactKnotVector<T>& ckv,
                                       const unsigned c_low,
                                       const unsigned c_high,
                                       const gsCompactKnotVector<T>& fkv,
                                       const unsigned f_low,
                                       const unsigned f_high,
                                       std::vector<T>& knots)
    {

        unsigned f_knot_mltpl; // finest knot multiplicity
        unsigned c_knot_mltpl; // coarse knot multiplicity

        T f_knot, c_knot;// finest and coarsest knot

        unsigned c_index = c_low;
        unsigned f_index = f_low;

        while (f_index <= f_high)
        {
            f_knot = fkv.uValue(f_index); // finest knot
            c_knot = ckv.uValue(c_index); // coarse knot

            f_knot_mltpl = fkv.u_multiplicityIndex(f_index);

            if (f_knot == c_knot)
            {
                c_knot_mltpl = ckv.u_multiplicityIndex(c_index);

                if (c_knot_mltpl < f_knot_mltpl)
                {
                    for (unsigned j = 0; j < f_knot_mltpl - c_knot_mltpl; ++j)
                        knots.push_back(f_knot);
                }

                ++f_index;
                ++c_index;
            }
            else // f_knot < c_knot
            {
                for (unsigned j = 0; j < f_knot_mltpl; ++j)
                    knots.push_back(f_knot);
                ++f_index;
            }
        }
    }

private:

    /// insert a domain into quadtree
    void insert_box(gsVector<unsigned,d> const & k1, 
                    gsVector<unsigned,d> const & k2, int lvl);
    
    void initialize_class(gsBasis<T> const&  tbasis, int nlevels);

    /// set all functions to active or passive- one by one
    void set_activ1(int level);

    /// set all functions to active or passive- recursive
    void set_activ1_rec(std::vector<int> sizes, std::vector<int> variables);



    ///returns a transfer matrix using the characteristic matrix of the old and new basis
    virtual gsMatrix<T> coarsening(std::vector<gsSortedVector<unsigned> > old,  std::vector<gsSortedVector<unsigned> > n, const gsSparseMatrix<T,RowMajor> & transfer) = 0;
    virtual gsMatrix<T> coarsening_direct(std::vector<gsSortedVector<unsigned> > old,std::vector<gsSortedVector<unsigned> > n,  const std::vector<gsSparseMatrix<T,RowMajor> >& transfer) = 0;

 //D
public:
    ///returns transfer matrix betweend the hirarchycal spline given by the characteristic matrix  "old" and this
    gsMatrix<T> transfer (std::vector<gsSortedVector<unsigned> > &old);
    /// create characteristic matrices for basis where "level" is the maximum level i.e. ignoring higher level refinements
    void setActiveToLvl(int level, std::vector<gsSortedVector<unsigned> >& x_matrix_lvl);


//    void local2globalIndex( gsVector<unsigned,d> const & index,
//                            unsigned lvl,
//                            gsVector<unsigned,d> & result
//        ) const;

//    void global2localIndex( gsVector<unsigned,d> const & index,
//                            unsigned lvl,
//                            gsVector<unsigned,d> & result
//        ) const;

}; // class gsHTensorBasis
 

} // namespace gismo

//////////////////////////////////////////////////
//////////////////////////////////////////////////

#ifndef GISMO_HEADERS_ONLY
#include GISMO_HPP_HEADER(gsHTensorBasis.hpp)
#endif
