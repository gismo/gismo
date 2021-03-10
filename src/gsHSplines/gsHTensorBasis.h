/** @file gsHTensorBasis.h

    @brief Provides definition of HTensorBasis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): G. Kiss, A. Mantzaflaris, J. Speh
*/

#pragma once

#include <gsCore/gsBasis.h>

#include <gsHSplines/gsHDomain.h>
#include <gsHSplines/gsHDomainIterator.h>
#include <gsHSplines/gsHDomainBoundaryIterator.h>

#include <gsCore/gsBoundary.h>

#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsNurbs/gsBSplineBasis.h> // for gsBasis::component(short_t)

#include <gsUtils/gsSortedVector.h>

namespace gismo
{

struct lvl_coef
{
    int pos; // flat index at grid of level \a lvl
    int unsigned lvl; // level
    real_t coef; // value of the coefficient (lvl,pos)
};

/**
 * @brief Class representing a (scalar) hierarchical tensor
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
 *    \ingroup HSplines
 */

template<short_t d, class T>
class GISMO_DEFAULT_VIS gsHTensorBasis: public gsBasis<T>
{
public:
    /// Shared pointer for gsHTensorBasis
    typedef memory::shared_ptr< gsHTensorBasis > Ptr;

    /// Unique pointer for gsHTensorBasis
    typedef memory::unique_ptr< gsHTensorBasis > uPtr;

    typedef gsHTensorBasis<d,T> Self_t;

    typedef T Scalar_t;

    typedef gsHDomain<d> hdomain_type;

    typedef typename hdomain_type::point point;

    typedef typename hdomain_type::box   box;

    typedef std::vector< box > boxHistory;

    typedef gsSortedVector< index_t > CMatrix; // charMatrix_

    typedef typename CMatrix::const_iterator cmatIterator;

    typedef gsKnotVector<T> KnotVectorType;

    typedef typename gsBSplineTraits<static_cast<short_t>(d-1),T>::Basis BoundaryBasisType;

    typedef typename gsBSplineTraits<d,T>::Basis tensorBasis;

    /// Dimension of the parameter domain
    static const short_t Dim = d;

public:

    /// Default empty constructor
    gsHTensorBasis()
    {
        initialize_class(tensorBasis());
        update_structure();
    }

    gsHTensorBasis( gsBasis<T> const&  tbasis)
    {
        initialize_class(tbasis);
        // Build the characteristic matrices
        update_structure();
    }

    gsHTensorBasis( gsTensorBSplineBasis<d,T> const&  tbasis,
                    const std::vector<index_t> & boxes)
    {
        initialize_class(tbasis);
        point i1;
        point i2;
        //i1.resize(d);
        //i2.resize(d);
        // Set all functions to active
        GISMO_ASSERT( (boxes.size()%(2*d+1))==0,
                      "The points did not define boxes properly. The basis was created without any domain structure.");

        for( size_t i = 0; i < (boxes.size()/(2*d+1)); i++)
        {
            for( short_t j = 0; j < d; j++)
            {
                i1[j] = boxes[(2*d+1)*i+j+1];
                i2[j] = boxes[(2*d+1)*i+j+d+1];
            }
            insert_box(i1,i2,boxes[i*(2*d+1)]);
        }

        // Build the characteristic matrices (note: call is non-vritual)
        update_structure();
    }

/**
   @brief gsHTensorBasis
   @param tbasis - tensor basis
   @param boxes - matrix containing boxes - each 2x2 submatrix
   contains the lover left and upper right corner of the box

   - the level where the box should be inserted is one higher as the
   level where it is completely contained
*/
    gsHTensorBasis( gsTensorBSplineBasis<d,T> const&  tbasis,
                    gsMatrix<T> const & boxes)
    {
        //assert(boxes.rows() == 2);    //can accept only 2D coordinates- Delete during the generalization of the lac to nD
        GISMO_ASSERT(boxes.rows() == d, "Points in boxes need to be of dimension d.");
        GISMO_ASSERT(boxes.cols()%2 == 0, "Each box needs two corners but you don't provied gsHTensorBasis constructor with them.");
        initialize_class(tbasis);

        point k1;
        point k2;

        for(index_t i = 0; i < boxes.cols()/2; i++)
        {
            for(short_t j = 0; j < d; j++)
            {
                k1[j] = this->m_bases.back()->knots(j).uFind(boxes(j,2*i)).uIndex();
                k2[j] = this->m_bases.back()->knots(j).uFind(boxes(j,2*i+1)).uIndex()+1;
            }
            int level = m_tree.query3(k1,k2,m_bases.size()-1);
            for(short_t j = 0; j < d; j++)
            {
                k1[j] = this->m_bases[level+1]->knots(j).uFind(boxes(j,2*i)).uIndex();
                k2[j] = this->m_bases[level+1]->knots(j).uFind(boxes(j,2*i+1)).uIndex()+1;
            }

            insert_box(k1,k2,level+1);

            // Build the characteristic matrices (note: call is non-vritual)
            update_structure();

        }
    }

/**
 * @brief gsHTensorBasis
 * @param tbasis
 * @param boxes - matrix containing boxes - eaxh 2x2 submatrix
 * contains the lover left and upper right corner of the box
 * @param levels
 */
    gsHTensorBasis( gsTensorBSplineBasis<d,T> const& tbasis,
                    gsMatrix<T> const & boxes,
                    const std::vector<index_t> & levels)
    {
        GISMO_ASSERT(boxes.rows() == d, "Points in boxes need to be of dimension d.");
        GISMO_ASSERT(boxes.cols()%2 == 0, "Each box needs two corners but you don't provied gsHTensorBasis constructor with them.");
        GISMO_ASSERT(unsigned (boxes.cols()/2) <= levels.size(), "We don't have enough levels for the boxes.");

        initialize_class(tbasis);

        gsVector<index_t,d> k1;
        gsVector<index_t,d> k2;

        const size_t mLevel = *std::max_element(levels.begin(), levels.end() );
        needLevel( mLevel );

        for(index_t i = 0; i < boxes.cols()/2; i++)
        {
            for(short_t j = 0; j < d; j++)
            {
                k1[j] = m_bases[levels[i]]->knots(j).uFind(boxes(j,2*i)).uIndex();
                k2[j] = m_bases[levels[i]]->knots(j).uFind(boxes(j,2*i+1)).uIndex()+1;
            }

            /* m_boxHistory.push_back( box(k1,k2,levels[i]) );  */
            this->m_tree.insertBox(k1,k2, levels[i]);

            // Build the characteristic matrices (note: call is non-vritual)
            update_structure();
        }
    }

    /// Copy constructor
    gsHTensorBasis( const gsHTensorBasis & o) : gsBasis<T>(o)
    {
        this->operator=(o);
    }

    gsHTensorBasis& operator=(const gsHTensorBasis & o)
    {
        if ( this != &o )
        {
            m_xmatrix_offset = o.m_xmatrix_offset;
            m_deg            = o.m_deg;
            m_tree           = o.m_tree;
            m_xmatrix        = o.m_xmatrix;

            freeAll( m_bases );
            m_bases.resize( o.m_bases.size() );
            cloneAll(o.m_bases.begin(), o.m_bases.end(), m_bases.begin());
        }
        return *this;
    }

#if EIGEN_HAS_RVALUE_REFERENCES
    gsHTensorBasis(gsHTensorBasis&& other)
    {
        this->operator=(other);
    }

    gsHTensorBasis & operator=(gsHTensorBasis&& other)
    {
        m_deg     = std::move(other.m_deg);
        freeAll( m_bases );
        m_bases   = std::move(other.m_bases);
        m_xmatrix = std::move(other.m_xmatrix);
        m_tree    = std::move(other.m_tree);
        m_xmatrix_offset = std::move(other.m_xmatrix_offset);
        return *this;
    }
#endif

    /// Destructor
    virtual ~gsHTensorBasis()
    {
        freeAll( m_bases );
    }

protected:

    // TO DO: remove these members after they are not used anymore
    std::vector<int> m_deg;

protected:

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
    mutable std::vector<tensorBasis*> m_bases;

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

public:
    // Needed since m_tree is 16B aligned
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
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
    std::vector<index_t> m_xmatrix_offset;

    //------------------------------------

public:
    const std::vector< CMatrix >& getXmatrix() const
    {
        return m_xmatrix;
    }

    /// \brief Returns the tensor B-spline space of all levels.
    /// \return
    ///
    const std::vector<tensorBasis*>& getBases() const
    {
        return m_bases;
    }

    /// Returns the dimension of the parameter space
    virtual short_t dim() const
    { return d; }

    /// Returns the number of breaks (distinct knot values) in
    /// direction \a k of level \a lvl
    int numBreaks(int lvl, int k) const
    {
        return m_tree.numBreaks(lvl,k);
    }

    /// Returns the number of knots in direction \a k of level \a lvl
    int numKnots(int lvl, int k) const
    {
        needLevel(lvl);

        return m_bases[lvl]->knots(k).size();
    }

    /// Returns the \a i-th knot in direction \a k at level \a lvl
    T knot(int lvl, int k, int i) const
    {
        needLevel(lvl);

        return m_bases[lvl]->component(k).knot(i);
        //return m_bases[lvl]->knot(k,i);
    }

    /// Returns the anchors points that represent the members of the
    /// basis
    virtual void anchors_into(gsMatrix<T> & result) const
    {
        result.resize( d, this->size()) ;
        unsigned k(0);

        gsVector<index_t, d> ind;
        for(size_t i = 0; i < m_xmatrix.size(); i++)
        {
            for( CMatrix::const_iterator it =
                     m_xmatrix[i].begin(); it != m_xmatrix[i].end(); it++)
            {
                ind = m_bases[i]->tensorIndex(*it);
                for ( short_t r = 0; r!=d; ++r )
                    result(r,k) = m_bases[i]->knots(r).greville( ind[r] );
                k++;
            }
        }
    }

    virtual void connectivity(const gsMatrix<T> & nodes, gsMesh<T> & mesh) const;
    void connectivity(const gsMatrix<T> & nodes, int level, gsMesh<T> & mesh) const;

    // Prints the characteristic matrices (ie. the indices of all basis
    // functions in the basis)
    void printCharMatrix(std::ostream &os = gsInfo) const
    {
        os<<"Characteristic matrix:\n";
        for(unsigned i = 0; i<= maxLevel(); i++)
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
        for(unsigned i = 0; i<= maxLevel(); i++)
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
            ",\nlevels="<< this->m_tree.getMaxInsLevel()+1 <<", size="
           << this->size()<<", tree_nodes="<< this->m_tree.size()
            //   << ", leaf_nodes="<< this->m_tree.leafSize();
            //const std::pair<int,int> paths  = this->m_tree.minMaxPath();
            //os << ", path lengths=("<<paths.first<<", "<<paths.second
           << ").\n";
        const gsMatrix<T> supp  = this->support();
        os << "Domain: ["<< supp.col(0).transpose()<< "]..["<<
            supp.col(1).transpose()<< "].\n";
        os <<"Size per level: ";
        for(unsigned i = 0; i<= this->m_tree.getMaxInsLevel(); i++)
            os << this->m_xmatrix[i].size()<< " ";
        os<<"\n";
    }

    // Look at gsBasis.h for the documentation of this function
    void active_into(const gsMatrix<T> & u, gsMatrix<index_t>& result) const;


    // Look at gsBasis.h for the documentation of this function
    gsMatrix<index_t> allBoundary( ) const;

    // Look at gsBasis.h for the documentation of this function
    virtual gsMatrix<index_t> boundaryOffset(boxSide const & s, index_t offset ) const;

    // Look at gsBasis.h for the documentation of this function
    // /// \todo impl. evalAllDers_into
    //void evalAllDers_into(const gsMatrix<T> & u, int n,
    //                      std::vector<gsMatrix<T> >& result) const;

    /// Returns a reference to m_tree
    const gsHDomain<d> & tree() const { return m_tree; }

    /// Returns a reference to m_tree
    gsHDomain<d> &       tree()       { return m_tree; }

    /// Cleans the basis, removing any inactive levels
    void makeCompressed();


    /// Returns the boundary basis for side s
    // GISMO_UPTR_FUNCTION_DEC(gsHTensorBasis<d,T>, boundaryBasis, boxSide const &)

    /// Returns a bounding box for the basis' domain
    gsMatrix<T> support() const;

    gsMatrix<T> support(const index_t & i) const;

    void elementSupport_into(const index_t i, gsMatrix<index_t, d, 2>& result) const
    {
        index_t lvl = levelOf(i);
        m_bases[lvl]->elementSupport_into(m_xmatrix[lvl][ i - m_xmatrix_offset[lvl] ],
                                          result);
    }

    gsMatrix<T> elementInSupportOf(index_t j) const
    {
        index_t lvl = levelOf(j);
        gsMatrix<index_t,d,2> sup;
        m_bases[lvl]->elementSupport_into(m_xmatrix[lvl][j-m_xmatrix_offset[lvl]], sup);
        std::pair<point,point> box =  m_tree.queryLevelCell(sup.col(0),sup.col(1),lvl);
        for ( short_t i = 0; i!=d; ++i) //get intersection
        {
            box.first[i]  = ( sup(i,0) >= box.first[i]  ? sup(i,0) : box.first[i] );
            box.second[i] = ( sup(i,1) <= box.second[i] ? sup(i,1) : box.second[i]);
        }
        sup.col(0) = (box.first+box.second)/2;
        sup.col(1) = sup.col(0).array() + 1;
        return m_bases[lvl]->elementDom(sup);
    }

    GISMO_UPTR_FUNCTION_PURE(gsHTensorBasis, clone)

    /// The number of basis functions in this basis
    index_t size() const;

    /// The number of nodes in the tree representation
    int treeSize() const
    {
        return m_tree.size();
    }

    /// The number of active basis functions at points \a u
    void numActive_into(const gsMatrix<T> & u, gsVector<index_t>& result) const;

    /// The 1-d basis for the i-th parameter component at the highest level
    virtual gsBSplineBasis<T> & component(short_t i)
    {
        return m_bases[ this->maxLevel() ]->component(i);
    }

    /// The 1-d basis for the i-th parameter component at the highest level
    virtual const gsBSplineBasis<T> & component(short_t i) const
    {
        return m_bases[ this->maxLevel() ]->component(i);
    }

    /// Returns the tensor basis member of level i
    tensorBasis & tensorLevel(index_t i) const
    {
        needLevel( i );
        return *this->m_bases[i];
    }

    // Refine the basis uniformly by inserting \a numKnots new knots on each knot span
    virtual void uniformRefine(int numKnots = 1, int mul=1);

    // Refine the basis uniformly and adjust the given matrix of coefficients accordingly
    //virtual void uniformRefine_withCoefs(gsMatrix<T>& coefs, int numKnots = 1, int mul = 1)

    // Refine the basis uniformly and produce a sparse matrix which
    // maps coarse coefficient vectors to refined ones
    //virtual void uniformRefine_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, int numKnots = 1, int mul = 1)

    // Refine the basis uniformly and adjust the given matrix of coefficients accordingly
    void uniformRefine_withCoefs(gsMatrix<T>& coefs, int numKnots = 1, int mul = 1);

    // Refine the basis and adjust the given matrix of coefficients accordingly
    void refine_withCoefs(gsMatrix<T> & coefs, gsMatrix<T> const & boxes);

    /** Refine the basis and adjust the given matrix of coefficients accordingly.
     * @param coefs is a matrix of coefficients as given, e.g., by gsTHBSpline<>::coefs();
     * @param boxes specify where to refine; each 5-tuple gives the level of the box,
     * then two indices (in the current level indexing) of the lower left corner and finally
     * two indices of the upper right corner, see gsHTensorBasis::refineElements() for details.
     */
    void refineElements_withCoefs   (gsMatrix<T> & coefs,std::vector<index_t> const & boxes);
    void refineElements_withTransfer(std::vector<index_t> const & boxes, gsSparseMatrix<T> &transfer);

    void refineElements_withCoefs2(gsMatrix<T> & coefs,std::vector<index_t> const & boxes);

    // see gsBasis for documentation
    void matchWith(const boundaryInterface & bi, const gsBasis<T> & other,
                   gsMatrix<index_t> & bndThis, gsMatrix<index_t> & bndOther) const;

    short_t maxDegree() const
    {
        short_t td = m_bases[0]->degree(0);
        // take maximum of coordinate bases degrees
        for (short_t k=1; k!=d; ++k)
            td = math::max(td, m_bases[0]->degree(k));
        return td;
    }

    short_t minDegree() const
    {
        short_t td = m_bases[0]->degree(0);
        // take maximum of coordinate bases degrees
        for (short_t k=1; k!=d; ++k)
            td = math::min(td, m_bases[0]->degree(k));
        return td;
    }

    /// @brief Reduces spline continuity at interior knots by \a i
    void reduceContinuity(int const & i = 1)
    {
        for (unsigned int lvl = 0; lvl <= maxLevel(); lvl++)
        {
            for (unsigned int dir = 0; dir < d; dir++)
            {
                // TODO check: max interior mult + i <= m_p+1

                // We iterate through unique knots, skipping the first and last knot
                // At level 0 we iterate through all unique knots,
                // At level >0 we iterate through all knots that are new, i.e. every other knot starting from 1
                for (gsKnotVector<>::uiterator it = m_bases[lvl]->knots(dir).ubegin() + 1; it < m_bases[lvl]->knots(dir).uend() - 1; it += (lvl == 0? 1 : 2))
                {
                    for(unsigned int j =lvl;j < m_bases.size();j++)
                        m_bases[j]->component(dir).insertKnot(*it,i);

                }
            }
        }
        update_structure();
    }

    /// @brief If the basis is a tensor product of (piecewise)
    /// polynomial bases, then this function returns the polynomial
    /// degree of the \a i-th component.
    virtual inline short_t degree(short_t i) const
    { return m_bases[0]->degree(i);}



    // S.K.
    /// @brief Returns the level(s) at point(s) in the parameter domain.
    ///
    /// \param[in] Pt gsMatrix of size <em>d</em> x <em>n</em>, where\n
    /// \em d is the dimension of the parameter domain and\n
    /// \em n is the number of evaluation points.\n
    /// Each column of \em Pts represents one evaluation point.
    /// \return levels gsMatrix of size <em>1</em> x <em>n</em>.\n
    /// <em>levels(0,i)</em> is the level of the point defined by the <em>i</em>-th column in \em Pts.
    index_t getLevelAtPoint(const  gsMatrix<T> & Pt ) const;

    // S.K.
    /// @brief Returns the level(s) and knot span(s) at point(s) in the parameter domain.
    ///
    /// \param[in] Pt gsMatrix of size <em>d</em> x <em>n</em>, where\n
    /// \em d is the dimension of the parameter domain and\n
    /// \em n is the number of evaluation points.\n
    /// Each column of \em Pts represents one evaluation point.
    /// \param[out] lvl gsVector of length \em n with the levels of the respective points.
    /// \param[out] loIdx  gsMatrix of size <em>d</em> x <em>n</em>.\n
    /// Each column contains
    /// the lower corner of the knot span containing <em>i</em>-th point. The corner is given
    /// in unique knot span indices of level lvl[i].
    void getLevelUniqueSpanAtPoints(const  gsMatrix<T> & Pt,
                                    gsVector<index_t> & lvl,
                                    gsMatrix<index_t> & loIdx ) const;

    /// Returns the level in which the indices are stored internally
    unsigned maxLevel() const
    {
        return m_tree.getMaxInsLevel();
    }

    /// Returns the level of the function indexed \a i (in continued indices)
    inline index_t levelOf(index_t i) const
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

    void degreeElevate(int const & i = 1, int const dir = -1);
    void degreeReduce(int const & i = 1, int const dir = -1);

    void degreeIncrease(int const & i= 1, int const dir = -1);
    void degreeDecrease(int const & i = 1, int const dir = -1);

    /** @brief Refine the basis to levels and in the areas defined by
     * \a boxes with an extension.
     *
     * \param[in] boxes gsMatrix of size \em d x \em n, where\n
     * \em n is the number of refinement boxes.\n
     * Every two consecutive columns specify the lower and upper corner of one refinement box
     * (See also documentation of refine() for the format of \em box)
     * \param[in] refExt is an integer specifying how many cells should also be
     * refined around the respective boxes.
     *
     */
    virtual void refine(gsMatrix<T> const & boxes, int refExt);

    std::vector<index_t> asElements(gsMatrix<T> const & boxes, int refExt = 0) const;

    // std::vector<index_t> asElements(gsMatrix<T> const & boxes, int refExt = 0) const;

    /** @brief Refine the basis to levels and in the areas defined by \a boxes.
     *
     * \param[in] boxes gsMatrix of size \em d x \em n, where\n
     * \em n is the number of refinement boxes.\n
     * Every two consecutive columns specify the lower and upper corner of one refinement box
     * (See also documentation of refine() for the format of \em box)
     */
    virtual void refine(gsMatrix<T> const & boxes);

    /**
     * @brief Insert the given boxes into the quadtree.
     *
     * Each box is defined by <em>2d+1</em> indices, where \em d is the dimension
     * of the parameter domain.
     * The first index defines the level in which the box should be inserted,
     * the next \a d indices the "coordinates" of the lower corner in the index space,
     * and the last \a d indices the "coordinates" of the upper corner.
     *
     * <b>Example:</b> Let <em>d=3</em> and
     *\f[ \mathsf{boxes} = [ L^1, \ell_x^1, \ell_y^1, \ell_z^1, u_x^1, u_y^1, u_z^1,
     L^2, \ell_x^2, \ell_y^2, \ell_z^2, u_x^2, u_y^2, u_z^2,
     L^3, \ell_x^3, \ell_y^3,
     \ldots ], \f]
     * then, the first box will be inserted in level \f$L^1\f$ and its
     * lower and upper corner will have the indices
     * \f$ (\ell_x^1, \ell_y^1, \ell_z^1)\f$ and \f$ (u_x^1, u_y^1, u_z^1) \f$
     * in the index space of level \f$L^1\f$, respectively.
     *
     *
     * @param boxes vector of size <em>N (2d+1)</em>, where\n
     * \em N is the number of boxes,\n
     * \em d is the dimension of the parameter domain.\n
     * See description above for details on the format.
     */
    virtual void refineElements(std::vector<index_t> const & boxes);

    /// Refines all the cells on the side \a side up to level \a lvl
    void refineSide(const boxSide side, index_t lvl);

    /// Refines the basis function with (hierarchical) index \a i
    void refineBasisFunction(const index_t i);

    // Look at gsBasis.h for the documentation of this function
    //virtual void uniformRefine(int numKnots = 1);

    typename gsBasis<T>::domainIter makeDomainIterator() const
    {
        return typename gsBasis<T>::domainIter(new gsHDomainIterator<T, d>(*this));
    }

    typename gsBasis<T>::domainIter makeDomainIterator(const boxSide & s) const
    {
        return ( s == boundary::none ?
                 typename gsBasis<T>::domainIter(new gsHDomainIterator<T, d>(*this)) :
                 typename gsBasis<T>::domainIter(new gsHDomainBoundaryIterator<T, d>(*this,s) )
            );
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
    index_t flatTensorIndexOf(const index_t i) const
    {

        const index_t level = this->levelOf(i);

        const index_t offset = this->m_xmatrix_offset[level];
        const index_t ind_in_level = this->m_xmatrix[level][i - offset];

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
    index_t flatTensorIndexOf(const index_t i, const index_t level) const
    {

        const index_t offset = this->m_xmatrix_offset[level];
        const index_t ind_in_level = this->m_xmatrix[level][i - offset];

        return ind_in_level;
    }

    /// @brief Gives polylines on the boundaries between different levels of the mesh.
    /// @param result variable where to write the polylines in the form
    /// < levels < polylines_in_one_level < one_polyline < one_segment (x1, y1, x2, y2) > > > > ,
    /// where <x1, y1, x2, y2 > are so that (x1, y1) <=LEX  (x2, y2)
    /// and where x1, y1, x2 and y2 are parameters (knots).
    /// @return bounding boxes of the polylines in the form
    /// < levels < polylines_in_one_level < x_ll, y_ll, x_ur, y_ur > > >, where "ur" stands for "upper right" and "ll" for "lower left".
    std::vector< std::vector< std::vector<index_t > > > domainBoundariesParams( std::vector< std::vector< std::vector< std::vector< T > > > >& result) const;

    /// @brief Gives polylines on the boundaries between different levels of the mesh.
    /// @param result variable where to write the polylines in the form
    /// < levels < polylines_in_one_level < one_polyline < one_segment (x1, y1, x2, y2) > > > >
    /// where <x1, y1, x2, y2 > are so that (x1, y1) <=LEX  (x2, y2)
    /// and where x1, y1, x2 and y2 are indices of the knots with respect to m_maxInsLevel.
    /// @return bounding boxes of the polylines in the form
    /// < levels < polylines_in_one_level < x_ll, y_ll, x_ur, y_ur > > >, where "ur" stands for "upper right" and "ll" for "lower left".
    std::vector< std::vector< std::vector<index_t > > > domainBoundariesIndices( std::vector< std::vector< std::vector< std::vector<index_t > > > >& result) const;
    // TO DO: use gsHDomainLeafIterator for a better implementation
    size_t numElements() const
    {
        gsHDomainIterator<T, d> domIter(*this);

        size_t numEl = 0;
        for (; domIter.good(); domIter.next())
        {
            numEl++;
        }

        return numEl;
    }
    using gsBasis<T>::numElements; //unhide

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
    int flatTensorIndexToHierachicalIndex(index_t index,const int level) const;

    /// Fills the vector actives with booleans, that determine if a function
    /// of the given level is active. The functions on the boundary are ordered
    /// in ascending patchindex order.
    ///
    /// \param[in] level : level of the boundary functions
    /// \param[in] s : boundary side
    /// \param[out] actives : the result, true if its active, false if not
    void activeBoundaryFunctionsOfLevel(const unsigned level,const boxSide & s,std::vector<bool>& actives) const;

    /// @brief Increases the multiplicity of a knot with the value \a knotValue in level \a lvl
    /// in direction \a dir by \a mult.
    /// If knotValue is not currently in the given knot vector its not added.
    ///
    /// \param[in] lvl : level
    /// \param[in] dir : direction
    /// \param[in] knotValue : value of the knot
    /// \param[in] mult : multiplicity

    virtual void increaseMultiplicity(index_t lvl, int dir, T knotValue, int mult = 1);

    /// @brief Increases the multiplicity of several knots with the value \a knotValue in level \a lvl
    /// in direction \a dir by \a mult.
    /// If knotValue is not currently in the given knot vector its not added.
    ///
    /// \param[in] lvl : level
    /// \param[in] dir : direction
    /// \param[in] knotValue : value of the knot
    /// \param[in] mult : multiplicity
    virtual void increaseMultiplicity(index_t lvl, int dir, const std::vector<T> & knotValue, int mult = 1);

protected:

    /// @brief Updates the basis structure (eg. charact. matrices, etc), to
    /// be called after any modifications.
    virtual void update_structure(); // to do: rename as updateCharMatrices

    /// @brief Makes sure that there are \a numLevels grids computed
    /// in the hierarachy
    void needLevel(int maxLevel) const;

    /// @brief Creates \a numLevels extra grids in the hierarchy
    void createMoreLevels(int numLevels) const;

    /// gets all the boxes along a slice in direction \a dir at parameter \a par.
    /// the boxes are given back in a std::vector<index_t> and are in the right format
    /// to be given to refineElements().
    void getBoxesAlongSlice( int dir, T par,std::vector<index_t>& boxes ) const;

private:

    /// \brief Inserts a domain into the basis
    void insert_box(point const & k1, point const & k2, int lvl);

    void initialize_class(gsBasis<T> const&  tbasis);

    /// \brief Returns the basis functions of \a level which have support on \a
    /// box, represented as an index box
    void functionOverlap(const point & boxLow, const point & boxUpp,
                         const int level, point & actLow, point & actUpp);

    // \brief Sets all functions of \a level to active or passive- one by one
    void set_activ1(int level);

    // \brief Computes the set of active basis functions in the basis
    void setActive();

    // \brief Computes the connectivity for a level on a mesh that has all vertices
    void addConnectivity(int level, gsMesh<T> & mesh) const;

    ///returns a transfer matrix using the characteristic matrix of the old and new basis
    virtual gsSparseMatrix<T> coarsening(const std::vector<CMatrix>& old,
                                         const std::vector<CMatrix>& n,
                                         const gsSparseMatrix<T,RowMajor> & transfer) const = 0;

    virtual gsSparseMatrix<T> coarsening_direct(const std::vector<gsSortedVector<index_t> >& old,
                                                const std::vector<gsSortedVector<index_t> >& n,
                                                const std::vector<gsSparseMatrix<T,RowMajor> >& transfer) const = 0;

    virtual gsSparseMatrix<T> coarsening_direct2(const std::vector<gsSortedVector<index_t> >& old,
                                                 const std::vector<gsSortedVector<index_t> >& n,
                                                 const std::vector<gsSparseMatrix<T,RowMajor> >& transfer) const = 0;

    /// \brief Implementation of the features common to domainBoundariesParams and domainBoundariesIndices. It takes both
    /// @param indices and @param params but fills in only one depending on @param indicesFlag (if true, then it returns indices).
    std::vector< std::vector< std::vector<index_t > > > domainBoundariesGeneric(std::vector< std::vector< std::vector< std::vector<index_t > > > >& indices,
                                                                                      std::vector< std::vector< std::vector< std::vector< T > > > >& params,
                                                                                      bool indicesFlag ) const;

    //D
public:
    /// \brief Returns transfer matrix betweend the hirarchycal spline given
    /// by the characteristic matrix "old" and this
    void transfer (const std::vector<gsSortedVector<index_t> > &old, gsSparseMatrix<T>& result);

    void transfer2 (const std::vector<gsSortedVector<index_t> > &old, gsSparseMatrix<T>& result);

    /// \brief Creates characteristic matrices for basis where "level" is the
    /// maximum level i.e. ignoring higher level refinements
    void setActiveToLvl(int level, std::vector<CMatrix>& x_matrix_lvl) const;


//    void local2globalIndex( gsVector<index_t,d> const & index,
//                            index_t lvl,
//                            gsVector<index_t,d> & result
//        ) const;

//    void global2localIndex( gsVector<index_t,d> const & index,
//                            index_t lvl,
//                            gsVector<index_t,d> & result
//        ) const;

}; // class gsHTensorBasis

// Next line disallows instantization of gsTensorBasis<0,T>
template<typename T> class gsHTensorBasis<0,T>
{using T::GISMO_ERROR_gsHTensorBasis_cannot_have_dimension_zero;};


} // namespace gismo

// ************************************************
// ************************************************

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsHTensorBasis.hpp)
#endif
