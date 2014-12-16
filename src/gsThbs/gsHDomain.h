/** @file gsHDomain.h

    @brief Provides declaration of the HDomain class.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): G. Kiss, A. Mantzaflaris
*/

# pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsThbs/gsHDomainLeafIter.h>
#include <gsThbs/gsAAPolyline.h>

#include <list>

namespace gismo {

template<typename T, unsigned d> class gsHDomainIterator;

template<class T> class gsSegment;

/**
\brief
Class with a <em>hierarchical domain structure</em> represented by a box
k-d-tree


The hierarchical domain structure represets a sequence of domains \f$\Omega^\ell  \subset \Omega\f$ which are nested in the sense that
\f[
\Omega = \Omega^0 \supset \Omega^1 \supset \Omega^2 \supset \ldots \supset \Omega^N \supset \Omega^{N+1} = \emptyset \f]
Each subdomain \f$\Omega^\ell\f$ is a (not necessarily connected) collection of axis-aligned elements/cells.

\remark In the context of HB-splines and THB-splines, these elements/cells are the knot spans of the tensor-product mesh at the respective level.

The information on the hierarchical domains is stored in a k-d-tree, where each leaf represents an axis-aligned box
\f$\omega \subseteq \Omega\f$, such that
\f$\omega \subseteq \Omega^\ell \land \omega \cap \Omega^{\ell+1} = \emptyset\f$ (i.e., each leaf of the tree can be assiciated with exactly one level of the hierarchy).

The implementation is, up to some technical differences, based on the technique described in the following publication

- G. Kiss, C. Giannelli, and B. Juettler.
Algorithms and data structures for truncated hierarchical B–splines. In M. Floater et al., editors, Mathematical Methods for Curves and Surfaces, volume 8177, pages 304–323.
Lecture Notes in Computer Science, 2014.

also available as a technical report

- G. Kiss, C. Giannelli, and B. Juettler.
Algorithms and Data Structures for
Truncated Hierarchical B–splines
DK-Report No. 2012-14, Doctoral Program Computational Mathematics: Numerical Anaylsis and Symbolic Computation, 2012.

Regarding the mentioned technical differences: A binary tree is used instead of a quad-tree (which was discussed in the above-mentioned publications). Also, the domains are not necessarily split at their middle, but according to the position of the domain of the next level.

Template parameters
\param d is the dimension
\param T is the box-index type
*/

template<unsigned d, class T = unsigned>
class gsHDomain
{
public:
    typedef kdnode<d,T> node;

    typedef typename node::point point; 

    typedef typename node::kdBox box; 

    typedef gsHDomainLeafIter<node,false> literator;

    typedef gsHDomainLeafIter<node,true> const_literator;

    template <class U, unsigned dd> friend class gsHDomainIterator;

private:
    struct query1_visitor;
    struct query2_visitor;
    struct query3_visitor;
    struct query4_visitor;

    struct numLeaves_visitor;
    struct numNodes_visitor;
    struct levelUp_visitor;

private:

    /// Pointer to the root node of the tree
    node * m_root;

    /// Keeps the highest upper indices
    point m_upperIndex;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /// The level of the box representation (global indices)
    unsigned m_indexLevel;

    /// Maximum level present in the tree
    unsigned m_maxInsLevel;

public:

    gsHDomain() : m_indexLevel(0)
    {
        m_root=NULL; 
        m_maxInsLevel = 0;
    }
    
    gsHDomain(point const & upp, unsigned index_level = 10 )
    { 
        init( upp, index_level); 
    }

/*
    /// Constructor using a given node 
    gsHDomain(const node & root_node, unsigned index_level_ , unsigned max_ins_level_ )
    {
        m_root          = new node(root_node) ;
        m_root->parent  = NULL                ;
        m_index_level = index_level_        ;
        max_ins_level = max_ins_level_      ;
    }
*/

    /// Copy constructor (makes a deep copy)
    gsHDomain( const gsHDomain & o) : 
        m_upperIndex(o.m_upperIndex),
        m_indexLevel(o.m_indexLevel),
        m_maxInsLevel(o.m_maxInsLevel)
    {
        m_root = new node(*o.m_root);
    }

    /// Assignment operator (makes a deep copy)
    gsHDomain& operator=( const gsHDomain & o)
    {
        if ( this == &o )
            return *this;
        
        m_root = new node(*o.m_root);

        m_upperIndex  = o.m_upperIndex;
        m_indexLevel = o.m_indexLevel;
        m_maxInsLevel = o.m_maxInsLevel;

        return *this;
    }

    /// Initialize the tree
    void init(point const & upp, unsigned index_level = 10) 
    { 
        m_indexLevel = index_level;
        m_maxInsLevel = 0;

        if ( m_root )
            delete m_root;

        for (unsigned i=0; i<d; ++i) 
            m_upperIndex[i] = (upp[i]<< m_indexLevel);

        m_root = new node(m_upperIndex);
    }

    /// Destructor deletes the whole tree
    ~gsHDomain() { delete m_root; }

    /// Clones the object
    gsHDomain * clone() const;

public:

    void computeFinestIndex( gsVector<unsigned,d> const & index,
                             unsigned lvl,
                             gsVector<unsigned,d> & result
        ) const;

    void computeLevelIndex( gsVector<unsigned,d> const & index,
                            unsigned lvl,
                            gsVector<unsigned,d> & result
        ) const;

    void local2globalIndex( gsVector<unsigned,d> const & index,
                            unsigned lvl,
                            gsVector<unsigned,d> & result
        ) const;

    void global2localIndex( gsVector<unsigned,d> const & index,
                            unsigned lvl,
                            gsVector<unsigned,d> & result
        ) const;

    const point & upperCorner() const
    {
        return m_upperIndex;
    }
    

    /* \brief The insert function which insert box
    defined by points \em lower and \em upper to level \em lvl.

    [\em lower, \em upper] are given by unique knot indices of level \em lvl

    \param lower the lower left corner of the box
    \param upper the upper right corner of the box
    \param *_node the current node
    \param lvl the desired level

    \remarks This function should only be called by the other query1().
    It will return an incorrect result, if the box corresponds to
    a leaf which is in a branch of the tree
    that cannot be reached from \em _node.
    */
    void insertBox (point const & lower, point const & upper,
                    node * _node, int lvl);

    /** \brief The insert function which insert box
    defined by points \em lower and \em upper to level \em lvl.

    [\em lower, \em upper] are given by unique knot indices of level \em lvl.

    \param lower the lower left corner of the box given in \em lvl representation
    \param upper the upper right corner of the box given in \em lvl representation
    \param lvl the desired level
    */
    void insertBox (point const & lower, point const & upper, int lvl)
    { insertBox(lower, upper, m_root, lvl); }

    /** \brief Sinks the box defined by points \em lower and \em upper
    to one level higher.

    [\em lower, \em upper] are given by unique knot indices of level \em lvl.

    \param lower the lower left corner of the box
    \param upper the upper right corner of the box
    \param lvl the level in which \a lower and \a upper are defined
    */
    void sinkBox (point const & lower, point const & upper, int lvl);

    /// Returns the internal coordinates of point \a point_idx of level \a lvl
    void internalIndex (point const & point_idx, int lvl, point & internal_idx)
    {
        for ( unsigned i = 0; i!=d; ++i )
            internal_idx[i] = point_idx[i] << (m_indexLevel-lvl) ;
    }


    /*
    \brief Returns true if the box defined by \em lower and \em upper
    is completely contained in level and
    does not overlap with any higher level.

    This query is called by the other query1()
    (the one without \em _node as input parameter) and starts
    checking at node \em _node of the k-d-tree.

     \param lower the lower left corner of the box
     \param upper the upper right corner of the box
     \param level the level to be checked against
     \param _node node of the k-d-tree where the search starts.

     \remarks This function should only be called by the other query1().
    It will return an incorrect result, if the box corresponds to
    a leaf which is in a branch of the tree
    that cannot be reached from \em _node.
    */
    bool query1 (point const & lower, point const & upper,
                 int level, node * _node ) const
        { return boxSearch< query1_visitor >(upper,lower,level,_node); }


    /** \brief Returns true if the box defined by \em lower and \em upper
    * is completely contained in \em level and
    * does not overlap with any higher level.
    *
    * More mathematically, returns
    * \em true, if
    * \f$\omega \subseteq \Omega^{\mathsf{level}} \land
    * \omega \cap
    * \Omega^{\mathsf{level}+1} = \emptyset \f$, where \f$\omega\f$
    * is the box  defined by \em lower and \em upper.
    *
    * \todo Specify input format of \em lower and \em upper.
    *
    * <b>Example:</b> in two dimensons, let\n
    * \f$\Omega^4 = (0.1,0.9)\times(0.2,0.7)\f$,\n
    * \f$\Omega^5 = (0.4,0.8)\times(0.3,0.6) \subset \Omega^4\f$,\n
    * \f$\Omega^6 = \emptyset \subset \Omega^5\f$.
    *
    * Testing the box \f$ (0.1,0.4)\times(0.3,0.5)\f$,
    * contained in \f$\Omega^4 \setminus \Omega^5\f$, returns:\n
    * query1( [0.1,0.4], [0.3,0.5], 4 ) = true\n
    * query1( [0.1,0.4], [0.3,0.5], 5 ) = false
    *
    * Testing the box \f$ (0.3,0.5)\times(0.4,0.5)\f$,
    * partly overlapping \f$\Omega^5\f$, returns:\n
    * query1( [0.3,0.4], [0.5,0.5], 4 ) = false\n
    * query1( [0.3,0.4], [0.5,0.5], 5 ) = false
    *
    * Testing the box \f$ (0.6,0.8)\times(0.4,0.5)
    * \subset \Omega^5\f$ returns:\n
    * query1( [0.6,0.4], [0.8,0.5], 4 ) = false\n
    * query1( [0.6,0.4], [0.8,0.5], 5 ) = true\n
    *
    * \param lower the lower left corner of the box
    * \param upper the upper right corner of the box
    * \param level the level to be checked against
    * \returns True, if
    * \f$\omega \subseteq \Omega^{\mathsf{level}} \land \omega
    * \cap \Omega^{\mathsf{level}+1} = \emptyset \f$,
    * where \f$\omega\f$ is the box defined by \em lower and \em upper.
    */
    bool query1 (point const & lower, point const & upper,
                 int level) const
        { return boxSearch< query1_visitor >(lower,upper,level,m_root); }

    /* \brief Returns true if the box defined by \em lower and \em upper
     is contained in a domain with a higher level than \em level
    \param lower lower left corner of the cube [k1,k2]
    \param upper upper right corner of the cube [k1,k2]
    \param level current level
    \param _node node of the k-d-tree where the search starts.

     \remarks This function should only be called by the other query2().
    It will return an incorrect result, if the box corresponds to
    a leaf which is in a branch of the tree
    that cannot be reached from \em _node.
    */
    bool query2 (point const & lower, point const & upper,
                 int level, node *_node ) const
    { return boxSearch< query2_visitor >(lower,upper,level,_node); }

    /** \brief Returns true if the box defined by \em lower and \em upper
     * is completely contained in a domain with a level different to \em level.
     *
     * \param lower the lower left corner of the box
     * \param upper the upper right corner of the box
     * \param level the level \f$ \ell_0 \f$ to be checked against
     * \returns True, if there exists a level \f$\ell \neq \ell_0 \f$, such that
     * \f$\omega \subseteq \Omega^\ell \land \omega \cap \Omega^{\ell+1} = \emptyset\f$,
     * where \f$\omega\f$ is the box defined by \em lower and \em upper.
     */
    bool query2 (point const & lower, point const & upper,
                 int level) const
    { return boxSearch< query2_visitor >(lower,upper,level,m_root); }

    // query3 is used if both query1 and query2 are false to decide if the
    // coresponding basis function is active or not.  it returns the
    // smallest level in which [k1,k2] is completely contained and not
    // completely overlaped by higher omega structure of the function is
    // the same as in case of query1 and query2 but insted of returning a
    // true or false value we remember the lowest level, which is returned
    // at the end.
    // TO DO "leaf" is misleading -- rename to node_
    int query3(point const & k1, point const & k2, 
               int level, node *_node ) const
    { return boxSearch< query3_visitor >(k1,k2,level,_node); }


    /** \brief Returns the lowest level \f$\ell\f$ s.t.
     * \f$\omega \subseteq \Omega^\ell \land \omega
     * \cap \Omega^{\ell+1} = \emptyset \f$.
     *
     * ...where \f$\omega\f$ is the box defined by \em lower
     * and \em upper.
     *
     * \param lower the lower left corner of the box
     * \param upper the upper right corner of the box
     * \param level specifies which level \em lower and \em upper refer to.
     *
     * \todo Check if this really is the use of \em level !
     */
    int query3(point const & lower, point const & upper,
               int level) const
    { return boxSearch< query3_visitor >(lower,upper,level,m_root); }


    // query4 returns the highest level with which box [k1, k2]
    // overlaps
    int query4(point const & lower, point const & upper,
               int level, node  *_node) const
    { return boxSearch< query4_visitor >(lower,upper,level,_node); }

    /** \brief Returns the highest level with which
     * the box defined by \em lower and \em upper
     * overlaps.
     *
     * \param lower the lower left corner of the box
     * \param upper the upper right corner of the box
     * \param level specifies which level \em lower and \em upper refer to.
     *
     * \todo Check if this is true!
     */
    int query4(point const & lower, point const & upper,
               int level) const
    { return boxSearch< query4_visitor >(lower,upper,level,m_root); }

    // to do: move to the hpp file do avoid need for instantization
    void incrementLevel()
    {
        m_maxInsLevel++;

        GISMO_ASSERT( m_maxInsLevel <= m_indexLevel,
        "Problem with indices, increase number of levels (to do).");

        leafSearch< levelUp_visitor >(); 
    }

    /// Multiply all coordinates by two
    void multiplyByTwo()
    {
        m_upperIndex *= 2;
        nodeSearch< liftCoordsOneLevel_visitor >(); 
    }

    // to do: move to the hpp file do avoid need for instantization
    void decrementLevel()
    {
        m_maxInsLevel--;
        leafSearch< levelDown_visitor >(); 
    }

    literator beginLeafIterator()
    {
        return literator(m_root, m_indexLevel);
    }

    const_literator beginLeafIterator() const
    {
        return const_literator(m_root, m_indexLevel);
    }

    void makeCompressed();
        
    /// Returns the number of nodes in the tree
    int size() const
    { return nodeSearch< numNodes_visitor >(); }

    /// Returns the number of distinct knots in direction \a k of level \a lvl
    int numBreaks(int lvl, int k) const
    {
        return (m_upperIndex[k] >> (m_indexLevel - lvl) );
    }

    /// Returns the number of leaves in the tree
    int leafSize() const 
    { return leafSearch< numLeaves_visitor >(); }

    void printLeaves() const
    { leafSearch< printLeaves_visitor >(); }

    /// Returns a list of boxes defined by left-bottom (b1) and
    /// right-top (b2) corners for the splitting in B-spline patches
    /// together with the corresponding levelOf
    /// \param b1 left bottom corners of boxes
    /// \param b2 right upper corners of boxs
    /// \param level corresponding level
    void getBoxes(gsMatrix<unsigned>& b1, 
                  gsMatrix<unsigned>& b2, 
                  gsVector<unsigned>& level) const;

    /// Returns a list of boxes defined by left-bottom (b1) and
    /// right-top (b2) corners for the splitting in B-spline patches
    /// together with the corresponding levelOf
    /// b1 and b2 are indexing in the corresponding level indices
    /// \param b1 left bottom corners of boxes
    /// \param b2 right upper corners of boxs
    /// \param level corresponding level
    void getBoxesInLevelIndex(gsMatrix<unsigned>& b1,
                  gsMatrix<unsigned>& b2,
                  gsVector<unsigned>& level) const;

    ///return a list of polylines- boundaries of each connected
    ///component for all levels in the parameter space
    std::vector< std::vector< std::vector< std::vector< unsigned int > > > > getPolylines() const;

    std::vector< std::vector< std::vector< unsigned int > > > getPolylinesSingleLevel(std::vector<gsVSegment<T> >& seg) const;


    inline unsigned getIndexLevel() const
    {
        return m_indexLevel;
    }

    inline unsigned getMaxInsLevel() const
    {
        return m_maxInsLevel;
    }



private:
    
    /// Returns true if the boxes defined by pairs [k1,k2] and [k3,k4]
    /// overlap
    /// \param k1 the lower left corner of the first box
    /// \param k2 the upper right corner of the first box
    /// \param k3 the lower left corner of the second box
    /// \param k4 the upper right corner of the second box
    static bool haveOverlap(box const & box1, box const & box2);

    /// Returns true if box1 is contained in box2
    static bool isContained(box const & box1, box const & box2);

    /// Returns two point- upper left and lower right corner of the
    /// overlaping part of [k1,k2] and [k3,k4]
    /// \param k1 the lower left corner of the first box
    /// \param k2 the upper right corner of the first box
    /// \param k3 the lower left corner of the second box
    /// \param k4 the upper right corner of the second box
    static std::pair<point,point> select_part(point const & k1, point const & k2,
                                              point const & k3, point const & k4);
    
    static void bisectBox(box const & original, 
                          int k, T coord,
                          box & leftBox,
                          box & rightBox );

    /// Sets the level of all descendents of the node <em>_node</em> to level
    /// lvl if the level of the node is smaller then \em lvl.
    /// Otherwise the level of node is not changed
    /// \param *_node the node whose descendents are to be set their levels
    /// \param lvl the new level for the descendents
    static void setLevel(node *_node, int lvl);

    static bool isDegenerate(box const & someBox);

    /// Adds \a nlevels new index levels in the tree
    void setIndexLevel(int nlevels) const
    {
        GISMO_NO_IMPLEMENTATION
    }

    /// Represents boxes of the tree in a big vector.
    /// \param[out] boxes each item corresponds to a box and is
    /// represented by a vector of unsigned ints: first the
    /// coordinates of the lower left and then the coordinates of
    /// upper right corner. Everything in terms of their level.
    void getBoxes_vec(std::vector<std::vector<unsigned int> >& boxes) const;

    ///connect the boxes returned from quadtree
    void connect_Boxes(std::vector<std::vector<unsigned int> > &boxes) const;
    void connect_Boxes_2(std::vector<std::vector<unsigned int> > &boxes) const;

    /// For each x-coordinate delete repeated parts of vertical segments.
    /// \a vert_seg_lists list, where for each x coordinate (sorted increasingly) a list of vertical segments
    /// with that coordinate is saved.
    void getRidOfOverlaps( std::list< std::list< gsVSegment<T> > >& vert_seg_lists ) const;

    /// Sweepline algorithm.
    void sweeplineConnectAndMerge( std::vector< std::vector< std::vector<unsigned int> > >& result,
                                   std::list< std::list< gsVSegment<T> > >& vert_seg_lists ) const;
    

private:

    /// Iterates on the leafs of the tree and applies \ visitor.  The
    /// visitor controls the return type and the update of the result
    /// type at every leaf node
    template<typename visitor>
    typename visitor::return_type
    boxSearch(point const & k1, point const & k2, 
              int level, node  *_node) const;

    /// Iterates on the leafs of the tree and applies \ visitor.  The
    /// visitor controls the operation to be performed
    template<typename visitor>
    typename visitor::return_type
    leafSearch() const;

    /// Iterates on all the nodes of the tree and applies \ visitor.
    /// The visitor controls the operation to be performed
    template<typename visitor>
    typename visitor::return_type
    nodeSearch() const;

    // Query 1
    struct query1_visitor
    {
        typedef bool return_type;
        
        // initialize result as true
        static const return_type init = true;
        
        static void visitLeaf(kdnode<d,T> * leafNode , int level, return_type & res)
        {
            // If we hit a leaf, then it overlaps qBox, so take
            // minimum with current result
            // (assumes that no degenerate leaves exist in the tree)
            //GISMO_ASSERT( !isDegenerate(*curNode->box), "Encountered an empty leaf");
            
            if ( (!isDegenerate(*leafNode->box)) && leafNode->level != level )
                // if (leafNode->level != level )
                res = false;
        }
    };
    
    // Query 2
    struct query2_visitor
    {
        typedef bool return_type;
        
        // initialize result as true
        static const return_type init = true;
        
        static void visitLeaf(kdnode<d,T> * leafNode , int level, return_type & res)
        {
            // If we hit a leaf, then it overlaps qBox, so 
            // we checj the level of this leaf
            // (assumes that no degenerate leaves exist in the tree)
            //GISMO_ASSERT( !isDegenerate(*curNode->box), "Encountered an empty leaf");
            
            if ( (!isDegenerate(*leafNode->box)) && leafNode->level <= level )
                // if (leafNode->level <= level )
                res = false;
        }
    };
    
    // Query 3
    struct query3_visitor
    {
        typedef int return_type;
        
        // initialize result as a max possible value, since we are looking
        // for a minimum
        static const return_type init = 1000000;
        
        static void visitLeaf(kdnode<d,T> * leafNode , int level, return_type & res)
        {
            // If we hit a leaf, then it overlaps qBox, so take
            // minimum with current result
            // (assumes that no degenerate leaves exist in the tree)
            //GISMO_ASSERT( !isDegenerate(*curNode->box), "Encountered an empty leaf");
            if ( (!isDegenerate(*leafNode->box)) && leafNode->level < res )
                res = leafNode->level;
        }
    };
    
    // Query 4
    struct query4_visitor
    {
        typedef int return_type;
        
        // initialize result as a minimum possible value, since we are
        // looking for a maximum
        static const return_type init = -1;
        
        static void visitLeaf(kdnode<d,T> * leafNode , int level, return_type & res)
        {
            // If we hit a leaf, then it overlaps qBox, so take
            // minimum with current result
            // (assumes that no degenerate leaves exist in the tree)
            //GISMO_ASSERT( !isDegenerate(*curNode->box), "Encountered an empty leaf");
            if ( (!isDegenerate(*leafNode->box)) && leafNode->level > res )
                res = leafNode->level;
        }
    };
    
    // Increases the level by 1 for all leaves
    struct levelUp_visitor
    {
        typedef int return_type;
        static const return_type init = 0;
        
        static void visitLeaf(kdnode<d,T> * leafNode, return_type & i)
        {
            GISMO_UNUSED(i);
            leafNode->level++;
        }
    };

    // Decreases the level by 1 for all leaves
    struct levelDown_visitor
    {
        typedef int return_type;
        static const return_type init = 0;
        
        static void visitLeaf(kdnode<d,T> * leafNode, return_type & i)
        {
            GISMO_UNUSED(i);
            leafNode->level--;
        }
    };

    /// Counts number of nodes in the tree
    struct numLeaves_visitor
    {
        typedef int return_type;
        static const return_type init = 0;
        
        static void visitLeaf(kdnode<d,T> * leafNode, return_type & i)
        {
            i++;
        }
    };

    /// Counts number of nodes in the tree
    struct numNodes_visitor
    {
        typedef int return_type;
        static const return_type init = 0;
        
        static void visitNode(kdnode<d,T> * _node, return_type & i)
        {
            i++;
        }
    };

    /// Multiplies everything by 2
    struct liftCoordsOneLevel_visitor
    {
        typedef int return_type;
        static const return_type init = 0;
        
        static void visitNode(kdnode<d,T> * leafNode, return_type & i)
        {
            leafNode->multiplyByTwo();
        }
    };

    /// Counts number of nodes in the tree
    struct printLeaves_visitor
    {
        typedef int return_type;
        static const return_type init = 0;
        
        static void visitLeaf(kdnode<d,T> * leafNode, return_type & i)
        {
            GISMO_UNUSED(i);
            gsInfo << *leafNode;
        }
    };

    /*
    struct toGlobalIndex 
    {
        toGlobalIndex(T lvl, T ilevel) : m_pow(ilevel-lvl) {}
        const T operator()(const T & x) const { return x << m_pow; }        
        T m_pow;
    };

    struct toLocalIndex 
    {
        toLocalIndex(T lvl, T ilevel) : m_pow(ilevel-lvl) {}
        const T operator()(const T & x) const { return x >> m_pow; }        
        T m_pow;
    };
    */

};


}// end namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsHDomain.hpp)
#endif

