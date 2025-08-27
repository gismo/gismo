#pragma once

#include <gsDomain/gsKdTree.h>
#include <gsHSplines/gsHTreeData.h>


namespace gismo { //todo: move to gsForwardDeclarations
template<short_t d, class Z = index_t> class gsHTree;
}

namespace
{
    // Query 1
    struct query1_visitor
    {
        typedef bool return_type;

        // initialize result as true
        static return_type init() {return true;}

        template<short_t d, class Z>
        static void visitLeaf(const gismo::gsKdTree<d,Z,gismo::gsHTreeData<d,Z> > * leafNode, int level, return_type & res)
        {
            if ( leafNode->nodeData().level() != level )
                res = false;
        }
    };

    // Query 2
    struct query2_visitor
    {
        typedef bool return_type;

        // initialize result as true
        static return_type init() {return true;}

        template<short_t d, class Z>
        static void visitLeaf(const gismo::gsKdTree<d,Z,gismo::gsHTreeData<d,Z> > * leafNode, int level, return_type & res)
        {
            if ( leafNode->nodeData().level() <= level )
                res = false;
        }
    };

    // Query 3
    struct query3_visitor
    {
        typedef int return_type;

        // initialize result as a max possible value, since we are looking
        // for a minimum
        static return_type init() {return 1000000;}

        template<short_t d, class Z>
        static void visitLeaf(const gismo::gsKdTree<d,Z,gismo::gsHTreeData<d,Z> > * leafNode, int , return_type & res)
        {
            if ( leafNode->nodeData().level() < res )
                res = leafNode->nodeData().level();
        }
    };

    // Query 4
    struct query4_visitor
    {
        typedef int return_type;

        // initialize result as a minimum possible value, since we are
        // looking for a maximum
        static return_type init() {return -1;}

        template<short_t d, class Z>
        static void visitLeaf(const gismo::gsKdTree<d,Z,gismo::gsHTreeData<d,Z> > * leafNode , int , return_type & res)
        {
            if ( leafNode->nodeData().level() > res )
                res = leafNode->nodeData().level();
        }
    };

    /// Decreases the level by 1 for all leaves
    struct levelDown_visitor
    {
        typedef int return_type;
        static return_type init() {return 0;}

        template<short_t d, class Z>
        static void visitLeaf(gismo::gsKdTree<d,Z,gismo::gsHTreeData<d,Z> > * leafNode, return_type &)
        {
            leafNode->nodeData().level()--;
        }
    };

    struct maxLevel_visitor
    {
        typedef int return_type;
        static return_type init() {return 0;}

        template<short_t d, class Z>
        static void visitLeaf(gismo::gsKdTree<d,Z,gismo::gsHTreeData<d,Z> > * leafNode, return_type &i)
        {
            if (leafNode->nodeData().level()>i) i=leafNode->nodeData().level();
        }
    };

} //namespace

namespace gismo {

// This is the root node of a gsHTree. It contains children nodes of type gsKdTree
template<short_t d, class Z>
class gsHTree : public gsKdTree<d,Z,gsHTreeData<d,Z> >
{
public:
    typedef gsKdTree<d,Z,gsHTreeData<d,Z> > node;
    typedef typename node::data_t data_t;

    typedef typename node::point_t point; // it's a gsVector<index_t,d>

    typedef typename data_t::kdBox box; // it's a gsAABB<d,unsigned>

    typedef typename node::literator        literator;
    typedef typename node::const_literator const_literator;

private:

    /// Keeps the highest upper indices (at level gsHTree::m_indexLevel)
    point m_upperIndex;
// #   define Eigen gsEigen
//     EIGEN_MAKE_ALIGNED_OPERATOR_NEW
// #   undef Eigen

    /// The level of the box representation (global indices)
    unsigned m_indexLevel;

    /// Maximum level present in the tree
    unsigned m_maxInsLevel;

    unsigned m_maxPath;

public:
    gsHTree() : m_indexLevel(0)
    {
        m_maxInsLevel = 0;
        m_maxPath = 0;
    }

    gsHTree(point const & upp)
    {
        m_maxInsLevel = 0;
        m_maxPath = 0;
        init(upp);
    }

    /// Copy constructor (makes a deep copy)
    gsHTree(const gsHTree & o) : node(o),
        m_upperIndex(o.m_upperIndex),
        m_indexLevel(o.m_indexLevel),
        m_maxInsLevel(o.m_maxInsLevel),
        m_maxPath(o.m_maxPath)
    { }

    /// Assignment operator (makes a deep copy)
    gsHTree& operator=( const gsHTree & o)
    {
        if ( this != &o )
        {
            node::operator=(o);
            m_upperIndex  = o.m_upperIndex;
            m_indexLevel  = o.m_indexLevel;
            m_maxInsLevel = o.m_maxInsLevel;
            m_maxPath    = o.m_maxPath;
        }
        return *this;
    }

/*
#if EIGEN_HAS_RVALUE_REFERENCES
    gsHTree(gsHTree&& o) :
    m_upperIndex(std::move(o.m_upperIndex)),
    m_indexLevel(o.m_indexLevel),
    m_maxInsLevel(o.m_maxInsLevel),
    m_maxPath(o.m_maxPath)
    {

    }

    gsHTree & operator=(gsHTree&& o)
    {
        m_upperIndex  = std::move(o.m_upperIndex);
        m_indexLevel  = o.m_indexLevel;
        m_maxInsLevel = o.m_maxInsLevel;
        m_maxPath     = o.m_maxPath;
        return *this;
    }
#endif
*/

    /// Initialize the tree
    void init(point const & upp, unsigned index_level)
    {
        m_indexLevel = index_level;
        m_maxInsLevel = 0;

        for (short_t i=0; i<d; ++i)
            m_upperIndex[i] = (upp[i]<< m_indexLevel);

        m_maxPath = 1;

        // To be improved... bad programming
        delete this->data;
        this->data = new gsHTreeData<d,Z>(m_upperIndex);
        this->axis = -1;
        delete this->parent;
        delete this->left;
        delete this->right;
        this->parent = this->left = this->right = nullptr;
    }

	/// Initialize the tree with computing the index_level.
    void init(point const & upp)
    {
        // Idea: find index_level so that for each i, upp[i] * (2 ^ index_level) does not overflow in Z.
        // See issue #672 for more details.

        // backwards compatibility
        Z oldMax = 13;

        Z numMax = std::numeric_limits<Z>::max();

        std::vector<Z> logUpps(d);
        for(short_t i=0; i<d; i++)
        {
            // prevent division by zero
            if(upp[i] == 1)
                logUpps[i] = oldMax;
            else
            {
                // floor of log_2 (numMax / upp[i]):
                logUpps[i] = math::floor( (math::log(numMax) - math::log(upp[i])) / math::log(2) );
            }
        }

        // If the computed number would be too big we take 13 as we used to.
        init(upp, std::min( *std::min_element(logUpps.begin(), logUpps.end()), oldMax) );
    }

    /// Destructor deletes the whole tree
    ~gsHTree() { }

    /// Clones the object
    gsHTree * clone() const;

public:

    /// Returns the number of distinct knots in direction \a k of level \a lvl
    int numBreaks(int lvl, int k) const
    {
        return m_upperIndex[k] >> (m_indexLevel - lvl);
    }

    inline unsigned getIndexLevel() const
    {
        return m_indexLevel;
    }

    inline unsigned getMaxInsLevel() const
    {
        return m_maxInsLevel;
    }


    void computeMaxInsLevel()
    {
        m_maxInsLevel = this->template leafSearch< maxLevel_visitor >();
    }

    inline void decrementLevel()
    {
        m_maxInsLevel--;
        this->template leafSearch< levelDown_visitor >();
    }

    /** \brief Returns the boxes which make up the hierarchical domain
    * and the respective levels.
    *
    * Returns a list of boxes defined by left-bottom (b1) and
    * right-top (b2) corners for the splitting in B-spline patches
    * together with the corresponding levelOf()
    *
    * The boxes \em b1 and \em b2 are given as matrices
    * of size <em>n</em> x <em>d</em>, where \em d is the dimension
    * of the domain, and where \em n is the number of boxes of
    * the gsHTree.\n
    *
    * The numbers in \em b1 and \em b2 are given as
    * unique knot indices of gsHTree::m_maxInsLevel
    *
    * \param[out] b1 <em>n</em> x <em>d</em>-matrix, left bottom corners of boxes
    * \param[out] b2 <em>n</em> x <em>d</em>-matrix, right upper corners of boxes
    * \param[out] level vector of length \em n, corresponding levels
    */
    // getBoxes-functions might get removed at some point of time.
    // Use iterators instead whenever possible.
    void getBoxes(gsMatrix<Z>& b1,
                  gsMatrix<Z>& b2,
                  gsVector<Z>& level) const;

    /// Returns a list of boxes defined by left-bottom (b1) and
    /// right-top (b2) corners for the splitting in B-spline patches
    /// together with the corresponding levelOf
    /// b1 and b2 are indexing in the corresponding level indices
    /// \param b1 left bottom corners of boxes
    /// \param b2 right upper corners of boxs
    /// \param level corresponding levels
    // getBoxes-functions might get removed at some point of time.
    // Use iterators instead whenever possible.
    void getBoxesInLevelIndex(gsMatrix<Z>& b1,
                              gsMatrix<Z>& b2,
                              gsVector<Z>& level) const;

    /// \brief connect the boxes returned from quadtree getBoxes_vec()
    ///
    /// If two neighbouring boxes could be represented by a single
    /// box (i.e., if the union of two axis-aligned boxes is again
    /// an axis-aligned box), then these two are merged into a single box.
    ///
    /// \param[in,out] boxes Format as of getBoxes_vec(), i.e., each box
    /// is represented as vector of size <em>2*d + 1</em> containing
    /// [ [lower corner],[upper corner], Level ], where the corners
    /// are defined by the knot indices on level gsHTree::m_maxInsLevel.
    void connect_Boxes(std::vector<std::vector<Z> > &boxes) const;

    /// Returns the level of the point \a p
    int levelOf(point const & p, int level) const
    { return this->knotSearch(p,level)->nodeData().level();}

    /// Return the upper corner of the tree in level 0
    const point upperCornerIndex() const
    {
        point ind = m_upperIndex;
        for (short_t i=0; i<d; ++i)
            ind[i] = (ind[i] >> m_indexLevel);
        return ind;
    }

    void insertBox(point const & k1, point const & k2, int lvl)
    {
    GISMO_ENSURE( lvl <= static_cast<int>(m_indexLevel), "Max index level reached..");

    // Make a box
    data_t iData(k1,k2,lvl);
    if( iData.aabb().isDegenerate() )
        return;

    // Represent box in the index level
    // iBox.first .unaryExpr(toGlobalIndex(lvl, m_index_level) );
    // iBox.second.unaryExpr(toGlobalIndex(lvl, m_index_level) );
    local2globalIndex( iData.lowerCorner(), static_cast<unsigned>(lvl), iData.lowerCorner() );
    local2globalIndex( iData.upperCorner(), static_cast<unsigned>(lvl), iData.upperCorner());

    // Ensure that the box is within the valid limits
    if ( ( iData.lowerCorner().array() >= m_upperIndex.array() ).any() )
    {
        gsWarn<<" Invalid box coordinate "<<  k1.transpose() <<" at level" <<lvl<<".\n";
        return;
    }

    // Initialize stack
    std::vector<node*> stack;
    stack.reserve( 2 * (m_maxPath + d) );
    stack.push_back(this);  //push(_node);

    node * curNode;
    while ( ! stack.empty() )
    {
        curNode = stack.back(); //top();
        stack.pop_back();       //pop();

        if ( curNode->isLeaf() ) // reached a leaf
        {
            // Since we reached a leaf, it should overlap with iBox

            // If this leaf is already in level lvl, then we have nothing to do
            if ( lvl <= curNode->nodeData().level() )
                continue;

            // Split the leaf (if possible)
            //node * newLeaf = curNode->adaptiveSplit(iBox);
            const unsigned h = 1 << (m_indexLevel - curNode->nodeData().level()) ;
            node * newLeaf = curNode->adaptiveAlignedSplit(iData, h);

            // If curNode is still a leaf, its domain is almost
            // contained in iBox
            if ( !newLeaf ) //  curNode->isLeaf()
            {
                // Increase level and recurse
                if ( ++curNode->nodeData().level() != lvl)
                    stack.push_back(curNode);
            }
            else // treat new child
            {
                stack.push_back(newLeaf);
            }
        }
        else // roll down the tree
        {
            if ( iData.upperCorner()[curNode->axis] <= curNode->pos)
                // iBox overlaps only left child of this split-node
                stack.push_back(curNode->left);
            else if  ( iData.lowerCorner()[curNode->axis] >= curNode->pos)
                // iBox overlaps only right child of this split-node
                stack.push_back(curNode->right);
            else
            {
                // iBox overlaps both children of this split-node
                stack.push_back(curNode->left );
                stack.push_back(curNode->right);
            }
        }
    }

    // Update maximum inserted level
    if ( static_cast<unsigned>(lvl) > m_maxInsLevel)
        m_maxInsLevel = lvl;
}

void sinkBox (point const & k1,
              point const & k2, int lvl)
{
    GISMO_ENSURE( m_maxInsLevel+1 <= m_indexLevel,
                  "Max index level might be reached..");

    // Make a box
    data_t iData(k1,k2,lvl);
    if( iData.aabb().isDegenerate() )
        return;

    // Represent box in the index level
    local2globalIndex( iData.lowerCorner() , static_cast<unsigned>(lvl), iData.lowerCorner() );
    local2globalIndex( iData.upperCorner(), static_cast<unsigned>(lvl), iData.upperCorner());

    // Ensure that the box is within the valid limits
    if ( ( iData.lowerCorner().array() >= m_upperIndex.array() ).any() )
    {
        //gsWarn<<" Invalid box coordinate "<<  k1.transpose() <<" at level" <<lvl<<".\n";
        return;
    }

    // Initialize stack
    std::stack<node*, std::vector<node*> > stack;
    stack.push(this);

    node * curNode;
    while ( ! stack.empty() )
    {
        curNode = stack.top();
        stack.pop();

        if ( curNode->isLeaf() ) // reached a leaf
        {
            // Since we reached a leaf, it should overlap with iBox.
            // Split the leaf (if possible)
            const unsigned h = 1 << (m_indexLevel - curNode->nodeData().level()) ;
            node * newLeaf = curNode->adaptiveAlignedSplit(iData, h);

            // If curNode is still a leaf, its domain is almost
            // contained in iBox
            if ( !newLeaf ) //  implies curNode was a leaf
            {
                // Increase level
                if ( ++curNode->nodeData().level() > static_cast<int>(m_maxInsLevel) )
                    m_maxInsLevel = curNode->nodeData().level();
            }
            else // treat new child
            {
                stack.push(newLeaf);
            }
        }
        else // walk down the tree
        {
            if ( iData.upperCorner()[curNode->axis] <= curNode->pos)
                // iBox overlaps only left child of this split-node
                stack.push(curNode->left);
            else if  ( iData.lowerCorner()[curNode->axis] >= curNode->pos)
                // iBox overlaps only right child of this split-node
                stack.push(curNode->right);
            else
            {
                // iBox overlaps both children of this split-node
                stack.push(curNode->left );
                stack.push(curNode->right);
            }
        }
    }
}

void clearBox ( point const & k1, point const & k2,
                int lvl) // CONSTRAINT: lvl is "minimum level"
{
    GISMO_ENSURE( lvl <= static_cast<int>(m_indexLevel), "Max index level reached..");

    // Make a box
    data_t iData(k1,k2,lvl);
    if( iData.aabb().isDegenerate() )
        return;

    // Represent box in the index level
    // iBox.first .unaryExpr(toGlobalIndex(lvl, m_index_level) );
    // iData.upperCorner().unaryExpr(toGlobalIndex(lvl, m_index_level) );
    local2globalIndex( iData.lowerCorner() , static_cast<unsigned>(lvl), iData.lowerCorner() );
    local2globalIndex( iData.upperCorner(), static_cast<unsigned>(lvl), iData.upperCorner());

    // Ensure that the box is within the valid limits
    if ( ( iData.lowerCorner().array() >= m_upperIndex.array() ).any() )
    {
        gsWarn<<" Invalid box coordinate "<<  k1.transpose() <<" at level" <<lvl<<".\n";
        return;
    }

    // Initialize stack
    std::vector<node*> stack;
    stack.reserve( 2 * (m_maxPath + d) );
    stack.push_back(this);  //push(m_root);

    node * curNode;
    while ( ! stack.empty() )
    {
        curNode = stack.back(); //top();
        stack.pop_back();       //pop();

        if ( curNode->isLeaf() ) // reached a leaf
        {
            // Since we reached a leaf, it should overlap with iBox

            // If this leaf is already in level lvl, then we have nothing to do
            if ( curNode->nodeData().level() <= lvl )
                continue;

            // Split the leaf (if possible)
            const unsigned h = 1 << (m_indexLevel - curNode->nodeData().level()) ;
            node * newLeaf = curNode->adaptiveAlignedSplit(iData, h);
            // If curNode is still a leaf, its domain is almost
            // contained in iBox
            if ( !newLeaf ) //  curNode->isLeaf()
            {
                // Decrease level and reccurse
                if ( --curNode->nodeData().level() != lvl)
                    stack.push_back(curNode);
            }
            else // treat new child
            {
                stack.push_back(newLeaf);
            }
        }
        else // roll down the tree
        {
            if ( iData.upperCorner()[curNode->axis] <= curNode->pos)
                // iBox overlaps only left child of this split-node
                stack.push_back(curNode->left);
            else if  ( iData.lowerCorner()[curNode->axis] >= curNode->pos)
                // iBox overlaps only right child of this split-node
                stack.push_back(curNode->right);
            else
            {
                // iBox overlaps both children of this split-node
                stack.push_back(curNode->left );
                stack.push_back(curNode->right);
            }
        }
    }

    m_maxInsLevel = this->template leafSearch< maxLevel_visitor >(); // compute again globally
}



    /// Iterates on the leafs of the tree and applies \ visitor.  The
    /// visitor controls the return type and the update of the result
    /// type at every leaf node
    template<typename visitor>
    typename visitor::return_type
    boxSearch(point k1, point k2, int level) const
    {
        local2globalIndex( k1, static_cast<unsigned>(level), k1 );
        local2globalIndex( k2, static_cast<unsigned>(level), k2 );
        // box aa(k1,k2, level);
        // return this->template rangeSearch< visitor>( aa, m_maxPath);
        return this->template rangeSearch< visitor>(k1,k2,level, m_maxPath);
    }

    const node * knotSearch(const point & p, int level) const
    {
        point pp;
        local2globalIndex(p, static_cast<unsigned>(level), pp);

        GISMO_ASSERT( ( pp.array() <= m_upperIndex.array() ).all(),
                      "pointSearch: Wrong input: "<< p.transpose()<<", level "<<level<<".\n" );
        return node::pointSearch(pp,level,m_maxPath);
    }

    int query1(point const & k1, point const & k2,
               int level) const
    {
        return boxSearch< query1_visitor >(k1,k2,level);
    }

    int query2(point const & k1, point const & k2,
               int level) const
    {
        return boxSearch< query2_visitor >(k1,k2,level);
    }

    int query3(point const & k1, point const & k2,
               int level) const
    {
        return boxSearch< query3_visitor >(k1,k2,level);
    }

    int query4(point const & k1, point const & k2,
               int level) const
    {
        return boxSearch< query4_visitor >(k1,k2,level);
    }

    /// Returns an cell/element box of a requested level
    // TODO: stop traverse as soon as it is found for the first time..
    struct get_cell_visitor
    {
        typedef std::pair<point,point> return_type;

        // initialize result
        static return_type init()
        {
            return std::make_pair(point::Zero(),point::Zero());
            //return return_type();//!does not properly initialize the points
        }

        static void visitLeaf(const gismo::gsKdTree<d,Z,gsHTreeData<d,Z> > * leafNode , int level, return_type & res)
        {
            if ( leafNode->nodeData().level() == level )
            {
                res.first  = leafNode->nodeData().lowerCorner();
                res.second = leafNode->nodeData().upperCorner();
            }
        }
    };

    std::pair<point, point>
    queryLevelCell(point const & lower, point const & upper,
                   int level) const
    {
        std::pair<point,point> tmp = boxSearch< get_cell_visitor >(lower,upper,level);
        global2localIndex(tmp.first,level,tmp.first);
        global2localIndex(tmp.second,level,tmp.second);
        return tmp;
    }

    inline void local2globalIndex(gsVector<Z, d> const & index,
                                  unsigned lvl,
                                  gsVector<Z, d> & result) const
    {
        for(short_t i = 0; i!=d; ++i)
            result[i] = index[i] << (m_indexLevel-lvl) ;
    }

    inline void global2localIndex(gsVector<Z, d> const & index,
                                  unsigned lvl,
                                  gsVector<Z, d> & result) const
    {
        for(short_t i = 0; i!=d; ++i)
            result[i] = index[i] >> (this->m_indexLevel-lvl) ;
    }

    void computeFinestIndex(gsVector<Z, d> const & index,
                            unsigned lvl,
                            gsVector<Z, d> & result) const
    {
        for(short_t i = 0; i!=d; ++i)
            result[i] = index[i] << (m_maxInsLevel-lvl) ;
    }

    void computeLevelIndex(gsVector<Z, d> const & index,
                           unsigned lvl,
                           gsVector<Z, d> & result) const
    {
        for(short_t i = 0; i!=d; ++i)
            result[i] = index[i] >> (m_maxInsLevel-lvl) ;
    }

    void getBoxes_vec(std::vector<std::vector<Z>>& boxes) const;
};

template<short_t d, class Z>
void gsHTree<d, Z>::getBoxes(gsMatrix<Z>& b1,
                             gsMatrix<Z>& b2,
                             gsVector<Z>& level) const
{
    std::vector<std::vector<Z> > boxes;

    // get all boxes in vector-format
    getBoxes_vec(boxes);

    // connect boxes which have the same levels and are
    // are aligned such that their union again is an
    // axis-aligned box.
    connect_Boxes(boxes);

    // write the result into b1, b2, and level
    b1.resize(boxes.size(),d);
    b2.resize(boxes.size(),d);
    level.resize(boxes.size());
    for(size_t i = 0; i < boxes.size(); i++)
    {
        for(short_t j = 0; j < d; j++)
        {
            b1(i,j) = boxes[i][j];
            b2(i,j) = boxes[i][j+d];
        }
        level[i] = boxes[i][2*d];
    }
}

template<short_t d, class Z>
void gsHTree<d, Z>::getBoxesInLevelIndex(gsMatrix<Z>& b1,
                                         gsMatrix<Z>& b2,
                                         gsVector<Z>& level) const
{
    std::vector<std::vector<Z> > boxes;
    getBoxes_vec(boxes);
    GISMO_ASSERT(d==2 || d==3, "Wrong dimension, should be 2 or 3.");
    //is this test really necessary? florian b.
    for(size_t i = 0; i < boxes.size(); i++)
    {
        if ((boxes[i][0]==boxes[i][d+0]) || (boxes[i][1]==boxes[i][1+d]))
        {
            boxes.erase(boxes.begin()+i);
            i--;
        }
        else if((d == 3) && (boxes[i][2]==boxes[i][d+2]))
        {
            boxes.erase(boxes.begin()+i);
            i--;
        }
    }
    gsVector<Z, d> lowerCorner;
    gsVector<Z, d> upperCorner;
    connect_Boxes(boxes);
    b1.resize(boxes.size(), d);
    b2.resize(boxes.size(), d);
    level.resize(boxes.size());
    for(size_t i = 0; i < boxes.size(); i++)
    {
        for(short_t j = 0; j < d; j++)
        {
            lowerCorner[j] = boxes[i][j];
            upperCorner[j] = boxes[i][j+d];
        }
        level[i] = boxes[i][2*d];
        computeLevelIndex(lowerCorner, level[i], lowerCorner);
        computeLevelIndex(upperCorner, level[i], upperCorner);
        b1.row(i) = lowerCorner;
        b2.row(i) = upperCorner;
    }
}

template<short_t d, class Z> void
gsHTree<d, Z>::getBoxes_vec(std::vector<std::vector<Z> >& boxes) const
{
    boxes.clear();

    std::stack<const node*, std::vector<const node*> > stack;
    //stack.reserve( 2 * m_maxPath );
    stack.push(this);
    const node * curNode;
    point lower, upper;
    while ( ! stack.empty() )
    {
        curNode = stack.top();
        stack.pop();

        if ( curNode->isLeaf() )
        {
            // We need to convert the indices to those of m_maxInsLevel
            // to be able to reconstruct the earlier results.
            const point & lowerGlob = curNode->nodeData().lowerCorner();
            const point & upperGlob = curNode->nodeData().upperCorner();
            unsigned int level = this->m_maxInsLevel;

            global2localIndex(lowerGlob,level,lower);
            global2localIndex(upperGlob,level,upper);

            boxes.push_back(std::vector<Z>());
            for(short_t i = 0; i < d; i++)
            {
                boxes.back().push_back(lower[i]);
            }
            for(short_t i = 0; i < d; i++)
            {
                boxes.back().push_back(upper[i]);
            }
            boxes.back().push_back( curNode->nodeData().level() );
        }
        else
        {
            stack.push(curNode->left) ;
            stack.push(curNode->right);
        }
    }
}

template<short_t d, class Z> void
gsHTree<d, Z>::connect_Boxes(std::vector<std::vector<Z>> &boxes) const
{
    bool change = true;
    while(change)
    {
        change =  false;
        size_t s = boxes.size();
        for(size_t i = 0; i < s; i++)
        {
            for(size_t j = i+1; j < s; j++)
            {
                if(boxes[i][2*d]==boxes[j][2*d]) // if( the levels are the same )
                {
                    unsigned nmCoordLo = 0;
                    unsigned nmCoordUp = 0;
                    unsigned nmCountUp = 0;
                    unsigned nmCountLo = 0; //...the "nm" is for non-matching

                    // Compare the lower and upper corners of the boxes
                    // coordinate-wise, and check if there are differences.
                    // If there are differences, count and store the coordinate
                    for(short_t k=0; k < d; k++)
                    {
                        if( boxes[i][k] != boxes[j][k] )
                        {
                            nmCountLo++;
                            nmCoordLo = k;
                        }

                        if( boxes[i][d+k] != boxes[j][d+k] )
                        {
                            nmCountUp++;
                            nmCoordUp = k;
                        }
                    }

                    // The boxes can only be merged if
                    // the lower and upper corners are the same,
                    // except in one coordinate direction.
                    if( nmCountLo == 1
                        && nmCountUp == 1
                        && nmCoordLo == nmCoordUp )
                    {

                        if( boxes[i][nmCoordLo] == boxes[j][d+nmCoordUp] )
                        {
                            // box i is "on top" of box j.
                            // It inherits the lower corner from box j:
                            boxes[i][nmCoordLo] = boxes[j][nmCoordLo];
                            boxes.erase( boxes.begin()+j );
                            s--;
                            j--;
                            change = true;
                        }

                        if( boxes[i][d+nmCoordUp] == boxes[i][nmCoordLo] )
                        {
                            // box i is "below" of box j.
                            // It inherits the upper corner from box j:
                            boxes[i][d+nmCoordUp] = boxes[j][d+nmCoordUp];
                            boxes.erase( boxes.begin()+j );
                            s--;
                            j--;
                            change = true;
                        }
                    }
                } // if boxes same level
            } // for j
        } // for i
    }
}

/*
template<short_t d, class Z>
void gsHTree<d, Z>::getBoxes(gsMatrix<Z>& b1, gsMatrix<Z>& b2, gsVector<Z>& level) const
{
    std::vector<std::vector<Z> > boxes;

    // get all boxes in vector-format
    getBoxes_vec(boxes);

    // connect boxes which have the same levels and are
    // are aligned such that their union again is an
    // axis-aligned box.
    connect_Boxes(boxes);

    // write the result into b1, b2, and level
    b1.resize(boxes.size(),d);
    b2.resize(boxes.size(),d);
    level.resize(boxes.size());
    for(size_t i = 0; i < boxes.size(); i++)
    {
        for(short_t j = 0; j < d; j++)
        {
            b1(i,j) = boxes[i][j];
            b2(i,j) = boxes[i][j+d];
        }
        level[i] = boxes[i][2*d];
    }
}
*/

}//namespace gismo
