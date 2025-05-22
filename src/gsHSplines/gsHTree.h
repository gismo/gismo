

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
        static void visitLeaf(const typename gismo::gsHTree<d, Z>::node * leafNode , int level, return_type & res)
        {
            if ( leafNode->level != level )
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
        static void visitLeaf(const typename gismo::gsHTree<d, Z>::node * leafNode , int level, return_type & res)
        {
            if ( leafNode->level <= level )
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
        static void visitLeaf(const typename gismo::gsHTree<d, Z>::node * leafNode , int , return_type & res)
        {
            if ( leafNode->level < res )
                res = leafNode->level;
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
        static void visitLeaf(const typename gismo::gsHTree<d, Z>::node * leafNode , int , return_type & res)
        {
            if ( leafNode->level > res )
                res = leafNode->level;
        }
    };

    /// Decreases the level by 1 for all leaves
    struct levelDown_visitor
    {
        typedef int return_type;
        static return_type init() {return 0;}

        template<short_t d, class Z>
        static void visitLeaf(typename gismo::gsHTree<d, Z>::node * leafNode, return_type &)
        {
            leafNode->nodeData().level()--;
        }
    };

    struct maxLevel_visitor
    {
        typedef int return_type;
        static return_type init() {return 0;}

        template<short_t d, class Z>
        static void visitLeaf(typename gismo::gsHTree<d, Z>::node * leafNode, return_type &i)
        {
            if (leafNode->level>i) i=leafNode->level;
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
    gsHTree( const gsHTree & o) :
        m_upperIndex(o.m_upperIndex),
        m_indexLevel(o.m_indexLevel),
        m_maxInsLevel(o.m_maxInsLevel),
        m_maxPath(o.m_maxPath)
    {
    }

    /// Assignment operator (makes a deep copy)
    gsHTree& operator=( const gsHTree & o)
    {
        if ( this == &o )
            return *this;

        m_upperIndex  = o.m_upperIndex;
        m_indexLevel  = o.m_indexLevel;
        m_maxInsLevel = o.m_maxInsLevel;
        m_maxPath    = o.m_maxPath;

        return *this;
    }

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

    /// Initialize the tree
    void init(point const & upp, unsigned index_level)
    {
        m_indexLevel = index_level;
        m_maxInsLevel = 0;

        for (short_t i=0; i<d; ++i)
            m_upperIndex[i] = (upp[i]<< m_indexLevel);

        m_maxPath = 1;
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

    inline void decrementLevel()
    {
        m_maxInsLevel--;
        this->template leafSearch< levelDown_visitor >();
    }

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
    box iBox(k1,k2);
    const unsigned h = m_indexLevel - lvl;
    if( isDegenerate(iBox) )
        return;

    // Represent box in the index level
    // iBox.first .unaryExpr(toGlobalIndex(lvl, m_index_level) );
    // iBox.second.unaryExpr(toGlobalIndex(lvl, m_index_level) );
    local2globalIndex( iBox.first , static_cast<unsigned>(lvl), iBox.first );
    local2globalIndex( iBox.second, static_cast<unsigned>(lvl), iBox.second);

    // Ensure that the box is within the valid limits
    if ( ( iBox.first.array() >= m_upperIndex.array() ).any() )
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
            if ( lvl <= curNode->level )
                continue;

            // Split the leaf (if possible)
            //node * newLeaf = curNode->adaptiveSplit(iBox);
            node * newLeaf = curNode->adaptiveAlignedSplit(iBox, h);

            // If curNode is still a leaf, its domain is almost
            // contained in iBox
            if ( !newLeaf ) //  curNode->isLeaf()
            {
                // Increase level and recurse
                if ( ++curNode->level != lvl)
                    stack.push_back(curNode);
            }
            else // treat new child
            {
                stack.push_back(newLeaf);
            }
        }
        else // roll down the tree
        {
            if ( iBox.second[curNode->axis] <= curNode->pos)
                // iBox overlaps only left child of this split-node
                stack.push_back(curNode->left);
            else if  ( iBox.first[curNode->axis] >= curNode->pos)
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
    box iBox(k1,k2);
    if( isDegenerate(iBox) )
        return;

    // Represent box in the index level
    local2globalIndex( iBox.first , static_cast<unsigned>(lvl), iBox.first );
    local2globalIndex( iBox.second, static_cast<unsigned>(lvl), iBox.second);

    // Ensure that the box is within the valid limits
    if ( ( iBox.first.array() >= m_upperIndex.array() ).any() )
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
            node * newLeaf = curNode->adaptiveAlignedSplit(iBox, m_indexLevel);

            // If curNode is still a leaf, its domain is almost
            // contained in iBox
            if ( !newLeaf ) //  implies curNode was a leaf
            {
                // Increase level
                if ( ++curNode->level > static_cast<int>(m_maxInsLevel) )
                    m_maxInsLevel = curNode->level;
            }
            else // treat new child
            {
                stack.push(newLeaf);
            }
        }
        else // walk down the tree
        {
            if ( iBox.second[curNode->axis] <= curNode->pos)
                // iBox overlaps only left child of this split-node
                stack.push(curNode->left);
            else if  ( iBox.first[curNode->axis] >= curNode->pos)
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
    box iBox(k1,k2);
    if( isDegenerate(iBox) )
        return;

    // Represent box in the index level
    // iBox.first .unaryExpr(toGlobalIndex(lvl, m_index_level) );
    // iBox.second.unaryExpr(toGlobalIndex(lvl, m_index_level) );
    local2globalIndex( iBox.first , static_cast<unsigned>(lvl), iBox.first );
    local2globalIndex( iBox.second, static_cast<unsigned>(lvl), iBox.second);

    // Ensure that the box is within the valid limits
    if ( ( iBox.first.array() >= m_upperIndex.array() ).any() )
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
            if ( curNode->level <= lvl )
                continue;

            // Split the leaf (if possible)
            node * newLeaf = curNode->adaptiveAlignedSplit(iBox, m_indexLevel);

            // If curNode is still a leaf, its domain is almost
            // contained in iBox
            if ( !newLeaf ) //  curNode->isLeaf()
            {
                // Decrease level and reccurse
                if ( --curNode->level != lvl)
                    stack.push_back(curNode);
            }
            else // treat new child
            {
                stack.push_back(newLeaf);
            }
        }
        else // roll down the tree
        {
            if ( iBox.second[curNode->axis] <= curNode->pos)
                // iBox overlaps only left child of this split-node
                stack.push_back(curNode->left);
            else if  ( iBox.first[curNode->axis] >= curNode->pos)
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

    node * knotSearch(const point & p, int level) const
    {
        point pp;
        local2globalIndex(p, static_cast<unsigned>(level), pp);
        
        GISMO_ASSERT( ( pp.array() <= m_upperIndex.array() ).all(),
                      "pointSearch: Wrong input: "<< p.transpose()<<", level "<<level<<".\n" );
        return node::pointSeach(pp,level,m_maxPath);
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

        static void visitLeaf(const node * leafNode , int level, return_type & res)
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

    
};

}//namespace gismo
