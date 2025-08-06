/** @file gsKdNode2.h

    @brief Provides declaration of the tree node.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Mantzaflaris
*/

# pragma once

#include <gsDomain/gsKdTreeLeafIterator.h>

namespace gismo {

/**
    @brief Struct representing a kd-tree node

    The nodes are of two types:
    - Split nodes
    - Leaf nodes

    Template parameters
    \param d is the dimension
    \param Z is the box-coordinate index type

    \ingroup HSplines
*/
template<short_t d, class Z, class leafData> // TODO: use template template params to avoid duplivation
struct gsKdTree
{
    // Defines the type of the box
    typedef typename leafData::point point_t;
    typedef leafData        data_t;
    typedef gsKdTree       node_t;

    // Tree iterators
    typedef gsKdTreeLeafIter<node_t,false> literator;
    typedef gsKdTreeLeafIter<node_t,true > const_literator;
    //typedef gsHDomainLeafIter<node,false> literator;
    ///typedef gsHDomainLeafIter<node,true> const_literator;

    /// axis in which the children of this node split the domain
    /// special value -1 denotes a leaf node
    int axis;

    /// Split coordinate (meaningfull only for split nodes)
    Z pos;

    /// Pointer to the parent node
    gsKdTree * parent;

    /// Pointer to the left child of this split node (if it is one)
    gsKdTree * left;

    /// Pointer to the right child of this split node (if it is one)
    gsKdTree * right;

    /// The data held in this leaf node (if leaf)
    /// box->first is the lower left corner of the box
    /// box->second is the upper right corner of the box
    leafData * data;

    /// Constructor (empty node)
    gsKdTree() : axis(-2), parent(0), left(0), right(0), data(nullptr)
    { }

    /// Constructor (root node)
    gsKdTree(const leafData & leaf) : axis(-1), parent(0), left(0) , right(0),
                                       data(new leafData(leaf))
    { }

    /// Recursively copies the whole subtree under \a o, and sets it's
    /// parent to \a parentNode
    gsKdTree(const gsKdTree & o, gsKdTree * parentNode = NULL) : axis(o.axis)
    {
        parent = parentNode;
        if ( axis == -1 )
        {
            GISMO_ASSERT( (o.left == 0) && (o.right == 0),
                          "Problem: leaf with children." );
            data  = new leafData(*o.data);
            left = right = nullptr;
        }
        else
        {
            GISMO_ASSERT( o.data == nullptr,
                          "Problem: split node with box." );
            pos   = o.pos;
            left  = new gsKdTree(*o.left , this);
            right = new gsKdTree(*o.right, this);
            data  = nullptr;
        }
    }

    /// Recursively deletes the whole subtree under this node
    ~gsKdTree()
    {
        // TODO: non-recursive
        if ( isLeaf() )
        {
            delete data;
        }
        else
        {
            delete left;
            delete right;
        }
    }

public: // Member functions related to the tree starting at \a this node

    /// Iterates on all the nodes of the tree and applies \ visitor.
    /// The visitor controls the operation to be performed
//    template<typename visitor>
//    typename visitor::return_type
//    nodeSearch() const;

    gsKdTree * pointSearch(const point_t & p) const;

    literator beginLeafIterator()
    { return literator(this); }

    const_literator beginLeafIterator() const
    { return const_literator(this); }

    void makeCompressed();

    /// Returns the number of nodes in the tree
    int size() const;

    /// Returns the number of leaves in the tree
    int leafSize() const;

    /// Prints out the leaves of the kd-tree
    void printLeaves() const;

private: // Structs related to tree operations

    /// Counts number of nodes in the tree
    struct numLeaves_visitor
    {
        typedef int return_type;
        static return_type init() {return 0;}

        static void visitLeaf(node_t * , return_type & i)
        {
            i++;
        }
    };

    /// Counts number of nodes in the tree
    struct numNodes_visitor
    {
        typedef int return_type;
        static return_type init() {return 0;}

        static void visitNode(const node_t * , return_type & i)
        {
            i++;
        }
    };

    /// Counts number of nodes in the tree
    struct printLeaves_visitor
    {
        typedef int return_type;
        static return_type init() {return 0;}

        static void visitLeaf(node_t * leafNode, return_type &)
        {
            gsInfo << *leafNode;
        }
    };


public: // Member functions related to \a this node

    // Data Accessors
    leafData & nodeData()
    {
        GISMO_ASSERT(data, "Asked for lowCorner at node without box data.");
        return *data;
    }

    const leafData & nodeData() const { return const_cast<gsKdTree*>(this)->nodeData(); }

    bool isLeaf() const { return axis == -1; }

    bool isRoot() const { return parent == NULL; }

    bool isTerminal() const
    { return (axis!=-1) && (left->axis==-1) && (right->axis==-1); }

    bool isLeftChild()  const { return parent!=NULL && this==parent->left; }

    bool isRightChild() const { return parent!=NULL && this==parent->right; }

    gsKdTree * sibling() const
    {
        GISMO_ASSERT( parent != 0, "Root does not have a sibling.");
        return (parent->left == this ? parent->right : parent->left );
    }

    void multiplyByTwo()
    {
        if ( isLeaf() )
        {
            data->multiplyByTwo();
        }
        else
        {
            pos *= 2;
        }
    }

    void divideByTwo()
    {
        if ( isLeaf() )
        {
            data->divideByTwo();
        }
        else
        {
            pos /= 2;
        }
    }

    /// Splits the node (i.e., two children are added)
    inline void split()
    {
        GISMO_ASSERT( (left == 0) && (right == 0),
                      "Can only split leaf nodes.");
        GISMO_ASSERT( axis > -1, "Split axis not prescribed.");

        // Make new left and right children
        left          = new gsKdTree;
        right         = new gsKdTree;
        // Set axis to -1 (since they are leaves)
        left ->axis   =
        right->axis   = -1;
        // Set parent to this node
        left ->parent =
        right->parent = this;
        // Set box and level
        left ->data   = data;
        right->data   = new leafData(*data);
        // Detach box from parent (is now at left child)
        data = NULL;

        leafData::split(*this);
        // substitutes
        //left ->box->second[axis] =
        //right->box->first [axis] = pos;
    }

    /// Merges terminal node (i.e., two children are joined)
    inline void merge()
    {
        GISMO_ASSERT( (left->isLeaf()) && (right->isLeaf()),
                      "Can only merge terminal nodes.");

        // merge the data of right into left
        leafData::mergeToLeft(*this);
        // substitutes:
        //box = left->box;
        //box->second[axis] = right->box->second[axis];
        //level = left->level;

        // Recover box
        axis  = - 1;
        data = left->data;
        left->data = NULL;

        // Delete children
        delete  left;
        left  = NULL;
        delete right;
        right = NULL;
    }


    /// Splits the node (i.e., two children are added)
    void split(int splitAxis, Z splitPos)
    {
        //GISMO_ASSERT( box->second[splitAxis] != splitPos, "Degenerate split " << box->second[splitAxis] <<" != "<<splitPos);
        //GISMO_ASSERT( box->first [splitAxis] != splitPos, "Degenerate split " << box->first[splitAxis]  <<" != "<<splitPos);
        axis = splitAxis;
        pos  = splitPos;
        split();
    }

    /// Splits the node in the middle (ie. two children are added)
    // TODO: remove
    void nextMidSplit()
    {
        leafData::nextMidSplit(*this);
        //axis = ( parent == 0 ? 0 : (parent->axis+1)%d );
        //pos  = box->first [axis] +
        //    (box->second[axis] - box->first[axis])/2 ;
        split(); // Can be degenerate
    }

    /// Splits the node in the middle (ie. two children are added)
    /// If non-degenerate split is impossible, then this is a no-op
    void anyMidSplit(unsigned h)
    {
        int doSplit;
        leafData::adaptivedSplit(h,*this);
        if (doSplit)
            split();
    }

    /// Splits the node adaptively (i.e., two children are added)
    /// according to \a insBox.  If non-degenerate split is impossible,
    /// then this is a no-op.
    /// Splitting is done on a coordinate of the current \a level (aligned)
    /// returns the child that intersects \a insBox or NULL (if no split)
    gsKdTree * adaptiveAlignedSplit(const leafData & insData, unsigned h)
    {
        int doSplit;
        leafData::adaptiveAlignedSplit(insData, h, *this, doSplit);
        if (doSplit)
        {
            split();
            return (doSplit==1 ? right : left );
        }
        return NULL;
    }

    /// Splits the node adaptively (i.e., two children are added)
    /// according to \a insBox. If non-degenerate split is impossible,
    /// then this is a no-op
    // TODO: remove
    gsKdTree * adaptiveSplit(const leafData & insData)
    {
        // assumption: insBox intersects box
        int doSplit;
        leafData::adaptiveAlignedSplit(insData, *this, doSplit);
        if (doSplit)
        {
            split();
            return (doSplit==1 ? right : left );
        }
        return NULL;
    }

    friend std::ostream & operator<<(std::ostream & os, const gsKdTree & n)
    {
        if ( n.isLeaf() )
            os << "Leaf node "<< *n.data <<"\n";
        else
            os << "Split node, axis= "<< n.axis <<", pos="<< n.pos <<"\n";
        return os;
    }

    template<typename visitor>
    typename visitor::return_type
    rangeSearch(point_t const & k1, point_t const & k2, //change k1,k2 to a gsTreeData ??
                int level, size_t maxPath = 16) const
    {
        GISMO_ASSERT( !(k1.array() >= k2.array()).any(), // !isDegenerate(k1,k2)
                      "rangeSearch: Wrong order of points defining the box (or empty box): "
                      << k1.transpose() <<", "<< k2.transpose() <<".\n" );

        typename visitor::return_type res = visitor::init();

        std::vector<const node_t*> stack;
        stack.reserve( 2 * maxPath );
        stack.push_back(this);

        const node_t * curNode;
        while ( ! stack.empty() )
        {
            curNode = stack.back(); //top();
            stack.pop_back();       //pop();

            if ( curNode->isLeaf() )
            {
                // Visit the leaf
                GISMO_ASSERT( curNode->nodeData().check(), "Encountered an invalid leaf");
                visitor::visitLeaf(curNode->nodeData(), level, res );
            }
            else // this is a split-node
            {
                if ( k2[curNode->axis] <= curNode->pos)
                    // qBox overlaps only left child of this split-node
                    stack.push_back(curNode->left); //push(curNode->left);
                else if  ( k1[curNode->axis] >= curNode->pos)
                    // qBox overlaps only right child of this split-node
                    stack.push_back(curNode->right); //push(curNode->right);
                else
                {
                    // qBox overlaps both children of this split-node
                    stack.push_back(curNode->left ); //push(curNode->left );
                    stack.push_back(curNode->right); //push(curNode->right);
                }
            }
        }

        return res;
    }

    const node_t * pointSearch(const point_t & p, int level, size_t maxPath = 16) const
    {
        std::vector<const node_t*> stack;
        stack.reserve( 2 * maxPath );
        stack.push_back(this);

        const node_t * curNode;
        while ( ! stack.empty() )
        {
            curNode = stack.back(); //top();
            stack.pop_back();       //pop();

            if ( curNode->isLeaf() )
            {
                // Point found at current node
                return curNode;
            }
            else // this is a split-node
            {
                if ( p[curNode->axis] < curNode->pos)
                    stack.push_back(curNode->left);
                else
                    stack.push_back(curNode->right);
            }
        }
        GISMO_ERROR("pointSearch: Error ("<< p.transpose()<<").\n" );
    }

/// Iterates on the leafs of the tree and applies \ visitor.  The
/// visitor controls the operation to be performed
template<typename visitor>
typename visitor::return_type
leafSearch()
{
    typename visitor::return_type i = visitor::init();

    node_t * curNode = this;

    while(true)
    {
        if ( !curNode->isLeaf() )
        {   //property: tree has no singles (only childs)
            curNode = curNode->left;
        }
        else
        {
            // Visit the leaf
            visitor::template visitLeaf<d,Z>(curNode->nodeData(), i);

            while (curNode->parent != NULL &&
                   curNode != curNode->parent->left)
                curNode = curNode->parent;

            if ( curNode->isRoot() )
                break;
            else
                curNode = curNode->parent->right;
        }
    }
    return i;
}

template<typename visitor>
typename visitor::return_type
nodeSearch() const
{
    typename visitor::return_type i = visitor::init();

    const node_t * curNode = this;

    while(true)
    {
        visitor::visitNode(curNode, i);

        if ( !curNode->isLeaf() )
        {   //property: tree has no singles
            curNode = curNode->left;
        }
        else
        {
            while (curNode->parent != NULL &&
                   curNode != curNode->parent->left)
                curNode = curNode->parent;

            if ( curNode->isRoot() )
                break;
            else
                curNode = curNode->parent->right;
        }
    }
    return i;
}

};

template<short_t d, class Z, class leafData>
inline int gsKdTree<d, Z, leafData>::size() const
{
    return nodeSearch< numNodes_visitor >();
}

template<short_t d, class Z, class leafData>
inline int gsKdTree<d, Z, leafData>::leafSize() const
{
    return leafSearch< numLeaves_visitor >();
}

template<short_t d, class Z, class leafData>
inline void gsKdTree<d, Z, leafData>::printLeaves() const
{
    leafSearch< printLeaves_visitor >();
}

template<short_t d, class Z, class leafData>
inline void gsKdTree<d, Z, leafData>::makeCompressed()
{
    std::stack<node_t*, std::vector<node_t*> > tstack;
    node_t * curNode;

    // First step: gather all terminal nodes
    std::stack<node_t*, std::vector<node_t*> > stack;
    stack.push(this);
    while ( ! stack.empty() )
    {
        curNode = stack.top();
        stack.pop();

        if ( curNode->isTerminal() )
        {
            // Remember this terminal node
            tstack.push(curNode);
        }
        else if ( ! curNode->isLeaf() ) // this is a non-terminal split-node
        {
                stack.push(curNode->left );
                stack.push(curNode->right);
        }
    }

    // Second step: reccursively merge siblings that have the same level
    while ( ! tstack.empty() )
    {
        curNode = tstack.top();
        tstack.pop();

        if (curNode->left->nodeData().level() == curNode->right->nodeData().level())
        {
            // Merge left and right
            curNode->merge();
            if ( !curNode->isRoot() &&
                  curNode->parent->isTerminal() )
                tstack.push(curNode->parent );
        }
    }

    // // Store the max path length
    // m_maxPath = minMaxPath().second;
}


}// namespace gismo
