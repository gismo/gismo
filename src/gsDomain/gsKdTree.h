/** @file gsKdTree.h

    @brief Provides a k-d tree (binary space partition by axis-aligned
    splits) and the base class of its leaf-data policies.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Mantzaflaris
*/

# pragma once

#include <cstddef>
#include <ostream>
#include <stack>
#include <vector>

#include <gsCore/gsForwardDeclarations.h>
#include <gsDomain/gsKdTreeLeafIterator.h>

namespace gismo {

/**
    @brief Base class of the leaf-data policies of gsKdTree.

    A policy is the type \a leafData of gsKdTree<d,Z,leafData>. It
    is stored in the leaves of the tree and defines how a leaf is split
    and merged. A policy derives from this class and implements a
    hook by declaring a member of the same name, which hides the
    default below by ordinary name lookup; \a operator<< is instead
    implemented by a friend taking the policy type, which overload
    resolution prefers over the default. The hooks that are not
    implemented by the policy fall back to the defaults below.

    The hooks used by gsKdTree are:
    - \a point: the type of a point of the tree (default gsVector<Z,d>).
    - copy construction \a leafData(const leafData&): used when a leaf is
      copied or split.
    - \a multiplyByTwo(), \a divideByTwo(): scale the coordinates held in
      the leaf data (no default: error).
    - static \a split(const Node & n): called after the leaf \a n was
      split. The policy must set the bounds of the left and right
      child along \a n.axis, so that they meet at \a n.pos (no default:
      error).
    - static \a mergeToLeft(const Node & n): called before the children of
      the terminal node \a n are merged. The policy must extend the
      data of the left child to cover the right child (no default: error).
    - static \a anyMidSplit(unsigned h, Node & n, int & doSplit): chooses
      axis and pos of the leaf \a n and sets \a doSplit to non-zero iff
      a non-degenerate split exists (no default: error).
    - static \a adaptiveAlignedSplit(const Data & ins, unsigned h,
      Node & n, int & doSplit): as anyMidSplit, according to the data
      \a ins; \a doSplit is 1 (resp. -1) if the part of \a ins that
      matters is in the right (resp. left) child, and 0 if no split is
      done (no default: error).
    - \a check() const: returns true iff the data are valid (default true).
    - static \a canMerge(const Node & n): true iff the children of the
      terminal node \a n may be merged by gsKdTree::makeCompressed(). The
      default compares the \a level() of the data of the two children,
      and so compiles only for policies that provide \a level().
    - \a operator<<(std::ostream&, const leafData&): prints the data
      (default prints a placeholder).

    The base class is empty and non-polymorphic, so it adds
    neither storage nor dispatch cost to the leaves.

    \tparam d the dimension
    \tparam Z the coordinate index type

    \ingroup HSplines
*/
template<short_t d, class Z>
struct gsKdTreeData
{
    /// The type of a point of the tree
    typedef gsVector<Z,d> point;

    /// Scales the coordinates of the leaf data by a factor of 2
    void multiplyByTwo()
    { GISMO_ERROR("gsKdTreeData: multiplyByTwo must be provided by the leaf-data policy"); }

    /// Divides the coordinates of the leaf data by a factor of 2
    void divideByTwo()
    { GISMO_ERROR("gsKdTreeData: divideByTwo must be provided by the leaf-data policy"); }

    /// Updates the data of the children of the freshly split node \a n
    template<class Node>
    static void split(const Node &)
    { GISMO_ERROR("gsKdTreeData: split must be provided by the leaf-data policy"); }

    /// Updates the data of the left child of the terminal node \a n,
    /// before merging
    template<class Node>
    static void mergeToLeft(const Node &)
    { GISMO_ERROR("gsKdTreeData: mergeToLeft must be provided by the leaf-data policy"); }

    /// Chooses axis and position of a middle split of \a n
    template<class Node>
    static void anyMidSplit(unsigned, Node &, int &)
    { GISMO_ERROR("gsKdTreeData: anyMidSplit must be provided by the leaf-data policy"); }

    /// Chooses axis and position of a split of \a n that is aligned
    /// to the data \a ins
    template<class Data, class Node>
    static void adaptiveAlignedSplit(const Data &, unsigned, Node &, int &)
    { GISMO_ERROR("gsKdTreeData: adaptiveAlignedSplit must be provided by the leaf-data policy"); }

    /// Returns true iff the data are valid
    bool check() const { return true; }

    /// Returns true iff the children of the terminal node \a n may be
    /// merged; requires the policy to provide \a level()
    template<class Node>
    static bool canMerge(const Node & n)
    { return n.left->nodeData().level() == n.right->nodeData().level(); }

    friend std::ostream & operator<<(std::ostream & os, const gsKdTreeData &)
    { return os << "leaf data\n"; }

protected:
    gsKdTreeData() = default;
    gsKdTreeData(const gsKdTreeData &) = default;
    gsKdTreeData & operator=(const gsKdTreeData &) = default;
    ~gsKdTreeData() = default;
};

/**
    @brief Struct representing a k-d tree node.

    The nodes are of two types:
    - Split nodes
    - Leaf nodes

    The tree is a pure spatial partition: it stores, for every split
    node, an axis and a position, and for every leaf the data of the
    policy \a leafData (see gsKdTreeData).

    The only structural invariant is the following. Each split node
    partitions its cell along \a axis at the integer coordinate \a pos.
    Within one split node, \a pos and the bounds of the two children's
    cells along \a axis are expressed in the same integer grid. The leaf-data
    policy defines what that grid is. gsKdTree itself never interprets
    coordinates. For example, gsHTreeData stores global indices on the
    finest grid, while gsCutTreeData stores per-leaf level-local
    indices. Code outside a policy specialisation must not read raw leaf
    data. A specialisation exposes it through its own typed accessors.

    The value of \a axis denotes the type of the node: -2 is an empty
    node (default-constructed or moved-from), -1 is a leaf, and a value
    >= 0 is the split axis of a split node.

    The members multiplyByTwo() and divideByTwo() scale \a pos on a split
    node and forward to the policy on a leaf. Applied to every node of a
    subtree, they refine or coarsen the grid of that subtree by a factor
    of 2 and preserve the invariant. divideByTwo() is exact only if every
    coordinate is even.

    Template parameters
    \param d is the dimension
    \param Z is the coordinate index type
    \param leafData is the leaf-data policy, derived from gsKdTreeData

    \ingroup HSplines
*/
template<short_t d, class Z, class leafData> // TODO: use template template params to avoid duplication
struct gsKdTree
{
    // Defines the type of the box
    typedef typename leafData::point point_t;
    typedef leafData        data_t;
    typedef gsKdTree       node_t;

    // Tree iterators
    typedef gsKdTreeLeafIter<node_t,false> literator;
    typedef gsKdTreeLeafIter<node_t,true > const_literator;

    /// axis in which the children of this node split the domain
    /// special value -1 denotes a leaf node, -2 an empty node
    int axis;

    /// Split coordinate (meaningful only for split nodes)
    Z pos;

    /// Pointer to the parent node
    gsKdTree * parent;

    /// Pointer to the left child of this split node (if it is one)
    gsKdTree * left;

    /// Pointer to the right child of this split node (if it is one)
    gsKdTree * right;

    /// The leaf data (policy-defined), null on split and empty nodes
    leafData * data;

    /// Constructor (empty node)
    gsKdTree() : axis(-2), pos(0), parent(nullptr), left(nullptr), right(nullptr), data(nullptr)
    { }

    /// Constructor (root node)
    explicit gsKdTree(const leafData & leaf) : axis(-1), pos(0), parent(0), left(0) , right(0),
                                       data(new leafData(leaf))
    { }

    /// Recursively copies the whole subtree under \a o, and sets it's
    /// parent to \a parentNode
    gsKdTree(const gsKdTree & o, gsKdTree * parentNode = NULL)
    {
        _deep_copy(o,parentNode);
    }

    /// Move constructor: takes over the subtree of \a o, which is left
    /// as an empty node. The new node is a root.
    gsKdTree(gsKdTree && o) noexcept
    : axis(o.axis), pos(o.pos), parent(nullptr),
      left(o.left), right(o.right), data(o.data)
    {
        if (left)  left->parent  = this;
        if (right) right->parent = this;
        o.axis = -2;
        o.left = o.right = nullptr;
        o.data = nullptr;
    }

    /// Assignment operator (makes a deep copy). \a o may lie inside the
    /// subtree of this node. The parent of this node is kept.
    gsKdTree& operator=( const gsKdTree & o)
    {
        if ( this != &o )
        {
            // The copy is complete before anything is freed
            gsKdTree tmp(o);
            *this = give(tmp);
        }
        return *this;
    }

    /// Move assignment: replaces the subtree of this node by the one of
    /// \a o, which is left as an empty node. The parent of this node is
    /// kept. \a o may lie inside the subtree of this node, but must not
    /// be an ancestor of this node.
    gsKdTree& operator=( gsKdTree && o) noexcept
    {
        if ( this != &o )
        {
            const int    oldAxis  = axis;
            gsKdTree *   oldLeft  = left;
            gsKdTree *   oldRight = right;
            leafData *   oldData  = data;

            axis  = o.axis;
            pos   = o.pos;
            left  = o.left;
            right = o.right;
            data  = o.data;
            if (left)  left->parent  = this;
            if (right) right->parent = this;

            o.axis = -2;
            o.left = o.right = nullptr;
            o.data = nullptr;

            if ( oldAxis == -1 )
            {
                delete oldData;
            }
            else
            {
                delete oldLeft;
                delete oldRight;
            }
        }
        return *this;
    }

    /// Recursively deletes the whole subtree under this node
    ~gsKdTree()
    {
        // TODO: non-recursive
        _free_payload();
    }

protected:

    /// Deletes the leaf data, or the children, of this node
    void _free_payload()
    {
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

    void _deep_copy(const gsKdTree & o, gsKdTree * parentNode)
    {
        parent = parentNode;
        axis   = o.axis;
        if ( axis == -2 )
        {
            pos  = o.pos;
            left = right = nullptr;
            data = nullptr;
        }
        else if ( axis == -1 )
        {
            GISMO_ASSERT( (o.left == 0) && (o.right == 0),
                          "Problem: leaf with children." );
            pos   = o.pos;
            data  = new leafData(*o.data);
            left  = right = nullptr;
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

public: // Member functions related to the tree starting at \a this node

    /// Iterator over the leaves; an empty node has no leaves
    literator beginLeafIterator()
    { return axis == -2 ? literator() : literator(this); }

    const_literator beginLeafIterator() const
    { return axis == -2 ? const_literator() : const_literator(this); }

    /// Merges, bottom-up, the terminal nodes whose children are
    /// mergeable according to leafData::canMerge
    void makeCompressed();

    /// Returns the number of nodes in the tree
    int size() const;

    /// Returns the number of leaves in the tree
    int leafSize() const;

    /// Prints out the leaves of the kd-tree
    void printLeaves() const;

    void printNodes() const;

private: // Structs related to tree operations

    /// Counts number of leaves in the tree
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

    /// Prints every leaf node to gsInfo
    struct printLeaves_visitor
    {
        typedef int return_type;
        static return_type init() {return 0;}

        static void visitLeaf(const node_t * leafNode, return_type &)
        {
            gsInfo << *leafNode;
        }
    };

    /// Prints every node to gsInfo
    struct printNodes_visitor
    {
        typedef int return_type;
        static return_type init() {return 0;}

        static void visitNode(const node_t * _Node, return_type &)
        {
            gsInfo << *_Node;
        }
    };


public: // Member functions related to \a this node

    // Data Accessors
    leafData & nodeData()
    {
        GISMO_ASSERT(data, "nodeData() called on a node without leaf data.");
        return *data;
    }

    const leafData & nodeData() const { return const_cast<gsKdTree*>(this)->nodeData(); }

    bool isLeaf() const { return axis == -1; }

    bool isRoot() const { return parent == NULL; }

    bool isTerminal() const
    { return (axis>=0) && (left->axis==-1) && (right->axis==-1); }

    bool isLeftChild()  const { return parent!=NULL && this==parent->left; }

    bool isRightChild() const { return parent!=NULL && this==parent->right; }

    gsKdTree * sibling() const
    {
        GISMO_ASSERT( parent != 0, "Root does not have a sibling.");
        return (parent->left == this ? parent->right : parent->left );
    }

    /// Multiplies by 2 the split position (split node) or the leaf data (leaf)
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

    /// Divides by 2 the split position (split node) or the leaf data (leaf)
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

        left          = new gsKdTree;
        right         = new gsKdTree;
        // Set axis to -1 (since they are leaves)
        left ->axis   =
        right->axis   = -1;
        left ->parent =
        right->parent = this;
        // Set box and level
        left ->data   = data;
        right->data   = new leafData(*data);
        // Detach box from parent (is now at left child)
        data = NULL;

        leafData::split(*this);
    }

    /// Merges terminal node (i.e., two children are joined)
    inline void merge()
    {
        GISMO_ASSERT( (left->isLeaf()) && (right->isLeaf()),
                      "Can only merge terminal nodes.");

        // merge the data of right into left
        leafData::mergeToLeft(*this);

        // Recover box
        axis  = - 1;
        data = left->data;
        left->data = NULL;

        delete  left;
        left  = NULL;
        delete right;
        right = NULL;
    }


    /// Splits the node (i.e., two children are added)
    void split(int splitAxis, Z splitPos)
    {
        axis = splitAxis;
        pos  = splitPos;
        split();
    }

    /// Splits the node in the middle (ie. two children are added)
    /// If non-degenerate split is impossible, then this is a no-op
    void anyMidSplit(unsigned h)
    {
        int doSplit;
        leafData::anyMidSplit(h, *this, doSplit);
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

    friend std::ostream & operator<<(std::ostream & os, const gsKdTree & n)
    {
        if ( n.isLeaf() )
        {
            os << "Leaf node "<< *n.data;
            os << "       -> axis= "<< n.axis <<"\n";
        }
        else
            os << "Split node, axis= "<< n.axis <<", pos="<< n.pos <<"\n";
        return os;
    }

    template<typename visitor>
    typename visitor::return_type
    rangeSearch(point_t const & k1, point_t const & k2,
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
            curNode = stack.back();
            stack.pop_back();

            if ( curNode->axis == -2 ) // empty node: no leaves
                continue;
            else if ( curNode->isLeaf() )
            {
                // Visit the leaf
                GISMO_ASSERT( curNode->nodeData().check(), "Encountered an invalid leaf");
                visitor::visitLeaf(curNode, level, res );
            }
            else // this is a split-node
            {
                if ( k2[curNode->axis] <= curNode->pos)
                    // qBox overlaps only left child of this split-node
                    stack.push_back(curNode->left);
                else if  ( k1[curNode->axis] >= curNode->pos)
                    // qBox overlaps only right child of this split-node
                    stack.push_back(curNode->right);
                else
                {
                    // qBox overlaps both children of this split-node
                    stack.push_back(curNode->left );
                    stack.push_back(curNode->right);
                }
            }
        }

        return res;
    }

    /// Returns the leaf containing \a p (a point on the tree's integer grid).
    /// A point on a split plane belongs to the right child. O(depth), no allocation.
    const node_t * pointSearch(const point_t & p) const
    {
        const node_t * curNode = this;
        while ( curNode->axis >= 0 )
            curNode = ( p[curNode->axis] < curNode->pos ) ? curNode->left : curNode->right;

        GISMO_ENSURE( curNode->isLeaf(), "pointSearch: empty tree ("<< p.transpose()<<")." );
        return curNode;
    }

    template<typename visitor>
    typename visitor::return_type
    leafSearch() const
    {
        return const_cast<node_t*>(this)->template leafSearch<visitor>();
    }

    /// Iterates on the leafs of the tree and applies \ visitor.  The
    /// visitor controls the operation to be performed
    template<typename visitor>
    typename visitor::return_type
    leafSearch()
    {
        typename visitor::return_type i = visitor::init();

        if ( axis == -2 ) // empty node: no leaves
            return i;

        node_t * curNode = this;

        while(true)
        {
            if ( curNode->axis >= 0 )
            {   //property: tree has no singles (only childs)
                curNode = curNode->left;
            }
            else
            {
                // Visit the leaf (empty nodes have no leaves)
                if ( curNode->isLeaf() )
                    visitor::visitLeaf(curNode, i);

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
        return const_cast<node_t*>(this)->template nodeSearch<visitor>();
    }

    template<typename visitor>
    typename visitor::return_type
    nodeSearch()
    {
        typename visitor::return_type i = visitor::init();

        if ( axis == -2 ) // empty node: no nodes to visit
            return i;

        node_t * curNode = this;

        while(true)
        {
            if ( curNode->axis != -2 ) // empty nodes are not visited
                visitor::visitNode(curNode, i);

            if ( curNode->axis >= 0 )
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
    const_cast<gsKdTree*>(this)->leafSearch< printLeaves_visitor >();
}

template<short_t d, class Z, class leafData>
inline void gsKdTree<d, Z, leafData>::printNodes() const
{
    this->nodeSearch< printNodes_visitor >();
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

        if ( curNode->axis == -2 ) // empty node: nothing to merge
            continue;
        else if ( curNode->isTerminal() )
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

    // Second step: recursively merge siblings that can be merged
    while ( ! tstack.empty() )
    {
        curNode = tstack.top();
        tstack.pop();

        if (leafData::canMerge(*curNode))
        {
            // Merge left and right
            curNode->merge();
            if ( !curNode->isRoot() &&
                  curNode->parent->isTerminal() )
                tstack.push(curNode->parent );
        }
    }
}


}// namespace gismo
