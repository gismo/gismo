/** @file gsKDTree.h

    @brief Provides declaration of kd-tree interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): M. Möller, A. Mantzaflaris, inspired by Junjie Dong

    The code of Junjie Dong can be found at
    https://github.com/junjiedong/KDTree
*/

#pragma once

//#include <stdexcept>
//#include <cmath>
//#include <vector>
//#include <unordered_map>
//#include <utility>
//#include <algorithm>

#include <gsUtils/gsBoundedPriorityQueue.h>

namespace gismo
{

template<typename Z = index_t, typename Data_t = int>
class gsKDTree
{
public:
    index_t axis;  /// split axis
    Z pos;         /// split position
    Data_t * data;  /// pointer to node data

    gsKDTree * parent, * left, * right; /// parent and children

    //typedef gsVector<Z, d> point;

    /// Constructor (empty node)
    gsKDTree() : axis(-2), data(0), parent(0), left(0), right(0) { }

    gsKDTree(index_t _axis, Z _position, Data_t * _data = nullptr)
    : axis(_axis), pos(_position), data(_data), parent(nullptr), left(nullptr), right(nullptr) {}

    /// Recursively copies the whole subtree under \a o, and sets it's
    /// parent to \a parentNode
    gsKDTree(const gsKDTree & o, gsKDTree * parentNode = NULL) : axis(o.axis)
    {
        this->operator=(o);
        parent = parentNode;
    }

    /// Recursively deletes the whole subtree under this node
    ~gsKDTree()
    {
        // TODO: non-recursive
        if ( hasData() ) delete data;

        if ( !isLeaf() ) 
        {
            delete left;
            delete right;
        }
    }


    gsKDTree & operator=(const gsKDTree& o)
    {
        parent = nullptr;
        if (this != &o)
        {
            if ( isLeaf() )
            {
                GISMO_ASSERT( (o.left == 0) && (o.right == 0), "Problem: leaf with children." );
                left = right = NULL;
            }
            else
            {
                GISMO_ASSERT(o.data == 0, "Problem: split node with data." );
                pos   = o.pos;
                left  = new gsKDTree(*o.left , this);
                right = new gsKDTree(*o.right, this);
            }
            
            if (o.hasData())
                data  = new Data_t(*o.data);
            else
                data   = nullptr;            
        }

        return *this;
    }

    
public:

    bool hasData() const { return (nullptr != data); }
    //Data_t & data() { return *data; }

    bool isLeaf() const { return axis == -1; }
    bool isRoot() const { return parent == NULL; }
    bool isTerminal() const { return (axis!=-1) && (left->axis==-1) && (right->axis==-1); }

    bool isLeftChild()  const { return parent!=NULL && this==parent->left; }
    bool isRightChild() const { return parent!=NULL && this==parent->right; }

    gsKDTree * sibling() const
    { 
        GISMO_ASSERT( parent != 0, "Root does not have a sibling.");
        return (parent->left == this ? parent->right : parent->left ); 
    }

public:

    /// Splits the node (i.e., two children are added)
    inline void split()
    {
        GISMO_ASSERT( (left == 0) && (right == 0),
                      "Can only split leaf nodes.");
        GISMO_ASSERT( axis > -1, "Split axis not prescribed.");

        // Make new left and right children
        left          = new gsKDTree();
        right         = new gsKDTree();
        // Set axis to -1 (since they are leaves)
        left ->axis   =
        right->axis   = -1;
        // Set parent to this node
        left ->parent = 
        right->parent = this;
        // Set data to both children
        left ->data    = data;    
        right->data    = new Data_t(*data);
        // Detach data from parent (is now at left child)
        data = nullptr;
        //left->data->setLeft (axis,pos);
        //left->data->setRight(axis,pos);
    }

    /// Splits the node (i.e., two children are added)
    void split(int splitAxis, Z splitPos)
    {
        axis = splitAxis;
        pos  = splitPos;
        split();
    }

    /// Merges terminal node (i.e., two children are joined)
    inline void merge()
    {
        GISMO_ASSERT( (left->isLeaf()) && (right->isLeaf()),
                      "Can only merge terminal nodes.");

        // Recover left data
        data = left->data;
        left->data = NULL;
        //data->setLeft(axis, right->data->getLeftPos(axis) );
        axis  = - 1;

        // Delete children
        delete  left;
        left  = NULL;
        delete right;
        right = NULL;
    }

    inline int numLeaves() const
    { return leafSearch< counter_visitor >(); }

    inline int printLeaves() const
    { return leafSearch< printLeaves_visitor >(); }

    inline int numNodes() const
    { return nodeSearch< counter_visitor >(); }

private:

    /// Counts visited nodes
    struct counter_visitor
    {
        typedef int return_type;
        static return_type init() {return 0;}

        static void visitNode(const gsKDTree * , return_type & i)
        {
            i++;
        }
    };
    
    /// Prints the nodes in the tree
    struct printLeaves_visitor
    {
        typedef int return_type;
        static return_type init() {return 0;}
        
        static void visitLeaf(const gsKDTree * leafNode, return_type &)
        {
            gsInfo << *leafNode;
        }
    };


public:

    template<typename visitor>
    typename visitor::return_type leafSearch() const
    {
        typename visitor::return_type i = visitor::init();

        const gsKDTree * curNode = this;

        while(true)
        {
            if ( !curNode->isLeaf() )
            {   //property: tree has no singles (only childs)
                curNode = curNode->left;
            }
            else
            {
                // Visit the leaf
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
    typename visitor::return_type nodeSearch() const
    {
        typename visitor::return_type i = visitor::init();
        
        const gsKDTree * curNode = this;

        while(true)
        {
            visitor::visitNode(curNode, i);
            //gsInfo << "curnode "<< curNode <<"\n";
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

    template<typename visitor, typename point>
    typename visitor::return_type rangeSearch(point const & k1, point const & k2) const
    {

        typename visitor::return_type res = visitor::init();
        
        std::vector<const gsKDTree*> stack;
        size_t m_maxPath = 20;
        stack.reserve( 2 * m_maxPath );
        stack.push_back(this);

        const gsKDTree * curNode = this;
        while ( ! stack.empty() )
        {
            curNode = stack.back(); //top();
            stack.pop_back();       //pop();
            
            if ( curNode->isLeaf() )
            {
                // Visit the leaf
                visitor::visitLeaf(curNode, res );
                // TODO
                // if (visitor::visitLeaf(curNode, res ) ) return res;
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

    template<typename point>
    gsKDTree * pointSearch(const point & pp) const
    {
        //TODO:
        // return this->rangeSearch(pp,pp);
        std::vector<const gsKDTree*> stack;
        int m_maxPath = 20;
        stack.reserve( 2 * m_maxPath );
        stack.push_back(this);

        const gsKDTree * curNode;
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
                if ( pp[curNode->axis] < curNode->pos)
                    stack.push_back(curNode->left);
                else
                    stack.push_back(curNode->right);
            }
        }

        GISMO_ERROR("pointSearch: Error ("<< pp.transpose()<<").\n" );
    }

    std::pair<int,int> minMaxPath() const
    {
        const gsKDTree * curNode = this;
        int min = 1000000000, max = -1, cur = 0;

        while(true)
        {
            if ( !curNode->isLeaf() )
            {   //property: tree has no singles
                curNode = curNode->left;
                cur++;
            }
            else
            {
                // Update min-max
                min = math::min(min,cur);
                max = math::max(max,cur);

                while (curNode->parent != NULL &&
                       curNode != curNode->parent->left)
                {
                    curNode = curNode->parent;
                    cur--;
                }

                if ( curNode->isRoot() )
                    break;
                else
                    curNode = curNode->parent->right;
            }
        }
        return std::make_pair(min,max);
    }

    friend std::ostream & operator<<(std::ostream & os, const gsKDTree & n)
    {
        if ( n.isLeaf() ) 
            os << "Leaf node. ";
        else
            os << (n.isTerminal() ? "Terminal" : "Split" ) << " node, axis= "
                   << n.axis <<", pos="<< n.pos <<". ";
        if ( n.hasData() ) os << "Node has data."; //<< *n.data;
        os <<" \n";
        return os;
    }

};

/**
   \brief A kd-tree in some number of dimensions \a d
*/
template<std::size_t d, typename Z = index_t, typename Data_t = int>
class gsKDTreeExample
{
    typedef gsKDTree<Z,Data_t> Node;
    Node * m_root; /// Root node of the gsKDTree

public:

    // Constructs an empty gsKDTree.
    gsKDTreeExample();
  
    // Frees up all the dynamically allocated resources
    ~gsKDTreeExample();

    // Frees up all the dynamically allocated resources
    void clear();
  
    // Deep-copies the contents of another gsKDTree into this one.
    gsKDTreeExample(const gsKDTreeExample& other);

    // Deep-copies the contents of another gsKDTree into this one.
    gsKDTreeExample& operator=(const gsKDTreeExample& other);
    
    // Returns the dimension of the data stored in this gsKDTree.
    static std::size_t dimension() { return d; }

    int numLeaves() const { return m_root->numLeaves(); };
    int numNodes () const { return m_root->numNodes(); };

    // Returns true if this gsKDTree is empty and false otherwise.
    bool empty() const { return 0;}

    /// Prints the object as a string.
    void print(std::ostream &os) const;

    /// Print (as string) operator
    friend std::ostream &operator<<(std::ostream &os, const gsKDTreeExample &obj)
    {
        obj.print(os);
        return os;
    }
  
}; // class gsKDTreeExample
  
} //namespace gismo
