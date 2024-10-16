/** @file gsKdNode.h

    @brief Provides declaration of the tree node.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Mantzaflaris
*/

# pragma once

#include <gsHSplines/gsAABB.h>

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
template<short_t d, class Z = index_t>
struct gsKdNode
{
    // Defines the type of the box
    typedef          gsAABB<d, Z> kdBox;
    typedef typename kdBox::point point;

    /// axis in which the children of this node split the domain
    /// special value -1 denotes a leaf node
    int axis; 

    /// Split coordinate (meaningfull only for split nodes)
    Z pos;

    /// level in which the box in the node is completely contained
    /// special value -1 denotes unknown level (in case of a split node)
    int level ;

    /// The box held in this leaf node (if leaf)
    /// box->first is the lower left corner of the box
    /// box->second is the upper right corner of the box
    kdBox * box;

    /// Pointer to the parent node
    gsKdNode * parent;

    /// Pointer to the left child of this split node (if it is one)
    gsKdNode * left;

    /// Pointer to the right child of this split node (if it is one)
    gsKdNode * right;

    /// Constructor (empty node)
    gsKdNode() : axis(-2) ,level(0), box(0),
                 parent(0), left(0), right(0)
    { }

    /// Constructor (root node)
    gsKdNode(point const & upp) : axis(-1) , level(0),
                                  parent(0), left(0) , right(0)
    { 
        // Initial box, upp is expected to be indexed in finest level
        box = new kdBox( point::Zero(), upp);
    }

    /// Constructor (root node)
    gsKdNode(point const & low, point const & upp) : axis(-1) , level(0),
                                                     parent(0), left(0) , right(0)
    {
        // Initial box, upp is expected to be indexed in finest level
        box = new kdBox( low, upp);
    }

    /// Constructor (leaf node)
    gsKdNode(kdBox const & bb) : axis(-1) , level(0),
                                 parent(0), left(0) , right(0)
    { 
        // ..
    }
    
    /// Recursively copies the whole subtree under \a o, and sets it's
    /// parent to \a parentNode
    gsKdNode(const gsKdNode & o, gsKdNode * parentNode = NULL) : axis(o.axis), level(o.level)
    {
        parent = parentNode;
        if ( axis == -1 )
        {
            GISMO_ASSERT( (o.left == 0) && (o.right == 0), 
                          "Problem: leaf with children." );
            box  = new kdBox(*o.box);
            left = right = NULL;
        }
        else
        {
            GISMO_ASSERT( o.box == 0, 
                          "Problem: split node with box." );
            pos   = o.pos;
            left  = new gsKdNode(*o.left , this);
            right = new gsKdNode(*o.right, this);
            box   = NULL;
        }
    }

    /// Recursively deletes the whole subtree under this node
    ~gsKdNode()
    {
        // TODO: non-recursive

        if ( isLeaf() ) 
        {
            delete box;
        }
        else
        {
            delete left;
            delete right;
        }
    }

    // Box Accessors
    const point & lowCorner() const 
    { 
        GISMO_ASSERT(box, "Asked for lowCorner at node without box data.");
        return box->first ; 
    }

    const point & uppCorner() const 
    { 
        GISMO_ASSERT(box, "Asked for uppCorner at node without box data.");
        return box->second; 
    }

    bool isLeaf() const { return axis == -1; }

    bool isRoot() const { return parent == NULL; }

    bool isTerminal() const 
    { return (axis!=-1) && (left->axis==-1) && (right->axis==-1); }

    bool isLeftChild()  const { return parent!=NULL && this==parent->left; }

    bool isRightChild() const { return parent!=NULL && this==parent->right; }

    bool isDegenerate() const
    { return (box->first.array() >= box->second.array()).any(); }

    gsKdNode * sibling() const
    { 
        GISMO_ASSERT( parent != 0, "Root does not have a sibling.");
        return (parent->left == this ? parent->right : parent->left ); 
    }

    void multiplyByTwo()
    {
        if ( isLeaf() )
        {
            box->first .array() *= 2;
            box->second.array() *= 2;
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
            box->first .array() /= 2;
            box->second.array() /= 2;
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
        left          = new gsKdNode;
        right         = new gsKdNode;
        // Set axis to -1 (since they are leaves)
        left ->axis   =
        right->axis   = -1;
        // Set parent to this node
        left ->parent = 
        right->parent = this;
        // Set level
        left ->level  = 
        right->level  = level;
        // Set box
        left ->box    = box;    
        right->box    = new kdBox(*box);
        // Detach box from parent (is now at left child)
        box = NULL;
        // Resize properly the box coordinates
        left ->box->second[axis] = 
        right->box->first [axis] = pos;
    }

    /// Merges terminal node (i.e., two children are joined)
    inline void merge()
    {
        GISMO_ASSERT( (left->isLeaf()) && (right->isLeaf()),
                      "Can only merge terminal nodes.");

        // Recover box
        box = left->box;
        left->box = NULL;
        box->second[axis] = right->box->second[axis];
        axis  = - 1;
        level = left->level;

        // Delete children
        delete  left;
        left  = NULL;
        delete right;
        right = NULL;
    }


    /// Splits the node (i.e., two children are added)
    void split(int splitAxis, Z splitPos)
    {
        GISMO_ASSERT( box->second[splitAxis] != splitPos, "Degenerate split " << box->second[splitAxis] <<" != "<<splitPos);
        GISMO_ASSERT( box->first [splitAxis] != splitPos, "Degenerate split " << box->first[splitAxis]  <<" != "<<splitPos);
        axis = splitAxis;
        pos  = splitPos;
        split();
    }

    /// Splits the node in the middle (ie. two children are added)
    // TODO: remove
    void nextMidSplit()
    {        
        axis = ( parent == 0 ? 0 : (parent->axis+1)%d );        
        pos  = box->first [axis] + 
            (box->second[axis] - box->first[axis])/2 ;
        split(); // Can be degenerate
    }

    /// Splits the node in the middle (ie. two children are added)
    /// If non-degenerate split is impossible, then this is a no-op
    void anyMidSplit(int index_level)
    {        
        const unsigned h = 1 << (index_level - level) ;
        const unsigned mask = ~(h - 1);
        for ( unsigned i = 0; i < d; ++i )
        {
            const unsigned c = 
                (box->first [i] + (box->second[i] - box->first[i])/2) & mask ;
            if ( c != box->first [i] ) // avoid degenerate split
            {
                split(i, c);
                return;
            }
        }
    }


    /// Splits the node adaptively (i.e., two children are added)
    /// according to \a insBox.  If non-degenerate split is impossible,
    /// then this is a no-op.
    /// Splitting is done on a coordinate of the current \a level (aligned)
    /// returns the child that intersects \a insBox or NULL (if no split)
    gsKdNode * adaptiveAlignedSplit(kdBox const & insBox, int index_level)
    {
        const unsigned h = 1 << (index_level - level) ;
        
        for (short_t i = 0; i < d; ++i)
        {
            const Z c1 = insBox. first[i] - insBox. first[i] % h; //floor
            const Z cc = insBox.second[i] % h;
            const Z c2 = insBox.second[i] + (cc ? h-cc : 0 ); // ceil

            if ( c1 > box->first[i] )
            {
                // right child intersects insBox
                split(i, c1 );
                return right;
            }
            else if ( c2 < box->second[i]  )
            {
                // left child intersects insBox
                split(i, c2 );
                return left;
            }
        }
        return NULL;
    }

    /// Splits the node adaptively (i.e., two children are added)
    /// according to \a insBox. If non-degenerate split is impossible,
    /// then this is a no-op
    // TODO: remove
    gsKdNode * adaptiveSplit(kdBox const & insBox)
    {
        // assumption: insBox intersects box
        for ( unsigned i = 0; i < d; ++i )
        {
            // TODO: strategy: Try to split as close to the middle as possible.
            if ( insBox.first[i] > box->first[i] )
            {
                axis = i;
                pos  = insBox.first[i];
                split();
                return right;
            }
            else if ( insBox.second[i] < box->second[i] )
            {
                axis = i;
                pos  = insBox.second[i];
                split();
                return left ;
            }
        }
        // insBox is equal to this->box or they do not overlap
        return NULL;
    }

    friend std::ostream & operator<<(std::ostream & os, const gsKdNode & n)
    {
        if ( n.isLeaf() ) 
        {
            os << "Leaf node ("<< n.box->first.transpose() <<"), ("
               << n.box->second.transpose() <<"). level="<<n.level<<" \n";
        }
        else
        {
            os << "Split node, axis= "<< n.axis <<", pos="<< n.pos <<"\n";
        }

        return os;
    }

};


}// namespace gismo
