/** @file gsHDomain.hpp

    @brief Provides implementation of the HDomain class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): G. Kiss, A. Mantzaflaris
*/


/*
// Note 1: stride(lvl): 1 << (m_index_level-lvl);
// Note 2:Translate index from lvl1 to lvl2
inline unsigned translateIndex( unsigned const & i,
                                unsigned const & lvl1,
                                unsigned const & lvl2)
{
    if ( lvl1< lvl2 )
        return i << (lvl2-lvl1);
    else // lvl1 >= lvl2
      {
    if ( lvl1< lvl2 )
        return i >> (lvl1-lvl2);
    else
        return i;
      }
}
*/

#include <gsHSplines/gsAAPolyline.h>
#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsBoundary.h>

#include <queue>

namespace
{
    // Query 1
    struct query1_visitor
    {
        typedef bool return_type;

        // initialize result as true
        static const return_type init = true;

        template<short_t d, class T >
        static void visitLeaf(gismo::kdnode<d,T> * leafNode , int level, return_type & res)
        {
            //if ( (!isDegenerate(*leafNode->box)) && leafNode->level != level )
            if ( leafNode->level != level )
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

        template<short_t d, class T >
        static void visitLeaf(gismo::kdnode<d,T> * leafNode , int level, return_type & res)
        {
            //if ( (!isDegenerate(*leafNode->box)) && leafNode->level <= level )
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
        static const return_type init = 1000000;

        template<short_t d, class T >
        static void visitLeaf(gismo::kdnode<d,T> * leafNode , int , return_type & res)
        {
            //if ( (!isDegenerate(*leafNode->box)) && leafNode->level < res )
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
        static const return_type init = -1;

        template<short_t d, class T >
        static void visitLeaf(gismo::kdnode<d,T> * leafNode , int , return_type & res)
        {
            //if ( (!isDegenerate(*leafNode->box)) && leafNode->level > res )
            if ( leafNode->level > res )
                res = leafNode->level;
        }
    };

} //namespace

namespace gismo {


template<short_t d, class T >
gsHDomain<d,T> * gsHDomain<d,T>::clone() const
{
    return new gsHDomain(*this);
}


template<short_t d, class T > inline bool
gsHDomain<d,T>::haveOverlap(box const & box1, box const & box2)
{
    return !( (box1.second.array() <= box1.first .array()).any() ||
              (box2.first .array() >= box1.second.array()).any() );

/* // Equivalent implelementation :

    for( unsigned i = 0; i < d; i++ )
    {
        if( (box2.second[i] <= box1.first[i] ) ||
            (box2.first[i]  >= box1.second[i]) )
            return false;
    }
    return true;
*/
}

template<short_t d, class T > inline bool
gsHDomain<d,T>::isContained(box const & box1, box const & box2)
{
    return !( (box1.first .array() < box2.first .array()).any() ||
              (box1.second.array() > box2.second.array()).any() ) ;

/* // Equivalent implelementation :

    for( unsigned i = 0; i < d; i++ )
    {
        if( (box1.first [i] <  box2.first [i] ) ||
            (box1.second[i] >  box2.second[i]) )
            return false;
    }
    return true;
*/
}

template<short_t d, class T > void
gsHDomain<d,T>::setLevel(node *_node, int lvl)
{
    // to do: non-recc.
    if ( _node->isLeaf() )
    {
        // if the level in the *_node is smaller the value of lvl we
        // modify it
        _node->level = lvl;
    }
    else
    {
        // we call the function set level for all existing children of
        // the node *_node
        setLevel(_node->left , lvl);
        setLevel(_node->right, lvl);
    }
}

template<short_t d, class T> inline bool
gsHDomain<d,T>::isDegenerate(box const & someBox)
{
    return (someBox.first.array() >= someBox.second.array()).any() ;

/* // Equivalent implelementation :
    for( unsigned i = 0; i < d; i++ )
        if( someBox.first[i] >= someBox.second[i] )
        {
            return true;
        }
    return false;
*/
}

template<short_t d, class T > void
gsHDomain<d,T>::insertBox ( point const & k1, point const & k2,
                            node *_node, int lvl) // CONSTRAINT: lvl is "minimum level"
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
    stack.push_back(_node);  //push(_node);

    node * curNode;
    while ( ! stack.empty() )
    {
        curNode = stack.back(); //top();
        stack.pop_back();       //pop();

/*
        if ( curNode->is_Node() ) // reached a leaf
        {
            if ( isDegenerate(*curNode->box) )
                continue;

            if ( isContained(*curNode->box, iBox) )
            {
                if ( lvl > curNode->level )
                    curNode->level = lvl;
            }
            else if ( haveOverlap(*curNode->box, iBox) )
            {
                curNode->nextMidSplit();
                stack.push_back(curNode);
            }
        }
        //else..
*/
        if ( curNode->isLeaf() ) // reached a leaf
        {
            // Since we reached a leaf, it should overlap with iBox

            // If this leaf is already in level lvl, then we have nothing to do
            if ( lvl <= curNode->level )
                continue;

            // Split the leaf (if possible)
            //node * newLeaf = curNode->adaptiveSplit(iBox);
            node * newLeaf = curNode->adaptiveAlignedSplit(iBox, m_indexLevel);

            // If curNode is still a leaf, its domain is almost
            // contained in iBox
            if ( !newLeaf ) //  curNode->isLeaf()
            {
                // Increase level and reccurse
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

template<short_t d, class T > void
gsHDomain<d,T>::sinkBox ( point const & k1,
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
    //stack.reserve( 2 * m_maxPath );
    stack.push(m_root);

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

template<short_t d, class T > void
gsHDomain<d,T>::makeCompressed()
{
    std::stack<node*, std::vector<node*> > tstack;
    node * curNode;

    // First step: gather all terminal nodes
    std::stack<node*, std::vector<node*> > stack;
    stack.push(m_root);
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

        if (curNode->left->level == curNode->right->level)
        {
            // Merge left and right
            curNode->merge();
            if ( !curNode->isRoot() &&
                  curNode->parent->isTerminal() )
                tstack.push(curNode->parent );
        }
    }

    // Store the max path length
    m_maxPath = minMaxPath().second;
}

template<short_t d, class T >
bool gsHDomain<d,T>::query1(point const & lower, point const & upper,
               int level, node  *_node) const
{ return boxSearch< query1_visitor >(upper,lower,level,_node); }

template<short_t d, class T >
bool gsHDomain<d,T>::query1(point const & lower, point const & upper,
               int level) const
{ return boxSearch< query1_visitor >(upper,lower,level,m_root); }

template<short_t d, class T >
bool gsHDomain<d,T>::query2(point const & lower, point const & upper,
               int level, node  *_node) const
{ return boxSearch< query2_visitor >(lower,upper,level,_node); }

template<short_t d, class T >
bool gsHDomain<d,T>::query2 (point const & lower, point const & upper,
                 int level) const
{ return boxSearch< query2_visitor >(lower,upper,level,m_root); }

template<short_t d, class T >
int gsHDomain<d,T>::query3(point const & lower, point const & upper,
               int level, node  *_node) const
{ return boxSearch< query3_visitor >(lower,upper,level,_node); }

template<short_t d, class T >
int gsHDomain<d,T>::query3(point const & lower, point const & upper,
               int level) const
{ return boxSearch< query3_visitor >(lower,upper,level,m_root); }

template<short_t d, class T >
int gsHDomain<d,T>::query4(point const & lower, point const & upper,
               int level, node  *_node) const
{ return boxSearch< query4_visitor >(lower,upper,level,_node); }

template<short_t d, class T >
int gsHDomain<d,T>::query4(point const & lower, point const & upper,
               int level) const
{ return boxSearch< query4_visitor >(lower,upper,level,m_root); }

template<short_t d, class T >
std::pair<typename gsHDomain<d,T>::point, typename gsHDomain<d,T>::point>
gsHDomain<d,T>::select_part(point const & k1, point const & k2,
                            point const & k3, point const & k4)
{
    // intersect boxes
    std::pair<point,point> result;

    for ( unsigned i = 0; i<d; ++i)
    {
        //find the lower left corner
        result.first[i]  = ( k1[i] >= k3[i] ? k1[i] : k3[i] );
        //find the upper right corner
        result.second[i] = ( k2[i] >= k4[i] ? k4[i] : k2[i] );
    }
    return result;
}

template<short_t d, class T > void
gsHDomain<d,T>::bisectBox(box const & original, int k, T coord,
                          box & leftBox, box & rightBox )
{
    GISMO_ASSERT( ! isDegenerate(original) , "Invalid box .");
    GISMO_ASSERT( (k>=0) && (k< static_cast<int>(d)) , "Invalid axis "<< k <<".");
    leftBox = rightBox = original;
    leftBox.second[k] = rightBox.first[k] = coord;
}


template<short_t d, class T>
template<typename visitor>
typename visitor::return_type
gsHDomain<d,T>::boxSearch(point const & k1, point const & k2,
                          int level, node  *_node ) const
{
    // Make a box
    box qBox(k1,k2);
    local2globalIndex( qBox.first , static_cast<unsigned>(level), qBox.first );
    local2globalIndex( qBox.second, static_cast<unsigned>(level), qBox.second);

    GISMO_ASSERT( !isDegenerate(qBox),
                  "boxSearch: Wrong order of points defining the box (or empty box): "
                  << qBox.first.transpose() <<", "<< qBox.second.transpose() <<".\n" );

    typename visitor::return_type res = visitor::init;

/*  // under construction
    node * curNode = m_root;

    while (true)
    {
        if ( curNode->isLeaf() )
        {
            visitor::visitLeaf(curNode, level, res );

            //curNode->isRightChild();
            while ( curNode->parent != NULL &&
                    curNode != curNode->parent->left )
            {
                curNode = curNode->parent;

                if ( curNode->isLeftChild() &&
                     qBox.second[curNode->axis] <= curNode->pos )
                    //curNode = curNode->parent;
                    curNode = curNode->parent->right;
            }

            if ( curNode->parent == NULL )
                break;
            else// Found a left child,follow right simbling
                curNode = curNode->parent->right;


            // while ( curNode->parent != NULL &&
            //         curNode != curNode->parent->left)
            //     curNode = curNode->parent;

            // if ( curNode->parent == NULL )
            //     break;
            // else
            //     curNode = curNode->parent;

            // if  ( qBox.second[curNode->axis] > curNode->pos )  // overlap right ?
            //     curNode = curNode->right;
            // else
            // {
            //     while ( curNode->parent != NULL &&
            //             curNode != curNode->parent->left)
            //         curNode = curNode->parent;
            //     if ( curNode->parent == NULL )
            //         break;
            //     else
            //         curNode = curNode->parent;

            //     if  ( qBox.second[curNode->axis] > curNode->pos )  // overlap right ?
            //         curNode = curNode->right;
            //     else
            //         break;
            // }

        }
        else
        {
            if ( qBox.first[curNode->axis] < curNode->pos ) // overlap left ?
                curNode = curNode->left;
            else if  ( qBox.second[curNode->axis] > curNode->pos) // overlap right ?
                curNode = curNode->right;
            else
                break;
        }
    }

// */

// /* // implementation with stack

    std::vector<node*> stack;
    stack.reserve( 2 * m_maxPath );
    stack.push_back(_node);  //push(_node);

    node * curNode;
    while ( ! stack.empty() )
    {
        curNode = stack.back(); //top();
        stack.pop_back();       //pop();

        if ( curNode->isLeaf() )
        {
            // Visit the leaf
            GISMO_ASSERT( !isDegenerate(*curNode->box), "Encountered an empty leaf");
            visitor::visitLeaf(curNode, level, res );
        }
        else // this is a split-node
        {
            if ( qBox.second[curNode->axis] <= curNode->pos)
                // qBox overlaps only left child of this split-node
                stack.push_back(curNode->left); //push(curNode->left);
            else if  ( qBox.first[curNode->axis] >= curNode->pos)
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
//*/
    return res;
}



template<short_t d, class T>
typename gsHDomain<d,T>::node *
gsHDomain<d,T>::pointSearch(const point & p, int level, node  *_node ) const
{
    point pp;
    local2globalIndex(p, static_cast<unsigned>(level), pp);

    GISMO_ASSERT( ( pp.array() <= m_upperIndex.array() ).all(),
        "pointSearch: Wrong input: "<< p.transpose()<<", level "<<level<<".\n" );

    std::vector<node*> stack;
    stack.reserve( 2 * m_maxPath );
    stack.push_back(_node);  //push(_node);

    node * curNode;
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
                stack.push_back(curNode->left); //push(curNode->left);
            else
                stack.push_back(curNode->right); //push(curNode->right);
        }
    }
    GISMO_ERROR("pointSearch: Error ("<< p.transpose()<<").\n" );
}



/*
   Function to traverse the tree nodes without recursion and without
   stack, when parents are not stored. Not thread-safe
//
template<short_t d, class T>
template<typename visitor>
typename visitor::return_type
gsHDomain<d,T>::nodeSearchMorris() const
{
    node * curNode, * preNode;

    typename visitor::return_type i = visitor::init;
    curNode = m_root;

    while(curNode != NULL)
    {
        if(curNode->left == NULL)
        {
            visitor::visitNode(curNode, i);
            curNode = curNode->right;
        }
        else
        {
            // Find the inorder predecessor of curNode
            preNode = curNode->left;
            while(preNode->right != NULL && preNode->right != curNode)
                preNode = preNode->right;

            // Make curNode as right child of its inorder predecessor
            if(preNode->right == NULL)
            {
                preNode->right = curNode;
                curNode = curNode->left;
            }
            else
            {   // Revert the changes made in "if" part to restore the original
                preNode->right = NULL;
                visitor::visitNode(curNode, i);
                curNode = curNode->right;
            } // end if(preNode->right == NULL)
        } // end if (curNode->left == NULL)
    } // end while

    return i;
}
*/

template<short_t d, class T>
template<typename visitor>
typename visitor::return_type
gsHDomain<d,T>::nodeSearch() const
{
    typename visitor::return_type i = visitor::init;

    node * curNode = m_root;

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

/*
// node search version with stack (if no parents available)
template<short_t d, class T>
template<typename visitor>
typename visitor::return_type
gsHDomain<d,T>::nodeSearch() const
{
    typename visitor::return_type i = visitor::init;
    std::stack<node*, std::vector<node*> > stack;
    stack.push(m_root);

    node * curNode;
    while ( ! stack.empty() )
    {
        curNode = stack.top();
        stack.pop();
        visitor::visitNode(curNode, i);

        if ( ! curNode->isLeaf() )
        {
                stack.push(curNode->left );
                stack.push(curNode->right);
        }
    }
    return i;
}
//*/

template<short_t d, class T>
template<typename visitor>
typename visitor::return_type
gsHDomain<d,T>::leafSearch() const
{
    typename visitor::return_type i = visitor::init;

    node * curNode = m_root;

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


/*
// leaf search version with stack (if no parents available)
template<short_t d, class T>
template<typename visitor>
typename visitor::return_type
gsHDomain<d,T>::leafSearch() const
{
    typename visitor::return_type i = visitor::init;
    std::stack<node*, std::vector<node*> > stack;
    stack.push(m_root);

    node * curNode;
    while ( ! stack.empty() )
    {
        curNode = stack.top();
        stack.pop();

        if ( curNode->isLeaf() )
        {
            // Visit the leaf
            visitor::visitLeaf(curNode, i);
        }
        else // this is a split-node
        {
                stack.push(curNode->left );
                stack.push(curNode->right);
        }
    }
    return i;
}
//*/

template<short_t d, class T >
std::pair<int,int>
gsHDomain<d,T>::minMaxPath() const
{
    node * curNode = m_root;
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


/*
template<short_t d, class T > int
gsHDomain<d,T>::query3Recur(box const & qBox, node *_node) const
{
    // Note: reccursive implementation of query3, qBox assumed in m_index_level indices.
    // Is kept here as an example of reccursive implementation
    GISMO_ASSERT(_node != NULL, "invalid node.");

    if( isDegenerate(qBox) )
        GISMO_ERROR("query3 says: Wrong order of points defining the box (or empty box)."
                    << qBox.first.transpose() <<", "<< qBox.second.transpose() <<".\n" );

    if ( _node->isLeaf() )
    {
        return _node->level;
    }
    else // Move down the tree
    {
        if ( qBox.second[_node->axis] <= _node->pos) //  qBox inside the left child
            return query3Recur(qBox, _node->left);
        else if  ( qBox.first[_node->axis] >= _node->pos) // qBox inside the right child
            return query3Recur(qBox, _node->right);
        else // qBox needs to be bisected into left and right
        {
            // Bisect qBox at (_node->axis, _node->pos)
            box leftBox, rightBox;
            bisectBox( qBox, _node->axis, _node->pos, leftBox, rightBox);
            // Call reccusively for both children
            const int levLeft  = query3Recur(leftBox , _node->left );
            const int levRight = query3Recur(rightBox, _node->right);
            // return the minimum
            return ( levLeft < levRight ? levLeft : levRight );
        }
    }
}
*/


template<short_t d, class T>
void gsHDomain<d,T>::getBoxes(gsMatrix<index_t>& b1, gsMatrix<index_t>& b2, gsVector<index_t>& level) const
{
    std::vector<std::vector<index_t> > boxes;

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
    for(size_t i = 0; i < boxes.size(); i++){
        for(short_t j = 0; j < d; j++){
            b1(i,j) = boxes[i][j];
            b2(i,j) = boxes[i][j+d];
        }
        level[i] = boxes[i][2*d];
    }
}


template<short_t d, class T>
void gsHDomain<d,T>::getBoxesOnSide(boundary::side s, gsMatrix<index_t>& b1, gsMatrix<index_t>& b2, gsVector<index_t>& level) const
{

    getBoxes( b1, b2, level);
    std::vector<index_t> onSide;

    unsigned remainder = (s-1) % 2;
    // remainder will be
    // 0 for sides 1 (west), 3 (south/down), or 5 (front)
    // 1 for sides 2 (east), 4 (north/up),   or 6 (back)
    unsigned quotient = ( (s-1) - remainder ) / 2;
    // quotient will be
    // 0 for east/west
    // 1 for south/north or down/up
    // 2 for front/back

    if( remainder == 0 ) // sides 1 (west), 3 (south/down), or 5 (front)
    {
        // if a box touches
        for( index_t i = 0; i < b1.rows(); i++)
            if( b1(i, quotient ) == 0 )
                onSide.push_back(i);
    }
    else // remainder == 1, sides east, north/up, or back
    {
        for( index_t i = 0; i < b1.rows(); i++)
        {
            // index of upper corner
            index_t B2( b2(i, quotient ) );
            // transform to index-level
            B2 = B2 << (m_indexLevel - m_maxInsLevel);

            // check if this index is on the boundary
            if( B2 == m_upperIndex[ quotient ] )
                onSide.push_back( i );
        }
    }

    // select only the boxes on side s:
    for( size_t i=0; i < onSide.size(); i++)
    {
        b1.row(i) = b1.row( onSide[i] );
        b2.row(i) = b2.row( onSide[i] );
        level[i]  = level( onSide[i] );
    }
    b1.conservativeResize( onSide.size(), b1.cols() );
    b2.conservativeResize( onSide.size(), b2.cols() );
    level.conservativeResize( onSide.size() );
}

template<short_t d, class T>
void gsHDomain<d,T>::getBoxesInLevelIndex(gsMatrix<index_t>& b1,
              gsMatrix<index_t>& b2,
              gsVector<index_t>& level) const{
    std::vector<std::vector<index_t> > boxes;
    getBoxes_vec(boxes);
    GISMO_ASSERT(d==2 || d==3,"Wrong dimension, should be 2 or 3.");
    //is this test really necessary? florian b.
    for(size_t i = 0; i < boxes.size(); i++){
        if ((boxes[i][0]==boxes[i][d+0]) || (boxes[i][1]==boxes[i][1+d]))
        {
            boxes.erase(boxes.begin()+i);
            i--;
        }
        else if( (d == 3) && ( boxes[i][2]==boxes[i][d+2] ) )
        {
            boxes.erase(boxes.begin()+i);
            i--;
        }
    }
    gsVector<index_t,d>lowerCorner;
    gsVector<index_t,d>upperCorner;
    connect_Boxes(boxes);
    b1.resize(boxes.size(),d);
    b2.resize(boxes.size(),d);
    level.resize(boxes.size());
    for(size_t i = 0; i < boxes.size(); i++)
    {
        for(short_t j = 0; j < d; j++)
        {
//            b1(i,j) = boxes[i][j];
//            b2(i,j) = boxes[i][j+d];
            lowerCorner[j] = boxes[i][j];
            upperCorner[j] = boxes[i][j+d];
        }
        level[i] = boxes[i][2*d];
        computeLevelIndex( lowerCorner , level[i], lowerCorner );
        computeLevelIndex( upperCorner , level[i], upperCorner );
        b1.row(i) = lowerCorner;
        b2.row(i) = upperCorner;
    }
}

// Old version that only works for 2D
// Should be replaced by the "new" connect_Boxes.
// Keeping the code for the moment, in order not to loose
// the old code before the new one is properly tested.
template<short_t d, class T> void
gsHDomain<d,T>::connect_Boxes2d(std::vector<std::vector<index_t> > &boxes) const
{
    GISMO_ASSERT( d == 2, "This one only works for 2D");
    bool change = true;
    while(change){
        change =  false;
        int s = boxes.size();
        for(int i = 0; i < s;i++)
        {
            for(int j = i+1; j < s; j++)
            {
                // if the levels are the same:
                if(boxes[i][4]==boxes[j][4])
                {
                    if( (boxes[i][0]==boxes[j][0]) && (boxes[i][2]==boxes[j][2]))
                    {
                        if(boxes[i][1]==boxes[j][3])
                        {
                            boxes[i][1] = boxes[j][1];
                            boxes.erase(boxes.begin()+j);
                            s--;
                            j--;
                            change =  true;
                        }
                        if(boxes[i][3]==boxes[j][1])
                        {
                            boxes[i][3] = boxes[j][3];
                            boxes.erase(boxes.begin()+j);
                            s--;
                            j--;
                            change =  true;
                        }
                    }
                    if( (boxes[i][1]==boxes[j][1]) && (boxes[i][3]==boxes[j][3])
)
                    {
                        if(boxes[i][0]==boxes[j][2])
                        {
                            boxes[i][0] = boxes[j][0];
                            boxes.erase(boxes.begin()+j);
                            s--;
                            j--;
                            change =  true;
                        }
                        if(boxes[i][2]==boxes[j][0])
                        {
                            boxes[i][2] = boxes[j][2];
                            boxes.erase(boxes.begin()+j);
                            s--;
                            j--;
                            change =  true;
                        }
                    }
                }
            }
        }
    }
    //gsDebug<<"in the connecting fucntion"<<boxes.size()<<std::endl;
}

template<short_t d, class T> void
gsHDomain<d,T>::connect_Boxes(std::vector<std::vector<index_t> > &boxes) const
{
    bool change = true;
    while(change)
    {
        change =  false;
        int s = boxes.size();
        for(int i = 0; i < s;i++)
        {
        for(int j = i+1; j < s; j++)
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
            for( unsigned k=0; k < d; k++)
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

            // The boxes can only be merged, if
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


template<short_t d, class T>
void gsHDomain<d,T>::connect_Boxes_2(std::vector<std::vector<index_t> > &boxes) const
{
    bool change = true;
    while(change){
        change =  false;
        int s = boxes.size();
        for(int i = 0; i < s;i++)
        {
            for(int j = i+1; j < s; j++)
            {
                if(boxes[i][4]==boxes[j][4]){
                    if( (boxes[i][0]==boxes[j][0]) && (boxes[i][2]==boxes[j][2]))
                    {
                        if(boxes[i][1]==boxes[j][3])
                        {
                            boxes[i][1] = boxes[j][1];
                            boxes.erase(boxes.begin()+j);
                            s--;
                            j--;
                            change =  true;
                        }
                        if(boxes[i][3]==boxes[j][1])
                        {
                            boxes[i][3] = boxes[j][3];
                            boxes.erase(boxes.begin()+j);
                            s--;
                            j--;
                            change =  true;
                        }
                    }
                }
            }
        }
    }
    change = true;
    while(change){
        change =  false;
        int s = boxes.size();
        for(int i = 0; i < s;i++)
        {
            for(int j = i+1; j < s; j++)
            {
                if( (boxes[i][1]==boxes[j][1]) && (boxes[i][3]==boxes[j][3]))
                {
                    if(boxes[i][0]==boxes[j][2])
                    {
                        boxes[i][0] = boxes[j][0];
                        boxes.erase(boxes.begin()+j);
                        s--;
                        j--;
                        change =  true;
                    }
                    if(boxes[i][2]==boxes[j][0])
                    {
                        boxes[i][2] = boxes[j][2];
                        boxes.erase(boxes.begin()+j);
                        s--;
                        j--;
                        change =  true;
                    }
                }
            }
        }
    }




    //gsDebug<<"in the connecting fucntion"<<boxes.size()<<std::endl;
}



template<short_t d, class T> void
gsHDomain<d,T>::getBoxes_vec(std::vector<std::vector<index_t> >& boxes) const
{
    boxes.clear();

    std::stack<node*, std::vector<node*> > stack;
    //stack.reserve( 2 * m_maxPath );
    stack.push(m_root);
    node * curNode;
    while ( ! stack.empty() )
    {
        curNode = stack.top();
        stack.pop();

        if ( curNode->isLeaf() )
        {
            // We need to convert the indices to those of m_maxInsLevel
            // to be able to reconstruct the earlier results.
            const point & lowerGlob = curNode->lowCorner();
            const point & upperGlob = curNode->uppCorner();
            unsigned int level = this->m_maxInsLevel;
            point lower;
            point upper;

            global2localIndex(lowerGlob,level,lower);
            global2localIndex(upperGlob,level,upper);

            boxes.push_back(std::vector<index_t>());
            for(unsigned i = 0; i < d; i++)
            {
                boxes.back().push_back(lower[i]);
            }
            for(unsigned i = 0; i < d; i++)
            {
                boxes.back().push_back(upper[i]);
            }
            boxes.back().push_back(curNode->level);
        }
        else
        {
            stack.push(curNode->left) ;
            stack.push(curNode->right);
        }
    }
}

/*
 * functions for returning the boudaries of domains
 */
template<short_t d, class T>
std::vector< std::vector<std::vector< std::vector< index_t > > > >
gsHDomain<d,T>::getPolylines() const
{
/*
 1: Get boxes from the quadtree.
 2: Return their borderlines.

 Return structure:
 < levels < polylines_in_one_level < one_polyline < one_segment (x1, y1, x2, y2) > > > > result
 note that <x1, y1, x2, y2 > are so that (x1, y1) <=LEX  (x2, y2)
*/
    std::vector<std::vector<index_t> > boxes;
    getBoxes_vec(boxes);// Returns all leaves.

    // Get rid of boxes that are not of full dimension.
    for( std::vector< std::vector< index_t> >::iterator it = boxes.begin(); it != boxes.end(); ++it )
    {
        if( ( (*it)[0] == (*it)[2] ) || ( (*it)[0] == (*it)[2] ) )
            it = boxes.erase(it);
    }

    std::vector<std::vector<gsVSegment<T> > > seg;
    seg.resize(m_maxInsLevel+1);

    // For each level prepare the vertical lines separately
    for (unsigned int i = 0; i < boxes.size() ; i ++)
    {
        seg[boxes[i][4]].push_back(gsVSegment<index_t>(boxes[i][0],boxes[i][1],boxes[i][3], false) );
        seg[boxes[i][4]].push_back(gsVSegment<index_t>(boxes[i][2],boxes[i][1],boxes[i][3], false) );
    }

    // Process vertical lines from each level separately
    std::vector< std::vector<std::vector< std::vector<index_t > > > > result;
    for(unsigned int i = 0; i < m_maxInsLevel+1; i++)
    {
       //result.push_back(getPoly(seg[i]));
        result.push_back( getPolylinesSingleLevel( seg[i] ));
    }
    //result.push_back(getPolylinesSingleLevel(seg[0]));

    return result;
}


template<short_t d, class T>
std::vector<std::vector< std::vector<index_t > > > gsHDomain<d,T>::getPolylinesSingleLevel(std::vector<gsVSegment<T> >& seg) const
{
    // For didactic purposes the interior of the function has been refined into two procedures
    // (overal length of the older version was intimidating and people would not like to read it then =)).

    // For each x coordinate there will be a list of vertical segments with this x coord.
    std::list< std::list< gsVSegment< T > > > vert_seg_lists;

    // For returning
    std::vector< std::vector< std::vector<index_t > > > result;

    // This magically sorts seg according to x value
    std::sort( seg.begin(), seg.end() );

    // Put stuff from seg into vert_seg_lists
    std::list< gsVSegment< index_t > > segs_x; // segments with the same particular x coord
    for( typename std::vector< gsVSegment< T > >::const_iterator it_seg = seg.begin(); it_seg != seg.end(); ++it_seg )
    {
        if( segs_x.empty() || (*it_seg).getX() == segs_x.front().getX() )
            segs_x.push_back( *it_seg );
        else
        {
            vert_seg_lists.push_back( segs_x );
            segs_x.erase( segs_x.begin(), segs_x.end() );
            segs_x.push_back( *it_seg );
        }
    // Possible simplification: if( !empty && != ){ verts_seg.push_back; erase }segs_x.push_back;
    }
    vert_seg_lists.push_back( segs_x );
    // So vert_seg_lists is now the list as specified before and segs_x we do not need any further.

    // Get rid of overlaps of the segments.
    getRidOfOverlaps( vert_seg_lists );

    // And now the key point: sweepline driven connect and merge procedure
    sweeplineConnectAndMerge( result, vert_seg_lists );
    return result;
}

template <short_t d, class T>
void gsHDomain<d,T>::getRidOfOverlaps( std::list< std::list< gsVSegment<T> > >& vert_seg_lists ) const
{
    bool need_to_erase = false;
    for( typename std::list< std::list< gsVSegment<T> > >::iterator it_x = vert_seg_lists.begin(); it_x != vert_seg_lists.end(); )
    {
        // For each x coordinate substract overlapping segments.
        for( typename std::list< gsVSegment<T> >::iterator it_slow = (*it_x).begin(); it_slow != (*it_x).end(); ++it_slow)
        {
            // Don't we need to start the next already at (*it_x).begin() ?
            for( typename std::list< gsVSegment<T> >::iterator it_quick = it_slow; it_quick != (*it_x).end(); ++it_quick)
            {
                if( it_quick != it_slow )
                {
                    need_to_erase = it_slow->cannotDeleteOverlap( *it_quick) || need_to_erase;
                    // We do not actually use need_to_erase as initially intended
                    // but this prevents warnings (cannot_delete_overlap returns bool).
                }
            }
        }

        // Get rid of segments of zero length
        //if( need_to_erase )
        //{
        for( typename std::list< gsVSegment<T> >::iterator it_slow = (*it_x).begin(); it_slow != (*it_x).end();)
        {
            if( (*it_slow).length() == 0 )
            {
                (*it_x).erase( it_slow++ );
            }
            else
                ++it_slow; // Here instead of the internal for.
        }
        //need_to_erase = false;
        //}

        // There might be some empty lists;
        if( (*it_x).empty() )
            vert_seg_lists.erase( it_x++ );
        else
            ++it_x; // Here instead of the outer for.
    }
}

template <short_t d, class T>
void gsHDomain<d,T>::sweeplineConnectAndMerge( std::vector< std::vector< std::vector<index_t> > >& result,
                                               std::list< std::list< gsVSegment<T> > >& vert_seg_lists ) const
{
    /*========================
     * Algorithm description:
     *========================
     * The segments in vert_seg_lists are vertical sides of the boxes from the quadtree after getting rid of overlapping ones.
     * The sweepline proceeds from left to right.
     * Events: new x-coordinate
     * State: list of active (or open) polylines (their part to the left from the sweepline is non-empty and does not form an open polyline yet).
     * In each event we:
     * 1] Try to append vertical segments to the active lines. Report those that get closed.
     * 2] Remaining segments become polylines.
     * 3] Try to merge active polylines. Report those that get closed.
     * Remark: in kissing vertices (point, where two vertical and two horizontal line segments meet) we are allowed to turn left or right
     * BUT we are not allowed to go straight there. In other words, the interior must remain on the same hand side (left or right) the whole time.
     */

    // vector of active polylines
    std::list< gsAAPolyline<T> > act_poly;
    std::queue< gsAAPolyline<T> > poly_queue;

    for( typename std::list< std::list< gsVSegment<T> > >::iterator it_x = vert_seg_lists.begin(); it_x != vert_seg_lists.end(); ++it_x )
    {
        // Try to extend polylines by adding segments to them.
        // The funny queue procedure is a trick to make sure that none of the polylines
        // would be extended too much (i.e., over a vertex, where it should be closed).
        for( typename std::list< gsAAPolyline<T> >::const_iterator it = act_poly.begin(); it != act_poly.end(); ++it )
            poly_queue.push( *it );

        act_poly.erase( act_poly.begin(), act_poly.end() );

        gsAAPolyline<index_t> curr_poly;
        while( !poly_queue.empty() )
        {
            curr_poly = poly_queue.front();
            poly_queue.pop();
            bool is_active = true;
            for( typename std::list< gsVSegment<T> >::iterator it_seg = (*it_x).begin(); it_seg != (*it_x).end(); )
            {
                if( curr_poly.canBeExtended( *it_seg) )
                {
                    if( curr_poly.almostClosed() )
                        result.push_back( curr_poly.writeParasolid () );
                    else
                        poly_queue.push( curr_poly ); // This is the prevention of wrong behaviour in kissing vertices (cf. Remark above).


                    (*it_x).erase( it_seg++ );
                    is_active = false;
                    break;
                }
                else
                    ++it_seg;
            }
            if( is_active )
            {
                act_poly.push_back( curr_poly );
            }
        }

        // The remaining segments become polylines.
        for( typename std::list<gsVSegment<T> >::const_iterator it = (*it_x).begin(); it != (*it_x).end(); ++it )
        {
            gsAAPolyline<T> new_polyline( *it );
            act_poly.push_back( new_polyline );
        }

        // Merge, what's possible.
        // Trick: go through act_poly from left to right. Check each of them with all to the right of it in act_poly.
        // If they can be merged, move it to the end. Those more to the left are already finalised.
        // Could be replaced by sorting them first according to y-coordinate (David's suggestion) and then connecting
        // but I ran into some implementation trouble.
        for( typename std::list< gsAAPolyline<T> >::iterator it_slow = act_poly.begin(); it_slow != act_poly.end(); ++it_slow )
        {
            for( typename std::list< gsAAPolyline<T> >::iterator it_quick = it_slow; it_quick != act_poly.end(); ++it_quick)
            {
                if(( it_slow != it_quick ) && ((*it_slow).mergeWith(*it_quick)))
                {
                    act_poly.erase( it_quick );
                    act_poly.push_back( *it_slow );
                    it_slow = act_poly.erase( it_slow );
                    it_quick = it_slow;
                }
            }
        }

        // Check for closed lines as a result of the merging.
        for( typename std::list< gsAAPolyline<T> >::iterator it = act_poly.begin(); it != act_poly.end(); )
        {
            if( (*it).almostClosed() )
            {
                result.push_back( (*it).writeParasolid() );
                act_poly.erase( it++ );
            }
            else
                ++it;

        }
    }
}
    // REMARK to remember: whenever 'discard qualifiers' error, it refers to something that isn't const or to volatile


template<short_t d, class T> inline void
 gsHDomain<d,T>::computeFinestIndex( gsVector<index_t,d> const & index,
                                     unsigned lvl,
                                     gsVector<index_t,d> & result ) const
{
    for ( short_t i = 0; i!=d; ++i )
        result[i] = index[i] << (m_maxInsLevel-lvl) ;
}

template<short_t d, class T> inline void
 gsHDomain<d,T>::computeLevelIndex( gsVector<index_t,d> const & index,
                                    unsigned lvl,
                                    gsVector<index_t,d> & result ) const
{
    for ( short_t i = 0; i!=d; ++i )
        result[i] = index[i] >> (m_maxInsLevel-lvl) ;
}

template<short_t d, class T> inline void
 gsHDomain<d,T>::local2globalIndex( gsVector<index_t,d> const & index,
                    unsigned lvl,
                    gsVector<index_t,d> & result
                    ) const
{
    for ( short_t i = 0; i!=d; ++i )
        result[i] = index[i] << (m_indexLevel-lvl) ;
}

template<short_t d, class T> inline void
gsHDomain<d,T>::global2localIndex( gsVector<index_t,d> const & index,
                                        unsigned lvl,
                                        gsVector<index_t,d> & result
    ) const
{
    for ( short_t i = 0; i!=d; ++i )
        result[i] = index[i] >> (this->m_indexLevel-lvl) ;
}

template<short_t d, class T> inline void
gsHDomain<d,T>::incrementLevel()
{
    m_maxInsLevel++;

    GISMO_ASSERT( m_maxInsLevel <= m_indexLevel,
                  "Problem with indices, increase number of levels (to do).");

    leafSearch< levelUp_visitor >();
}

template<short_t d, class T> inline void
gsHDomain<d,T>::multiplyByTwo()
    {
        m_upperIndex *= 2;
        nodeSearch< liftCoordsOneLevel_visitor >();
    }

template<short_t d, class T> inline void
gsHDomain<d,T>::decrementLevel()
    {
        m_maxInsLevel--;
        leafSearch< levelDown_visitor >();
    }

template<short_t d, class T> inline int
gsHDomain<d,T>::size() const
    {
        return nodeSearch< numNodes_visitor >();
    }

template<short_t d, class T> inline int
gsHDomain<d,T>::leafSize() const
    { return leafSearch< numLeaves_visitor >(); }

template<short_t d, class T> inline void
gsHDomain<d,T>::printLeaves() const
    { leafSearch< printLeaves_visitor >(); }

}// end namespace gismo
