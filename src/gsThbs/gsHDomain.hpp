/** @file gsHDomain.hpp

    @brief Provides implementation of the HDomain class.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): G. Kiss, A. Mantzaflaris
*/

#include <gsThbs/gsAAPolyline.h>
#include <gsCore/gsLinearAlgebra.h>

#include <queue>

namespace gismo {


template<unsigned d, class T >
gsHDomain<d,T> * gsHDomain<d,T>::clone() const
{
    return new gsHDomain(*this);
}


template<unsigned d, class T > inline bool
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

template<unsigned d, class T > inline bool
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

template<unsigned d, class T > void
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

template<unsigned d, class T> inline bool 
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

template<unsigned d, class T > void
gsHDomain<d,T>::insertBox ( point const & k1, point const & k2,
                             node *_node, int lvl)
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
    std::stack<node*> stack;
    stack.push(_node); //start from root
    
    node * curNode;
    while ( ! stack.empty() )
    {
        curNode = stack.top();
        stack.pop();
        
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
                stack.push(curNode);
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
                    stack.push(curNode);
            }
            else // treat new child
            {
                stack.push(newLeaf);
            }
        }
        else // roll down the tree
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

    // Update maximum inserted level
    if ( static_cast<unsigned>(lvl) > m_maxInsLevel)
        m_maxInsLevel = lvl;
}

template<unsigned d, class T > void
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
    std::stack<node*> stack;
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

template<unsigned d, class T > void
gsHDomain<d,T>::makeCompressed()
{
    std::stack<node*> tstack;
    node * curNode;

    // First step: gather all terminal nodes
    std::stack<node*> stack;
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
}


template<unsigned d, class T >
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

template<unsigned d, class T > void
gsHDomain<d,T>::bisectBox(box const & original, int k, T coord,
                          box & leftBox, box & rightBox )
{
    GISMO_ASSERT( ! isDegenerate(original) , "Invalid box .");
    GISMO_ASSERT( (k>=0) && (k< static_cast<int>(d)) , "Invalid axis "<< k <<".");
    leftBox = rightBox = original;
    leftBox.second[k] = rightBox.first[k] = coord; 
}


template<unsigned d, class T> 
template<typename visitor>
typename visitor::return_type
gsHDomain<d,T>::boxSearch(point const & k1, point const & k2, 
                          int level, node  *_node ) const
{
    // Make a box
    box qBox(k1,k2);
    local2globalIndex( qBox.first , static_cast<unsigned>(level), qBox.first );
    local2globalIndex( qBox.second, static_cast<unsigned>(level), qBox.second);

    if( isDegenerate(qBox) )
        GISMO_ERROR("query3 says: Wrong order of points defining the box (or empty box): "
                    << qBox.first.transpose() <<", "<< qBox.second.transpose() <<".\n" );

    std::stack<node*> stack;
    stack.push(_node); 
    // initialize result
    typename visitor::return_type res = visitor::init;

    node * curNode;
    while ( ! stack.empty() )
    {
        curNode = stack.top();
        stack.pop();
        
        if ( curNode->isLeaf() )
        {
            // Visit the leaf
            visitor::visitLeaf(curNode, level, res );
        }
        else // this is a split-node
        {
            if ( qBox.second[curNode->axis] <= curNode->pos)
                // qBox overlaps only left child of this split-node
                stack.push(curNode->left);
            else if  ( qBox.first[curNode->axis] >= curNode->pos)
                // qBox overlaps only right child of this split-node
                stack.push(curNode->right);
            else
            {   
                // qBox overlaps both children of this split-node 
                stack.push(curNode->left );
                stack.push(curNode->right);
            }
        }
    }   

    return res;
}

template<unsigned d, class T> 
template<typename visitor>
typename visitor::return_type
gsHDomain<d,T>::leafSearch() const
{
    typename visitor::return_type i = visitor::init;
    std::stack<node*> stack;
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

template<unsigned d, class T> 
template<typename visitor> 
typename visitor::return_type
gsHDomain<d,T>::nodeSearch() const
{
    typename visitor::return_type i = visitor::init;
    std::stack<node*> stack;
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


/*
template<unsigned d, class T > int
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


template<unsigned d, class T>
void gsHDomain<d,T>::getBoxes(gsMatrix<unsigned>& b1, gsMatrix<unsigned>& b2, gsVector<unsigned>& level) const
{
    std::vector<std::vector<unsigned int> > boxes;
    getBoxes_vec(boxes);
    for(unsigned int i = 0; i < boxes.size(); i++){
        if ((boxes[i][0]==boxes[i][2]) || (boxes[i][1]==boxes[i][3])){
            boxes.erase(boxes.begin()+i);
            i--;
        }
    }

    connect_Boxes(boxes);
    b1.resize(boxes.size(),d);
    b2.resize(boxes.size(),d);
    level.resize(boxes.size());
    for(std::size_t i = 0; i < boxes.size(); i++){
        for(unsigned j = 0; j < d; j++){
            b1(i,j) = boxes[i][j];
            b2(i,j) = boxes[i][j+d];
        }
        level[i] = boxes[i][2*d];
    }

    /*ad hoc example- no oscillation fillet- keep for some time
     *int nboxes = 36;
    level.resize(nboxes);
    //gsMatrix<unsigned> b1, b2;
    b1.resize(nboxes,2);
    b2.resize(nboxes,2);

    level[0] = 1; b1(0,0) =  0; b1(0,1) =  24; b2(0,0) = 2; b2(0,1) = 28;
    level[1] = 1; b1(1,0) =  0; b1(1,1) =  10; b2(1,0) = 10; b2(1,1) = 24;
    level[2] = 1; b1(2,0) =  0; b1(2,1) =  6; b2(2,0) = 14; b2(2,1) = 10;
    level[3] = 1; b1(3,0) =  0; b1(3,1) =  0; b2(3,0) = 40; b2(3,1) = 6;
    level[4] = 1; b1(4,0) =  24; b1(4,1) =  6; b2(4,0) = 40; b2(4,1) = 14;
    level[5] = 1; b1(5,0) =  28; b1(5,1) =  14; b2(5,0) = 40; b2(5,1) = 20;
    level[6] = 1; b1(6,0) =  30; b1(6,1) =  20; b2(6,0) = 40; b2(6,1) = 22;
    level[7] = 1; b1(7,0) =  38; b1(7,1) =  22; b2(7,0) = 40; b2(7,1) = 26;
    //level 2
    level[8] = 2; b1(8,0) =  0; b1(8,1) =  35; b2(8,0) = 17; b2(8,1) = 40;
    level[9] = 2; b1(9,0) =  0; b1(9,1) =  33; b2(9,0) = 14; b2(9,1) = 35;
    level[10] = 2; b1(10,0) =  0; b1(10,1) =  28; b2(10,0) = 13; b2(10,1) = 33;
    level[11] = 2; b1(11,0) =  2; b1(11,1) =  25; b2(11,0) = 13; b2(11,1) = 28;
    level[12] = 2; b1(12,0) =  2; b1(12,1) =  24; b2(12,0) = 15; b2(12,1) = 25;
    level[13] = 2; b1(13,0) =  10; b1(13,1) =  22; b2(13,0) = 15; b2(13,1) = 24;
    level[14] = 2; b1(14,0) =  10; b1(14,1) =  20; b2(14,0) = 30; b2(14,1) = 22;
    level[15] = 2; b1(15,0) =  10; b1(15,1) =  14; b2(15,0) = 28; b2(15,1) = 20;
    level[16] = 2; b1(16,0) =  10; b1(16,1) =  10; b2(16,0) = 24; b2(16,1) = 14;
    level[17] = 2; b1(17,0) =  14; b1(17,1) =  6; b2(17,0) = 24; b2(17,1) = 10;
    level[18] = 2; b1(18,0) =  23; b1(18,1) =  22; b2(18,0) = 38; b2(18,1) = 26;
    level[19] = 2; b1(19,0) =  25; b1(19,1) =  26; b2(19,0) = 40; b2(19,1) = 30;
    level[20] = 2; b1(20,0) =  26; b1(20,1) =  30; b2(20,0) = 40; b2(20,1) = 35;
    level[21] = 2; b1(21,0) =  26; b1(21,1) =  35; b2(21,0) = 32; b2(21,1) = 36;
    level[22] = 2; b1(22,0) =  25; b1(22,1) =  36; b2(22,0) = 32; b2(22,1) = 37;
    level[23] = 2; b1(23,0) =  25; b1(23,1) =  37; b2(23,0) = 28; b2(23,1) = 39;
    level[24] = 2; b1(24,0) =  24; b1(24,1) =  39; b2(24,0) = 28; b2(24,1) = 40;
    level[25] = 2; b1(25,0) =  37; b1(25,1) =  35; b2(25,0) = 40; b2(25,1) = 40;
    //level 3
    level[26] = 3; b1(26,0) =  17; b1(26,1) =  39; b2(26,0) = 24; b2(26,1) = 40;
    level[27] = 3; b1(27,0) =  17; b1(27,1) =  36; b2(27,0) = 25; b2(27,1) = 39;
    level[28] = 3; b1(28,0) =  17; b1(28,1) =  35; b2(28,0) = 26; b2(28,1) = 36;
    level[29] = 3; b1(29,0) =  14; b1(29,1) =  33; b2(29,0) = 26; b2(29,1) = 35;
    level[30] = 3; b1(30,0) =  13; b1(30,1) =  30; b2(30,0) = 26; b2(30,1) = 33;
    level[31] = 3; b1(31,0) =  13; b1(31,1) =  26; b2(31,0) = 25; b2(31,1) = 30;
    level[32] = 3; b1(32,0) =  13; b1(32,1) =  25; b2(32,0) = 23; b2(32,1) = 30;
    level[33] = 3; b1(33,0) =  15; b1(33,1) =  22; b2(33,0) = 23; b2(33,1) = 25;
    level[34] = 3; b1(34,0) =  28; b1(34,1) =  37; b2(34,0) = 37; b2(34,1) = 40;
    level[35] = 3; b1(35,0) =  32; b1(35,1) =  35; b2(35,0) = 37; b2(35,1) = 37;*/
}
template<unsigned d, class T>
void gsHDomain<d,T>::getBoxesInLevelIndex(gsMatrix<unsigned>& b1,
              gsMatrix<unsigned>& b2,
              gsVector<unsigned>& level) const{
    std::vector<std::vector<unsigned int> > boxes;
    getBoxes_vec(boxes);
    for(unsigned int i = 0; i < boxes.size(); i++){
        if ((boxes[i][0]==boxes[i][2]) || (boxes[i][1]==boxes[i][3])){
            boxes.erase(boxes.begin()+i);
            i--;
        }
    }
    gsVector<unsigned,d>temp;
    gsVector<unsigned,d>temp1;
    connect_Boxes(boxes);
    b1.resize(boxes.size(),d);
    b2.resize(boxes.size(),d);
    level.resize(boxes.size());
    for(std::size_t i = 0; i < boxes.size(); i++){
        for(unsigned j = 0; j < d; j++){
//            b1(i,j) = boxes[i][j];
//            b2(i,j) = boxes[i][j+d];
            temp[j] = boxes[i][j];
            temp1[j] = boxes[i][j+d];
        }
        level[i] = boxes[i][2*d];
        global2localIndex( temp , level[i], temp );
        global2localIndex( temp1 , level[i], temp1 );
        b1.row(i) = temp;
        b2.row(i) = temp1;
    }
}


template<unsigned d, class T> void 
gsHDomain<d,T>::connect_Boxes(std::vector<std::vector<unsigned int> > &boxes) const
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
    //std::cout<<"in the connecting fucntion"<<boxes.size()<<std::endl;
}


template<unsigned d, class T>
void gsHDomain<d,T>::connect_Boxes_2(std::vector<std::vector<unsigned int> > &boxes) const
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




    //std::cout<<"in the connecting fucntion"<<boxes.size()<<std::endl;
}



template<unsigned d, class T> void
gsHDomain<d,T>::getBoxes_vec(std::vector<std::vector<unsigned int> >& boxes) const
{
    boxes.clear();

    std::stack<node*> stack;
    stack.push(m_root);    
    node * curNode;
    while ( ! stack.empty() )
    {
        curNode = stack.top();
        stack.pop();
        
        if ( curNode->isLeaf() )
        {
            const point & lower = curNode->lowCorner();
            const point & upper = curNode->uppCorner();

            boxes.push_back(std::vector<unsigned int>());
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


/////////////////////functions for returning the boudaries of domains////////////////////////
template<unsigned d, class T>
std::vector< std::vector<std::vector< std::vector< unsigned int > > > >
gsHDomain<d,T>::getPolylines() const
{
/*
 1: Get boxes from the quadtree.
 2: Return their borderlines.

 Return structure:
 < levels < polylines_in_one_level < one_polyline < one_segment (x1, y1, x2, y2) > > > > result
 note that <x1, y1, x2, y2 > are so that (x1, y1) <=LEX  (x2, y2)
*/
    std::vector<std::vector<unsigned int> > boxes;
    getBoxes_vec(boxes);//boxes are in the highest level indices, returns all leaves
    /*for(unsigned int i = 0; i < boxes.size(); i++){
        if ((boxes[i][0]==boxes[i][2]) || (boxes[i][1]==boxes[i][3])){
            boxes.erase(boxes.begin()+i);
            i--;
        }
    }*/

    // Get rid of boxes that are not of full dimension.
    for( std::vector< std::vector< unsigned int> >::iterator it = boxes.begin(); it != boxes.end(); ++it )
    {
        if( ( (*it)[0] == (*it)[2] ) || ( (*it)[0] == (*it)[2] ) )
            it = boxes.erase(it);
    }

    std::vector<std::vector<gsVSegment<T> > > seg;
    seg.resize(m_maxInsLevel+1);

    // For each level prepare the vertical lines separately
    for (unsigned int i = 0; i < boxes.size() ; i ++)
    {
        seg[boxes[i][4]].push_back(gsVSegment<unsigned int>(boxes[i][0],boxes[i][1],boxes[i][3], false) );
        seg[boxes[i][4]].push_back(gsVSegment<unsigned int>(boxes[i][2],boxes[i][1],boxes[i][3], false) );
    }

    // Process vertical lines from each level separately
    std::vector< std::vector<std::vector< std::vector<unsigned int > > > > result;
    for(unsigned int i = 0; i < m_maxInsLevel+1; i++)
    {
       //result.push_back(getPoly(seg[i]));
        result.push_back( getPolylinesSingleLevel( seg[i] ));
    }
    //result.push_back(getPolylinesSingleLevel(seg[0]));

    return result;
}



/*
template<unsigned d, class T>
int gsHDomain<d,T>::closes(gsSegment<T>& seg,
                            std::vector<gsSegment<T>  >& active_spans)const
{
    // Tells, whether the segment seg closes a polyline represented by an active span (the difference between y coordinates of its endpoints).
    // If so, returns its index.
    for(unsigned int i = 0; i < active_spans.size();i++){
        if(active_spans[i].m_y==seg.m_y){
            return i;
        }
    }

    return -1;
}

template<unsigned d, class T>
int gsHDomain<d,T>::appends(gsSegment<T>& seg,
                            std::vector<gsSegment<T>  >& active_spans, bool& front)const
{
  ///find index
    for(unsigned int i = 0; i < active_spans.size();i++){
        if( (seg.m_y.first==active_spans[i].m_y.second) || (seg.m_y.second==active_spans[i].m_y.second)){
            front = false;
            return i;
        }
        if( (seg.m_y.first==active_spans[i].m_y.first) || (seg.m_y.second==active_spans[i].m_y.first) ){
            front = true;
            return i;
        }
    }

    return -1;

}
*/





template<unsigned d, class T>
std::vector<std::vector< std::vector<unsigned int > > > gsHDomain<d,T>::getPolylinesSingleLevel(std::vector<gsVSegment<T> >& seg) const
{
    // For didactic purposes the interior of the function has been refined into two procedures
    // (overal length of the older version was intimidating and people would not like to read it then =)).

    // For each x coordinate there will be a list of vertical segments with this x coord.
    std::list< std::list< gsVSegment< T > > > vert_seg_lists;

    // For returning
    std::vector< std::vector< std::vector<unsigned int > > > result;

    // This magically sorts seg according to x value
    std::sort( seg.begin(), seg.end() );

    // Put stuff from seg into vert_seg_lists
    std::list< gsVSegment< unsigned int > > segs_x; // segments with the same particular x coord
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

template <unsigned d, class T>
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

template <unsigned d, class T>
void gsHDomain<d,T>::sweeplineConnectAndMerge( std::vector< std::vector< std::vector<unsigned int> > >& result, std::list< std::list< gsVSegment<T> > >& vert_seg_lists ) const
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

        gsAAPolyline<unsigned int> curr_poly;
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
                        result.push_back( curr_poly.writeParasolidUnsigned () );
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
                result.push_back( (*it).writeParasolidUnsigned() );
                act_poly.erase( it++ );
            }
            else
                ++it;

        }
    }
}
    // REMARK to remember: whenever 'discard qualifiers' error, it refers to something that isn't const or to volatile


template<unsigned d, class T> inline void
 gsHDomain<d,T>::computeFinestIndex( gsVector<unsigned,d> const & index,
                                     unsigned lvl,
                                     gsVector<unsigned,d> & result ) const
{
    for ( unsigned i = 0; i!=d; ++i )
        result[i] = index[i] << (m_maxInsLevel-lvl) ;
}

template<unsigned d, class T> inline void
 gsHDomain<d,T>::computeLevelIndex( gsVector<unsigned,d> const & index,
                                    unsigned lvl,
                                    gsVector<unsigned,d> & result ) const
{
    for ( unsigned i = 0; i!=d; ++i )
        result[i] = index[i] >> (m_maxInsLevel-lvl) ;
}

template<unsigned d, class T> inline void
 gsHDomain<d,T>::local2globalIndex( gsVector<unsigned,d> const & index,
                    unsigned lvl,
                    gsVector<unsigned,d> & result
                    ) const
{
    for ( unsigned i = 0; i!=d; ++i )
        result[i] = index[i] << (m_indexLevel-lvl) ;
}

template<unsigned d, class T> inline void
 gsHDomain<d,T>::global2localIndex( gsVector<unsigned,d> const & index,
                                        unsigned lvl,
                                        gsVector<unsigned,d> & result
    ) const
{
    for ( unsigned i = 0; i!=d; ++i )
        result[i] = index[i] >> (this->m_indexLevel-lvl) ;
}
}// end namespace gismo
