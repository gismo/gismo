
#pragma once

#include <gsHSplines/gsAABB.h>

namespace gismo {


template<short_t d, class Z>
class gsHTreeData
{
public:

    typedef          gsAABB<d, Z> kdBox;
    typedef typename kdBox::point point;

private:
    
    int m_level;
    kdBox box;

public:

    bool check() const
    { return (box.first.array() >= box.second.array()).any(); }
    
    static void split(const gsKdTree<d,Z,gsHTreeData> & node)
    {
        node.left->nodeData().upperCorner()[node.axis] = 
            node.right->nodeData().lowerCorner()[node.axis] = node.pos;
    }
    
    // Merges the data of right child into left child of the node
    //template<short_t d, class Z = index_t>
    static void mergeToLeft(const gsKdTree<d,Z,gsHTreeData> & node)
    {
        const kdBox * lbox = &node.left ->nodeData().aabb();
        const kdBox * rbox = &node.right->nodeData().aabb();
        lbox->second[node.axis] = rbox->second[node.axis];
        //rbox->first[node.axis] = lbox->first[axis];
    }

    void multiplyByTwo()
    {
        box.first .array() *= 2;
        box.second.array() *= 2;   
    }

    void divideByTwo()
    {
        box.first .array() /= 2;
        box.second.array() /= 2;   
    }

    static void nextMidSplit(gsKdTree<d,Z,gsHTreeData> &  node)
    {
        node.axis = ( node.parent == 0 ? 0 : (node.parent->axis+1)%d );
        node.pos  = node.box->first[node.axis] +
            (node.box->second[node.axis] - node.box->first[node.axis])/2 ;
    }

    static void anyMidSplit(unsigned h, gsKdTree<d,Z,gsHTreeData> &  node,
                            int & doSplit)
    {
        const unsigned mask = ~(h - 1);
        for ( unsigned i = 0; i < d; ++i )
        {
            const unsigned c = 
                (node.box.first [i] + (node.box.second[i] - node.box.first[i])/2) & mask ;
            if ( c != node.box.first [i] ) // avoid degenerate split
            {
                node.axis = i;
                node.pos  = c;
                doSplit   = 1;
                return;
            }
        }
        doSplit = 0;
    }

    static void adaptiveAlignedSplit(const gsHTreeData & insData, unsigned h,
                                     gsKdTree<d,Z,gsHTreeData> &  node,
                                     int & doSplit)
    {
        for (short_t i = 0; i < d; ++i)
        {
            const Z c1 = insData.lowerCorner()[i] - insData.lowerCorner()[i] % h; //floor
            const Z cc = insData.upperCorner()[i] % h;
            const Z c2 = insData.upperCorner()[i] + (cc ? h-cc : 0 ); // ceil

            if ( c1 > node.nodeData().lowerCorner()[i] )
            {
                // right child intersects insBox
                node.axis = i;
                node.pos = c1;
                doSplit  = 1;
                return;
            }
            else if ( c2 <node.nodeData().upperCorner()[i]  )
            {
                // left child intersects insBox
                node.axis = i;
                node.pos = c2;
                doSplit  = -1;
                return;
            }
        }
        doSplit  = 0;

    }

    static void adaptiveSplit(gsHTreeData & insData,
                              const gsKdTree<d,Z,gsHTreeData> &  node,
                              int & doSplit)
    {
        for (short_t i = 0; i < d; ++i)
        {
            if ( insData.box->first[i] > node.data->box->first[i] )
            {
                // right child intersects insBox
                node.axis = i;
                node.pos = insData.box->first[i];
                doSplit  = 1;
                return;
            }
            else if ( insData.box->second[i] <node.data->box->second[i]  )
            {
                // left child intersects insBox
                node.axis = i;
                node.pos = insData.box->second[i];
                doSplit  = -1;
                return;
            }
        }
        doSplit  = 0;
    }

    int level() const {return m_level;}
    int & level() {return m_level;}

    const point & lowerCorner() const { return box.first; }
          point & lowerCorner()       { return box.first; }

    const point & upperCorner() const { return box.second; }
          point & upperCorner()       { return box.second; }

    const kdBox & aabb() const { return box; }
    kdBox & aabb() { return box; }

    bool isAligned(unsigned index_level) const
    {
        const unsigned h = 1 << (index_level - m_level);
        
        for ( index_t i = 0; i!=box.first.size(); ++i )
        {
            if (box.second[i] % h != 0 ||
                box.first[i]  % h != 0 )
                return false;
        }
        return true;
    }

};


/// CRTP ??

}

