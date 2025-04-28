
//template<typename T>
class gsHTreeData
{
    typedef          gsAABB<d, Z> kdBox;
    typedef typename kdBox::point point;

    int level;
    kdBox box;

public:

    static void split(const gsKdNode2<d,Z,gsHtreeData> & node)
    {
        node.left->box->second[node.axis] = 
        node.right->box->first [node.axis] = node.pos;
    }
    
    // Merges the data of right child into left child of the node
    //template<short_t d, class Z = index_t>
    static void mergeToLeft(const gsKdNode2<d,Z,gsHtreeData> & node)
    {
        kdBox * lbox = node.left->data->box;
        kdBox * rbox = node.right->data->box;
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

    static void nextMidSplit(gsKdNode2<d,Z,gsHtreeData> &  node)
    {
        node.axis = ( node.parent == 0 ? 0 : (node.parent->axis+1)%d );
        node.pos  = node.box->first[node.axis] +  (node.box->second[axis] - node.box->first[axis])/2 ;
    }

    static void anyMidSplit(unsigned h, gsKdNode2<d,Z,gsHtreeData> &  node)
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

    static void adaptiveAlignedSplit(gsHtreeData & insData, unsigned h,
                                     const gsKdNode2<d,Z,gsHtreeData> &  node,
                                     int & doSplit)
    {
        for (short_t i = 0; i < d; ++i)
        {
            const Z c1 = insBox. first[i] - insBox. first[i] % h; //floor
            const Z cc = insBox.second[i] % h;
            const Z c2 = insBox.second[i] + (cc ? h-cc : 0 ); // ceil

            if ( c1 > node.data->box->first[i] )
            {
                // right child intersects insBox
                node.axis = i;
                node.pos = c1;
                doSplit  = 1;
                return;
            }
            else if ( c2 <node.data->box->second[i]  )
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

    static void adaptiveSplit(gsHtreeData & insData,
                              const gsKdNode2<d,Z,gsHtreeData> &  node,
                              int & doSplit)
    {
        for (short_t i = 0; i < d; ++i)
        {
            if ( insBox.box->first[i] > node.data->box->first[i] )
            {
                // right child intersects insBox
                node.axis = i;
                node.pos = insBox.box->first[i];
                doSplit  = 1;
                return;
            }
            else if ( insBox.box->second[i] <node.data->box->second[i]  )
            {
                // left child intersects insBox
                node.axis = i;
                node.pos = insBox.box->second[i];
                doSplit  = -1;
                return;
            }
        }
        doSplit  = 0;
    }

    point lowerCorner() const
    {
        point result = curNode->box->first;
        const int lvl = curNode->level;

        //result = result.array() / (1>> (m_index_level-lvl)) ;
        for ( index_t i = 0; i!= result.size(); ++i )
            result[i] = result[i] >> (m_index_level-lvl) ;

        return result;
    }

    point upperCorner() const
    {
        point result = curNode->box->second;
        const int lvl = curNode->level;

        for ( index_t i = 0; i!=result.size(); ++i )
            result[i] = result[i] >> (m_index_level-lvl) ;

        return result;
    }

    bool isAligned() const
    {
        const unsigned h = 1 << (m_index_level - curNode->level);
        
        for ( index_t i = 0; i!=curNode->box->first.size(); ++i )
        {
            if (curNode->box->second[i] % h != 0 ||
                curNode->box->first[i]  % h != 0 )
                return false;
        }
        return true;
    }

};

/// CRTP ??
