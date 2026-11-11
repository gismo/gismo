/** @file gsHTree.hpp

    @brief Provides the implementation of the class gsHTree.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): G. Kiss, A. Mantzaflaris
*/

#pragma once

namespace gismo {

template<short_t d, class Z>
gsHTree<d,Z> gsHTree<d,Z>::merge(const gsHTree<d,Z>& tree1, const gsHTree<d,Z>& tree2)
{
    GISMO_ASSERT( (tree1.m_upperIndex.array() == tree2.m_upperIndex.array()).all(),
        "gsHTree::merge: trees must have the same domain." );
    GISMO_ASSERT( tree1.m_indexLevel == tree2.m_indexLevel,
        "gsHTree::merge: trees must have the same index level." );

    gsHTree<d,Z> result(tree1);

    for (const_literator it = tree2.beginLeafIterator(); it.good(); it.next())
        if (it.level() > 0)
            result.insertBox(it.lowerCorner(), it.upperCorner(), it.level());

    result.makeCompressed();
    return result;
}

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
    stack.push(this);
    const node * curNode;
    point lower, upper;
    while ( ! stack.empty() )
    {
        curNode = stack.top();
        stack.pop();

        if ( curNode->axis == -2 ) // empty node: no leaves
            continue;

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

template<short_t d, class Z>
void gsHTree<d, Z>::getBoxesOnSide(boundary::side s, gsMatrix<Z>& b1, gsMatrix<Z>& b2, gsVector<Z>& level) const
{

    getBoxes( b1, b2, level);
    std::vector<Z> onSide;

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
        for(index_t i = 0; i < b1.rows(); i++)
            if( b1(i, quotient ) == 0 )
                onSide.push_back(i);
    }
    else // remainder == 1, sides east, north/up, or back
    {
        for(index_t i = 0; i < b1.rows(); i++)
        {
            // index of upper corner
            Z B2( b2(i, quotient ) );
            // transform to index-level
            B2 = B2 << (m_indexLevel - m_maxInsLevel);

            // check if this index is on the boundary
            if( B2 == m_upperIndex[ quotient ] )
                onSide.push_back( i );
        }
    }

    // select only the boxes on side s:
    for(size_t i=0; i < onSide.size(); i++)
    {
        b1.row(i) = b1.row( onSide[i] );
        b2.row(i) = b2.row( onSide[i] );
        level[i]  = level( onSide[i] );
    }
    b1.conservativeResize( onSide.size(), b1.cols() );
    b2.conservativeResize( onSide.size(), b2.cols() );
    level.conservativeResize( onSide.size() );
}

template<short_t d, class Z>
std::vector<std::vector<std::vector< std::vector<Z> > > >
gsHTree<d, Z>::getPolylines() const
{
/*
 1: Get boxes from the quadtree.
 2: Return their borderlines.

 Return structure:
 < levels < polylines_in_one_level < one_polyline < one_segment (x1, y1, x2, y2) > > > > result
 note that <x1, y1, x2, y2 > are so that (x1, y1) <=LEX  (x2, y2)
*/
    std::vector<std::vector<Z> > boxes;
    getBoxes_vec(boxes);// Returns all leaves.

    // Get rid of boxes that are not of full dimension.
    for(auto it = boxes.begin(); it != boxes.end(); ++it)
    {
        if( ( (*it)[0] == (*it)[2] ) || ( (*it)[0] == (*it)[2] ) )
            it = boxes.erase(it);
    }

    std::vector<std::vector<gsVSegment<Z> > > seg;
    seg.resize(m_maxInsLevel+1);

    // For each level prepare the vertical lines separately.
    for(size_t i = 0; i < boxes.size() ; i ++)
    {
        seg[boxes[i][4]].push_back(gsVSegment<Z>(boxes[i][0], boxes[i][1], boxes[i][3], false));
        seg[boxes[i][4]].push_back(gsVSegment<Z>(boxes[i][2], boxes[i][1], boxes[i][3], false));
    }

    // Process vertical lines from each level separately.
    std::vector< std::vector<std::vector< std::vector<Z> > > > result;
    for(unsigned int i = 0; i < m_maxInsLevel+1; i++)
    {
        result.push_back( getPolylinesSingleLevel( seg[i] ));
    }

    return result;
}


template<short_t d, class Z>
std::vector<std::vector< std::vector<Z> > > gsHTree<d, Z>::getPolylinesSingleLevel(std::vector<gsVSegment<Z> >& seg) const
{

    // For each x coordinate there will be a list of vertical segments with this x coord.
    std::list< std::list< gsVSegment<Z> > > vert_seg_lists;

    // For returning
    std::vector< std::vector< std::vector<Z> > > result;

    // This magically sorts seg according to x value
    std::sort( seg.begin(), seg.end() );

    // Put stuff from seg into vert_seg_lists
    std::list< gsVSegment<Z> > segs_x; // segments with the same particular x coord
    for(auto it_seg = seg.begin(); it_seg != seg.end(); ++it_seg )
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

template <short_t d, class Z>
void gsHTree<d, Z>::getRidOfOverlaps( std::list<std::list<gsVSegment<Z> > >& vert_seg_lists ) const
{
    bool need_to_erase = false;
    for(auto it_x = vert_seg_lists.begin(); it_x != vert_seg_lists.end(); )
    {
        // For each x coordinate substract overlapping segments.
        for(auto it_slow = (*it_x).begin(); it_slow != (*it_x).end(); ++it_slow)
        {
            for(auto it_quick = it_slow; it_quick != (*it_x).end(); ++it_quick)
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
        for(auto it_slow = (*it_x).begin(); it_slow != (*it_x).end(); )
        {
            if( (*it_slow).length() == 0 )
            {
                (*it_x).erase( it_slow++ );
            }
            else
                ++it_slow; // Here instead of the internal for.
        }

        // There might be some empty lists;
        if( (*it_x).empty() )
            vert_seg_lists.erase( it_x++ );
        else
            ++it_x; // Here instead of the outer for.
    }
}

template <short_t d, class Z>
void gsHTree<d, Z>::sweeplineConnectAndMerge( std::vector< std::vector< std::vector<Z> > >& result,
                                               std::list< std::list< gsVSegment<Z> > >& vert_seg_lists ) const
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
    std::list< gsAAPolyline<Z> > act_poly;
    std::queue< gsAAPolyline<Z> > poly_queue;

    for(auto it_x = vert_seg_lists.begin(); it_x != vert_seg_lists.end(); ++it_x )
    {
        // Try to extend polylines by adding segments to them.
        // The funny queue procedure is a trick to make sure that none of the polylines
        // would be extended too much (i.e., over a vertex, where it should be closed).
        for(auto it = act_poly.begin(); it != act_poly.end(); ++it )
            poly_queue.push( *it );

        act_poly.erase( act_poly.begin(), act_poly.end() );

        gsAAPolyline<Z> curr_poly;
        while( !poly_queue.empty() )
        {
            curr_poly = poly_queue.front();
            poly_queue.pop();
            bool is_active = true;
            for(auto it_seg = (*it_x).begin(); it_seg != (*it_x).end(); )
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
        for(auto it = (*it_x).begin(); it != (*it_x).end(); ++it )
        {
            gsAAPolyline<Z> new_polyline( *it );
            act_poly.push_back( new_polyline );
        }

        // Merge what's possible.
        // Trick: go through act_poly from left to right. Check each of them with all to the right of it in act_poly.
        // If they can be merged, move it to the end. Those more to the left are already finalised.
        for(auto it_slow = act_poly.begin(); it_slow != act_poly.end(); ++it_slow )
        {
            for(auto it_quick = it_slow; it_quick != act_poly.end(); ++it_quick)
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
        for(auto it = act_poly.begin(); it != act_poly.end(); )
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


} // namespace gismo
