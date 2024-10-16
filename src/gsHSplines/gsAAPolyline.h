/** @file gsAAPolyline.h

    @brief Class for handling an axis-aligned polyline in 2D.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): D. Mokris, G. Kiss
*/

#pragma once

#include <gsHSplines/gsVSegment.h>

namespace gismo
{

template<class T>
class gsAAPolyline
{
    /** \brief Axis-aligned polyline in the plane.  Represented in
     * terms of vertices of vertical segments forming it. Attempts to
     * get rid of redundant vertices where reasonable.
     */

public:

    gsAAPolyline()
	{
		m_closed = false;
	}

    gsAAPolyline(gsVSegment<T> VSeg)
    {
        std::vector<T> newPoint(2);
        newPoint[0] = VSeg.getX();
        newPoint[1] = VSeg.getYDown();
        m_vertices.push_back( newPoint );
        newPoint[1] = VSeg.getYUp();
        m_vertices.push_back( newPoint );
        m_closed = false;
    }

    /// Returns the polyline in a format suitable for further treatment by PARASOLID
    std::vector< std::vector< index_t > > writeParasolid();

    /// Returns the polyline in a format suitable for further treatment by PARASOLID
    std::vector< std::vector<unsigned int > > writeParasolidUnsigned();

    /// Tries to add the vertical segment (\a x, \a y0, \a x, \a y1) to polyline's end or beginning.
    inline bool addVerticalSegment( int x, int y0, int y1 )
	{
		return ((pushBack(x, y0, y1)) ||
				(pushBack(x, y1, y0)) ||
				(pushFront(x, y0, y1)) ||
				(pushFront(x, y1, y0)));
	}

    /// Tries to add the vertical segment \a vert_seg to polyline's end or beginning.
    inline bool canBeExtended( gsVSegment<T> vert_seg )
    {
        return addVerticalSegment( vert_seg.getX(), vert_seg.getYDown(), vert_seg.getYUp() );
    }

    /// Could be made closed by adding a horizontal segment?
	/// If the last and first segment could be merged, it tries so.
    bool almostClosed()
    {
        if ( frontY() == backY() )
        {
            // It could have happened that the vertices forming the last and first line are colinear.
            if( frontX() == backX() ) // WARNING: Assumes the last end first line to be vertical.
            {
                m_vertices.pop_back();
                m_vertices.front()[1] = backY();
                m_vertices.pop_back();
            }
            return true;
        }
        return false;
    }

    /// Tries to merge with \a other_poly
    bool mergeWith( gsAAPolyline<T>& other_poly );

    inline std::list< std::vector<T> > getVertices()
	{
		return m_vertices;
	}

    inline T frontX()
	{
		return m_vertices.front()[0];
	}

    inline T backX()
    {
		return m_vertices.back()[0];
	}

    inline T frontY()
	{
		return m_vertices.front()[1];
	}

    inline T backY()
    {
		return m_vertices.back()[1];
	}

protected:

	/// Vertices of of the polyline, each represented as a vector of two coordinates.
    std::list< std::vector<T> > m_vertices;

	/// Whether the polyline is a closed curve.
    bool m_closed;

    /// Tries to add vertical segment (x,y0,x,y1) at the end of the polyline
    bool pushBack( T x, T y0, T y1 );

    /// Tries to insert vertical segment (x,y0,x,y1) to the beginning of the polyline.
    bool pushFront( T x, T y0, T y1 );
};

template <class T>
bool gsAAPolyline<T>::pushBack( T x, T y0, T y1 )
{
    if( m_vertices.empty() || (y0 == backY() ))
    {
        if( m_vertices.empty() || x != backX() ) // new segment
        {
            std::vector< T > newVertex (2);

            newVertex[0] = x;
            newVertex[1] = y0;
            m_vertices.push_back( newVertex );

            newVertex[1] = y1; // x coord is the same
            m_vertices.push_back( newVertex );
        }
        else // prolongation of the last segment
            m_vertices.back()[1] = y1;

        m_closed = almostClosed(); // This also looks, if the last and first segment can be merged.
        return true;
    }

    else // appending would not be axis aligned
        return false;
}

template <class T>
bool gsAAPolyline<T>::pushFront( T x, T y0, T y1 ) // Symmetric case to push_back
{
    if( m_vertices.empty() || y0 == frontY() )
    {
        if( m_vertices.empty() || x != frontX() ) // Adding a segment with a different x coordinate
        {
            std::vector< T > newVertex (2);

            newVertex[0] = x;
            newVertex[1] = y0;
            m_vertices.push_front( newVertex );

            newVertex[1] = y1;
            m_vertices.push_front( newVertex );
        }
        else
            m_vertices.front()[1] = y1; // Adding a colinear point is just a prolongation.

        m_closed = almostClosed();
        return true;
    }
    else
        return false;
}

template<class T>
std::vector< std::vector< index_t > > gsAAPolyline<T>::writeParasolid()
{
    std::vector< std::vector< index_t > > result;

    if( m_vertices.size() < 4 )
    {
        gsWarn << "Function write_parasolid() says: only " << m_vertices.size() << " vertices.\n";
        return result;
    }

    if( !almostClosed() )
    {
        gsWarn << "Function write_parasolid() says: your curve should be closed but it isn't.";
        return result;
    }

    std::vector< index_t > item(4);

    // it will go through m_vertices; it_fwd will be pointing always to it's next % m_vertices.size().
    auto it_fwd = m_vertices.begin();
    ++it_fwd;
    for(auto it = m_vertices.begin(); it != m_vertices.end(); ++it, ++it_fwd)
    {
        if( it_fwd == m_vertices.end() ) // cyclic list
            it_fwd = m_vertices.begin();

        // The lines are expected oriented from left to right, from bottom up.
        item[0] = std::min( (*it)[0], (*it_fwd)[0] );
        item[1] = std::min( (*it)[1], (*it_fwd)[1] );
        item[2] = std::max( (*it)[0], (*it_fwd)[0] );
        item[3] = std::max( (*it)[1], (*it_fwd)[1] );

        result.push_back( item );
    }
    return result;
}


template<class T>
std::vector< std::vector< unsigned int > > gsAAPolyline<T>::writeParasolidUnsigned()
{
    std::vector< std::vector<unsigned  int > > result;

    if( m_vertices.size() < 4 )
    {
        gsWarn << "write_parasolid() says: only " << m_vertices.size() << " vertices.\n";
        return result;
    }

    if( !almostClosed() )
    {
        gsWarn << "write_parasolid() says: your curve should be closed but it isn't.";
        return result;
    }

    std::vector< unsigned int > item(4,0);

    // it will go through m_vertices; it_fwd will be pointing always to it's next % m_vertices.size().
    auto it_fwd = m_vertices.begin();
    ++it_fwd;
    for(auto it = m_vertices.begin(); it != m_vertices.end(); ++it, ++it_fwd)
    {
        if( it_fwd == m_vertices.end() ) // cyclic list
            it_fwd = m_vertices.begin();

        // The lines are expected oriented from left to right, from bottom up.
        item[0] = std::min( (*it)[0], (*it_fwd)[0] );
        item[1] = std::min( (*it)[1], (*it_fwd)[1] );
        item[2] = std::max( (*it)[0], (*it_fwd)[0] );
        item[3] = std::max( (*it)[1], (*it_fwd)[1] );

        result.push_back( item );
    }
    return result;
}

template <class T>
bool gsAAPolyline<T>::mergeWith( gsAAPolyline<T>& otherPoly )
{
    // We assume both polylines to be ending with a vertical segment (!).
    std::list< std::vector<T> > newVertices;

    // adding in front of m_vertices
    if( otherPoly.backY() == frontY() )
    {
        newVertices = otherPoly.getVertices();
        if( newVertices.back() == m_vertices.front() )
        {
            newVertices.pop_back();
            m_vertices.pop_front();
        }
        m_vertices.splice( m_vertices.begin(), newVertices );
        m_closed = almostClosed();
        return true;
    }
    // adding behind vertices in reverse order
    else if( otherPoly.backY() == backY() )
    {
        newVertices = otherPoly.getVertices();
        newVertices.reverse();
        if( newVertices.front() == m_vertices.back() )
        {
            newVertices.pop_front();
            m_vertices.pop_back();
        }
        m_vertices.splice( m_vertices.end(), newVertices );
        m_closed = almostClosed();
        return true;
    }
    // adding in front of m_vertices in reverse order
    else if( otherPoly.frontY() == frontY() )
    {
        newVertices = otherPoly.getVertices();
        newVertices.reverse();
        if( newVertices.back() == m_vertices.front() )
        {
            newVertices.pop_back();
            m_vertices.pop_front();
        }
        m_vertices.splice( m_vertices.begin(), newVertices );
        m_closed = almostClosed();
        return true;
    }
    // adding behind m_vertices
    else if( otherPoly.frontY() == backY() )
    {
        newVertices = otherPoly.getVertices();
        if( newVertices.front() == m_vertices.back() )
        {
            newVertices.pop_front();
            m_vertices.pop_back();
        }
        m_vertices.splice( m_vertices.end(), newVertices );
        m_closed = almostClosed();
        return true;
    }
    // none of the previous was possible
    else
        return false;
}


} // namespace gismo
