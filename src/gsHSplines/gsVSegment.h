/** @file gsVSegment.h

    @brief Helper class for gsAAPolyline.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Authors: D. Mokris, G. Kiss
*/

#pragma once

namespace gismo
{
	
/** \brief
 *  Class for representing a vertical line segment in 2D. Helper for the class gsAAPolyline.
 *
 *  \ingroup HSplines
*/

template<class T>
class gsVSegment
{
public:
    /// Default empty constructor.
    //gsVSegment() { }

    /// Constructor with a line segment.
    gsVSegment(int x, int bottomY, int upperY, bool l)
    {
        m_x = x;
        m_y = std::pair<int, int> (bottomY, upperY);
        left = l;
    }

    /// Default destructor.
    //~gsVSegment() { }

    /// Test the equality with the segment \em other.
    /// \param other
    bool isClosing(const gsVSegment& other) const;

    /// Test conection to the segment \em other.
    /// \param other
    bool isConnected(const gsVSegment& other) const;

    T length() const;

    /// Returns true if the segments (i.e., this and \param other) are disjoint;
	/// otherwise gets rid of their overlap and returns false.
    bool cannotDeleteOverlap(gsVSegment<T> & other);

	// getters

    inline T getX() const
	{
		return m_x;
	}

    inline T getYUp() const
	{
		return m_y.second;
	}

    inline T getYDown() const
	{
		return m_y.first;
	}

    inline std::pair<T,T> getY() const
	{
		return m_y;
	}

	// setters

    inline void setX( const T& new_x )
	{
		m_x = new_x;
	}

    inline void setYUp( const T& new_y_up )
	{
		m_y.second = new_y_up;
	}

    inline void setYDown( const T& new_y_down )
	{
		m_y.first = new_y_down;
	}

    inline void setY( const std::pair<T,T>& new_m_y )
	{
		m_y = new_m_y;
	}

    bool operator< (const gsVSegment<T>& other) const
	{
        if(m_x == other.getX())
        {
            return m_y.first < other.getYDown();
        }
		else
        {
            return m_x < other.getX();
        }
    }

    bool operator== ( const gsVSegment<T>& other ) const
	{
        if( (m_x == other.getX() ) && (m_y.first == other.getYDown() ) && (m_y.second == other.getYUp() ) )
            return true;
        else
            return false;
    }

protected:

	/// The x-coordinate.
    T m_x;

	/// Lower and upper y--coordinates.
    std::pair<T, T> m_y;

	/// Whether the segment is on the left hand-side of the enclosed domain.
    bool left;
};

template<class T>
bool gsVSegment<T>::isClosing(const gsVSegment& other ) const
{
    return ( ( other.getYDown() == m_y.first) && (( other.getYUp() == m_y.second)) );
}

template<class T>
bool gsVSegment<T>::isConnected(const gsVSegment& other ) const
{
    return ( ( other.getYUp() == m_y.first) || (( other.getYDown() == m_y.second)) );
}

template<class T>
T gsVSegment<T>::length() const
{
    return (m_y.second - m_y.first);
}

template<class T>
bool gsVSegment<T>::cannotDeleteOverlap(gsVSegment<T> & other)
{
    gsVSegment<T> tempa = *this;
    gsVSegment<T> tempb = other;
    std::vector<T> v(4);

    v[0] = m_y.first;
    v[1] = m_y.second;
    v[2] = other.m_y.first;
    v[3] = other.m_y.second;

    // This is the key trick:
    std::sort(v.begin(), v.end());

    m_y.first = v[0];
    m_y.second = v[1];
    other.m_y.first = v[2];
    other.m_y.second = v[3];

    // iff nothing changed (or just swapped), return true
    if( ((tempa==*this)&&(tempb==other)) || ((tempa==other)&&(tempb==*this)))
        return true;
    else
        return false;
}

/// Print (as string) operator
template<class T>
std::ostream &operator<<(std::ostream &os, const gsVSegment<T> a)
{
    os<<"x:"<<a.getX() <<" , y: ["<< a.getYDown() <<" ,"<< a.getYUp() <<"] ";
    return os;
}

} // namespace gismo
