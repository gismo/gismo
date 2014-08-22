
#pragma once

#include <ostream>

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsGeometry.h>

namespace gismo
{

/** \brief
    A truncated hierarchical B-Spline function, in \em d dimensions.

    This is the geometry type associated with gsTHBSplineBasis.

    R^d -> R

    \tparam T is the coefficient type

    \ingroup geometry
*/

template<class T>
class gsSegment
{
private:
T m_x;
std::pair<T, T> m_y;
bool left;
public:
    /// Default empty constructor
    gsSegment() { };

    /// constructor with a line segment
    gsSegment(int x, int b_y, int u_y, bool l)
    {
        m_x = x;
        m_y = std::pair<int, int> (b_y, u_y);
        left = l;
    };

    ///default destructor
    ~gsSegment() { }; //destructor

    /// test the equality with segment a
    bool is_closing(const gsSegment& other)const;

    /// test conection to segment a
    bool is_connected(const gsSegment& other)const;
    T length()const;

    /// returns true, if the segments are disjoint; otherwise gets rid of their overlap and returns false.
    bool cannotDeleteOverlap(gsSegment<T> & other);

    inline T get_m_x() const { return m_x; }
    inline T get_m_y_up() const { return m_y.second; }
    inline T get_m_y_down() const { return m_y.first;}
    inline std::pair<T,T> get_m_y() const { return m_y; }

    inline void set_m_x( const T& new_x ){ m_x = new_x; }
    inline void set_m_y_up( const T& new_y_up ){ m_y.second = new_y_up; }
    inline void set_m_y_down( const T& new_y_down ){ m_y.first = new_y_down; }
    inline void set_m_y( const std::pair<T,T>& new_m_y ){ m_y = new_m_y; }

    bool operator< ( const gsSegment<T>& other )const{
        if(m_x == other.get_m_x() )
        {
            return m_y.first < other.get_m_y_down() ;
        }else
        {
            return m_x < other.get_m_x();
        }
    }

    bool operator== ( const gsSegment<T>& other )const{
        if( (m_x == other.get_m_x() ) && (m_y.first == other.get_m_y_down() ) && (m_y.second == other.get_m_y_up() ) )
            return true;
        else
            return false;
    }

};

template<class T>
bool gsSegment<T>::is_closing(const gsSegment& other )const{
    if ( ( other.get_m_y_down() == m_y.first) && (( other.get_m_y_up() == m_y.second)) )
        return true;
    else
        return false;
}

template<class T>
bool gsSegment<T>::is_connected(const gsSegment& other )const{
    if ( ( other.get_m_y_up() == m_y.first) || (( other.get_m_y_down() == m_y.second)) )
        return true;
    else
        return false;
}
template<class T>
T gsSegment<T>::length()const{
    return (m_y.second - m_y.first);
}

template<class T>
bool gsSegment<T>::cannotDeleteOverlap(gsSegment<T> & other)
{
    gsSegment<T> tempa = *this;
    gsSegment<T> tempb = other;
    std::vector<T> v(4);
    //v.resize(4);

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
std::ostream &operator<<(std::ostream &os, const gsSegment<T> a)
{
    os<<"x:"<<a.get_m_x() <<" , y: ["<< a.get_m_y_down() <<" ,"<< a.get_m_y_up() <<"] ";
    return os;
}
}; // namespace gismo


//#ifndef GISMO_HEADERS_ONLY
//#include GISMO_HPP_HEADER(gsTHBSpline.hpp)
//#endif
