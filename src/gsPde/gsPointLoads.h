/** @file gsPointLoads.h

    @brief Provides a simple container for point loads on multi-patch domains

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

namespace gismo
{



/** @brief 
    Struct defining a point together with a scalar or vector load.
    
    \ingroup Pde
*/
template<class T>
struct point_load
{
    point_load(const gsVector<T> & _point, 
               const T             _value, 
               int _patch = 0, 
               bool _parametric = true)
    :
    patch(_patch), value(_value), point(1), parametric(_parametric)
    { 
        point[0] = _value;
    }

    point_load(const gsVector<T> & _point, 
               const gsVector<T> & _value, 
               int _patch = 0, 
               bool _parametric = true)
    :
    patch(_patch), value(_value), point(_point), parametric(_parametric)
    { }

    int patch;

    gsVector<T> value;

    gsVector<T> point;

    bool parametric;
};



/** @brief Class containing a set of points on a multi-patch
    isogeometric domain, together with boundary conditions.
    
    \ingroup Pde
*/
template<class T>
class gsPointLoads
{
public:
    typedef point_load<T> pLoad;
    typedef typename std::vector<pLoad> plContainer;

    typedef typename std::vector<pLoad>::iterator iterator;

    typedef typename std::vector<pLoad>::const_iterator const_iterator;

public:

    /// Prints the object as a string.
    std::ostream & print(std::ostream &os) const
    { 
        os << "gsPointLoads: "<<m_pointLoads.size()<<"\n";
        return os; 
    }

public:

    /// Default empty constructor
    gsPointLoads() 
    { }

    ~gsPointLoads() // Destructor
    { }
    
public:

    void clear()
    {
        m_pointLoads.clear();
    }

    inline pLoad   operator [] (size_t i) const { return m_pointLoads[i]; }
    inline pLoad & operator [] (size_t i) { return m_pointLoads[i]; }

    void addLoad(const gsVector<T> & _point, 
                 const gsVector<T> & _value, 
                 int _patch = 0, 
                 bool _parametric = true)
    {
        m_pointLoads.push_back( pLoad(_point,_value,_patch,_parametric) );
    }

    void addLoad(const gsVector<T> & _point, 
                 const T             _value, 
                 int _patch = 0, 
                 bool _parametric = true)
    {
        m_pointLoads.push_back( pLoad(_point,_value,_patch,_parametric) );
    }

    size_t numLoads() const { return  m_pointLoads.size(); }

private:

    plContainer  m_pointLoads; ///< List of Point loads

}; // class gsPointLoads

/// Print (as string)
template<class T>
std::ostream &operator<<(std::ostream &os, const gsPointLoads<T>& pls)
{return pls.print(os); }

} // namespace gismo

