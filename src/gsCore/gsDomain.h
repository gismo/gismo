/** @file gsDomain.h

    @brief Abstracgt Base class representing a domain. i.e. a
    collection of elements (triangles, rectangles, cubes, simplices.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>

namespace gismo
{

/** 
    @brief Class representing a domain. i.e. a collection of
    elements (triangles, rectangles, cubes, simplices.

    \warning  This interface is under development and is not used yet...
    
    \ingroup Core
*/
template<class T>
class gsDomain
{
public:

    virtual ~gsDomain() { }

#if EIGEN_HAS_RVALUE_REFERENCES && EIGEN_GNUC_AT_MOST(4,7) && !EIGEN_COMP_PGI
    // defaulted declaration required at least in Gcc 4.7.2
    gsDomain() = default;
    gsDomain(const gsDomain&) = default;
    gsDomain(gsDomain&&) = default;
    gsDomain & operator=(const gsDomain&) = default;
    gsDomain & operator=(gsDomain&&) = default;
#endif

public:

    /// dimension of the domain
    virtual short_t dim() const
    { gsWarn << "gsDomain: dimension() not defined at "<< *this << "\n"; return 0; }

    /// Returns a bounding box for the domain
    /// eg. This coincides to the domain in case of tensor-product domains
    virtual gsMatrix<T> boundingBox()
    {GISMO_NO_IMPLEMENTATION}

    /// Returns a list of elements
    virtual gsMatrix<T> elements()
    {GISMO_NO_IMPLEMENTATION}

    /// Returns the mesh..
    virtual gsMatrix<T> mesh()
    {GISMO_NO_IMPLEMENTATION }

    virtual T minMeshSize()
    {GISMO_NO_IMPLEMENTATION}

    /// Returns the breaks..
    virtual std::vector<T> breaks() const
    {GISMO_NO_IMPLEMENTATION}
    
    /// Clone function. Used to make a copy of the (derived) geometry
    virtual gsDomain* clone() const = 0;

    /// Prints the object as a string.
    virtual std::ostream &print(std::ostream &os) const = 0;

    /*
      Member functions that may be implemented or not in the derived class
    */

    virtual std::vector<T> unique() const 
    { gsWarn<<"gsDomain: unique() was not defined at"<< *this <<"\n"; return std::vector<T>() ;}

    void merge(gsDomain<T> * other ) 
    { gsWarn<<"gsDomain: merge(..) was not defined in "<< *this <<"\n"; return;}

}; // class gsDomain

/// Print (as string) operator to be used by all derived classes
template<class T>
std::ostream &operator<<(std::ostream &os, const gsDomain<T>& b)
{ return b.print(os); }


} // namespace gismo
