/** @file gsStdVectorRef.h

    @brief Small wrapper for std::vector

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

namespace gismo
{

/** 
    @brief Simple wrapper class for a vector of objects.

    The list casts to object by returning the first object in the list
    
    \ingroup Core
*/
template<class obj>
class gsStdVectorRef
{
public:

    /// Constructor from a vector of objs
    inline gsStdVectorRef(const std::vector<obj> & refVec) :
        m_ref(refVec)
    { }


    /// Accessor
    inline const obj & operator[](size_t i) const
    { return m_ref[i]; }

    /// Cast to std::vector
    inline operator const std::vector<obj> & () const
    { return m_ref; }

    /// Cast to by returning the first element
    inline operator const obj &() const
    { return m_ref.front(); }

    /// Front
    inline const obj & front() const
    { return m_ref.front(); }

    /// Back
    inline const obj & back () const
    { return m_ref.back(); }

    /// Size
    inline size_t size () const
    { return m_ref.size(); }

private:
    const std::vector<obj> & m_ref;

private:
    gsStdVectorRef();
    gsStdVectorRef(const gsStdVectorRef &);
};


} // namespace gismo
