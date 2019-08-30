/** @file gsBasisRefs.h

    @brief Small class holding references to a list of gsBasis

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

namespace gismo
{

/** @brief Simple class to hold a list of gsBasis which discretize a
 * PDE system on a given patch.
 *
 * \ingroup Core
*/
template <class T>
class gsBasisRefs
{
private:
    typedef typename std::vector<gsMultiBasis<T> >::const_iterator mbIterator;

public:

    /// Constructor from a gsMultiBasis \a bases and a patch index \a k
    inline gsBasisRefs(std::vector<gsMultiBasis<T> > & bases, const size_t k)
    {
        GISMO_ASSERT(bases.size()>0, "Cannot construct empty list of gsBasis.");
        m_refs.reserve(bases.size());
        for(mbIterator  it = bases.begin(); it != bases.end(); ++it)
            m_refs.push_back( &(*it)[k] );
    }

    /// Accessor for a certain gsBasis
    inline const gsBasis<T> & operator[](size_t i) const
    { return *m_refs[i]; }

    /// Cast to gsBasis to emulate scalar problem
    inline operator const gsBasis<T> &() const 
    { return *m_refs.front(); }

    /// Front
    inline const gsBasis<T> & front() const 
    { return *m_refs.front(); }

    /// Back
    inline const gsBasis<T> & back () const 
    { return *m_refs.back(); }

    /// Size
    inline size_t size () const
    { return m_refs.size(); }

    /// Return the parametric dimension of the bases (assumed to be the same for all the bases)
    inline short_t dim() const { return m_refs.front()->dim();}

private:
    std::vector<const gsBasis<T>*> m_refs;

private:
    gsBasisRefs();
    gsBasisRefs(const gsBasisRefs &);
};

} // namespace gismo
