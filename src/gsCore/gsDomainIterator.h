/** @file gsDomainIterator.h

    @brief Provides declaration of DomainIterator abstract interface.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither
*/

#pragma once

//#include <gsCore/gsBasis.h> // todo: remove
#include <gsCore/gsDomain.h>
#include <gsCore/gsBoundary.h>

namespace gismo
{

/** The gsDomainIterator \n
 * \brief Class which enables iteration over all elements of a parameter domain.
 *
 *
 * It also includes some additional functionality which is typically used
 * when assembling the system matrix for numerically solving a PDE.
 *
 * - <b>Iteration through the elements:</b>\n
 * The function next() jumps to the "next" element and sets up the quadrature
 * nodes and weights on that element.
 * The specific implementation
 * of how to select the next element depends on the structure of the underlying mesh.\n
 * The function good() indicates whether there still is a "next" element to be found.
 *
 *
 * Note that the features of the gsDomainIterator strongly depend on the underlying basis.
 * Hence the gsBasis is given as an input argument to the constructor.
 *
 * An example of the typical use of gsDomainIterator (remark: replace
 * the constructor by the constructor of the actually used derived
 * class):
 *
 * \verbatim
     gsDomainIterator domIter( basis );         // constructor

     for (; domIter.good(); domIter.next() )    // loop over all elements
     {
         // Your source code using
         domIter.centerPoint();
         domIter.lowerCorner();
         domIter.upperCorner();

     }
     \endverbatim
     
     \ingroup Core
 *
 *
 *
 */



template <class T>
class gsDomainIterator
{
public:
    /// Shared pointer for gsDomainIterator
    typedef memory::shared_ptr< gsDomainIterator > Ptr;
    /// Unique pointer for gsDomainIterator
    typedef memory::unique_ptr< gsDomainIterator > uPtr;

public:

    gsDomainIterator( ) : m_basis(NULL), m_isGood( true ) { }

    /// \brief Constructor using a basis 
    gsDomainIterator( const gsBasis<T>& basisParam, const boxSide & s = boundary::none)
        : center( gsVector<T>::Zero(basisParam.dim()) ), m_basis( &basisParam ), 
          m_isGood( true ), m_side(s)
    { }

    virtual ~gsDomainIterator() { }

public:

    /** @brief Proceeds to the next element.
     *
     * The function returns true if there are still elements remaining that have not been treated.\n
     * For the typical usage of this function, see the example in the
     * documentation of gsDomainIterator.
     */
    virtual bool next() = 0;

    /// \brief Proceeds to the next element (skipping \p increment elements).
    virtual bool next(index_t increment) = 0;

    /// Resets the iterator so that it points to the first element
    virtual void reset()
    {
        GISMO_NO_IMPLEMENTATION
    }

public:
    /// Is the iterator still pointing to a valid element?
    bool good() const   { return m_isGood; }

    /// Return dimension of the elements
    short_t dim() const   { return center.size(); }

    /// Updates \a other with and adjacent element
    /// \todo upgrade to return adjacent range instead
    virtual void adjacent( const gsVector<bool> & ,
                           gsDomainIterator & )
    {
        GISMO_NO_IMPLEMENTATION
    }

    /// \brief Returns the center of the current element.
    ///
    /// The current element is a <em>d</em>-dimensional hypercube.
    /// The coordinates of its upper corner is returned as a gsVector of length \a d.\n
    /// \n
    /// E.g., if the current two-dimensional element is defined by <em>[a,b]x[c,d]</em>, then <em>[b,d]</em> is returned (see also lowerCorner()).
    const gsVector<T>& centerPoint () const
    { return center; }

    /// \brief Returns the lower corner of the current element.
    ///
    /// The current element is a <em>d</em>-dimensional hypercube.
    /// The coordinates of its lower corner is returned as a gsVector of length \a d.\n
    /// \n
    /// E.g., if the current two-dimensional element is defined by <em>[a,b]x[c,d]</em>, then <em>[a,c]</em> is returned (see also upperCorner()).
    virtual const gsVector<T>& lowerCorner() const
    {
        GISMO_NO_IMPLEMENTATION
    }

    /// \brief Returns the upper corner of the current element.
    ///
    /// The current element is a <em>d</em>-dimensional hypercube.
    /// The coordinates of its upper corner is returned as a gsVector of length \a d.\n
    /// \n
    /// E.g., if the current two-dimensional element is defined by <em>[a,b]x[c,d]</em>, then <em>[b,d]</em> is returned (see also lowerCorner()).
    virtual const gsVector<T>& upperCorner() const
    {
        GISMO_NO_IMPLEMENTATION
    }

    /// \brief Returns the perpendicular cell size of boundary iterator.
    ///
    /// Only works for boundary iterators. Returns the length from
    /// the boundary side to the parallel side not on the boundary.
    virtual const T getPerpendicularCellSize() const
    {
        GISMO_NO_IMPLEMENTATION
    }

    /// Return the diagonal of the element
    T getCellSize() const
    {
        return (upperCorner() - lowerCorner()).norm();
    }
    
    /// Return the length of the smallest edge of the element
    T getMinCellLength() const
    {
        return (upperCorner() - lowerCorner()).minCoeff();
    }

    /// Return the length of the largest edge of the element
    T getMaxCellLength() const
    {
        return (upperCorner() - lowerCorner()).maxCoeff();
    }

    /// Return the volume of the element
    T volume() const
    { return (upperCorner() - lowerCorner()).prod(); }

    /// Returns the number of elements.
    virtual size_t numElements() const
    {
        //\todo Remove this implementation. Probably using a shallow
        //copy, "reset" and "next" would do this better.

        // Buggy, and probably a terrible implementation,
        // but needed and therefore can be useful
        // sometimes.
        typename gsBasis<T>::domainIter domIter = m_basis->makeDomainIterator(m_side);

        size_t numEl = 0;
        for (; domIter->good(); domIter->next(), numEl++){}

        return numEl;
    }


    inline boxSide side() const {return m_side;}

public:

    /// Coordinates of a central point in the element (in the parameter domain).
    gsVector<T> center;

protected:
    /// The basis on which the domain iterator is defined.
    const gsBasis<T> * m_basis;

    /// Flag indicating whether the domain iterator is "good". If it
    /// is "good", the iterator can continue to the next element.
    bool m_isGood;

    boxSide m_side;

private:
    // disable copying
    gsDomainIterator( const gsDomainIterator& );
    gsDomainIterator& operator= ( const gsDomainIterator& );
}; // class gsDomainIterator


/// Print (as string) operator to be used by all derived classes
//template<class T>
//std::ostream &operator<<(std::ostream &os, const gsDomainIterator<T>& b)
//{return b.print(os); }


} // namespace gismo
