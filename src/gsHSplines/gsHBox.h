/** @file gsHBox.h

    @brief Provides definition of the gsHBox

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): G. Kiss, A. Mantzaflaris, J. Speh
*/

#pragma once

#include <gsHSplines/gsHDomainIterator.h>
#include <gsHSplines/gsHTensorBasis.h>

namespace gismo
{

template<short_t d, class T> class gsHBoxContainer;

template<int d, class T>
class gsHBox
{
public:
    typedef gsVector<index_t,d> point;
    typedef typename Eigen::aligned_allocator<gsHBox<d,T>> aalloc;

    typedef typename std::list<gsHBox<d,T>,typename gsHBox<d,T>::aalloc>    Container;
    typedef typename std::vector<gsHBox<d,T>,typename gsHBox<d,T>::aalloc>  SortedContainer;
    typedef typename std::vector<Container>                                 HContainer;
    typedef typename Container::iterator        Iterator;
    typedef typename Container::const_iterator  cIterator;
    typedef typename HContainer::iterator       HIterator;
    typedef typename HContainer::const_iterator cHIterator;


public:
    ~gsHBox() {  }

    gsHBox() { }

    gsHBox(const gsHDomainIterator<T,d> * domHIt);

    gsHBox(const point & low,const point & upp, index_t level);

    gsHBox(const point & low,const point & upp, index_t level, const gsHTensorBasis<d,T> * basis);

    gsHBox(const gsAabb<d,index_t> & box);

    gsHBox(const gsAabb<d,index_t> & box, const gsHTensorBasis<d,T> * basis);

    gsHBox(const std::vector<index_t> & indices);

    /// Copy constructor (makes deep copy)
    gsHBox( const gsHBox<d,T> & other );

    /// Move constructor
    gsHBox( gsHBox<d,T> && other );

    /// Assignment operator
    gsHBox<d,T> & operator= ( const gsHBox<d,T> & other );

    /// Move assignment operator
    gsHBox<d,T> & operator= ( gsHBox<d,T> && other );

    /**
     * @brief      Checks if the \a other cell is contained in \a this cell
     *
     * @param      other  The other cell
     *
     * @return     True if the specified other is contained in the cell, False otherwise.
     */
    bool isContained(const gsHBox<d,T> & other) const;

    /**
     * @brief      Checks if the \a other cell is contains in \a this cell
     *
     * @param      other  The other cell
     *
     * @return     True if the specified other is contains in the cell, False otherwise.
     */
    bool contains(const gsHBox<d,T> & other) const;

    bool isSame(const gsHBox<d,T> & other) const;

    bool isActive();

    const gsMatrix<T> & getCoordinates();

    const gsMatrix<T> & getCenter();

    gsVector<T,d> lowerCorner() const;
    gsVector<T,d> upperCorner() const;

    point lowerIndex() const;
    point upperIndex() const;

    index_t level() const;

    gsHBox<d,T> getParent() const;

    gsHBox<d,T> getAncestor(index_t k) const;

    bool isActiveOnLvl(index_t lvl) const;

    HContainer getSupportExtension();

    HContainer getMultiLevelSupportExtension(index_t k);

    Container getHneighborhood(index_t m);

    Container getTneighborhood(index_t m);

    HContainer           toContainer();
    gsHBoxContainer<d,T> toHBoxContainer();

    std::ostream& print( std::ostream& os ) const;

protected:
    void _computeCoordinates();

    void _computeIndices();
    gsAabb<d,index_t> _computeIndices(const gsMatrix<T> & coords, index_t level);
    gsAabb<d,index_t> _computeIndices(const gsMatrix<T> & coords, const gsMatrix<T> & center);
    gsAabb<d,index_t> _computeIndices(const gsMatrix<T> & coords);

    gsAabb<d,index_t> _elevateBox(const gsAabb<d,index_t> & box) const;

protected:
    gsAabb<d,index_t> m_indices;

    gsMatrix<T> m_coords;
    gsMatrix<T> m_center;
    const gsHTensorBasis<d,T> * m_basis;


public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
}; // class gsHBox


template<int d, class T>
std::ostream& operator<<( std::ostream& os, const gsHBox<d,T>& b )
{
    return b.print( os );
}

} // namespace gismo

// *****************************************************************
#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsHBox.hpp)
#else
#ifdef gsHBox_EXPORT
#include GISMO_HPP_HEADER(gsHBox.hpp)
#undef  EXTERN_CLASS_TEMPLATE
#define EXTERN_CLASS_TEMPLATE CLASS_TEMPLATE_INST
#endif
namespace gismo
{
EXTERN_CLASS_TEMPLATE gsHBox<1,real_t>;
EXTERN_CLASS_TEMPLATE gsHBox<2,real_t>;
EXTERN_CLASS_TEMPLATE gsHBox<3,real_t>;
EXTERN_CLASS_TEMPLATE gsHBox<4,real_t>;
}
#endif
// *****************************************************************
