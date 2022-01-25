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

template<int d, class T>
class gsHBox
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

public:
    typedef gsVector<index_t,d> point;
    typedef typename Eigen::aligned_allocator<gsHBox<d,T>> aalloc;

    typedef typename std::vector<index_t>                                   RefBox;
    typedef typename std::list<gsHBox<d,T>,typename gsHBox<d,T>::aalloc>    Container;
    typedef typename std::vector<gsHBox<d,T>,typename gsHBox<d,T>::aalloc>  SortedContainer;
    typedef typename std::vector<Container>                                 HContainer;
    typedef typename Container::iterator        Iterator;
    typedef typename Container::const_iterator  cIterator;
    typedef typename HContainer::iterator       HIterator;
    typedef typename HContainer::const_iterator cHIterator;


public:

    gsHBox() { }

    gsHBox(const gsHDomainIterator<T,d> * domHIt);

    gsHBox(const typename gsHBox<d,T>::point & low,const typename gsHBox<d,T>::point & upp, index_t level, const gsHTensorBasis<d,T> * basis);

    gsHBox(const gsAabb<d,index_t> & box, const gsHTensorBasis<d,T> * basis);

    gsHBox(const std::vector<index_t> & indices, const gsHTensorBasis<d,T> * basis);

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

    /**
     * @brief      Determines whether the \a this box is the same as the \a other box
     *
     * @param[in]  other  The other box
     *
     * @return     True if the boxes are the same, False otherwise.
     */
    bool isSame(const gsHBox<d,T> & other) const;

    /**
     * @brief      Determines if the box is active on its current level.
     *
     * @return     True if active, False otherwise.
     */
    bool isActive();

    /**
     * @brief      Gets the coordinates of the box (first column lower corner, second column higher corner).
     *
     * @return     The coordinates of the box.
     */
    const gsMatrix<T> & getCoordinates();

    /**
     * @brief      Gets the center of the box.
     *
     * @return     The center of the box.
     */
    const gsMatrix<T> & getCenter();

    /**
     * @brief      Gets the lower corner of the box
     *
     * @return     The lower corner
     */
    gsVector<T,d> lowerCorner() const;

    /**
     * @brief      Gets the upper corner of the box
     *
     * @return     The upper corner
     */
    gsVector<T,d> upperCorner() const;

    /**
     * @brief      Gets the lower index of the box
     *
     * @return     The lower index
     */
    const point & lowerIndex() const;
    /**
     * @brief      Gets the upper index of the box
     *
     * @return     The upper index
     */
    const point & upperIndex() const;

    /**
     * @brief      Gets the level of the object
     *
     * @return     The level of the object
     */
    index_t level() const;

    /**
     * @brief      Gets the parent of the object.
     *
     * @return     The parent of the object.
     */
    gsHBox<d,T> getParent() const;

    /**
     * @brief      Gets the ancestor of the object on level \a k.
     *
     * @param[in]  k     The reference level
     *
     * @return     The ancestor.
     */
    gsHBox<d,T> getAncestor(index_t k) const;

    /**
     * @brief      Gets the support extension.
     *
     * @return     The support extension.
     */
    Container getSupportExtension();

    /**
     * @brief      Gets the multi-level support extension.
     *
     * @param[in]  k     The reference level
     *
     * @return     The multi-level support extension.
     */
    Container getMultiLevelSupportExtension(index_t k);

    /**
     * @brief      Gets the H-neighborhood.
     *
     * @param[in]  m     The jump parameter
     *
     * @return     The H-neighborhood.
     */
    Container getHneighborhood(index_t m);

    /**
     * @brief      Gets the T-neighborhood.
     *
     * @param[in]  m     The jump parameter
     *
     * @return     The T-neighborhood.
     */
    Container getTneighborhood(index_t m);

    /**
     * @brief      Returns a hierarchical container representation of the object.
     *
     * @return     Hierarchical container representation of the object.
     */
    HContainer           toContainer();

    /**
     * @brief      Prints the object
     *
     */
    std::ostream& print( std::ostream& os ) const;

    /**
     * @brief      Returns a box representation of the object
     *
     * @return     Box representation of the object.
     */
    RefBox toBox() const;

    /**
     * @brief      Returns a box representation of the object on the lower level (needed for refinement).
     *
     * @return     Refinement box representation of the object.
     */
    RefBox toRefBox() const;

    // Helper functions
    HContainer           boxUnion(const HContainer & container1, const HContainer & container2) const;

    const gsHTensorBasis<d,T> & basis() { return *m_basis; }

    /**
     * @brief      Returns unit boxes representation of the object.
     *
     * @return     Unit boxes representation of the object.
     */
    Container toUnitBoxes() const;

    /**
     * @brief      Checks if the box has positive indices
     *
     * @return     True if the box has positive indices, false otherwise
     */
    bool good() const;

    void clean(Container & container) const;


protected:
    void _computeCoordinates();

    void _computeIndices();
    gsAabb<d,index_t> _computeIndices(const gsMatrix<T> & coords, index_t level);
    gsAabb<d,index_t> _computeIndices(const gsMatrix<T> & coords, const gsMatrix<T> & center);
    gsAabb<d,index_t> _computeIndices(const gsMatrix<T> & coords);

    gsAabb<d,index_t> _elevateBox(const gsAabb<d,index_t> & box) const;

    Container  _getParents(typename gsHBox<d,T>::Container  & container) const;
    HContainer _getParents(typename gsHBox<d,T>::HContainer & container) const;


    Container _boxUnion(const Container & container1, const Container & container2) const;
    Container _makeUnique(const Container & container) const;

protected:
    gsAabb<d,index_t> m_indices;

    gsMatrix<T> m_coords;
    gsMatrix<T> m_center;
    const gsHTensorBasis<d,T> * m_basis;

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
