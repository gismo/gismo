/** @file gsHBox.h

    @brief Provides definition of the gsHBox

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): G. Kiss, A. Mantzaflaris, J. Speh
*/

#pragma once

#include <gsHSplines/gsHBox.h>

namespace gismo
{

template<short_t d, class T>
class gsHBoxContainer
{
public:
    // std::list does not provide .at(k) but it provides iterators
    typedef typename gsHBox<d,T>::RefBox            RefBox;
    typedef typename gsHBox<d,T>::Container         Container;
    typedef typename gsHBox<d,T>::SortedContainer   SortedContainer;
    typedef typename gsHBox<d,T>::HContainer        HContainer;
    typedef typename gsHBox<d,T>::Iterator          Iterator;
    typedef typename gsHBox<d,T>::cIterator         cIterator;
    typedef typename gsHBox<d,T>::HIterator         HIterator;
    typedef typename gsHBox<d,T>::cHIterator        cHIterator;

public:
    gsHBoxContainer();

    gsHBoxContainer(const gsHBox<d,T> & box  );
    gsHBoxContainer(const Container   & boxes);
    gsHBoxContainer(const HContainer  & boxes);

    /// Returns the size of the container on \a level
    size_t size(index_t level)  const { return m_boxes[level].size(); }
    /// Returns the number of levels stored in the container
    size_t nLevels()            const { return m_boxes.size(); }
    /// Returns the total number of boxes
    size_t totalSize()          const;

    /// Checks if the hierarchical container is correctly defined
    bool check() { return _check(this->boxes()); };

    /// Adds a single box
    void add(const gsHBox<d,T>          & box  );
    /// Adds boxes stored in a container
    void add(const Container            & boxes);
    /// Adds boxes stored in a hierarchical container
    void add(const HContainer           & boxes);
    /// Adds boxes stored in a \ref gsHBoxContainer
    void add(const gsHBoxContainer<d,T> & boxes);

    /// Prints the container
    std::ostream& print( std::ostream& os ) const;


    /**
     * @brief      Takes the union of two hierarchical containers
     *
     * @param[in]  container1  The container 1
     * @param[in]  container2  The container 2
     *
     * @return     The hierarchical container with the union.
     */
    HContainer           boxUnion(const HContainer & container1, const HContainer & container2) const;
    /**
     * @brief      Takes the union of \a this and the \a other \ref gsHBoxContainer
     *
     * @param[in]  other  The other gsHBoxContainer
     *
     * @return     The gsHBoxContainer with the union.
     */
    gsHBoxContainer<d,T> boxUnion(const gsHBoxContainer<d,T> & other) const;
    /**
     * @brief      Takes the union of two \ref gsHBoxContainer
     *
     * @param[in]  container1  The container 1
     * @param[in]  container2  The container 2
     *
     * @return     The gsHBoxContainer with the union.
     */
    gsHBoxContainer<d,T> boxUnion(const gsHBoxContainer<d,T> & container1, const gsHBoxContainer<d,T> & container2) const;

    /// Removes duplicate boxes
    void makeUnique();

    /// Returns the actives on \a level
    Container &  getActivesOnLevel(index_t lvl);
    /// Returns the actives on \a level
    const Container & getActivesOnLevel(index_t lvl) const;
    /// Gives a hierarchical container with all the parents of the boxes stored in \a this
    HContainer getParents() const;


    /**
     * @brief      Marks H-recursively
     *
     * @param      marked  The marked boxes
     * @param[in]  lvl     The level
     * @param[in]  m       The jump parameter
     *
     * @return     The resulting hierarchical container.
     */
    HContainer  markHrecursive(HContainer & marked, index_t lvl, index_t m) const;
    /// Applies \ref markHrecursive on \a this
    void        markHrecursive(index_t lvl, index_t m);
    /**
     * @brief      Marks T-recursively
     *
     * @param      marked  The marked boxes
     * @param[in]  lvl     The level
     * @param[in]  m       The jump parameter
     *
     * @return     The resulting hierarchical container.
     */
    HContainer  markTrecursive(HContainer & marked, index_t lvl, index_t m) const;
    /// Applies \ref markTrecursive on \a this
    void        markTrecursive(index_t lvl, index_t m);

    /**
     * @brief      Performs H-admissible refinement
     *
     * @param      marked  The marked boxes
     * @param[in]  m       The jump parameter
     *
     * @return     The resulting hierarchical container.
     */
    void        markHadmissible(HContainer & marked, index_t m) const;
    /// Applies \ref markHadmissible on \a this
    void        markHadmissible(index_t m);
    /**
     * @brief      Performs T-admissible refinement
     *
     * @param      marked  The marked boxes
     * @param[in]  m       The jump parameter
     *
     * @return     The resulting hierarchical container.
     */
    void        markTadmissible(HContainer & marked, index_t m) const;
    /// Applies \ref markTadmissible on \a this
    void        markTadmissible(index_t m);

    /// Returns a heirarchical container with the boxes stored in the container
    HContainer & boxes() { return m_boxes; }
    /// Returns a heirarchical container with the boxes stored in the container
    const HContainer & boxes() const { return m_boxes; }

    /// Returns the maximum level in the container
    index_t maxLevel() {return m_boxes.size()-1; }

    /**
     * @brief      Returns boxes representation of the object.
     *
     * @return     Boxes representation of the object.
     */
    RefBox toBoxes()    const;

    /**
     * @brief      Returns refinement box representation of the object.
     *
     * @return     Refinement box representation of the object.
     */
    RefBox toRefBoxes() const;

    /**
     * @brief      Returns coarsening box representation of the object.
     *
     * @return     Coarsening box representation of the object.
     */
    RefBox toCrsBoxes() const;

    /**
     * @brief      Returns box coordinate represenation of the object
     *
     * @return     Box coordinate representation of the object
     */
    gsMatrix<T> toCoords();

    /**
     * @brief      Transforms the boxes in \a container as unit boxes
     *
     * @param[in]  container  A hierarchical container of boxes
     *
     * @return     A hierarchical container containing the unit boxes.
     */
    HContainer toUnitBoxes(const HContainer & container) const;

    /**
     * @brief      Returns the boxes inside as unit boxes
     *
     * @return     A hierarchical container containing the unit boxes.
     */
    HContainer toUnitBoxes() const;


    /**
     * @brief      Transforms/splits the boxes inside the container to unit boxes
     */
    void makeUnitBoxes();


protected:
    /// Helper to take the box union
    Container _boxUnion(const Container & container1, const Container & container2) const;

    /// Constructs a new level
    void _makeLevel(index_t lvl);

    /// Checks the container
    bool _check(const HContainer & boxes);


protected:
    HContainer m_boxes;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

}; // class gsHBoxContainer

template<short_t d, class T>
std::ostream& operator<<( std::ostream& os, const gsHBoxContainer<d,T>& b )
{
    return b.print( os );
}

} // namespace gismo

// *****************************************************************
#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsHBoxContainer.hpp)
#else
#ifdef gsHBoxContainer_EXPORT
#include GISMO_HPP_HEADER(gsHBoxContainer.hpp)
#undef  EXTERN_CLASS_TEMPLATE
#define EXTERN_CLASS_TEMPLATE CLASS_TEMPLATE_INST
#endif
namespace gismo
{
EXTERN_CLASS_TEMPLATE gsHBoxContainer<1,real_t>;
EXTERN_CLASS_TEMPLATE gsHBoxContainer<2,real_t>;
EXTERN_CLASS_TEMPLATE gsHBoxContainer<3,real_t>;
EXTERN_CLASS_TEMPLATE gsHBoxContainer<4,real_t>;
}
#endif
// *****************************************************************
