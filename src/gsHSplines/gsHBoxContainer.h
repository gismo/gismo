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
    typedef typename std::list<gsHBox<d,T>,typename gsHBox<d,T>::aalloc>    Container;
    typedef typename std::vector<gsHBox<d,T>,typename gsHBox<d,T>::aalloc>  SortedContainer;
    typedef typename std::vector<Container>                                 HContainer;
    typedef typename Container::iterator        Iterator;
    typedef typename Container::const_iterator  cIterator;
    typedef typename HContainer::iterator       HIterator;
    typedef typename HContainer::const_iterator cHIterator;
    // typedef typename Container::const_iterator citerator;

public:
    gsHBoxContainer();

    gsHBoxContainer(const gsHBox<d,T> & box  );
    gsHBoxContainer(const Container   & boxes);
    gsHBoxContainer(const HContainer  & boxes);

    index_t size() const { return m_boxes.size(); }

    // iterator begin() const {return m_boxes.begin();}
    // iterator end()   const {return m_boxes.end();}

    bool check() { return _check(this->boxes()); };

    void add(const gsHBox<d,T>          & box  );
    void add(const Container            & boxes);
    void add(const HContainer           & boxes);
    void add(const gsHBoxContainer<d,T> & boxes);

    std::ostream& print( std::ostream& os ) const;

    gsHBoxContainer<d,T> boxUnion(const gsHBoxContainer<d,T> & other) const;
    void makeUnique();

    Container &  getActivesOnLevel(index_t lvl);
    const Container & getActivesOnLevel(index_t lvl) const;
    HContainer getParents() const;

    gsHBoxContainer<d,T> markHrecursive(index_t lvl, index_t m) const;
    gsHBoxContainer<d,T> markTrecursive(index_t lvl, index_t m) const;


    HContainer & boxes() { return m_boxes; }
    const HContainer & boxes() const { return m_boxes; }

    // getBoundingBox()

protected:
    Container _boxUnion(const Container & container1, const Container & container2) const;

    void _makeLevel(index_t lvl);

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
