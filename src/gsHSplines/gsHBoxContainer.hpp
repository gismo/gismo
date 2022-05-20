/** @file gsHTensorBasis.hpp

    @brief Provides implementation of HTensorBasis common operations.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): G. Kiss, A. Mantzaflaris, J. Speh
*/

#pragma once

#include <gsHSplines/gsHBoxUtils.h>


namespace gismo
{

template <short_t d, class T>
gsHBoxContainer<d, T>::gsHBoxContainer()
{ }

template <short_t d, class T>
gsHBoxContainer<d, T>::gsHBoxContainer(const gsHBox<d,T> & box)
{
    this->_makeLevel(box.level());
    m_boxes[box.level()].push_back(box);
}

/**
 * @brief      Takes a Container (which can have boxes with different levels)
 *
 * @param[in]  box   The box
 *
 * @tparam     d     { description }
 * @tparam     T     { description }
 */
template <short_t d, class T>
gsHBoxContainer<d, T>::gsHBoxContainer(const Container & boxes)
{
    for (cIterator it=boxes.begin(); it!=boxes.end(); it++)
    {
        this->_makeLevel(it->level());
        m_boxes[it->level()].push_back(*it);
    }

}

template <short_t d, class T>
gsHBoxContainer<d, T>::gsHBoxContainer(const HContainer & boxes)
{
    if (_check(boxes)) // if the container is well-defined
    {
        m_boxes = boxes;
    }
    else
    {
        for (cHIterator hit = boxes.begin(); hit!=boxes.end(); hit++)
        {
            for (cIterator it = hit->begin(); it!=hit->end(); it++)
            {
                this->_makeLevel(it->level());
                m_boxes[it->level()].push_back(*it);
            }
        }
    }
}

template <short_t d, class T>
size_t gsHBoxContainer<d, T>::totalSize() const
{
    std::function<size_t(const size_t &, const Container &)> size_sum = [](const size_t & sum, const Container & a)
    {
        return sum+a.size();
    };
    return std::accumulate(m_boxes.begin(), m_boxes.end(), 0, size_sum);
}

/**
 * @brief      Checks if the level of the boxes in the container matches the container level
 *
 * @tparam     d     { description }
 * @tparam     T     { description }
 */
template <short_t d, class T>
bool gsHBoxContainer<d, T>::_check(const HContainer & boxes)
{
    bool result = true;
    for (size_t l = 0; l!=boxes.size(); l++)
        for (cIterator it = boxes[l].begin(); it!=boxes[l].end(); it++)
            result &= ( (size_t) it->level()==l);
    return result;
}


template <short_t d, class T>
void gsHBoxContainer<d, T>::add(const gsHBox<d,T> & box)
{
    this->_makeLevel(box.level());
    m_boxes[box.level()].push_back(box);
}

template <short_t d, class T>
void gsHBoxContainer<d, T>::add(const Container & boxes)
{
    for (cIterator it = boxes.begin(); it!=boxes.end(); it++)
        this->add(*it);
}

template <short_t d, class T>
void gsHBoxContainer<d, T>::add(const HContainer & boxes)
{
    for (cHIterator hit = boxes.begin(); hit!=boxes.end(); hit++)
        this->add(*hit);
}

template <short_t d, class T>
void gsHBoxContainer<d, T>::add(const gsHBoxContainer<d,T> & boxes)
{
    this->add(boxes.boxes());
}

template <short_t d, class T>
gsHBoxContainer<d,T> gsHBoxContainer<d, T>::boxUnion(const gsHBoxContainer<d,T> & other) const
{
    return boxUnion(*this,other);
}

template <short_t d, class T>
gsHBoxContainer<d,T> gsHBoxContainer<d, T>::boxUnion(const gsHBoxContainer<d,T> & container1, const gsHBoxContainer<d,T> & container2) const
{
    HContainer result;
    HContainer region1(container1.boxes());
    HContainer region2(container2.boxes());

    index_t lmax = std::max(region1.size(),region2.size());
    region1.resize(lmax);
    region2.resize(lmax);
    result.resize(lmax);

    for (index_t l = 0; l!=lmax; l++)
        result[l] = _boxUnion(region1[l],region2[l]);

    return gsHBoxContainer<d,T>(result);
}

template <short_t d, class T>
typename gsHBoxContainer<d,T>::HContainer gsHBoxContainer<d, T>::boxUnion(const HContainer & container1, const HContainer & container2) const
{
    HContainer result, region1, region2;

    region1 = container1;
    region2 = container2;

    index_t lmax = std::max(region1.size(),region2.size());
    region1.resize(lmax);
    region2.resize(lmax);
    result.resize(lmax);

    for (index_t l = 0; l!=lmax; l++)
        result[l] = _boxUnion(region1[l],region2[l]);

    return result;
}

template <short_t d, class T>
typename gsHBoxContainer<d, T>::Container gsHBoxContainer<d, T>::_boxUnion(const Container & container1, const Container & container2) const
{
    // SortedContainer sortedResult;

    // SortedContainer scontainer1 = gsHBoxSort<d,T>(container1);
    // SortedContainer scontainer2 = gsHBoxSort<d,T>(container2);

    // sortedResult.reserve(scontainer1.size() + scontainer2.size());
    // if (scontainer1.size()!=0 && scontainer2.size()!=0)
    // {
    //     std::set_union( scontainer1.begin(),scontainer1.end(),
    //                     scontainer2.begin(),scontainer2.end(),
    //                     std::inserter(sortedResult,sortedResult.begin()),
    //                     gsHBoxCompare<d,T>());
    // }
    // else if (scontainer1.size()!=0 && container2.size()==0)
    //     sortedResult.insert(sortedResult.end(),scontainer1.begin(),scontainer1.end());
    // else if (scontainer1.size()==0 && container2.size()!=0)
    //     sortedResult.insert(sortedResult.end(),scontainer2.begin(),scontainer2.end());
    // else    { /* Do nothing */ }

    Container result = gsHBoxUnion<d,T>(container1,container2);

    return result;
}

// template <short_t d, class T>
// void gsHBoxContainer<d, T>::makeUnique()
// {
//     auto comp = [](auto a, auto b)
//                     {
//                         return a.isSame(b);
//                     };

//     std::unique( this->begin(),this->end(),comp);
// }


// template <short_t d, class T>
// typename gsHBoxContainer<d, T>::Container & gsHBoxContainer<d, T>::getActives() const
// {
//     HContainer boxes(m_boxes.size());
//     for (index_t lvl = 0; lvl != m_boxes.size(); lvl++)
//         boxes[lvl] = this->getActivesOnLevel(lvl);

//     return boxes;
// }

// template <short_t d, class T>
// void gsHBoxContainer<d, T>::filterActives()
// {
//     m_boxes = getActives();
// }

// template <short_t d, class T>
// typename gsHBoxContainer<d, T>::Container & gsHBoxContainer<d, T>::getActivesOnLevel(index_t lvl) const
// {
//     Container boxes;
//     for (cIterator it=m_boxes[lvl].begin(); it!=m_boxes[lvl].end(); it++)
//         if (it->isActive())
//             boxes.push_back(*it);

//     return boxes;
// }


template <short_t d, class T>
typename gsHBoxContainer<d, T>::Container & gsHBoxContainer<d, T>::getActivesOnLevel(index_t lvl)
{
    return m_boxes[lvl];
}

template <short_t d, class T>
const typename gsHBoxContainer<d, T>::Container & gsHBoxContainer<d, T>::getActivesOnLevel(index_t lvl) const
{
    return m_boxes[lvl];
}

template <short_t d, class T>
typename gsHBoxContainer<d, T>::HContainer gsHBoxContainer<d, T>::getParents() const
{
    HContainer result;
    result.resize(this->boxes().size()-1);

    // Handle level 0 separately (operation cannot be performed on level 0)
    GISMO_ASSERT(m_boxes[0].size()==0,"Boxes at level 0 cannot have a parent. Did something go wrong? You can run check() to see if the boxes are allocated coorectly");

    HIterator resIt = result.begin();
    for (cHIterator hit = std::next(m_boxes.begin()); hit!=m_boxes.end(); hit++)
        for (cIterator it=hit->begin(); it!=hit->end(); it++)
            resIt->push_back(*it);

    return result;
}

template <short_t d, class T>
typename gsHBoxContainer<d, T>::HContainer gsHBoxContainer<d, T>::markTrecursive(HContainer & marked, index_t lvl, index_t m) const
{
    Container marked_l = marked[lvl];
    Container marked_k;

    gsHBoxContainer<d,T> neighbors;
    for (Iterator it = marked_l.begin(); it!=marked_l.end(); it++)
        neighbors.add(it->getTneighborhood(m));

    index_t k = lvl - m + 1;
    if (neighbors.boxes().size()!=0)
    {
        marked_k = marked[k];
        gsHBoxContainer<d,T> boxunion = boxUnion(neighbors,gsHBoxContainer<d,T>(marked_k));
        marked[k] = boxunion.getActivesOnLevel(k);
        marked = this->markTrecursive(marked,k,m);
    }
    return marked;
}

template <short_t d, class T>
void gsHBoxContainer<d, T>::markTadmissible(HContainer & marked, index_t m) const
{
    HContainer unitBoxes = this->toUnitBoxes(marked);
    for (size_t l = 0; l!=unitBoxes.size(); l++)
        this->markTrecursive(unitBoxes,l,m);
    marked = unitBoxes;
}

template <short_t d, class T>
void gsHBoxContainer<d, T>::markTadmissible(index_t m)
{
    this->markTadmissible(this->boxes(),m);
}

template <short_t d, class T>
void gsHBoxContainer<d, T>::markTrecursive(index_t lvl, index_t m)
{
    m_boxes = this->markTrecursive(m_boxes,lvl,m);
}

template <short_t d, class T>
typename gsHBoxContainer<d, T>::HContainer gsHBoxContainer<d, T>::markHrecursive(HContainer & marked, index_t lvl, index_t m) const
{
    Container marked_l = marked[lvl];
    Container marked_k;

    gsHBoxContainer<d,T> neighbors;
    for (Iterator it = marked_l.begin(); it!=marked_l.end(); it++)
        neighbors.add(it->getHneighborhood(m));

    index_t k = lvl - m + 1;
    if (neighbors.boxes().size()!=0)
    {
        marked_k = marked[k];
        gsHBoxContainer<d,T> boxunion = boxUnion(neighbors,gsHBoxContainer<d,T>(marked_k));
        marked[k] = boxunion.getActivesOnLevel(k);
        marked = this->markHrecursive(marked,k,m);
    }
    return marked;
}

template <short_t d, class T>
void gsHBoxContainer<d, T>::markHrecursive(index_t lvl, index_t m)
{
    m_boxes = this->markHrecursive(m_boxes,lvl,m);
}

template <short_t d, class T>
void gsHBoxContainer<d, T>::markHadmissible(HContainer & marked, index_t m) const
{
    HContainer unitBoxes = this->toUnitBoxes(marked);
    for (size_t l = 0; l!=unitBoxes.size(); l++)
        this->markHrecursive(unitBoxes,l,m);
    marked = unitBoxes;
}

template <short_t d, class T>
void gsHBoxContainer<d, T>::markHadmissible(index_t m)
{
    this->markHadmissible(this->boxes(),m);
}

template <short_t d, class T>
void gsHBoxContainer<d, T>::_makeLevel( index_t lvl )
{
    if (m_boxes.size() < static_cast<unsigned>(lvl + 1))
        m_boxes.resize(lvl+1);
}

template <short_t d, class T>
std::ostream& gsHBoxContainer<d, T>::print( std::ostream& os ) const
{
    for (cHIterator hit = m_boxes.begin(); hit!=m_boxes.end(); hit++)
        for (cIterator it = hit->begin(); it!=hit->end(); it++)
            os<<*it<<"\n";
    return os;
}

template <short_t d, class T>
typename gsHBoxContainer<d,T>::RefBox gsHBoxContainer<d, T>::toBoxes() const
{
    size_t N = this->totalSize();
    RefBox result;
    result.reserve(( N * (2*d+1) ));
    RefBox box;
    for (cHIterator hit = m_boxes.begin(); hit!=m_boxes.end(); hit++)
        for (cIterator it = hit->begin(); it!=hit->end(); it++)
        {
            box = it->toBox();
            for (typename RefBox::const_iterator boxIt = box.begin(); boxIt != box.end(); boxIt++)
                result.push_back(*boxIt);
        }

    return result;
}

template <short_t d, class T>
typename gsHBoxContainer<d,T>::RefBox gsHBoxContainer<d, T>::toRefBoxes() const
{
    size_t N = this->totalSize();
    RefBox result;
    result.reserve(( N * (2*d+1) ));
    RefBox box;
    for (cHIterator hit = m_boxes.begin(); hit!=m_boxes.end(); hit++)
        for (cIterator it = hit->begin(); it!=hit->end(); it++)
        {
            box = it->toRefBox();
            for (typename RefBox::const_iterator boxIt = box.begin(); boxIt != box.end(); boxIt++)
                result.push_back(*boxIt);
        }

    return result;
}

template <short_t d, class T>
typename gsHBoxContainer<d,T>::RefBox gsHBoxContainer<d, T>::toCrsBoxes() const
{
    size_t N = this->totalSize();
    RefBox result;
    result.reserve(( N * (2*d+1) ));
    RefBox box;
    for (cHIterator hit = m_boxes.begin(); hit!=m_boxes.end(); hit++)
        for (cIterator it = hit->begin(); it!=hit->end(); it++)
        {
            box = it->toCrsBox();
            for (typename RefBox::const_iterator boxIt = box.begin(); boxIt != box.end(); boxIt++)
                result.push_back(*boxIt);
        }

    return result;
}

template <short_t d, class T>
gsMatrix<T> gsHBoxContainer<d, T>::toCoords()
{
    size_t N = this->totalSize();
    gsMatrix<T> boxes(d,2*N);
    index_t boxCount = 0;
    for (HIterator hit = m_boxes.begin(); hit!=m_boxes.end(); hit++)
        for (Iterator it = hit->begin(); it!=hit->end(); it++)
        {
            boxes.block(0,2*boxCount,d,2) = it->getCoordinates();
            boxCount++;
        }

    return boxes;
}

template <short_t d, class T>
typename gsHBoxContainer<d,T>::HContainer gsHBoxContainer<d, T>::toUnitBoxes(const HContainer & container) const
{
    HContainer result(container.size());
    HIterator  resIt  = result.begin();
    Container  boxes;

    for (cHIterator hit = container.begin(); hit!=container.end(); hit++, resIt++)
        for (cIterator it = hit->begin(); it!=hit->end(); it++)
        {
            boxes = it->toUnitBoxes();
            for (cIterator boxIt = boxes.begin(); boxIt != boxes.end(); boxIt++)
                resIt->push_back(*boxIt);
        }

    return result;
}

template <short_t d, class T>
typename gsHBoxContainer<d,T>::HContainer gsHBoxContainer<d, T>::toUnitBoxes() const
{
    return this->toUnitBoxes(this->m_boxes);
}

template <short_t d, class T>
void gsHBoxContainer<d, T>::makeUnitBoxes()
{
    this->m_boxes = this->toUnitBoxes();
}

} // namespace gismo
