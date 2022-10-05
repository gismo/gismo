/** @file gsHTensorBasis.hpp

    @brief Provides implementation of HTensorBasis common operations.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): G. Kiss, A. Mantzaflaris, J. Speh
*/

#pragma once

#include <gsIO/gsXml.h>

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

template <short_t d, class T>
void gsHBoxContainer<d, T>::clear()
{
    for (HIterator hit = m_boxes.begin(); hit!=m_boxes.end(); hit++)
        hit->clear();
    m_boxes.clear();
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
gsHBoxContainer<d,T> gsHBoxContainer<d, T>::patch(const index_t patchID) const
{
    gsHBoxContainer<d,T> result;
    for (cHIterator hit = m_boxes.begin(); hit!=m_boxes.end(); hit++)
        for (cIterator it = hit->begin(); it!=hit->end(); it++)
            if (it->patch()==patchID || it->patch()==-1)
                result.add(*it);
    return result;
}

// template <short_t d, class T>
// gsHBoxContainer<d,T> gsHBoxContainer<d, T>::_boxUnion(const gsHBoxContainer<d,T> & other) const
// {
//     return _boxUnion(*this,other);
// }

// template <short_t d, class T>
// gsHBoxContainer<d,T> gsHBoxContainer<d, T>::_boxUnion(const gsHBoxContainer<d,T> & container1, const gsHBoxContainer<d,T> & container2) const
// {
//     HContainer result;
//     HContainer region1(container1.boxes());
//     HContainer region2(container2.boxes());

//     index_t lmax = std::max(region1.size(),region2.size());
//     region1.resize(lmax);
//     region2.resize(lmax);
//     result.resize(lmax);

//     for (index_t l = 0; l!=lmax; l++)
//         result[l] = __boxUnion(region1[l],region2[l]);

//     return gsHBoxContainer<d,T>(result);
// }

// template <short_t d, class T>
// typename gsHBoxContainer<d,T>::HContainer gsHBoxContainer<d, T>::_boxUnion(const HContainer & container1, const HContainer & container2) const
// {
//     HContainer result, region1, region2;

//     region1 = container1;
//     region2 = container2;

//     index_t lmax = std::max(region1.size(),region2.size());
//     region1.resize(lmax);
//     region2.resize(lmax);
//     result.resize(lmax);

//     for (index_t l = 0; l!=lmax; l++)
//         result[l] = __boxUnion(region1[l],region2[l]);

//     return result;
// }

// template <short_t d, class T>
// typename gsHBoxContainer<d, T>::Container gsHBoxContainer<d, T>::__boxUnion(const Container & container1, const Container & container2) const
// {
//     // SortedContainer sortedResult;

//     // SortedContainer scontainer1 = gsHBoxUtils<d,T>::Sort(container1);
//     // SortedContainer scontainer2 = gsHBoxUtils<d,T>::Sort(container2);

//     // sortedResult.reserve(scontainer1.size() + scontainer2.size());
//     // if (scontainer1.size()!=0 && scontainer2.size()!=0)
//     // {
//     //     std::set_union( scontainer1.begin(),scontainer1.end(),
//     //                     scontainer2.begin(),scontainer2.end(),
//     //                     std::inserter(sortedResult,sortedResult.begin()),
//     //                     gsHBoxCompare<d,T>());
//     // }
//     // else if (scontainer1.size()!=0 && container2.size()==0)
//     //     sortedResult.insert(sortedResult.end(),scontainer1.begin(),scontainer1.end());
//     // else if (scontainer1.size()==0 && container2.size()!=0)
//     //     sortedResult.insert(sortedResult.end(),scontainer2.begin(),scontainer2.end());
//     // else    { /* Do nothing */ }

//     Container result = gsHBoxUtils<d,T>::Union(container1,container2);

//     return result;
// }

// template <short_t d, class T>
// void gsHBoxContainer<d, T>::makeUnique()
// {
//     auto comp = [](const gsHBox<d,T> & a, const gsHBox<d,T> & b)
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
typename gsHBoxContainer<d, T>::Container gsHBoxContainer<d, T>::getParents() const
{
    Container result;
    // result.reserve(this->boxes().size()-1);

    // Handle level 0 separately (operation cannot be performed on level 0)
    GISMO_ASSERT(m_boxes[0].size()==0,"Boxes at level 0 cannot have a parent. Did something go wrong? You can run check() to see if the boxes are allocated coorectly");

    for (cHIterator hit = std::next(m_boxes.begin()); hit!=m_boxes.end(); hit++)
        for (cIterator it=hit->begin(); it!=hit->end(); it++)
            result.push_back(it->getParent());

    return result;
}

template <short_t d, class T>
typename gsHBoxContainer<d, T>::Container gsHBoxContainer<d, T>::getChildren() const
{
    Container result, children;

    for (cHIterator hit = m_boxes.begin(); hit!=m_boxes.end(); hit++)
        for (cIterator it=hit->begin(); it!=hit->end(); it++)
        {
            children = it->getChildren();
            for (cIterator cit=children.begin(); cit!=children.end(); cit++)
                result.push_back(*cit);
        }

    return result;
}

template <short_t d, class T>
void gsHBoxContainer<d, T>::markTadmissible(index_t m)
{
    m_boxes = gsHBoxUtils<d,T>::markTadmissible(m_boxes,m);
}

template <short_t d, class T>
void gsHBoxContainer<d, T>::markHadmissible(index_t m)
{
    m_boxes = gsHBoxUtils<d,T>::markHadmissible(m_boxes,m);
}

template <short_t d, class T>
void gsHBoxContainer<d, T>::markAdmissible(index_t m)
{
    m_boxes = gsHBoxUtils<d,T>::markAdmissible(m_boxes,m);
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
typename gsHBoxContainer<d,T>::RefBox gsHBoxContainer<d, T>::toBoxes(const index_t patchID) const
{
    size_t N = this->totalSize();
    RefBox result;
    result.reserve(( N * (2*d+1) ));
    RefBox box;
    gsHBoxContainer<d,T> patchBoxes = this->patch(patchID);
    for (cHIterator hit = patchBoxes.begin(); hit!=patchBoxes.end(); hit++)
        for (cIterator it = hit->begin(); it!=hit->end(); it++)
        {
            box = it->toBox();
            for (typename RefBox::const_iterator boxIt = box.begin(); boxIt != box.end(); boxIt++)
                result.push_back(*boxIt);
        }

    return result;
}

template <short_t d, class T>
typename gsHBoxContainer<d,T>::RefBox gsHBoxContainer<d, T>::toRefBoxes(const index_t patchID) const
{
    gsHBoxContainer<d,T> patchBoxes = this->patch(patchID);
    size_t N = patchBoxes.totalSize();
    RefBox result;
    result.reserve(( N * (2*d+1) ));
    RefBox box;
    for (cHIterator hit = patchBoxes.begin(); hit!=patchBoxes.end(); hit++)
        for (cIterator it = hit->begin(); it!=hit->end(); it++)
        {
            box = it->toRefBox();
            for (typename RefBox::const_iterator boxIt = box.begin(); boxIt != box.end(); boxIt++)
                result.push_back(*boxIt);
        }

    return result;
}

template <short_t d, class T>
typename gsHBoxContainer<d,T>::RefBox gsHBoxContainer<d, T>::toCrsBoxes(const index_t patchID) const
{
    gsHBoxContainer<d,T> patchBoxes = this->patch(patchID);
    size_t N = patchBoxes.totalSize();
    RefBox result;
    result.reserve(( N * (2*d+1) ));
    RefBox box;
    for (cHIterator hit = patchBoxes.begin(); hit!=patchBoxes.end(); hit++)
        for (cIterator it = hit->begin(); it!=hit->end(); it++)
        {
            box = it->toCrsBox();
            for (typename RefBox::const_iterator boxIt = box.begin(); boxIt != box.end(); boxIt++)
                result.push_back(*boxIt);
        }

    return result;
}

template <short_t d, class T>
gsMatrix<T> gsHBoxContainer<d, T>::toCoords(const index_t patchID) const
{
    gsHBoxContainer<d,T> patchBoxes = this->patch(patchID);
    size_t N = patchBoxes.totalSize();
    gsMatrix<T> boxes(d,2*N);
    index_t boxCount = 0;
    for (HIterator hit = patchBoxes.begin(); hit!=patchBoxes.end(); hit++)
        for (Iterator it = hit->begin(); it!=hit->end(); it++)
        {
            boxes.block(0,2*boxCount,d,2) = it->getCoordinates();
            boxCount++;
        }

    return boxes;
}

template <short_t d, class T>
void gsHBoxContainer<d, T>::makeUnitBoxes()
{
    this->m_boxes = this->toUnitHBoxes();
}

namespace internal
{

/// @brief Get a FunctionsExpr from XML data
template<short_t d, class T>
class gsXml< gsHBoxContainer<d,T> >
{
private:
    gsXml() { }
    typedef gsHBoxContainer<d,T> Object;
public:
    GSXML_COMMON_FUNCTIONS(Object);
    static std::string tag ()  { return "HBoxContainer"; }
    static std::string type () { return "HBoxContainer"+std::to_string(d); }

    GSXML_GET_POINTER(Object);

    static void get_into (gsXmlNode * node, Object & obj)
    {
        gsXmlNode * boxNode;
        for (boxNode = node->first_node("HBox");
             boxNode; boxNode = boxNode->next_sibling("HBox"))
        {
            gsHBox<d,T> * box = gsXml<gsHBox<d,T> >::get(boxNode);
            obj.add(*box);
        }
    }

    static gsXmlNode * put (const Object & obj,
                            gsXmlTree & data )
    {
        gsXmlNode * container = makeNode("HBoxContainer", data);
        container->append_attribute( makeAttribute("type",internal::gsXml<Object>::type().c_str(), data) );
        container->append_attribute(makeAttribute("size", obj.totalSize(), data));

        for (typename Object::cHIterator hit = obj.cbegin(); hit!=obj.cend(); hit++)
            for (typename Object::cIterator it = hit->begin(); it!=hit->end(); it++)
            {
                gsXmlNode * box = gsXml< gsHBox<d,T> >::put(*it,data);
                container->append_node(box);
            }
        return container;
    }
};

} // internal

} // namespace gismo
