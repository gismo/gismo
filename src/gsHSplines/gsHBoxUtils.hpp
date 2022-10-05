/** @file

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once
#include <gsHSplines/gsHBox.h>
#include <gsHSplines/gsHBSplineBasis.h>
#include <gsHSplines/gsTHBSplineBasis.h>

namespace gismo {

template <short_t d, class T>
typename gsHBox<d, T>::SortedContainer gsHBoxUtils<d,T>::Sort(const Container & container)
{
    SortedContainer scontainer(container.begin(), container.end());

    // First sort (otherwise unique is wrong)
    std::sort(scontainer.begin(),scontainer.end(),gsHBoxCompare<d,T>());

    return scontainer;
}

template <short_t d, class T>
typename gsHBox<d, T>::Container gsHBoxUtils<d,T>::Unique(const Container & container)
{
    SortedContainer scontainer = gsHBoxUtils<d,T>::Sort(container);

    // Get unique entries
    typename SortedContainer::iterator it = std::unique(scontainer.begin(),scontainer.end(),gsHBoxEqual<d,T>());
    scontainer.resize(std::distance(scontainer.begin(), it));
    Container result(scontainer.begin(),scontainer.end());
    return result;
}

template <short_t d, class T>
typename gsHBox<d, T>::HContainer gsHBoxUtils<d,T>::Unique(const HContainer & container)
{
    HContainer result(container.size());
    for (size_t k=0; k!=container.size(); k++)
        result[k] = gsHBoxUtils<d,T>::Unique(container[k]);
    return result;
}

template <short_t d, class T>
gsHBoxContainer<d, T> gsHBoxUtils<d,T>::Unique(const gsHBoxContainer<d,T> & container)
{
    HContainer result = gsHBoxUtils<d,T>::Unique(container.boxes());
    return gsHBoxContainer<d,T>(result);
}

template <short_t d, class T>
typename gsHBox<d, T>::Container gsHBoxUtils<d,T>::Union(const Container & container1, const Container & container2)
{
    SortedContainer sortedResult;

    SortedContainer scontainer1 = gsHBoxUtils<d,T>::Sort(container1);
    SortedContainer scontainer2 = gsHBoxUtils<d,T>::Sort(container2);

    sortedResult.reserve(scontainer1.size() + scontainer2.size());
    if (scontainer1.size()!=0 && scontainer2.size()!=0)
    {
        std::set_union( scontainer1.begin(),scontainer1.end(),
                        scontainer2.begin(),scontainer2.end(),
                        std::inserter(sortedResult,sortedResult.begin()),
                        gsHBoxCompare<d,T>());
    }
    else if (scontainer1.size()!=0 && container2.size()==0)
        sortedResult.insert(sortedResult.end(),scontainer1.begin(),scontainer1.end());
    else if (scontainer1.size()==0 && container2.size()!=0)
        sortedResult.insert(sortedResult.end(),scontainer2.begin(),scontainer2.end());
    else    { /* Do nothing */ }

    Container result(sortedResult.begin(),sortedResult.end());

    return result;
}

template <short_t d, class T>
typename gsHBox<d, T>::HContainer gsHBoxUtils<d,T>::Union(const HContainer & container1, const HContainer & container2)
{
    HContainer result, region1, region2;

    region1 = container1;
    region2 = container2;

    index_t lmax = std::max(region1.size(),region2.size());
    region1.resize(lmax);
    region2.resize(lmax);
    result.resize(lmax);

    for (index_t l = 0; l!=lmax; l++)
        result[l] = gsHBoxUtils<d,T>::Union(region1[l],region2[l]);

    return result;
}

template <short_t d, class T>
gsHBoxContainer<d, T> gsHBoxUtils<d,T>::Union(const gsHBoxContainer<d,T> & container1, const gsHBoxContainer<d,T> & container2)
{
    HContainer result;
    HContainer region1(container1.boxes());
    HContainer region2(container2.boxes());

    index_t lmax = std::max(region1.size(),region2.size());
    region1.resize(lmax);
    region2.resize(lmax);
    result.resize(lmax);

    for (index_t l = 0; l!=lmax; l++)
        result[l] = gsHBoxUtils<d,T>::Union(region1[l],region2[l]);

    return gsHBoxContainer<d,T>(result);
}

// Takes the difference container1 \ container2
template <short_t d, class T>
typename gsHBox<d, T>::Container gsHBoxUtils<d,T>::Difference(const Container & container1, const Container & container2)
{
    SortedContainer sortedResult;

    SortedContainer scontainer1 = gsHBoxUtils<d,T>::Sort(container1);
    SortedContainer scontainer2 = gsHBoxUtils<d,T>::Sort(container2);

    sortedResult.reserve(scontainer1.size());
    // if (scontainer1.size()!=0 && scontainer2.size()!=0)
    // {
        std::set_difference(scontainer1.begin(),scontainer1.end(),
                            scontainer2.begin(),scontainer2.end(),
                            std::inserter(sortedResult,sortedResult.begin()),
                            gsHBoxCompare<d,T>());
    // }
    // else if (scontainer1.size()!=0 && container2.size()==0)
    //     sortedResult.insert(sortedResult.end(),scontainer1.begin(),scontainer1.end());
    // else if (scontainer1.size()==0 && container2.size()!=0)
    //     sortedResult.insert(sortedResult.end(),scontainer2.begin(),scontainer2.end());
    // else    { /* Do nothing */ }

    Container result(sortedResult.begin(),sortedResult.end());

    // Container result;

    // bool intersects;
    // for (cIterator it1=container1.begin(); it1!=container1.end(); it1++)
    //  for (cIterator it2=container2.begin(); it2!=container2.end(); it2++)
    //  {
    //      intersects = it1->contains(*it2) || it2->contains(*it1);
    //      if (!intersects)
    //          result.push_back(*it1);
    //      else
    //      {}
    //  }

    return result;
}

template <short_t d, class T>
typename gsHBox<d, T>::HContainer gsHBoxUtils<d,T>::Difference(const HContainer & container1, const HContainer & container2)
{
    HContainer result, region1, region2;

    region1 = container1;
    region2 = container2;

    index_t lmax = std::max(region1.size(),region2.size());
    region1.resize(lmax);
    region2.resize(lmax);
    result.resize(lmax);

    for (index_t l = 0; l!=lmax; l++)
        result[l] = gsHBoxUtils<d,T>::Difference(region1[l],region2[l]);

    return result;
}

template <short_t d, class T>
gsHBoxContainer<d, T> gsHBoxUtils<d,T>::Difference(const gsHBoxContainer<d,T> & container1, const gsHBoxContainer<d,T> & container2)
{
    HContainer result;
    HContainer region1(container1.boxes());
    HContainer region2(container2.boxes());

    index_t lmax = std::max(region1.size(),region2.size());
    region1.resize(lmax);
    region2.resize(lmax);
    result.resize(lmax);

    for (index_t l = 0; l!=lmax; l++)
        result[l] = gsHBoxUtils<d,T>::Difference(region1[l],region2[l]);

    return gsHBoxContainer<d,T>(result);
}

/**
 * @brief      Performs an intersection. Keeps the smallest boxes in overlapping regions, can also intersect partial boxes
 *
 * @param[in]  container1  The container 1
 * @param[in]  container2  The container 2
 *
 * @tparam     d           { description }
 * @tparam     T           { description }
 *
 * @return     { description_of_the_return_value }
 */
template <short_t d, class T>
typename gsHBox<d, T>::Container gsHBoxUtils<d,T>::Intersection(const Container & container1, const Container & container2)
{
    Container result;

    bool b1in2, b2in1;
    for (cIterator it1=container1.begin(); it1!=container1.end(); it1++)
        for (cIterator it2=container2.begin(); it2!=container2.end(); it2++)
        {
            b1in2 = it1->contains(*it2);
            b2in1 = it2->contains(*it1);

            if (b1in2 && b2in1)
                result.push_back(*it1);
            else if (b1in2 && !b2in1)
                result.push_back(*it2);
            else if (!b1in2 && b2in1)
                result.push_back(*it1);
            else
            {}
                // continue;
                // GISMO_ERROR("Something went wrong");
        }

    return result;
}

template <short_t d, class T>
typename gsHBox<d, T>::Container gsHBoxUtils<d,T>::ContainedIntersection(const Container & container1, const Container & container2)
{
    Container result;

    bool b1in2;
    for (cIterator it1=container1.begin(); it1!=container1.end(); it1++)
        for (cIterator it2=container2.begin(); it2!=container2.end(); it2++)
        {
            b1in2 = it1->contains(*it2);

            if (b1in2)
                result.push_back(*it2);
            else
            {}
        }

    return result;
}


/**
 * @brief      Performs an intersection; only keeps the boxes that are EXACTLY the same (also level is the same)
 *
 * @param[in]  container1  The container 1
 * @param[in]  container2  The container 2
 *
 * @tparam     d           { description }
 * @tparam     T           { description }
 *
 * @return     { description_of_the_return_value }
 */
template <short_t d, class T>
typename gsHBox<d, T>::Container gsHBoxUtils<d,T>::ExactIntersection(const Container & container1, const Container & container2)
{
    SortedContainer sortedResult;

    SortedContainer scontainer1 = gsHBoxUtils<d,T>::Sort(container1);
    SortedContainer scontainer2 = gsHBoxUtils<d,T>::Sort(container2);

    sortedResult.reserve(scontainer1.size());
    // if (scontainer1.size()!=0 && scontainer2.size()!=0)
    // {
        std::set_intersection(  scontainer1.begin(),scontainer1.end(),
                                scontainer2.begin(),scontainer2.end(),
                                std::inserter(sortedResult,sortedResult.begin()),
                                gsHBoxCompare<d,T>());
    // }
    // else if (scontainer1.size()!=0 && container2.size()==0)
    //     sortedResult.insert(sortedResult.end(),scontainer1.begin(),scontainer1.end());
    // else if (scontainer1.size()==0 && container2.size()!=0)
    //     sortedResult.insert(sortedResult.end(),scontainer2.begin(),scontainer2.end());
    // else    { /* Do nothing */ }

    Container result(sortedResult.begin(),sortedResult.end());
    return result;
}

// template <short_t d, class T>
// typename gsHBox<d, T>::Container gsHBoxUtils<d,T>::gsHBoxIntersection(const gsHBox<d,T> & box1, const gsHBox<d,T> & box2)
// {
//  Container result, container1, container2;

//  container1 = box1.toUnitBoxes();
//  container2 = box2.toUnitBoxes();

//  bool b1in2, b2in1;
//  for (cIterator it1=container1.begin(); it1!=container1.end(); it1++)
//      for (cIterator it2=container2.begin(); it2!=container2.end(); it2++)
//      {
//          bool b1in2 = it1->contains(*it2);
//          bool b2in1 = it2->contains(*it1);

//          if (b1in2 && b2in1)
//              result.push_back(*it1);
//          else if (b1in2 && !b2in1)
//              result.push_back(*it1);
//          else if (!b1in2 && b2in1)
//              result.push_back(*it2);
//          else
//              GISMO_ERROR("Something went wrong");
//      }

//  return result;
// }

template <short_t d, class T>
typename gsHBox<d,T>::Container gsHBoxUtils<d,T>::HContainer2Container( const HContainer & container )
{
    Container result;
    for (cHIterator hit = container.begin(); hit!=container.end(); hit++)
        for (cIterator it=hit->begin(); it!=hit->end(); it++)
            result.push_back(*it);

    return result;
}

template <short_t d, class T>
typename gsHBox<d,T>::HContainer gsHBoxUtils<d,T>::Container2HContainer( const Container & container )
{
    HContainer result(1);
    for (cIterator it=container.begin(); it!=container.end(); it++)
    {
        if (result.size() < static_cast<unsigned>(it->level() + 1))
            result.resize(it->level()+1);

        result[it->level()].push_back(*it);
    }
    return result;
}

template <short_t d, class T>
typename gsHBox<d,T>::Container gsHBoxUtils<d, T>::toContainer(const HContainer & container)
{
    Container result;

    for (cHIterator hit = container.begin(); hit!=container.end(); hit++)
        for (cIterator it = hit->begin(); it!=hit->end(); it++)
        {
            result.push_back(*it);
        }

    return result;
}

template <short_t d, class T>
typename gsHBox<d,T>::Container gsHBoxUtils<d, T>::toUnitBoxes(const HContainer & container)
{
    Container result;
    Container  boxes;

    for (cHIterator hit = container.begin(); hit!=container.end(); hit++)
        for (cIterator it = hit->begin(); it!=hit->end(); it++)
        {
            boxes = it->toUnitBoxes();
            for (cIterator boxIt = boxes.begin(); boxIt != boxes.end(); boxIt++)
                result.push_back(*boxIt);
        }

    return result;
}

template <short_t d, class T>
typename gsHBox<d,T>::HContainer gsHBoxUtils<d, T>::toUnitHBoxes(const HContainer & container)
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
typename gsHBoxUtils<d,T>::HContainer gsHBoxUtils<d, T>::markTadmissible(const HContainer & marked, index_t m)
{
    return gsHBoxUtils<d,T>::markAdmissible<Neighborhood::T>(marked,m);
}

template <short_t d, class T>
typename gsHBoxUtils<d,T>::HContainer gsHBoxUtils<d, T>::markTadmissible(const gsHBox<d,T> & marked, index_t m)
{
    return gsHBoxUtils<d,T>::markAdmissible<Neighborhood::T>(marked,m);
}


template <short_t d, class T>
typename gsHBoxUtils<d,T>::HContainer gsHBoxUtils<d, T>::markHadmissible(const HContainer & marked, index_t m)
{
    return gsHBoxUtils<d,T>::markAdmissible<Neighborhood::H>(marked,m);
}

template <short_t d, class T>
typename gsHBoxUtils<d,T>::HContainer gsHBoxUtils<d, T>::markHadmissible(const gsHBox<d,T> & marked, index_t m)
{
    return gsHBoxUtils<d,T>::markAdmissible<Neighborhood::H>(marked,m);
}


template <short_t d, class T>
template<enum Neighborhood _mode>
typename std::enable_if<_mode==Neighborhood::T || _mode==Neighborhood::H, typename gsHBoxUtils<d, T>::HContainer>::type
gsHBoxUtils<d, T>::_markRecursive(const HContainer & marked, index_t lvl, index_t m)
{
    HContainer marked_copy = marked;
    Container marked_l = marked[lvl];
    Container marked_k;

    gsHBoxContainer<d,T> neighbors;
    for (Iterator it = marked_l.begin(); it!=marked_l.end(); it++)
    {
        neighbors.add(it->template getNeighborhood<_mode>(m));
    }

    index_t k = lvl - m + 1;
    if (neighbors.boxes().size()!=0)
    {
        marked_k = marked_copy[k];
        gsHBoxContainer<d,T> boxUnion = gsHBoxUtils<d,T>::Union(neighbors,gsHBoxContainer<d,T>(marked_k));
        marked_copy[k] = boxUnion.getActivesOnLevel(k);
        marked_copy = gsHBoxUtils<d,T>::_markRecursive<_mode>(marked_copy,k,m);
    }
    return marked_copy;
}

template <short_t d, class T>
template<enum Neighborhood _mode>
typename std::enable_if<_mode!=Neighborhood::T && _mode!=Neighborhood::H, typename gsHBoxUtils<d, T>::HContainer>::type
gsHBoxUtils<d, T>::_markRecursive(const HContainer & marked, index_t lvl, index_t m)
{
    GISMO_ERROR("Mode must be T or H");
}

template <short_t d, class T>
template<enum Neighborhood _mode>
typename gsHBoxUtils<d,T>::HContainer gsHBoxUtils<d, T>::markAdmissible(const HContainer & marked, index_t m)
{
    HContainer unitBoxes = gsHBoxUtils<d,T>::toUnitHBoxes(marked);
    for (size_t l = 0; l!=unitBoxes.size(); l++)
        unitBoxes = gsHBoxUtils<d,T>::_markRecursive<_mode>(unitBoxes,l,m);

    unitBoxes = gsHBoxUtils<d,T>::Unique(unitBoxes);
    return unitBoxes;
}

template <short_t d, class T>
template<enum Neighborhood _mode>
typename gsHBoxUtils<d,T>::HContainer gsHBoxUtils<d, T>::markAdmissible(const gsHBox<d,T> & marked, index_t m)
{
    HContainer unitBoxes = gsHBoxUtils<d,T>::Container2HContainer(marked.toUnitBoxes());
    return gsHBoxUtils<d,T>::markAdmissible<_mode>(unitBoxes,m);
}

// template <short_t d, class T>
// typename gsHBoxUtils<d,T>::HContainer gsHBoxUtils<d, T>::markAdmissible(const gsHBox<d,T> & marked, index_t m)
// {
//     HContainer result;
//     if (dynamic_cast<const gsTHBSplineBasis<d,T>*>(marked.basis()))
//         result = gsHBoxUtils<d,T>::getNeighborhood<Neighborhood::T>(m);
//     else if (dynamic_cast<const gsHBSplineBasis<d,T>*>(marked.basis()))
//         result = gsHBoxUtils<d,T>::getNeighborhood<Neighborhood::T>(m);
//     else
//         GISMO_ERROR("Basis type should be gsTHBSplineBasis or gsHBSplineBasis");
//     return result;
// }

template <short_t d, class T>
typename gsHBoxUtils<d,T>::HContainer gsHBoxUtils<d, T>::markAdmissible(const HContainer & marked, index_t m)
{
    HContainer result;
    if (dynamic_cast<const gsTHBSplineBasis<d,T>*>(&marked.front().front().basis()))
        result = gsHBoxUtils<d,T>::markAdmissible<Neighborhood::T>(marked,m);
    else if (dynamic_cast<const gsHBSplineBasis<d,T>*>(&marked.front().front().basis()))
        result = gsHBoxUtils<d,T>::markAdmissible<Neighborhood::H>(marked,m);
    else
        GISMO_ERROR("Basis type should be gsTHBSplineBasis or gsHBSplineBasis");
    return result;
}

template <short_t d, class T>
typename gsHBoxUtils<d,T>::HContainer gsHBoxUtils<d, T>::markAdmissible(const gsHBox<d,T> & marked, index_t m)
{
    HContainer result;
    if (dynamic_cast<const gsTHBSplineBasis<d,T>*>(&marked.basis()))
        result = gsHBoxUtils<d,T>::markAdmissible<Neighborhood::T>(marked,m);
    else if (dynamic_cast<const gsHBSplineBasis<d,T>*>(&marked.basis()))
        result = gsHBoxUtils<d,T>::markAdmissible<Neighborhood::H>(marked,m);
    else
        GISMO_ERROR("Basis type should be gsTHBSplineBasis or gsHBSplineBasis");
    return result;
}

template <short_t d, class T>
bool gsHBoxUtils<d, T>::allActive(const Container & elements)
{
    bool check = true;
    for (cIterator it=elements.begin(); it!=elements.end() && check; it++)
    {
        check &= it->isActive();
    }
    return check;
}

template <short_t d, class T>
bool gsHBoxUtils<d, T>::allActive(const HContainer & elements)
{
    bool check = true;
    for (cHIterator hit = elements.begin(); hit!=elements.end() && check; hit++)
        check &= allActive(*hit);
    return check;
}


template <short_t d, class T>
bool gsHBoxCompare<d,T>::operator()(const gsHBox<d,T> & a, const gsHBox<d,T> & b) const
{
    return
     (a.patch() < b.patch())
    ||
    ((a.patch() == b.patch()) &&
     (a.level() < b.level())     )
    ||
    ((a.patch() == b.patch()) &&
     (a.level() == b.level()) &&
     std::lexicographical_compare(  a.lowerIndex().begin(), a.lowerIndex().end(),
                                    b.lowerIndex().begin(), b.lowerIndex().end())   )
    ||
    ((a.patch() == b.patch()) &&
     (a.level() == b.level()) &&
     (a.lowerIndex() == b.lowerIndex()) &&
     std::lexicographical_compare(  a.upperIndex().begin(), a.upperIndex().end(),
                                    b.upperIndex().begin(), b.upperIndex().end())    );
};

template <short_t d, class T>
bool gsHBoxEqual<d,T>::operator()(const gsHBox<d,T> & a, const gsHBox<d,T> & b) const
{
    return a.isSame(b);
};

// template <short_t d, class T>
// bool gsHBoxOverlaps<d,T>::operator()(const gsHBox<d,T> & a, const gsHBox<d,T> & b) const
// {
//     bool res = true;
//     if (a.level() <= b.level())
//     {
//      // Check if the indices of the ancestor are contained. If yes, the indices of the original box should be contained
//      for (index_t i=0; i!=d && res; i++)
//      {
//          res |= a.lowerIndex().at(i) >= b.getAncestor(a.level()).lowerIndex().at(i);
//          res |= a.upperIndex().at(i) <= b.getAncestor(a.level()).upperIndex().at(i);
//      }
//     }
//     else //a.level() > b.level()
//     {
//      // Check if the indices of the ancestor are contained. If yes, the indices of the original box should be contained
//      for (index_t i=0; i!=d && res; i++)
//      {
//          res |= b.lowerIndex().at(i) >= a.getAncestor(b.level()).lowerIndex().at(i);
//          res |= b.upperIndex().at(i) <= a.getAncestor(b.level()).upperIndex().at(i);
//      }
//     }

//     return res;
// };

template <short_t d, class T>
bool gsHBoxContains<d,T>::operator()(const gsHBox<d,T> & a, const gsHBox<d,T> & b) const
{
    bool res = true;
    res &= a.level() <= b.level();
    if (res)
    {
        // Check if the indices of the ancestor are contained. If yes, the indices of the original box should be contained
        for (index_t i=0; i!=d && res; i++)
        {
            res &= a.lowerIndex().at(i) >= b.getAncestor(a.level()).lowerIndex().at(i);
            res &= a.upperIndex().at(i) <= b.getAncestor(a.level()).upperIndex().at(i);
        }
    }
    return res;
};

// template <short_t d, class T>
// bool gsHBoxIsContained<d,T>::operator()(const gsHBox<d,T> & a, const gsHBox<d,T> & b) const
// {
//     bool res = true;
//     res &= a.level() <= b.level();
//     if (res)
//     {
//         for (index_t i=0; i!=d && res; i++)
//         {
//             res &= a.lowerIndex().at(i) >= b.getAncestor(a.level()).lowerIndex().at(i);
//             res &= a.upperIndex().at(i) <= b.getAncestor(a.level()).upperIndex().at(i);
//         }
//     }
//     return res;
// };

template <short_t d, class T>
bool gsHBoxIsContained<d,T>::operator()(const gsHBox<d,T> & a, const gsHBox<d,T> & b) const
{
    bool res = true;
    res &= a.level() >= b.level();
    if (res)
    {
        // Check if the indices of the ancestor are contained. If yes, the indices of the original box should be contained
        for (index_t i=0; i!=d && res; i++)
        {
            res &= a.getAncestor(b.level()).lowerIndex().at(i) >= b.lowerIndex().at(i);
            res &= a.getAncestor(b.level()).upperIndex().at(i) <= b.upperIndex().at(i);
        }
    }
    return res;
};


} // gismo
