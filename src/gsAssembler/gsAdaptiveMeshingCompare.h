/*/** @file gsAdaptiveMeshingCompare.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (TU Delft, 2019-)
*/

#pragma once


#include <iostream>
#include <map>
#include <gsHSplines/gsHBoxUtils.h>

namespace gismo
{


/**
 * @brief      Base class for performing checks on \ref gsHBox objects
 *
 * @tparam     d     { description }
 * @tparam     T     { description }
 */
template <short_t d, class T>
class gsHBoxCheck
{
public:
    virtual ~gsHBoxCheck() {};

    virtual bool check(const gsHBox<d,T> & box) const = 0;
};

/**
 * @brief      Checks if the level of a \ref gsHBox is bigger than a minimum level
 *
 * @tparam     d     parametric dimension
 * @tparam     T     real type
 */
template <short_t d, class T>
class gsMinLvlCompare : public gsHBoxCheck<d,T>
{
public:
    explicit gsMinLvlCompare(index_t minlevel = 0)
    :
    m_minLevel(minlevel)
    {}

    /// Checks the box
    bool check(const gsHBox<d,T> & box) const { return box.level() > m_minLevel; }

protected:
    index_t m_minLevel;
};

/**
 * @brief      Checks if the level of a \ref gsHBox is smaller than a maximum level
 *
 * @tparam     d     parametric dimension
 * @tparam     T     real type
 */
template <short_t d, class T>
class gsMaxLvlCompare : public gsHBoxCheck<d,T>
{
public:
    explicit gsMaxLvlCompare(index_t maxLevel)
    :
    m_maxLevel(maxLevel)
    {}

    bool check(const gsHBox<d,T> & box) const { return box.level() < m_maxLevel; }

protected:
    index_t m_maxLevel;
};

/**
 * @brief      Checks if the error of a \ref gsHBox is smaller than a threshold
 *
 * @tparam     d     parametric dimension
 * @tparam     T     real type
 */
template <short_t d, class T>
class gsSmallerErrCompare : public gsHBoxCheck<d,T>
{
public:
    gsSmallerErrCompare(const T & threshold)
    :
    m_threshold(threshold)
    {}

    bool check(const gsHBox<d,T> & box) const { return box.error() < m_threshold; }

protected:
    T m_threshold;
};

/**
 * @brief      Checks if the error of a \ref gsHBox is larger than a threshold
 *
 * @tparam     d     parametric dimension
 * @tparam     T     real type
 */
template <short_t d, class T>
class gsLargerErrCompare : public gsHBoxCheck<d,T>
{
public:
    gsLargerErrCompare(const T & threshold)
    :
    m_threshold(threshold)
    {}

    bool check(const gsHBox<d,T> & box) const { return box.error() > m_threshold; }

protected:
    T m_threshold;
};

/**
 * @brief      Checks if the coarsening neighborhood of a box is empty and if it overlaps with a refinement mask
 *
 * Checks if the coarsening neighborhood of a box is empty and if it overlaps with a refinement mask.
 * If so, the box can be coarsened admissibly.
 *
 * @tparam     d     parametric dimension
 * @tparam     T     real type
 */
template <short_t d, class T>
class gsOverlapCompare : public gsHBoxCheck<d,T>
{
public:
    /**
     * @brief      Construct a gsOverlapCompare
     *
     * @param[in]  markedRef  Container of elements marked for refinement
     * @param[in]  m          Jump parameter
     */
    gsOverlapCompare(const gsHBoxContainer<d,T> & markedRef, index_t m) //, patchHContainer & markedCrs
    :
    m_m(m)
    {
        gsHBoxContainer<d,T> tmp(markedRef);
        m_markedRefChildren = gsHBoxUtils<d,T>::Unique(tmp.getChildren());
    }

    /**
     * @brief      Checks if the coarsening neighborhood of the parent of \a box is clean
     *
     * The result depends on \a box only through its parent and on the mesh, so it
     * is cached per parent and a box and its siblings share one evaluation.
     * Precondition: the mesh (basis) does not change while the predicate is
     * alive; construct a new predicate after refining or coarsening. The cache is
     * not thread-safe.
     */
    bool check(const gsHBox<d,T> & box) const
    {
        // A level-0 box has no parent, so it cannot be coarsened; guard here rather
        // than relying on a preceding gsMinLvlCompare in the predicate list.
        if (box.level() == 0) return false;

        gsHBox<d,T> parent = box.getParent();
        typename Cache::const_iterator hit = m_cache.find(parent);
        if (hit != m_cache.end()) return hit->second;
        return m_cache.insert(typename Cache::value_type(parent,_check(parent))).first->second;
    }

protected:
    // Checks the coarsening extension (closely related to the coarsening neighborhood) of \a parent, i.e. the parent of the box that will be elevated.
    bool _check(gsHBox<d,T> parent) const
    {
        // 1) Check if the coarsening neighborhood is empty
        typename gsHBox<d,T>::Container Cextension = parent.getCextension(m_m);
        Cextension = gsHBoxUtils<d,T>::Unique(Cextension);

        for (typename gsHBox<d,T>::Iterator it = Cextension.begin(); it != Cextension.end(); it++)
        {
            it->computeCenter();
            // the level is even larger (i.e. even higher decendant); then it is not clean
            if (it->levelInCenter()>=it->level()) return false;
        }

        // 2) Now we check if the parents of any of the cells in the extensions overlap with the marked cells. If so, it would cause a problem with the coarsening.
        // Equivalent to an empty gsHBoxUtils::ContainedIntersection(Cextension,m_markedRefChildren), but stops at the first overlap. O(|Cextension| |m_markedRefChildren|) worst case.
        gsHBoxContains<d,T> contains;
        for (typename gsHBox<d,T>::Iterator it = Cextension.begin(); it != Cextension.end(); it++)
            for (typename gsHBox<d,T>::cIterator ch = m_markedRefChildren.begin(); ch != m_markedRefChildren.end(); ch++)
                if (contains(*it,*ch)) return false;

        return true;
    }

    typedef std::map<gsHBox<d,T>,bool,gsHBoxCompare<d,T> > Cache;

    typename gsHBox<d,T>::Container m_markedRefChildren;
    index_t m_m;
    mutable Cache m_cache; ///< result of _check per parent box
};


} // namespace gismo
