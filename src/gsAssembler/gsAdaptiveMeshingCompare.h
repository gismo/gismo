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
    gsMinLvlCompare(const T & minlevel = 0)
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
    gsMaxLvlCompare(const T & maxLevel)
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

    bool check(const gsHBox<d,T> & box) const
    {
        // We are going to check if the coarsening extension (closely related to the coarsening neighborhood) of the parent of \a box (since it will be elevated) fulfills the conditions of an empty coarsening neighborhood, as well as the condition of an empty coasening neighborhood provided that there is no element that will be refined herein.
        bool clean = true;

        // if (m_m>=2) // admissiblity part
        // {
            // 1) Check if the coarsening neighborhood is empty
            gsHBox<d,T> parent = box.getParent();
            typename gsHBox<d,T>::Container Cextension = parent.getCextension(m_m);
            Cextension = gsHBoxUtils<d,T>::Unique(Cextension);

            for (typename gsHBox<d,T>::Iterator it = Cextension.begin(); it != Cextension.end() && clean; it++)
            {
                it->computeCenter();
                clean &=
                        // the level is even larger (i.e. even higher decendant); then it is not clean
                            (!((it->levelInCenter()>=it->level()))
                        )
                        ;
            }

            if (!clean) return clean;
        // }

        // 2) Now we check if the parents of any of the cells in the extensions overlap with the marked cells. If so, it would cause a problem with the coarsening.
        typename gsHBox<d,T>::Container  intersection = gsHBoxUtils<d,T>::ContainedIntersection(Cextension,m_markedRefChildren);
        clean = intersection.size()==0;
        return clean;
    }

protected:
    typename gsHBox<d,T>::Container m_markedRefChildren;
    index_t m_m;
};


} // namespace gismo
