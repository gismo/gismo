/** @file gsAdaptiveRefUtils.h

    @brief Provides class for adaptive refinement.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (TU Delft, 2019-)
*/

#pragma once


#include <iostream>
#include <gsAssembler/gsAdaptiveRefUtils.h>
#include <gsIO/gsOptionList.h>
#include <gsCore/gsMultiPatch.h>
#include <gsCore/gsMultiBasis.h>
#include <gsHSplines/gsHBox.h>
#include <gsHSplines/gsHBoxContainer.h>

namespace gismo
{

// enum MarkingStrategy
// {
//     GARU=1,
//     PUCA=2,
//     BULK=3
// };

template <class T>
class gsAdaptiveMeshing
{
public:
    typedef typename std::vector<gsHBoxContainer<2,T>> patchHContainer;

public:

    gsAdaptiveMeshing(gsFunctionSet<T> & input)
    :
    m_input(&input)
    {
        defaultOptions();
        getOptions();
    }

    gsOptionList & options() {return m_options;}

    void defaultOptions();

    void getOptions();

    void mark(const std::vector<T> & errors) { mark(errors,m_maxLvl); }
    void mark(const std::vector<T> & errors, index_t maxLvl);

    void refine(const patchHContainer & markedRef);
    void unrefine(const patchHContainer & markedCrs);
    void adapt(const patchHContainer & markedRef, const patchHContainer & markedCrs);

    void refine(const std::vector<bool> & markedRef) { refine(_toContainer(markedRef)); }
    void unrefine(const std::vector<bool> & markedCrs) { unrefine(_toContainer(markedCrs)); };
    void adapt(const std::vector<bool> & markedRef,const std::vector<bool> & markedCrs) { adapt(_toContainer(markedRef),_toContainer(markedCrs)); }

    void refine() { refine(m_markedRef); }
    void unrefine() { unrefine(m_markedRef); };
    void adapt() { adapt(m_markedRef,m_markedCrs); }


    void flatten(const index_t level);
    void flatten() { flatten(m_maxLvl); } ;

    void unrefineThreshold(const index_t level);
    void unrefineThreshold(){ unrefineThreshold(m_maxLvl); };

private:
    void _refineMarkedElements(     const patchHContainer & container,
                                    index_t refExtension = 0);

    void _unrefineMarkedElements(   const patchHContainer & container,
                                    index_t refExtension = 0);

    void _processMarkedElements(const patchHContainer & elRefined,
                                const patchHContainer & elCoarsened,
                                index_t refExtension = 0,
                                index_t crsExtension = 0);

    void _flattenElementsToLevel(   const index_t level);

    void _unrefineElementsThreshold(const index_t level);

    void _markElements( const std::vector<T> & elError, int refCriterion, T refParameter, index_t maxLevel, patchHContainer & container, bool coarsen=false);
    void _markFraction( const std::vector<T> & elError, T refParameter, index_t maxLevel, std::vector<index_t> & elLevels, std::vector<bool> & elMarked, bool coarsen=false);
    void _markPercentage( const std::vector<T> & elError, T refParameter, index_t maxLevel, std::vector<index_t> & elLevels, std::vector<bool> & elMarked, bool coarsen=false);
    void _markThreshold( const std::vector<T> & elError, T refParameter, index_t maxLevel, std::vector<index_t> & elLevels, std::vector<bool> & elMarked, bool coarsen=false);

    void _markLevelThreshold( index_t level, std::vector<bool> & elMarked);
    void _getElLevels( std::vector<index_t> & elLevels);

    patchHContainer _toContainer(const std::vector<bool> & bools);


protected:
    // M & m_basis;
    gsFunctionSet<T> * m_input;
    // const gsMultiPatch<T> & m_patches;
    gsOptionList m_options;

    T               m_crsParam, m_refParam;
    MarkingStrategy m_crsRule, m_refRule;
    index_t         m_crsExt, m_refExt;
    index_t         m_maxLvl;

    index_t m_m;

    bool            m_admissible;

    patchHContainer m_markedRef, m_markedCrs;

};


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsAdaptiveMeshing.hpp)
#endif