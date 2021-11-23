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

    void mark(gsFunctionSet<T> * input, const std::vector<T> & errors, index_t level);
    void mark(const std::vector<T> & errors) { mark(m_input,errors,m_maxLvl); }
    void mark(const std::vector<T> & errors, index_t maxLvl) { mark(m_input,errors,maxLvl); };

    void refine(const std::vector<bool> & markedRef);

    void unrefine(const std::vector<bool> & markedCrs);

    void adapt(const std::vector<bool> & markedRef,const std::vector<bool> & markedCrs);

    void refine(){ refine(m_markedRef); }

    void unrefine() { unrefine(m_markedRef); }

    void adapt() { adapt(m_markedRef,m_markedCrs); }

    void flatten(const index_t level);
    void flatten() { flatten(m_maxLvl); } ;

    void unrefineThreshold(const index_t level);
    void unrefineThreshold(){ unrefineThreshold(m_maxLvl); };

private:
    void _refineMarkedElements( gsFunctionSet<T> * bases,
                                const std::vector<bool> & elMarked,
                                index_t refExtension = 0);

    void _unrefineMarkedElements(gsFunctionSet<T> * bases,
                                    const std::vector<bool> & elMarked,
                                    index_t refExtension = 0);

    void _processMarkedElements(gsFunctionSet<T> * bases,
                                const std::vector<bool> & elRefined,
                                const std::vector<bool> & elCoarsened,
                                index_t refExtension = 0,
                                index_t crsExtension = 0);

    void _flattenElementsToLevel(  gsFunctionSet<T> * bases,
                            const index_t level);

    void _unrefineElementsThreshold(  gsFunctionSet<T> * bases,
                            const index_t level);

    void _markElements( gsFunctionSet<T> * input, const std::vector<T> & elError, int refCriterion, T refParameter, index_t maxLevel, std::vector<bool> & elMarked, bool coarsen=false);
    void _markFraction( const std::vector<T> & elError, T refParameter, index_t maxLevel, std::vector<index_t> & elLevels, std::vector<bool> & elMarked, bool coarsen=false);
    void _markPercentage( const std::vector<T> & elError, T refParameter, index_t maxLevel, std::vector<index_t> & elLevels, std::vector<bool> & elMarked, bool coarsen=false);
    void _markThreshold( const std::vector<T> & elError, T refParameter, index_t maxLevel, std::vector<index_t> & elLevels, std::vector<bool> & elMarked, bool coarsen=false);

    void _markLevelThreshold( gsFunctionSet<T> * input, index_t level, std::vector<bool> & elMarked);
    void _getElLevels( gsFunctionSet<T> * input, std::vector<index_t> & elLevels);

protected:
    // M & m_basis;
    gsFunctionSet<T> * m_input;
    // const gsMultiPatch<T> & m_patches;
    gsOptionList m_options;

    T               m_crsParam, m_refParam;
    MarkingStrategy m_crsRule, m_refRule;
    index_t         m_crsExt, m_refExt;
    index_t         m_maxLvl;

    bool            m_admissible;


    std::vector<bool> m_markedRef, m_markedCrs;

};


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsAdaptiveMeshing.hpp)
#endif