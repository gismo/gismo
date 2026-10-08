/** @file gsHElementMarker.hpp

    @brief Provides a class for marking hierarchical elements in a mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#pragma once

#include <gsHSplines/gsHElementMarker.h>

namespace gismo
{

    template <short_t d, class T>
    bool gsHElementMarker<d,T>::CompareElementErrorPair::operator()(const std::pair<element_t, error_t> & a, const std::pair<element_t, error_t> & b) const
    {
        CompareElement comp;
        return a.second < b.second || (a.second == b.second && comp(a.first,b.first));
    };

    template <short_t d, class T>
    gsHElementMarker<d,T>::gsHElementMarker(const gsBasis<T> & basis, gsOptionList options)
    :
    gsHElementMarker(dynamic_cast<const gsHTensorBasis<d,T> &>(basis), give(options))
    {
    }

    template <short_t d, class T>
    gsHElementMarker<d,T>::gsHElementMarker(const gsHTensorBasis<d,T> & basis, gsOptionList options)
    :
    m_basis(basis),
    m_helper(basis),
    m_options(options)
    {
    }

    template <short_t d, class T>
    gsOptionList & gsHElementMarker<d,T>::options()
    {
        return m_options;
    }

    template <short_t d, class T>
    gsOptionList gsHElementMarker<d,T>::defaultOptions()
    {
        gsOptionList options;
        options.addInt("CoarsenRule","Rule used for coarsening: 1=GARU, 2=PUCA, 3=BULK.",1);
        options.addInt("RefineRule","Rule used for refinement: 1=GARU, 2=PUCA, 3=BULK.",1);
        options.addReal("CoarsenParam","Parameter used for coarsening",0.1);
        options.addInt("CoarsenGroupRule","Which children of a sibling group must have a small error for the group to be coarsened "
                       "(only with CoarsenRule=1, GARU): 0=any child (one child with error <= CoarsenParam*max suffices), "
                       "1=all children (max child error <= CoarsenParam*max), "
                       "2=summed (sqrt(sum of squared child errors) <= CoarsenParam*max; for norm-type indicators this is the error "
                       "of the current solution over the parent's area, not a prediction of the error after coarsening). "
                       "Rules 1 and 2 require all 2^d children to be active leaves and return all 2^d children of a candidate group; "
                       "rule 0 does not check activity.",0);
        options.addReal("RefineParam","Parameter used for refinement",0.1);
        options.addInt("MaxLevel","Maximum refinement level",3);
        options.addInt("Admissibility","Admissibility region, 0=T-admissibility (default), 1=H-admissibility",0);
        options.addSwitch("Admissible","Mark the admissible region",true);
        options.addInt("Jump","Jump parameter m",2);
        options.addSwitch("Extension", "Extend the refinement box of each marked element by floor(p/2) "
                          "finer-level spans on both sides; elements added by the admissible closure "
                          "are refined without extension", true);
        // options.addInt("Verbose","Verbosity level",0);
        return options;
    }

    template <short_t d, class T>
    void gsHElementMarker<d,T>::setErrors(const std::vector<error_t> & errors)
    {
        GISMO_ASSERT(errors.size() == m_basis.numElements(),
                        "The size of the errors vector must match the number of elements in the basis domain.");

        m_elementErrors.resize(m_basis.numElements());
        element_t element;
        level_t elementLevel;
        for (auto it = m_basis.domain()->beginAll(); it!= m_basis.domain()->endAll(); ++it)
        {
            // Create a pair of element and its associated error
            elementLevel = static_cast<const gsHDomainIterator<T,d> *>(it.get())->getLevel();
            element = m_helper.toElement(it.lowerCorner(), it.upperCorner(), elementLevel, it.patch());
            // m_elementErrors[elem] = it.id();
            m_elementErrors[it.id()] = std::make_pair(element, errors[it.id()]);
        }
        // Sort the element errors based on the error value in non-decreasing order
        std::stable_sort(m_elementErrors.begin(), m_elementErrors.end(),CompareElementErrorPair());
    }

    template <short_t d, class T>
    typename gsHElementMarker<d,T>::HElementContainer gsHElementMarker<d,T>::markRef() const

    {
        m_closureAdded.clear();
        m_lastRef.clear();
        HElementContainer result;
        switch (m_options.askInt("RefineRule",1))
        {
            case 1: // GARU
                result = _markRef_threshold();
                break;
            case 2: // PUCA
                result = _markRef_percentage();
                break;
            case 3: // BULK
                result = _markRef_fraction();
                break;
            default:
                GISMO_ERROR("Unknown refinement rule.");
        }
        m_lastRef = result;
        return result;
    }

    template <short_t d, class T>
    typename gsHElementMarker<d,T>::HElementContainer gsHElementMarker<d,T>::markCrs(const HElementContainer refined) const
    {
        // The extended refinement boxes of the seeds refine parts of the
        // neighbouring same-level cells; the admissible coarsening checks (and
        // the sibling exclusion) must treat those cells as refined too.
        const RefPlan plan = _refPlan(refined);
        HElementContainer refinedExt = plan.refined;
        if (m_options.askSwitch("Extension", true) && !plan.seeds.empty())
        {
            const HElementContainer cells = m_helper.extensionCells(plan.seeds);
            refinedExt.insert(cells.begin(), cells.end());
        }
        const index_t groupRule = m_options.askInt("CoarsenGroupRule",0);
        GISMO_ENSURE(groupRule >= 0 && groupRule <= 2,
                     "gsHElementMarker::markCrs: unknown CoarsenGroupRule " << groupRule << ", expected 0, 1 or 2");
        GISMO_ENSURE(groupRule == 0 || m_options.askInt("CoarsenRule",1) == 1,
                     "CoarsenGroupRule=" << groupRule << " requires CoarsenRule=1 (GARU); PUCA/BULK select by sorted position, where a group rule is undefined");
        HElementContainer result;
        switch (m_options.askInt("CoarsenRule",1))
        {
            case 1: // GARU
                result = _markCrs_threshold(refinedExt);
                break;
            case 2: // PUCA
                result = _markCrs_percentage(refinedExt);
                break;
            case 3: // BULK
                result = _markCrs_fraction(refinedExt);
                break;
            default:
                GISMO_ERROR("Unknown coarsening rule.");
        }
        if (refined.empty() || result.empty())
            return result;

        // The refinement boxes of the plan (toRefBoxes(refined)) reach into
        // cells that are not in \a refined. Coarsening a candidate (to level crsBox[0] = level-1) conflicts with
        // a refinement box (target level refBox[0]) exactly when the boxes
        // overlap and the coarsening target lies below the refinement target:
        // then the region would end below the level the refinement asks for.
        // If crsBox[0] >= refBox[0] both requests are met, so the candidate is
        // kept. Boxes are compared on the finest grid involved; touching boxes
        // do not overlap. Complexity O(|result| * |plan.refined|).
        const std::vector<index_t> refBoxes = _planBoxes(plan);
        const size_t stride = 2*d + 1;
        HElementContainer kept;
        for (const auto & elem : result)
        {
            const std::vector<index_t> crsBox = m_helper.toCrsBox(elem);
            bool overlap = false;
            for (size_t b = 0; b < refBoxes.size() && !overlap; b += stride)
            {
                const level_t lc = static_cast<level_t>(crsBox[0]);
                const level_t lr = static_cast<level_t>(refBoxes[b]);
                if (lc >= lr)
                    continue;
                const level_t lf = std::max(lc, lr);
                const index_t sc = index_t(1) << (lf - lc);
                const index_t sr = index_t(1) << (lf - lr);
                overlap = true;
                for (short_t i = 0; i < d && overlap; ++i)
                {
                    const index_t lo = std::max(crsBox[1+i]*sc, refBoxes[b+1+i]*sr);
                    const index_t hi = std::min(crsBox[1+d+i]*sc, refBoxes[b+1+d+i]*sr);
                    overlap = lo < hi;
                }
            }
            if (!overlap)
                kept.insert(elem);
        }
        return kept;
    }

    template <short_t d, class T>
    const gsHElementHelper<d,T> & gsHElementMarker<d,T>::helper()
    {
        return m_helper;
    }

    template <short_t d, class T>
    std::vector<index_t> gsHElementMarker<d,T>::toRefBoxes(const HElementContainer & elements) const
    {
        return _planBoxes(_refPlan(elements));
    }

    template <short_t d, class T>
    typename gsHElementMarker<d,T>::RefPlan gsHElementMarker<d,T>::_refPlan(const HElementContainer & elements) const
    {
        RefPlan plan;
        if (elements.empty())
            return plan;

        if (!m_options.askSwitch("Extension", true) || !m_options.askSwitch("Admissible", true))
        {
            plan.seeds = elements;
            plan.refined = elements;
        }
        else if (elements == m_lastRef)
        {
            // The container built by the last markRef(): the closure is
            // already contained in it, and m_closureAdded tells which of
            // its elements are not seeds.
            plan.refined = elements;
            for (const auto & elem : elements)
                (m_closureAdded.count(elem) ? plan.plain : plan.seeds).insert(elem);
        }
        else
        {
            plan.seeds = elements;
            plan.refined = _extendedClosure(elements);
            for (const auto & elem : plan.refined)
                if (!elements.count(elem))
                    plan.plain.insert(elem);
        }
        return plan;
    }

    template <short_t d, class T>
    std::vector<index_t> gsHElementMarker<d,T>::_planBoxes(const RefPlan & plan) const
    {
        const bool extension = m_options.askSwitch("Extension", true);
        if (!extension)
            return m_helper.toRefBoxes(plan.refined, false);

        // Boxes follow the (set) order of plan.refined; only seeds are extended.
        std::vector<index_t> refBoxes;
        refBoxes.reserve(plan.refined.size() * (2 * d + 1));
        for (const auto & elem : plan.refined)
        {
            std::vector<index_t> refBox = m_helper.toRefBox(elem, plan.seeds.count(elem) != 0);
            refBoxes.insert(refBoxes.end(), refBox.begin(), refBox.end());
        }
        return refBoxes;
    }

    template <short_t d, class T>
    std::vector<index_t> gsHElementMarker<d,T>::toCrsBoxes(const HElementContainer & elements) const
    {
        return m_helper.toCrsBoxes(elements);
    }

    template <short_t d, class T>
    typename gsHElementMarker<d,T>::HElementContainer gsHElementMarker<d,T>::_extendedClosure(const HElementContainer & seeds) const
    {
        const level_t jump = m_options.askInt("Jump",2);
        if (!m_options.askSwitch("Extension", true))
            return m_helper.markAdmissible(seeds, jump);

        // The extension refines parts of the neighbouring cells too, so
        // the closure is taken over the seeds AND those cells. The cells
        // themselves are only partly refined (through the extended boxes of the
        // seeds) and are dropped again, unless the closure of the seeds
        // alone needs them refined in full.
        const HElementContainer cells = m_helper.extensionCells(seeds);
        HElementContainer input = seeds;
        input.insert(cells.begin(), cells.end());
        HElementContainer result = m_helper.markAdmissible(input, jump);
        const HElementContainer plain = m_helper.markAdmissible(seeds, jump);
        for (const auto & cell : cells)
            if (!plain.count(cell))
                result.erase(cell);
        return result;
    }

    template <short_t d, class T>
    typename gsHElementMarker<d,T>::HElementContainer gsHElementMarker<d,T>::_markRef_admissible(const HElementContainer & refined) const
    {
        HElementContainer result = _extendedClosure(refined);
        m_closureAdded.clear();
        if (m_options.askSwitch("Extension", true))
            for (const auto & elem : result)
                if (!refined.count(elem))
                    m_closureAdded.insert(elem);
        return result;
    }

    template <short_t d, class T>
    typename gsHElementMarker<d,T>::HElementContainer gsHElementMarker<d,T>::_markCrs_admissible(const HElementContainer & coarsened, const HElementContainer & refined) const
    {
        HElementContainer result = coarsened;
        HElementContainer toErase;
        // Loop over the marked elements
        for (const auto & elem : result)
        {
            bool erase = false;
            // Compute the coarsening extension
            element_t parent = m_helper.getParent(elem);
            HElementContainer coarseningExtension = m_helper.getCextension(parent, m_options.askInt("Jump",2));
            // For all elements in the coarsening extension, check if the level is the same as in the basis or finer
            for (const auto & coarseningElem : coarseningExtension)
            {
                if (m_helper.levelInBasis(coarseningElem) >= coarseningElem.level())
                {
                    // If the level is the same or finer, the coarsening element is not admissible
                    erase = true;
                    break;
                }

                // Admissible coarsening needs N_c U N_rc = empty: the coarsening
                // neighborhood (Carraturo et al., CMAME 348 (2019), Def. 3.5) and
                // its counterpart over the elements marked for refinement (Verhelst
                // et al., Eng. Comput. 40 (2024), Eq. (44)), at level l for m = 2.
                for ( const auto & refElem : refined )
                {
                    if (refElem.level() >= elem.level() &&
                        (m_helper.contains(coarseningElem, refElem) ||
                         m_helper.contains(refElem, coarseningElem)))
                    {
                        // If the coarsening element contains a refinement element, the coarsening is not admissible
                        erase = true;
                        break;
                    }
                    if (erase)
                        break;
                }

            }

            if (erase)
            {
                // If the element is not admissible, erase it from the result
                toErase.insert(elem);
                continue;
            }
        }

        // Erase the non-admissible elements from the result
        for (const auto & elem : toErase)
        {
            result.erase(elem);
        }
        return result;
    }

    template <short_t d, class T>
    typename gsHElementMarker<d,T>::HElementContainer gsHElementMarker<d,T>::_markRef_threshold() const
    {
        HElementContainer result;
        T threshold = m_options.askReal("RefineParam",0.1) * m_elementErrors.back().second;
        for (typename std::vector<std::pair<element_t, error_t>>::const_reverse_iterator it = m_elementErrors.rbegin(); it != m_elementErrors.rend(); ++it)
        {
            // If the error is below the threshold, stop the iteration
            if (it->second < threshold)
                break;

            // If the level of the element is larger than the maximum level, skip it
            if (it->first.level() >= (level_t)m_options.askInt("MaxLevel",3))
                continue;

            // Add the element to the container
            result.insert(it->first);
        }
        return (m_options.askSwitch("Admissible",true) ? _markRef_admissible(result) : result);
    }

    template <short_t d, class T>
    typename gsHElementMarker<d,T>::HElementContainer gsHElementMarker<d,T>::_markCrs_threshold(const HElementContainer & refined) const
    {
        HElementContainer result;
        T threshold = m_options.askReal("CoarsenParam",0.1) * m_elementErrors.back().second;
        const index_t groupRule = m_options.askInt("CoarsenGroupRule",0);
        if (groupRule != 0)
        {
            // m_elementErrors holds exactly the active leaves, so membership in errMap is the activity test
            std::map<element_t, error_t, typename element_t::Compare> errMap(m_elementErrors.begin(), m_elementErrors.end());
            HElementContainer decided; // parents already evaluated
            for (typename std::vector<std::pair<element_t, error_t>>::const_iterator it = m_elementErrors.begin(); it != m_elementErrors.end(); ++it)
            {
                // Every child of a candidate group has error <= its statistic <= threshold (the max, resp. the
                // root of the sum of squares, bounds each child error), so all children lie in the ascending
                // prefix before the first error above the threshold and the group is reached through its smallest child.
                if (it->second > threshold)
                    break;
                if (it->first.level() == 0)
                    continue;
                if (!decided.insert(m_helper.getParent(it->first)).second)
                    continue;

                const HElementContainer siblings = m_helper.getSiblings(it->first,false);
                if (!refined.empty() &&
                    std::any_of(siblings.begin(), siblings.end(),[&refined](const element_t & elem) { return refined.find(elem) != refined.end(); }))
                    continue;

                error_t stat = 0;
                bool active = true;
                for (typename HElementContainer::const_iterator s = siblings.begin(); s != siblings.end(); ++s)
                {
                    typename std::map<element_t, error_t, typename element_t::Compare>::const_iterator f = errMap.find(*s);
                    if (f == errMap.end())
                    {
                        active = false;
                        break;
                    }
                    stat = (groupRule == 1 ? (std::max)(stat, f->second) : stat + f->second * f->second);
                }
                if (!active)
                    continue;
                if (groupRule == 2)
                    stat = math::sqrt(stat);

                if (stat <= threshold)
                    result.insert(siblings.begin(), siblings.end());
            }
            return (m_options.askSwitch("Admissible",true) ? _markCrs_admissible(result,refined) : result);
        }

        for (typename std::vector<std::pair<element_t, error_t>>::const_iterator it = m_elementErrors.begin(); it != m_elementErrors.end(); ++it)
        {
            // If the error is above the threshold, stop the iteration
            if (it->second > threshold)
                break;

            // If the level of the element is zero, skip it
            if (it->first.level() == 0)
                continue;

            if (!refined.empty())
            {
                // If the element or its siblings are in the refined set, skip it
                HElementContainer siblings = m_helper.getSiblings(it->first,false);
                if (std::any_of(siblings.begin(), siblings.end(),[&refined](const element_t & elem) { return refined.find(elem) != refined.end(); }))
                    continue;
                // If any of the siblings is not active, skip it
                if (std::any_of(siblings.begin(), siblings.end(),[&it](const element_t & elem) { return it->first.level() != elem.level(); }))
                    continue;
            }

            // Add the element to the container
            result.insert(it->first);
        }
        return (m_options.askSwitch("Admissible",true) ? _markCrs_admissible(result,refined) : result);
    }

    template <short_t d, class T>
    typename gsHElementMarker<d,T>::HElementContainer gsHElementMarker<d,T>::_markRef_percentage() const
    {
        HElementContainer result;
        T percentage = m_options.askReal("RefineParam",0.1);
        index_t numElements = m_elementErrors.size();
        index_t numToMark = cast<T,index_t>(math::floor(percentage * numElements));
        index_t numMarked = 0;
        for (typename std::vector<std::pair<element_t, error_t>>::const_reverse_iterator it = m_elementErrors.rbegin(); it != m_elementErrors.rend(); ++it, ++numMarked)
        {
            // If we have marked enough elements, stop the iteration
            if (numMarked >= numToMark)
                break;

            // If the level of the element is larger than the maximum level, skip it
            if (it->first.level() >= (level_t)m_options.askInt("MaxLevel",3))
                continue;

            // Add the element to the container
            result.insert(it->first);
        }
        return (m_options.askSwitch("Admissible",true) ? _markRef_admissible(result) : result);
    }

    template <short_t d, class T>
    typename gsHElementMarker<d,T>::HElementContainer gsHElementMarker<d,T>::_markCrs_percentage(const HElementContainer & refined) const
    {
        HElementContainer result;
        T percentage = m_options.askReal("CoarsenParam",0.1);
        index_t numElements = m_elementErrors.size();
        index_t numToMark = cast<T,index_t>(math::floor(percentage * numElements));
        index_t numMarked = 0;
        for (typename std::vector<std::pair<element_t, error_t>>::const_iterator it = m_elementErrors.begin(); it != m_elementErrors.end(); ++it, ++numMarked)
        {
            // If we have marked enough elements, stop the iteration
            if (numMarked >= numToMark)
                break;

            // If the level of the element is zero, skip it
            if (it->first.level() == 0)
                continue;

            // If the element or its siblings are in the refined set, skip it
            if (!refined.empty())
            {
                // If the element or its siblings are in the refined set, skip it
                HElementContainer siblings = m_helper.getSiblings(it->first,false);
                if (std::any_of(siblings.begin(), siblings.end(),[&refined](const element_t & elem) { return refined.find(elem) != refined.end(); }))
                    continue;
                // If any of the siblings is not active, skip it
                if (std::any_of(siblings.begin(), siblings.end(),[&it](const element_t & elem) { return it->first.level() < elem.level(); }))
                    continue;
            }

            // Add the element to the container
            result.insert(it->first);
        }
        return (m_options.askSwitch("Admissible",true) ? _markCrs_admissible(result,refined) : result);
    }

    template <short_t d, class T>
    typename gsHElementMarker<d,T>::HElementContainer gsHElementMarker<d,T>::_markRef_fraction() const
    {
        HElementContainer result;
        // Compute the total error
        T cummulErrMarked = T(0);
        T totalError = std::accumulate(m_elementErrors.begin(), m_elementErrors.end(), T(0),[](T sum, const std::pair<element_t, error_t> & elem) { return sum + elem.second; });
        T errorMarkSum = m_options.askReal("RefineParam",0.1) * totalError;
        for (typename std::vector<std::pair<element_t, error_t>>::const_reverse_iterator it = m_elementErrors.rbegin(); it != m_elementErrors.rend(); ++it)
        {
            // If the cumulative error exceeds the threshold, stop the iteration
            if (cummulErrMarked >= errorMarkSum)
                break;

            // If the level of the element is larger than the maximum level, skip it
            if (it->first.level() >= (level_t)m_options.askInt("MaxLevel",3))
                continue;

            // Add the element to the container
            result.insert(it->first);
            cummulErrMarked += it->second;
        }
        return (m_options.askSwitch("Admissible",true) ? _markRef_admissible(result) : result);
    }

    template <short_t d, class T>
    typename gsHElementMarker<d,T>::HElementContainer gsHElementMarker<d,T>::_markCrs_fraction(const HElementContainer & refined) const
    {
        HElementContainer result;
        // Compute the total error
        T cummulErrMarked = T(0);
        T totalError = std::accumulate(m_elementErrors.begin(), m_elementErrors.end(), T(0),[](T sum, const std::pair<element_t, error_t> & elem) { return sum + elem.second; });
        T errorMarkSum = m_options.askReal("CoarsenParam",0.1) * totalError;
        for (typename std::vector<std::pair<element_t, error_t>>::const_iterator it = m_elementErrors.begin(); it != m_elementErrors.end(); ++it)
        {
            // If the cumulative error exceeds the threshold, stop the iteration
            if (cummulErrMarked >= errorMarkSum)
                break;

            // If the level of the element is zero, skip it
            if (it->first.level() == 0)
                continue;

            // If there are no refinement elements provided, we don't check siblings
            if (!refined.empty())
            {
                // If the element or its siblings are in the refined set, skip it
                HElementContainer siblings = m_helper.getSiblings(it->first,false);
                if (std::any_of(siblings.begin(), siblings.end(),[&refined](const element_t & elem) { return refined.find(elem) != refined.end(); }))
                    continue;
                // If any of the siblings is not active, skip it
                if (std::any_of(siblings.begin(), siblings.end(),[&it](const element_t & elem) { return it->first.level() < elem.level(); }))
                    continue;
            }

            // Add the element to the container
            result.insert(it->first);
            cummulErrMarked += it->second;
        }
        return (m_options.askSwitch("Admissible",true) ? _markCrs_admissible(result,refined) : result);
    }

    template <short_t d, class T>
    std::ostream & gsHElementMarker<d,T>::print(std::ostream & os) const
    {
        // Print m_elementErrors
        os << "Element errors:\n";
        for (const auto & elemErr : m_elementErrors)
        {
            os << "Element: " << elemErr.first << ", Error: " << elemErr.second << "\n";
        }
        // Print the options
        os << "Options:\n";
        os << m_options<<"\n";
        // Print the basis
        os << "Basis:\n";
        os << m_basis;
        return os;
    }

}; //namespace gismo
