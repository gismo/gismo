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
        options.addReal("RefineParam","Parameter used for refinement",0.1);
        options.addInt("MaxLevel","Maximum refinement level",3);
        options.addInt("Admissibility","Admissibility region, 0=T-admissibility (default), 1=H-admissibility",0);
        options.addSwitch("Admissible","Mark the admissible region",true);
        options.addInt("Jump","Jump parameter m",2);
        options.addSwitch("Absolute","(For GARU marking) Compute threshold based on solution values without error scaling",false); // Computed threshold based on the actual values of the solution --> true if marking is done with absolute error, false if relative error is used
        options.addSwitch("Extension", "Extend marked elements regions", true);
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
        switch (m_options.askInt("RefineRule",1))
        {
            case 1: // GARU
                return _markRef_threshold();
            case 2: // PUCA
                return _markRef_percentage();
            case 3: // BULK
                return _markRef_fraction();
            default:
                GISMO_ERROR("Unknown refinement rule.");
        }
    }

    template <short_t d, class T>
    typename gsHElementMarker<d,T>::HElementContainer gsHElementMarker<d,T>::markCrs(const HElementContainer refined) const
    {
        switch (m_options.askInt("RefineRule",1))
        {
            case 1: // GARU
                return _markCrs_threshold(refined);
            case 2: // PUCA
                return _markCrs_percentage(refined);
            case 3: // BULK
                return _markCrs_fraction(refined);
            default:
                GISMO_ERROR("Unknown refinement rule.");
        }
    }

    template <short_t d, class T>
    const gsHElementHelper<d,T> & gsHElementMarker<d,T>::helper()
    {
        return m_helper;
    }

    template <short_t d, class T>
    std::vector<index_t> gsHElementMarker<d,T>::toRefBoxes(const HElementContainer & elements) const
    {
        return m_helper.toRefBoxes(elements, m_options.askSwitch("Extension", true));
    }

    template <short_t d, class T>
    std::vector<index_t> gsHElementMarker<d,T>::toCrsBoxes(const HElementContainer & elements) const
    {
        return m_helper.toCrsBoxes(elements);
    }

    template <short_t d, class T>
    typename gsHElementMarker<d,T>::HElementContainer gsHElementMarker<d,T>::_markRef_admissible(const HElementContainer & refined) const
    {
        HElementContainer result;
        // Mark the admissible region
        result = m_helper.markAdmissible(refined, m_options.askInt("Jump",2));
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

                // If the parents of any cells in the extension overlap with marked refinement cells, coarsening causes a problem
                for ( const auto & refElem : refined )
                {
                    if (m_helper.contains(coarseningElem, refElem))
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
        T maxError = m_options.getSwitch("Absolute") ? 1 : m_elementErrors.back().second;
        T threshold = m_options.askReal("RefineParam",0.1) * maxError;
        for (typename std::vector<std::pair<element_t, error_t>>::const_reverse_iterator it = m_elementErrors.rbegin(); it != m_elementErrors.rend(); ++it)
        {
            // If the error is below the threshold, stop the iteration
            if (it->second < threshold)
                break;

            // If the level of the element is larger than the maximum level, skip it
            if (it->first.level() >= (level_t)m_options.askInt("MaxLevel",10))
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
        T maxError = m_options.getSwitch("Absolute") ? 1 : m_elementErrors.back().second;
        T threshold = m_options.askReal("CoarsenParam",0.1) * maxError;
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
        index_t numToMark = static_cast<index_t>(math::floor(percentage * numElements));
        index_t numMarked = 0;
        for (typename std::vector<std::pair<element_t, error_t>>::const_reverse_iterator it = m_elementErrors.rbegin(); it != m_elementErrors.rend(); ++it, ++numMarked)
        {
            // If we have marked enough elements, stop the iteration
            if (numMarked >= numToMark)
                break;

            // If the level of the element is larger than the maximum level, skip it
            if (it->first.level() >= (level_t)m_options.askInt("MaxLevel",10))
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
        index_t numToMark = static_cast<index_t>(math::floor(percentage * numElements));
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
            if (it->first.level() >= (level_t)m_options.askInt("MaxLevel",10))
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
