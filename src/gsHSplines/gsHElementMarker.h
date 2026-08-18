/** @file gsHElementMarker.h

    @brief Provides a class for marking hierarchical elements in a mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#pragma once

#include <gsIO/gsOptionList.h>
#include <gsHSplines/gsHElement.h>
#include <gsHSplines/gsHElementHelper.h>
#include <gsHSplines/gsHTensorBasis.h>

namespace gismo
{

template <short_t d, class T>
class gsHElementMarker
/**
 * @brief Class for marking elements in hierarchical tensor bases for refinement or coarsening based on error indicators.
 *
 * This class provides functionality to mark elements in a hierarchical tensor basis for either
 * refinement or coarsening based on various criteria applied to error indicators associated with
 * the elements. The marking strategies include threshold-based, percentage-based, and fraction-based
 * approaches, with optional enforcement of admissibility constraints.
 *
 * @tparam d Dimension of the parameter domain
 * @tparam T Coefficient type
 */
{

public:
    typedef          gsHElement<d,T>                                    element_t;
    typedef typename gsHElement<d,T>::level_t                           level_t;
    typedef typename gsHElement<d,T>::Compare                           CompareElement;
    typedef typename std::set<element_t,typename element_t::Compare>    HElementContainer;
    typedef          T                                                  error_t;
    typedef typename std::vector<std::pair<element_t, error_t>>         HElementErrorContainer;

public:
    struct CompareElementErrorPair
    {
        bool operator()(const std::pair<element_t, error_t> & a, const std::pair<element_t, error_t> & b) const;
    };

public:

    /// @brief  Constructor
    /// @param basis The basis to use for element marking
    /// @param options Options for the marker, including rules for refinement and coarsening
    gsHElementMarker(const gsBasis<T> & basis, gsOptionList options = defaultOptions());

    /// @brief  Constructor
    /// @param basis The basis to use for element marking
    /// @param options Options for the marker, including rules for refinement and coarsening
    gsHElementMarker(const gsHTensorBasis<d,T> & basis, gsOptionList options = defaultOptions());

    /// @brief  Accessor for options
    /// @return Reference to the options used by the marker
    /// @example
    ///     marker.options().setInt("RefineRule", 1);
    gsOptionList & options();

    /// @brief  Default options for the marker
    /// @return A gsOptionList containing the default options for the marker
    /// @example
    ///     marker.options() = gsHElementMarker<d,T>::defaultOptions();
    static gsOptionList defaultOptions();

    /// @brief  Set the errors associated with the elements in the basis domain
    /// @param errors A vector of error indicators corresponding to the elements in the basis domain
    /// @note The size of the errors vector must match the number of elements in the basis domain.
    ///       The order of the errors must match the order of elements in the basis domain (obtained from domain iteration).
    void setErrors(const std::vector<error_t> & errors);

    /// @brief  Mark elements for refinement or coarsening based on the specified rules
    /// @return A container of elements marked for refinement or coarsening
    /// @note If the "Admissible" option is set to true, only elements that are admissible will be marked for refinement.
    /// @note The marking strategy is determined by the "RefineRule" option in the options list.
    ///       The available rules are:
    ///       - 1: GARU (greatest appearing eRror utilization)
    ///       - 2: PUCA (percentile-utilizing cutoff ascertainment)
    ///       - 3: BULK ("Doerfler-marking")
    /// @example
    ///     auto markedElements = marker.markRef();
    HElementContainer markRef() const;

    /// @brief  Mark elements for coarsening based on the specified rules
    /// @param refined A container of elements that have already been refined. If this is provided, elements marked for refinement will not be considered for coarsening.
    /// @return A container of elements marked for coarsening
    /// @note If the "Admissible" option is set to true, only elements that are admissible will be marked for coarsening.
    /// @note The marking strategy is determined by the "CoarsenRule" option in the options list.
    ///       The available rules are:
    ///       - 1: GARU (greatest appearing eRror utilization)
    ///       - 2: PUCA (percentile-utilizing cutoff ascertainment)
    ///       - 3: BULK ("Doerfler-marking")
    /// @note Behaviour change: until this fix, markCrs() dispatched on the
    ///       "RefineRule" option instead of "CoarsenRule". Callers that relied
    ///       on the old behaviour (for example by installing the coarsening
    ///       rule into "RefineRule" before calling markCrs()) will now select a
    ///       different coarsening strategy whenever CoarsenRule != RefineRule.
    ///       Results are unchanged when the two options are equal.
    /// @example
    ///     auto markedElements = marker.markCrs(refinedElements);
    HElementContainer markCrs(const HElementContainer refined = {}) const;

    /// @brief  Accessor for the helper class used for element operations
    /// @return A reference to the gsHElementHelper instance associated with this marker
    /// @note This helper class provides various utility functions for working with hierarchical elements.
    /// @example
    ///     auto markedElements = marker.markRef();
    ///     gsMatrix<real_t> boxes = marker.helper().toBoxes(markedElements);
    const gsHElementHelper<d,T> & helper();

    /// @brief Convert a container of elements to a vector of refinement box indices
    /// @param elements A container of elements to convert
    /// @return A vector of indices representing the refinement boxes corresponding to the elements
    /// @note The indices are based on the hierarchical structure of the elements in the basis.
    /// @example
    ///     std::vector<index_t> refBoxes = marker.toRefBoxes(markedElements);
    ///     basis.refineElements(refBoxes);
    std::vector<index_t> toRefBoxes(const HElementContainer & elements) const;

    /// @brief Convert a container of elements to a vector of coarsening box indices
    /// @param elements A container of elements to convert
    /// @return A vector of indices representing the coarsening boxes corresponding to the elements
    /// @note The indices are based on the hierarchical structure of the elements in the basis.
    /// @example
    ///     std::vector<index_t> crsBoxes = marker.toCrsBoxes(elements);
    ///     basis.unrefineElements(crsBoxes);
    std::vector<index_t> toCrsBoxes(const HElementContainer & elements) const;

private:

    /// @brief Applies admissibility marking to a container of already marked elements
    /// @param refined A container of elements that have already been marked for refinement
    /// @return A container of all elements that need to be refined (including the original marked elements)
    HElementContainer _markRef_admissible(const HElementContainer & refined) const;

    /// @brief Eliminates elements from a container based on admissibility criteria
    /// @param refined A container of elements that have already been marked for refinement
    /// @param coarsened A container of elements that have already been marked for coarsening
    /// @return A subset of \a coarsened that can be coarsened without violating admissibility
    HElementContainer _markCrs_admissible(const HElementContainer & refined, const HElementContainer & coarsened) const;

    /// @brief Marks elements for refinement based on a threshold criterion
    /// @return A container of elements marked for refinement based on a threshold applied to the error
    HElementContainer _markRef_threshold() const;

    /// @brief Marks elements for coarsening based on a threshold criterion
    /// @param refined A container of elements that have already been marked for refinement
    /// @return A container of elements marked for coarsening based on a threshold applied to
    HElementContainer _markCrs_threshold(const HElementContainer & refined) const;

    /// @brief Marks elements for refinement based on a percentage criterion
    /// @return A container of elements marked for refinement based on a percentage of the error
    HElementContainer _markRef_percentage() const;

    /// @brief Marks elements for coarsening based on a percentage criterion
    /// @param refined A container of elements that have already been marked for refinement
    /// @return A container of elements marked for coarsening based on a percentage of the
    HElementContainer _markCrs_percentage(const HElementContainer & refined) const;

    // NOTE: This does not take the extra error contributions due to admissibility into account.
    /// @brief Marks elements for refinement based on a fraction of the error
    /// @return A container of elements marked for refinement based on a fraction of the error
    HElementContainer _markRef_fraction() const;

    /// @brief Marks elements for coarsening based on a fraction of the error
    /// @param refined A container of elements that have already been marked for refinement
    /// @return A container of elements marked for coarsening based on a fraction of the
    HElementContainer _markCrs_fraction(const HElementContainer & refined) const;


public:
    /// @brief Print the marker information
    /// @param os The output stream to print to
    /// @return The output stream after printing the marker information
    std::ostream & print(std::ostream & os) const;

protected:
    const gsHTensorBasis<d,T> & m_basis; ///< The basis of the elements.
    gsHElementHelper<d,T> m_helper; ///< Helper for element operations.
    HElementErrorContainer m_elementErrors; ///< Container for elements and their associated errors.

    gsOptionList m_options; ///< Options for the marker.
};

template<short_t d, class T>
std::ostream& operator<<( std::ostream& os, const gsHElementMarker<d,T>& b )
{
    return b.print( os );
}


}; //namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsHElementMarker.hpp)
#endif