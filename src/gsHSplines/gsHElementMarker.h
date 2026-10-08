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
    /// @note The returned container is stored; toRefBoxes() and markCrs() recognise it
    ///       and treat its closure-added elements differently from its marked elements.
    /// @note Preconditions: the "Admissible" guarantee (the refined mesh stays admissible) holds
    ///       for meshes produced by this marker, i.e. markRef() -> toRefBoxes() ->
    ///       refineElements() starting from a tensor-product mesh, optionally with markCrs()
    ///       coarsening. A mesh refined externally (refineElements() on hand-picked boxes, or
    ///       read from file) can be admissible and still give a non-admissible result.
    /// @example
    ///     auto markedElements = marker.markRef();
    HElementContainer markRef() const;

    /// @brief  Mark elements for coarsening based on the specified rules
    /// @param refined A container of elements that have already been refined. If this is provided, elements marked for refinement will not be considered for coarsening.
    /// @return A container of elements marked for coarsening
    /// @note If the "Admissible" option is set to true, only elements that are admissible will be marked for coarsening.
    ///       The cells that count as refined in that check, and in the exclusion of siblings of
    ///       refined elements, are the elements that toRefBoxes(\a refined) refines
    ///       (see toRefBoxes for stored and foreign sets) plus the extension ring of the
    ///       elements that toRefBoxes extends (the seeds); closure-added elements add no ring.
    /// @note The overlap filter compares each returned element's coarsening box (toCrsBoxes) with
    ///       exactly the boxes toRefBoxes(\a refined) returns, also when "Extension" is off
    ///       (then the boxes are unextended). Where they overlap, the coarsening target level is
    ///       >= the refinement target level, so refineElements(toRefBoxes) followed by
    ///       unrefineElements(toCrsBoxes) leaves every refined region at >= its target level.
    /// @note Call it before the basis is modified; \a refined is interpreted against the basis
    ///       as it is at call time.
    /// @note The marking strategy is determined by the "CoarsenRule" option in the options list.
    ///       The available rules are:
    ///       - 1: GARU (greatest appearing eRror utilization)
    ///       - 2: PUCA (percentile-utilizing cutoff ascertainment)
    ///       - 3: BULK ("Doerfler-marking")
    /// @note The "CoarsenGroupRule" option (default 0) decides, for GARU only, which children of a
    ///       sibling group (the 2^d children of one parent, which are coarsened together) must have
    ///       a small error. With Y = "CoarsenParam" and max the largest element error:
    ///       - 0: any child, a child with error <= Y*max suffices;
    ///       - 1: all children, max_i err(c_i) <= Y*max;
    ///       - 2: summed, sqrt(sum_i err(c_i)^2) <= Y*max. For a norm-type indicator (L2, H1,
    ///         energy) this is the error of the current solution over the parent's area; it is
    ///         not a prediction of the error after coarsening, which over the parent is typically larger (no local guarantee).
    ///
    ///       Rules 1 and 2 only consider groups whose 2^d children are all active leaves, and
    ///       return all 2^d children of a candidate group. A nonzero value together with a
    ///       "CoarsenRule" other than 1, or a value outside {0,1,2}, throws.
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
    /// @note Which elements are extended by floor(p/2) finer-level spans depends on \a elements:
    ///       - the container returned by the last markRef() call ("Extension" on): its marked
    ///         elements are extended, the elements added by its admissible closure are not;
    ///       - any other container ("Extension" and "Admissible" on): its extended admissible
    ///         closure is computed with every element of \a elements a marked element, and that
    ///         closure is refined; the elements of \a elements are extended, the others are not.
    ///         The marker state is not changed;
    ///       - "Admissible" off ("Extension" on): every element is extended;
    ///       - "Extension" off: no element is extended, \a elements is refined as given.
    /// @note Call it before the basis is modified: the closure of a foreign container is computed
    ///       against the basis as it is at call time.
    /// @note The "Admissible" guarantee has the preconditions stated at markRef(): meshes produced by
    ///       this marker; an externally refined mesh can yield a non-admissible result.
    /// @note Complexity: a comparison of \a elements with the stored container, O(n); for a foreign
    ///       container additionally two markAdmissible() calls, the cost of one markRef().
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

    /// @brief Refinement plan of a container of elements, see _refPlan.
    struct RefPlan
    {
        HElementContainer seeds;   ///< Marked elements; refined with the floor(p/2) extension when "Extension" is on.
        HElementContainer plain;   ///< Elements refined without extension.
        HElementContainer refined; ///< All refined elements, seeds and plain.
    };

    /// @brief Maps a container of elements to the elements that toRefBoxes() refines, and the marked elements among them (extended when "Extension" is on).
    /// @param elements The container passed to toRefBoxes() / markCrs()
    /// @return The plan; no member is modified. See toRefBoxes for the four cases.
    RefPlan _refPlan(const HElementContainer & elements) const;

    /// @brief Refinement boxes of a plan, in the order of plan.refined; plan.seeds are extended when "Extension" is on, no element otherwise.
    std::vector<index_t> _planBoxes(const RefPlan & plan) const;

    /// @brief Extended admissible closure of a set of marked elements; modifies no member.
    /// @param seeds Elements marked for refinement, refined with their extension
    /// @return All elements that need to be refined (including \a seeds)
    HElementContainer _extendedClosure(const HElementContainer & seeds) const;

    /// @brief Applies admissibility marking to a container of already marked elements
    /// @param refined A container of elements that have already been marked for refinement
    /// @return A container of all elements that need to be refined (including the original marked elements)
    /// @note With "Extension" on, records the elements that the closure added in m_closureAdded; with it off, m_closureAdded is left empty.
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
    /// @return A container of elements marked for coarsening, with threshold CoarsenParam times the
    ///         largest element error. Which elements are returned depends on "CoarsenGroupRule":
    ///         - 0: every child (level > 0) whose own error is <= threshold, each on its own and not its
    ///           whole sibling group, so a group contributes between 1 and 2^d children; if \a refined is
    ///           non-empty, children with a sibling in \a refined are skipped. The activity of the
    ///           siblings is not checked.
    ///         - 1: the largest child error is <= threshold.
    ///         - 2: sqrt(sum of squared child errors) <= threshold. For a norm-type indicator this is
    ///           the error of the current solution over the parent's area, not a prediction of the
    ///           error after coarsening.
    ///
    ///         Rules 1 and 2 require all 2^d children to be active leaves and return all 2^d children
    ///         of a candidate group. Complexity of rules 1 and 2: O(n log n) to build the lookup
    ///         map plus O(k 2^d log n) for the k elements of the sub-threshold prefix.
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

    /// Elements added by the admissible closure in the last markRef() with
    /// "Extension"; toRefBoxes() refines them without extension when given the
    /// container markRef() returned, since the closure only accounted for the
    /// extended boxes of the marked elements.
    mutable HElementContainer m_closureAdded;

    /// The container returned by the last markRef() call; toRefBoxes() and
    /// markCrs() recognise it (and only it) as already closed under admissibility,
    /// provided markRef() ran with the same "Extension" and "Admissible" options
    /// as are set when toRefBoxes() or markCrs() is called.
    mutable HElementContainer m_lastRef;
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