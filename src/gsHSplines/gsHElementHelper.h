/** @file gsHElementHelper.h

    @brief Provides a class for creating and managing hierarchical elements in a mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#pragma once

#include <gsHSplines/gsElement.h>
#include <gsHSplines/gsElementHelper.h>
#include <gsHSplines/gsHTensorBasis.h>

namespace gismo
{

template <short_t d, class T>
class gsHElementHelper
{
public:
    typedef          gsHElement<d,T>                                    element_t;
    typedef typename gsHElement<d,T>::box_t                             box_t;
    typedef typename gsHElement<d,T>::point_t                           point_t;
    typedef typename gsHElement<d,T>::level_t                           level_t;
    typedef typename std::set<element_t,typename element_t::Compare>    HElementContainer;
    typedef typename gsElementHelper<d,T>::ElementContainer             ElementContainer;
    typedef typename gsHTensorBasis<d,T>::tensorBasis                   TensorBasis_t;

public:

    /// @brief Constructor
    /// @param basis The basis to use for element creation
    explicit gsHElementHelper(const gsBasis<T> & basis);


    /// @brief Constructor
    /// @param basis The basis to use for element creation
    explicit gsHElementHelper(const gsHTensorBasis<d,T> & basis);

public:

    /// @brief  Convert a box defined by parametric lower and upper corners to an element
    /// @param low  Lower corner of the box in parameter coordinates
    /// @param upp  Upper corner of the box in parameter coordinates
    /// @param patch  Patch ID of the element, if applicable
    /// @return An element representing the box in the mesh
    element_t toElement(const gsVector<T,d> & low, const gsVector<T,d> & upp, level_t level, index_t patch = 0) const;

    /// @brief  Convert the basis domain to a set of elements
    /// @param patch  Patch ID of the elements, if applicable
    /// @return A set of elements representing the basis domain
    HElementContainer toElements(index_t patch = 0) const;

    /// @brief Convert a set of elements to a set of hierarchical elements
    /// @param elements The set of elements to convert
    /// @param level The level of the elements in the hierarchy
    /// @return A set of hierarchical elements corresponding to the input elements
    HElementContainer toElements(const ElementContainer & elements, level_t level) const;

    /// @brief Convert an element to a box defined by its lower and upper corners in parameter coordinates
    /// @param element The element to convert
    /// @return A pair containing the element in parameter coordinates and its level
    std::pair<gsMatrix<T>,level_t> toBoxAndLevel(const element_t & element) const;

    /// @brief Convert an element to a box defined by its lower and upper corners in parameter coordinates
    /// @param element The element to convert
    /// @return A box representing the element in parameter coordinates
    gsMatrix<T> toBox(const element_t & element) const;

    /// @brief Convert a set of elements to a set of boxes and their levels
    /// @param elements The set of elements to convert
    /// @return A pair containing the elements in parameter coordinates and their levels
    std::pair<gsMatrix<T>,gsVector<level_t>> toBoxesAndLevels(const HElementContainer & elements) const;

    /// @brief Convert a set of elements to a set of boxes
    /// @param elements The set of elements to convert
    /// @return A matrix where each two columns represent a box defined by its lower and upper corners
    ///         in parameter coordinates
    gsMatrix<T> toBoxes(const HElementContainer & elements) const;

    /// @brief Get the level of the basis at the center of the element
    /// @param element The element for which to get the level
    /// @return The level of the basis at the center of the element
    level_t levelInBasis(const element_t & element) const;

    /// @brief Check if the element is active in the basis
    /// @param element The element to check
    /// @return True if the element is active in the basis, false otherwise
    /// @note An element is considered active if its level in the basis matches its level in the hierarchy
    bool isActive(const element_t & element) const;

    /// Get the parent of the element
    /// @param element The element for which the parent is requested.
    /// @return The parent element.
    /// @note The parent is the element with the same box, but one level lower.
    ///       If the element is at level 0, the parent is undefined.
    element_t getParent(const element_t & element) const;

    /// Get the ancestor of the element \a jump levels up in the hierarchy
    /// @param element The element for which the ancestor is requested.
    /// @param jump The number of levels to jump up in the hierarchy.
    /// @return The ancestor element at the specified jump level.
    element_t getAncestor(const element_t & element, level_t jump) const;

    /// Get the children of the element
    /// @param element The element for which the children are requested.
    /// @return A vector of children elements.
    /// @note The children are the elements with the same box, but one level higher.
    HElementContainer getChildren(const element_t & element) const;

    /// Get the descendants of the element \a jump levels down in the hierarchy
    /// @param element The element for which the descendants are requested.
    /// @param jump The number of levels to jump down in the hierarchy.
    /// @return A vector of descendant elements at the specified jump level.
    HElementContainer getDescendants(const element_t & element, level_t jump) const;

    /// Get the siblings of the element
    /// @param element The element for which the siblings are requested.
    /// @param erase_original If true, the original element is removed from the siblings.
    /// @return A vector of sibling elements.
    HElementContainer getSiblings(const element_t & element, bool erase_original = true) const;

    /// Get the support extension of the element
    /// @param element The element for which the support extension is requested.
    /// @return The support extension of the element as a box in parameter coordinates.
    box_t getSupportExtension(const element_t & element) const;

    /// Get the support extension of the element at a specific jump level
    /// @param element The element for which the support extension is requested.
    /// @param jump The number of levels to jump up in the hierarchy.
    /// @return The support extension of the element at the specified jump level as a box in
    box_t getMultiLevelSupportExtension(const element_t & element, level_t jump) const;

    /// Get the TNeighborhood of the element at a specific level
    /// @param element The element for which the TNeighborhood is requested.
    /// @param m The jump/offset parameter indicating how many levels to go back in the hierarchy.
    /// @return A set of elements representing the HNeighborhood of the element at the specified jump level
    HElementContainer getHNeighborhood(const element_t & element, level_t m) const;
    /// Get the TNeighborhood of the element at a specific jump level
    /// @param element The element for which the TNeighborhood is requested.
    /// @param m The jump/offset parameter indicating how many levels to go back in the hierarchy.
    /// @return A set of elements representing the TNeighborhood of the element at the specified jump level
    HElementContainer getTNeighborhood(const element_t & element, level_t m) const;
    /// Get the neighborhood of the element at a specific jump level
    /// @param element The element for which the neighborhood is requested.
    /// @param m The jump/offset parameter indicating how many levels to go back in the hierarchy.
    /// @return A set of elements representing the neighborhood of the element at the specified jump level
    /// @note If the basis is a tensor basis, the TNeighborhood is returned.
    ///       If the basis is a HBSpline basis, the HNeighborhood is returned.
    ///       If the basis is neither, an error is thrown.
    HElementContainer getNeighborhood(const element_t & element, level_t m) const;

    /// Get the coarsening extension of the element at a specific level
    /// @param element The element for which the coarsening extension is requested.
    /// @param m The level of the coarsening extension.
    /// @return A set of elements representing the coarsening extension of the element at the
    HElementContainer getCextension(const element_t & element, level_t m) const;

    /// Get the coarsening neighborhood of the element at a specific level
    /// @param element The element for which the coarsening neighborhood is requested.
    /// @param m The level of the coarsening neighborhood.
    /// @return A set of elements representing the coarsening neighborhood of the element at the
    ///         specified level
    HElementContainer getCneighborhood(const element_t & element, level_t m) const;

    /// Mark elements recursively based on their level and neighborhood
    /// @param elements The set of elements to mark.
    /// @param level The level at which to start marking.
    /// @param m The level of the neighborhood to consider.
    /// @return A set of marked elements.
    HElementContainer markRecursive(const HElementContainer & elements, level_t level, level_t m) const;

    /// Mark elements recursively based on their level and neighborhood
    /// @param marked The set of elements to mark.
    /// @param m The level of the neighborhood to consider.
    /// @return A set of marked elements.
    HElementContainer markAdmissible(const HElementContainer & marked, level_t m) const;

    /// @brief Leaf elements of the hierarchical mesh that violate admissibility of class m.
    ///
    /// Let \f$\Omega^0\supseteq\dots\supseteq\Omega^{N-1}\f$ be the hierarchy of
    /// domains. A leaf element \f$Q\f$ is a cell of level \f$k\f$ with
    /// \f$Q\subseteq\Omega^k\f$ and \f$Q\not\subseteq\Omega^{k+1}\f$. A basis function
    /// *acts on* \f$Q\f$ if it is not identically zero on \f$Q\f$. The mesh is
    /// - T-admissible of class m if, on every leaf, the *truncated* basis functions of
    ///   the THB-spline basis acting on the leaf belong to at most m successive levels;
    /// - H-admissible of class m if the same holds for the (untruncated) basis functions
    ///   of the HB-spline basis.
    ///
    /// With \f$\ell_{\min}(Q)\f$ and \f$\ell_{\max}(Q)\f$ the smallest and largest level
    /// of the functions acting on \f$Q\f$, the leaf is returned iff
    /// \f$\ell_{\max}-\ell_{\min}+1>m\f$.
    ///
    /// Since the support of a truncated function is contained in the support of the
    /// function it is truncated from, every T-violation is also an H-violation and
    /// H-admissibility implies T-admissibility on the same mesh.
    ///
    /// The refinement closure (markAdmissible) guarantees a strictly admissible mesh,
    /// which is a sufficient condition. This function tests the definition itself, so it
    /// may accept meshes that the closure would never construct.
    ///
    /// @param m Class of admissibility, m >= 1.
    /// @param type Kind of admissibility: -1 selects automatically by the basis class
    ///        (T for gsTHBSplineBasis, H for gsHBSplineBasis), 0 = T-admissibility
    ///        (only valid for a gsTHBSplineBasis, otherwise GISMO_ERROR), 1 = H-admissibility
    ///        (valid for both basis classes). Matches the marker's "Admissibility" option.
    /// @return The violating leaf elements, each at its own level, patch 0.
    ///
    /// Complexity: O(N_el * L * (p+1)^d * log n_l), with N_el the number of leaves, L the
    /// number of levels and n_l the number of active functions per level (binary search
    /// in the sorted active set); for T plus a merge scan of the sparse presentation of
    /// every truncated function, O(nnz) per truncated function per leaf. Memory is
    /// O(size of the returned set).
    ///
    /// References: A. Buffa, C. Giannelli, Adaptive isogeometric methods with
    /// hierarchical splines: error estimator and convergence, M3AS 26(1):1-25 (2016);
    /// C. Bracco, C. Giannelli, R. Vazquez, Refinement algorithms for adaptive
    /// isogeometric methods with hierarchical splines, Axioms 7(3):43 (2018).
    HElementContainer getNonAdmissibleElements(level_t m, index_t type = -1) const;

    /// @brief Check admissibility of class m of the hierarchical mesh.
    ///
    /// Same criterion as getNonAdmissibleElements, but stops at the first violating leaf.
    /// @param m Class of admissibility, m >= 1.
    /// @param type -1 (automatic), 0 (T-admissibility) or 1 (H-admissibility), see
    ///        getNonAdmissibleElements.
    /// @return True iff no leaf element violates admissibility of class m.
    bool isAdmissible(level_t m, index_t type = -1) const;

    /// Convert an element to refinement box indices
    /// @param element The element to convert.
    /// @param targetLevel The target level for the refinement box.
    /// @param extension If true, the box is extended by floor(p/2) target-level
    ///        knot spans on BOTH sides in every direction (clamped to the
    ///        domain), p being the degree in that direction.
    ///
    /// An interior extended box spans 2 + 2*floor(p/2) >= p + 1 target-level
    /// spans per direction, so it contains the full support of a target-level
    /// B-spline and the refinement of an isolated element always activates new
    /// THB functions. At the boundary the clamped box is shorter, but the open
    /// knot vector's boundary functions have shorter supports. Boxes snap to
    /// whole cells of the leaf they hit, one level at a time
    /// (gsHTree::insertBox), so the region raised to the target level is a
    /// one-cell ring of the element's level for p = 2..5 and a two-cell ring
    /// for p = 6, 7 (no extension for p = 1); a coarser neighbour cell is
    /// first raised to the element's level in full.
    /// The extended box reaches into neighbouring cells; the admissible
    /// closure has to see them, see extensionCells.
    /// @note The box alone, without a closure over the elements and their
    ///       extensionCells, gives no admissibility guarantee. With such a
    ///       closure (markAdmissible), the refined mesh stays admissible under
    ///       the preconditions of gsHElementMarker::markRef(): the mesh was
    ///       produced by gsHElementMarker. On an externally refined mesh, even
    ///       an admissible one, the result can be non-admissible (check with
    ///       isAdmissible()).
    /// @return A vector of indices representing the refinement box at the target level.
    std::vector<index_t> toRefBox(const element_t & element, level_t targetLevel, bool extension = true) const;
    std::vector<index_t> toRefBox(const element_t & element, bool extension = true) const;

    /// Same-level cells reached by the extension of toRefBox of each element,
    /// without the elements themselves.
    /// The cells need not be active. Feeding them to markAdmissible together
    /// with the elements makes the closure cover the whole extended region,
    /// because the closure only uses the level and position of its input cells.
    /// @note The closure of the elements plus their extensionCells under
    ///       markAdmissible keeps the refined mesh admissible under the
    ///       preconditions of gsHElementMarker::markRef(): the mesh was produced
    ///       by gsHElementMarker. On an externally refined mesh, even an
    ///       admissible one, the extended closure can give a non-admissible
    ///       mesh (check with isAdmissible()).
    /// @param elements The elements whose extensions are collected.
    /// @return The set of same-level cells covered partly by the extended boxes.
    HElementContainer extensionCells(const HElementContainer & elements) const;

    /// Convert a set of elements to refinement box indices
    /// @param elements The set of elements to convert.
    /// @param extension If true, elements are extended; see toRefBox.
    /// @return A vector of indices representing the refinement boxes of the elements.
    std::vector<index_t> toRefBoxes(const HElementContainer & elements, bool extension = true) const;

    /// Convert an element to coarsen box indices
    /// @param element The element to convert.
    /// @param targetLevel The target level for the refinement box.
    /// @return A vector of indices representing the coarsen box at the target level.
    std::vector<index_t> toCrsBox(const element_t & element, level_t targetLevel) const;
    std::vector<index_t> toCrsBox(const element_t & element) const;

    /// Convert a set of elements to refinement box indices
    /// @param elements The set of elements to convert.
    /// @return A vector of indices representing the refinement boxes of the elements.
    std::vector<index_t> toCrsBoxes(const HElementContainer & elements) const;

    /**
     * @brief Explodes a box into its elements at a given level
     * This function takes a box and a level, and returns a container of elements
     * that cover the box at the specified level. The elements are created by
     * splitting the box into smaller boxes that correspond to the basis functions
     * at the given level.
     * @param box The box to be exploded.
     * @param level The level at which to explode the box.
     * @param patch The patch ID to assign to the exploded elements (default is 0).
     * @return A container of elements that cover the box at the specified level.
     */
    HElementContainer explode(const point_t & low, const point_t & upp, level_t level, index_t patch = 0) const;

    /// @brief Explode a box defined by its lower and upper corners into a set of elements
    /// @param box  Box of the element, defined by its lower and upper corner in element indices
    /// @param level  Level of the elements in the hierarchy
    /// @param patch  Patch ID of the elements, if applicable
    /// @return A set of elements that cover the box
    HElementContainer explode(const box_t & box, level_t level, index_t patch = 0) const;

    /// @brief Check if element1 contains element2
    /// @param element1 The first element to check.
    /// @param element2 The second element to check.
    /// @return True if element1 contains element2, false otherwise.
    bool contains(const element_t & element1, const element_t & element2) const;

protected:
    /// Resolve the admissibility type (-1, 0, 1) to 0 (T) or 1 (H); raises GISMO_ERROR
    /// when the type is not available for the basis class.
    index_t admissibilityType(index_t type) const;

    /// True iff the functions of the given type acting on the leaf of level \a k with
    /// centre \a c span more than \a m successive levels.
    bool leafViolates(const gsMatrix<T> & c, level_t k, level_t m, index_t type) const;

    const gsHTensorBasis<d,T> & m_basis; ///< The basis of the elements.

};

}; //namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsHElementHelper.hpp)
#endif
