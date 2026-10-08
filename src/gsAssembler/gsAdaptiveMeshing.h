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
#include <gsAssembler/gsAdaptiveMeshingCompare.h>
#include <gsIO/gsOptionList.h>
#include <gsCore/gsMultiPatch.h>
#include <gsCore/gsMultiBasis.h>
#include <gsHSplines/gsHBox.h>
#include <gsHSplines/gsHBoxContainer.h>
#include <gsHSplines/gsHBoxUtils.h>

namespace gismo
{

/**
 * @brief      Provides adaptive meshing routines.
 *
 * Provided element errors, this class performs marking,
 * refinement and coarsening of a provided basis. The class
 * uses the \ref gsHBox and \ref gsHBoxContainer classes
 * to ensure admissible meshing.
 *
 * \deprecated Use \ref gsHElementMarker (header gsHSplines/gsHElementMarker.h,
 * not included by gismo.h) instead: gsHElementMarker::markRef() marks for
 * refinement, gsHElementMarker::toRefBoxes() converts the marked elements to
 * refinement boxes for gsHTensorBasis::refineElements, and
 * gsHElementMarker::markCrs(refined) marks for coarsening, taking the refinement
 * into account. Not yet available in gsHElementMarker: the PBULK marking rule
 * (RefineRule / CoarsenRule = 4), multi-patch bases (gsHElementMarker works on a
 * single basis), and RefineExtension > 0 (the marker only has the switch
 * Extension, i.e. the floor(p/2) extension of RefineExtension = 0).
 *
 * @tparam     T     { description }
 */
template <short_t _dim, class T>
class gsAdaptiveMeshing
{
public:
    typedef          gsHBox<_dim,T>                              HBox;
    typedef          gsHBox<_dim,T>*                              HBox_ptr;
    typedef          gsHBoxContainer<_dim,T>                     HBoxContainer;
    typedef typename HBox::SortedContainer                      boxContainer;
    typedef          std::map<gsHBox<_dim,T>,index_t,gsHBoxCompare<_dim,T>>  indexMapType;
    typedef          std::map<index_t,gsHBox<_dim,T>*>                 boxMapType;
    typedef          gsHBoxUtils<_dim,T>               HBoxUtils;

public:

    GISMO_DEPRECATED gsAdaptiveMeshing();

    GISMO_DEPRECATED gsAdaptiveMeshing(gsFunctionSet<T> & input);

    // ~gsAdaptiveMeshing();

    gsOptionList & options() {return m_options;}

    void defaultOptions();

    void getOptions();

    void rebuild();

    void container_into(const std::vector<T> & elError, HBoxContainer & result);

    /**
     * @brief      Marks elements for refinement.
     *
     * With \c Admissible the marked set contains the admissible closure of
     * every seed (element selected by the marking rule) over the seed and the
     * same-level cells that the extended refinement box of the seed reaches
     * (see \ref gsHBox::extensionCells). The boxes that the closure adds to
     * the seeds are remembered, see \ref refine.
     *
     * @param[in]  elError   The element errors
     * @param[out] elMarked  The marked elements
     */
    void markRef_into(const std::vector<T> & elError, HBoxContainer & elMarked);

    /**
     * @brief      Marks elements for coarsening, taking the refinement into account.
     *
     * No returned element has a parent cell that overlaps, with positive
     * volume, a box that \ref refine(\a markedRef) applies; this holds in
     * admissible and non-admissible mode. The admissibility predicates see
     * \a markedRef together with the cells reached by the extended
     * refinement boxes. The guarantee requires that the state of the last
     * \ref markRef_into (closure-added boxes) and the mesh are unchanged
     * between this call and \ref refine.
     *
     * @param[in]  elError     The element errors
     * @param[in]  markedRef   The elements marked for refinement
     * @param[out] elMarked    The elements marked for coarsening
     */
    void markCrs_into(const std::vector<T> & elError, const HBoxContainer & markedRef, HBoxContainer & elMarked);
    void markCrs_into(const std::vector<T> & elError, HBoxContainer & elMarked);

    void markRef(const std::vector<T> & errors);
    void markCrs(const std::vector<T> & errors);

    /**
     * @brief      Refines the elements in \a markedRef.
     *
     * With \c RefineExtension = 0, the refinement box of every element is
     * extended by floor(p/2) finer-level spans on both sides (clamped to the
     * domain), except for the boxes that the admissible closure of the last
     * \ref markRef_into added to the seeds: those are refined without
     * extension. An isolated seed therefore always activates new
     * functions; a box that was first marked by the closure of another seed
     * stays unextended.
     * With \c RefineExtension = k > 0, every box is extended by k cells of its
     * own level, which the admissible closure does not account for.
     *
     * @param[in]  markedRef  The marked elements
     *
     * @return     True if there was anything to refine
     */
    bool refine(const HBoxContainer & markedRef);
    bool unrefine(const HBoxContainer & markedCrs);

    bool refine(const std::vector<bool> & markedRef) { return refine(_toContainer(markedRef)); }
    bool unrefine(const std::vector<bool> & markedCrs) { return unrefine(_toContainer(markedCrs)); };

    bool refine() { return refine(m_markedRef); }
    bool unrefine() { return unrefine(m_markedRef); };

    bool refineAll();
    bool unrefineAll();

    // void flatten(const index_t level);
    // void flatten() { flatten(m_maxLvl); } ;

    // void unrefineThreshold(const index_t level);
    // void unrefineThreshold(){ unrefineThreshold(m_maxLvl); };

    index_t numBlocked() const;
    index_t numElements() const;

    void assignErrors(const std::vector<T> & elError);
    T blockedError() const;
    T nonBlockedError() const;

private:
    void _makeMap(const gsFunctionSet<T> * input, typename gsAdaptiveMeshing<_dim,T>::indexMapType & indexMap, typename gsAdaptiveMeshing<_dim,T>::boxMapType & boxMap);

    void _assignErrors(boxMapType & container, const std::vector<T> & elError);


    void _refineMarkedElements(     const HBoxContainer & container,
                                    index_t refExtension = 0,
                                    bool extension = true);

    /// Refinement boxes (RefBox format) of patch \a pn that \ref refine applies for \a markedRef
    std::vector<index_t> _refBoxes(const HBoxContainer & markedRef, index_t pn) const;

    /// Removes the elements of \a markedCrs whose parent cell overlaps a refinement box of \a markedRef
    /// whose target level exceeds the coarsening target level (the level of the parent cell). Hence, wherever the
    /// coarsening region of a kept element overlaps a refinement box, the coarsening target level is at least
    /// the refinement target level.
    /// Complexity O(|markedCrs|*|refinement boxes|).
    HBoxContainer _dropRefinementOverlap(const HBoxContainer & markedRef, const HBoxContainer & markedCrs) const;

    void _unrefineMarkedElements(   const HBoxContainer & container,
                                    index_t refExtension = 0);

    // void _flattenElementsToLevel(   const index_t level);

    // void _unrefineElementsThreshold(const index_t level);

    std::vector<index_t> _sortPermutation( const boxMapType & container);
    std::vector<index_t> _sortPermutationProjectedRef( const boxMapType & container);
    std::vector<index_t> _sortPermutationProjectedCrs( const boxMapType & container);
    // void _sortPermutated( const std::vector<index_t> & permutation, boxContainer & container);

    void _crsPredicates_into( std::vector<gsHBoxCheck<_dim,T>*> & predicates);
    void _crsPredicates_into(const HBoxContainer & markedRef, std::vector<gsHBoxCheck<_dim,T>*> & predicates);
    void _refPredicates_into( std::vector<gsHBoxCheck<_dim,T>*> & predicates);

    template<bool _coarsen,bool _admissible>
    void _markElements(  const std::vector<T> & elError, const index_t refCriterion, const std::vector<gsHBoxCheck<_dim,T>*> & predicates, HBoxContainer & elMarked) const;

    template<bool _coarsen,bool _admissible>
    void _markFraction( const boxMapType & elements, const std::vector<gsHBoxCheck<_dim,T>*> predicates, HBoxContainer & elMarked) const
    {
        _markFraction_impl<_coarsen,_admissible>(elements,predicates,elMarked);
    }
    template<bool _coarsen,bool _admissible>
    typename std::enable_if< _coarsen &&  _admissible, void>::type
    _markFraction_impl( const boxMapType & elements, const std::vector<gsHBoxCheck<_dim,T>*> predicates, HBoxContainer & elMarked) const;

    template<bool _coarsen,bool _admissible>
    typename std::enable_if< _coarsen && !_admissible, void>::type
    _markFraction_impl( const boxMapType & elements, const std::vector<gsHBoxCheck<_dim,T>*> predicates, HBoxContainer & elMarked) const;

    template<bool _coarsen,bool _admissible>
    typename std::enable_if<!_coarsen &&  _admissible, void>::type
    _markFraction_impl( const boxMapType & elements, const std::vector<gsHBoxCheck<_dim,T>*> predicates, HBoxContainer & elMarked) const;

    template<bool _coarsen,bool _admissible>
    typename std::enable_if<!_coarsen && !_admissible, void>::type
    _markFraction_impl( const boxMapType & elements, const std::vector<gsHBoxCheck<_dim,T>*> predicates, HBoxContainer & elMarked) const;

    template<bool _coarsen,bool _admissible>
    void _markProjectedFraction( const boxMapType & elements, const std::vector<gsHBoxCheck<_dim,T>*> predicates, HBoxContainer & elMarked) const
    {
        _markProjectedFraction_impl<_coarsen,_admissible>(elements,predicates,elMarked);
    }
    template<bool _coarsen,bool _admissible>
    typename std::enable_if< _coarsen &&  _admissible, void>::type
    _markProjectedFraction_impl( const boxMapType & elements, const std::vector<gsHBoxCheck<_dim,T>*> predicates, HBoxContainer & elMarked) const;

    template<bool _coarsen,bool _admissible>
    typename std::enable_if< _coarsen && !_admissible, void>::type
    _markProjectedFraction_impl( const boxMapType & elements, const std::vector<gsHBoxCheck<_dim,T>*> predicates, HBoxContainer & elMarked) const;

    template<bool _coarsen,bool _admissible>
    typename std::enable_if<!_coarsen &&  _admissible, void>::type
    _markProjectedFraction_impl( const boxMapType & elements, const std::vector<gsHBoxCheck<_dim,T>*> predicates, HBoxContainer & elMarked) const;

    template<bool _coarsen,bool _admissible>
    typename std::enable_if<!_coarsen && !_admissible, void>::type
    _markProjectedFraction_impl( const boxMapType & elements, const std::vector<gsHBoxCheck<_dim,T>*> predicates, HBoxContainer & elMarked) const;


    template<bool _coarsen,bool _admissible>
    void _markPercentage( const boxMapType & elements, const std::vector<gsHBoxCheck<_dim,T>*> predicates, HBoxContainer & elMarked) const
    {
        _markPercentage_impl<_coarsen,_admissible>(elements,predicates,elMarked);
    }
    template<bool _coarsen,bool _admissible>
    typename std::enable_if< _coarsen &&  _admissible, void>::type
    _markPercentage_impl( const boxMapType & elements, const std::vector<gsHBoxCheck<_dim,T>*> predicates, HBoxContainer & elMarked) const;

    template<bool _coarsen,bool _admissible>
    typename std::enable_if< _coarsen && !_admissible, void>::type
    _markPercentage_impl( const boxMapType & elements, const std::vector<gsHBoxCheck<_dim,T>*> predicates, HBoxContainer & elMarked) const;

    template<bool _coarsen,bool _admissible>
    typename std::enable_if<!_coarsen &&  _admissible, void>::type
    _markPercentage_impl( const boxMapType & elements, const std::vector<gsHBoxCheck<_dim,T>*> predicates, HBoxContainer & elMarked) const;

    template<bool _coarsen,bool _admissible>
    typename std::enable_if<!_coarsen && !_admissible, void>::type
    _markPercentage_impl( const boxMapType & elements, const std::vector<gsHBoxCheck<_dim,T>*> predicates, HBoxContainer & elMarked) const;

    template<bool _coarsen,bool _admissible>
    void _markThreshold( const boxMapType & elements, const std::vector<gsHBoxCheck<_dim,T>*> predicates, HBoxContainer & elMarked) const
    {
        _markThreshold_impl<_coarsen,_admissible>(elements,predicates,elMarked);
    }
    template<bool _coarsen,bool _admissible>
    typename std::enable_if< _coarsen &&  _admissible, void>::type
    _markThreshold_impl( const boxMapType & elements, const std::vector<gsHBoxCheck<_dim,T>*> predicates, HBoxContainer & elMarked) const;

    template<bool _coarsen,bool _admissible>
    typename std::enable_if< _coarsen && !_admissible, void>::type
    _markThreshold_impl( const boxMapType & elements, const std::vector<gsHBoxCheck<_dim,T>*> predicates, HBoxContainer & elMarked) const;

    template<bool _coarsen,bool _admissible>
    typename std::enable_if<!_coarsen &&  _admissible, void>::type
    _markThreshold_impl( const boxMapType & elements, const std::vector<gsHBoxCheck<_dim,T>*> predicates, HBoxContainer & elMarked) const;

    template<bool _coarsen,bool _admissible>
    typename std::enable_if<!_coarsen && !_admissible, void>::type
    _markThreshold_impl( const boxMapType & elements, const std::vector<gsHBoxCheck<_dim,T>*> predicates, HBoxContainer & elMarked) const;

    bool _checkBox  ( const          HBox            & box  , const std::vector<gsHBoxCheck<_dim,T>*> predicates) const;
    bool _checkBoxes( const typename HBox::Container & boxes, const std::vector<gsHBoxCheck<_dim,T>*> predicates) const;

    T _totalError(const boxMapType & elements);

    T _maxError(  const boxMapType & elements);

    void _addAndMark(          HBox            & box  , HBoxContainer & elMarked) const;
    void _addAndMark( typename HBox::Container & boxes, HBoxContainer & elMarked) const;

    void _setContainerProperties( typename HBox::Container & boxes ) const;

    HBox * _boxPtr(const HBox & box) const;

    typename gsAdaptiveMeshing<_dim,T>::HBoxContainer _toContainer( const std::vector<bool> & bools) const;

protected:
    // M & m_basis;
    gsFunctionSet<T> * m_input;
    // const gsMultiPatch<T> & m_patches;
    gsOptionList m_options;

    T               m_crsParam, m_crsParamExtra, m_refParam, m_refParamExtra;
    MarkingStrategy m_crsRule, m_refRule;
    index_t         m_crsExt, m_refExt;
    index_t         m_maxLvl;

    index_t         m_alpha, m_beta;

    index_t m_m;

    bool            m_admissible;

    index_t         m_verbose;

    HBoxContainer m_markedRef, m_markedCrs;
    /// Seeds (elements selected by the marking rule) of the last admissible \ref markRef_into
    mutable HBoxContainer m_refSeeds;
    /// Elements of the last \ref markRef_into that the admissible closure added to the seeds
    HBoxContainer m_closureAdded;
    // m_boxes is a container that does not contain patch IDs

    indexMapType m_indices;
    boxMapType   m_boxes;

    T m_totalError, m_maxError, m_uniformRefError, m_uniformCrsError;

    std::vector<index_t> m_refPermutation, m_crsPermutation;

    // std::map<index_t,std::shared_ptr<gsHBox<_dim,T>>> m_toindices;
    // std::map<std::shared_ptr<gsHBox<_dim,T>>,index_t> m_fromindices;


};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsAdaptiveMeshing.hpp)
#else
#ifdef gsAdaptiveMeshing_EXPORT
#include GISMO_HPP_HEADER(gsAdaptiveMeshing.hpp)
#undef  EXTERN_CLASS_TEMPLATE
#define EXTERN_CLASS_TEMPLATE CLASS_TEMPLATE_INST
#endif
namespace gismo
{
    EXTERN_CLASS_TEMPLATE gsAdaptiveMeshing<1,real_t>;
    EXTERN_CLASS_TEMPLATE gsAdaptiveMeshing<2,real_t>;
    EXTERN_CLASS_TEMPLATE gsAdaptiveMeshing<3,real_t>;
    EXTERN_CLASS_TEMPLATE gsAdaptiveMeshing<4,real_t>;
    // EXTERN_CLASS_TEMPLATE gsHBoxCheck<1,real_t>;
    // EXTERN_CLASS_TEMPLATE gsHBoxCheck<2,real_t>;
    // EXTERN_CLASS_TEMPLATE gsHBoxCheck<3,real_t>;
    // EXTERN_CLASS_TEMPLATE gsHBoxCheck<4,real_t>;

    // EXTERN_CLASS_TEMPLATE gsLvlCompare<1,real_t>;
    // EXTERN_CLASS_TEMPLATE gsLvlCompare<2,real_t>;
    // EXTERN_CLASS_TEMPLATE gsLvlCompare<3,real_t>;
    // EXTERN_CLASS_TEMPLATE gsLvlCompare<4,real_t>;

    // EXTERN_CLASS_TEMPLATE gsSmallerErrCompare<1,real_t>;
    // EXTERN_CLASS_TEMPLATE gsSmallerErrCompare<2,real_t>;
    // EXTERN_CLASS_TEMPLATE gsSmallerErrCompare<3,real_t>;
    // EXTERN_CLASS_TEMPLATE gsSmallerErrCompare<4,real_t>;

    // EXTERN_CLASS_TEMPLATE gsLargerErrCompare<1,real_t>;
    // EXTERN_CLASS_TEMPLATE gsLargerErrCompare<2,real_t>;
    // EXTERN_CLASS_TEMPLATE gsLargerErrCompare<3,real_t>;
    // EXTERN_CLASS_TEMPLATE gsLargerErrCompare<4,real_t>;

    // EXTERN_CLASS_TEMPLATE gsOverlapCompare<1,real_t>;
    // EXTERN_CLASS_TEMPLATE gsOverlapCompare<2,real_t>;
    // EXTERN_CLASS_TEMPLATE gsOverlapCompare<3,real_t>;
    // EXTERN_CLASS_TEMPLATE gsOverlapCompare<4,real_t>;
}
#endif